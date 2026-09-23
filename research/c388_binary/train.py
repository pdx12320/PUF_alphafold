#!/usr/bin/env python3
"""Binary C388 search using verified v4 CP preprocessing; never rewrites v4."""
import argparse
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

CLASSES = ['decrease', 'unchanged']
TARGETS = {'25_p9_gns', '30_p9_nps', '32_p9_ntq',
           'p8_r6_gve', 'p4_r5_sne_plus_p7_r5_sne'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binary_labels(d, task, a, b=0):
    import numpy as np
    if task != 'C388':
        raise ValueError('This runner is C388-only')
    z = d.C388_delta_pp.to_numpy(float)
    # Unavailable outcomes stay unavailable; never silently become unchanged.
    return np.where(np.isfinite(z), (z >= -a).astype(int), -1)


def install_binary(core):
    core.names = lambda task: CLASSES if task == 'C388' else (_ for _ in ()).throw(ValueError(task))
    core.labels = binary_labels


def family(c):
    if 'interface' in c['family']:
        return 'RNA_protein_interface'
    return 'all_PR' if c['retain'] == 'all' else 'high_variance_PR'


def save(out, name, obj):
    text = json.dumps(obj, indent=2, allow_nan=False,
                      default=lambda v: v.item() if hasattr(v, 'item') else str(v))
    p = out / name
    if p.exists():
        raise FileExistsError(p)
    if name.endswith('.gz'):
        p.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        p.write_text(text + '\n')


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--v4', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True, help='Must not exist')
    p.add_argument('--mode', choices=['audit', 'pilot', 'fixed-five'], default='audit')
    args = p.parse_args()
    provenance = json.loads((Path(__file__).resolve().parents[1] /
                            'fixed_test_20260923/model2/C388_protocol.json').read_text())
    for rel, expected in provenance['source_sha256'].items():
        if digest(args.v4 / rel) != expected:
            raise ValueError(f'Unexpected archived source: {rel}')
    if args.output.exists():
        raise FileExistsError('Use a new output directory')
    sys.path.insert(0, str(args.v4 / 'src'))
    core = importlib.import_module('core')
    install_binary(core)
    train = importlib.import_module('train')  # imports the binary functions
    import numpy as np
    import pandas as pd
    import sklearn
    import scipy
    import joblib

    ds = core.Dataset('C388', 'combined_PR')
    d = ds.d
    if len(d) != 81 or d.construct_id.nunique() != 81:
        raise ValueError('Expected frozen 81-construct cohort')
    if not np.isfinite(d.C388_delta_pp).all():
        raise ValueError('Missing experimental outcome')
    if not TARGETS.issubset(set(d.construct_id)):
        raise ValueError('Missing fixed-test construct')
    held = set(d.loc[d.construct_id.isin(TARGETS), 'sequence_identity'])
    te = np.flatnonzero(d.sequence_identity.isin(held))
    tr = np.flatnonzero(~d.sequence_identity.isin(held))
    if set(d.iloc[tr].construct_id) != set(provenance['train']) or len(te) != 5:
        raise ValueError('Fixed split differs from repository provenance')
    observed = d.iloc[te].copy()
    outcome_cols = [c for c in d if c.endswith(('_fraction', '_delta_pp'))]
    d.loc[te, outcome_cols] = np.nan
    cfg = core.CFG
    # A fixed, moderate grid drawn from the archived protocol. No TRM sequence inputs.
    reps = [c for c in cfg['representations']
            if c['family'] in ['full_PR', 'PR_reference_interface', 'PR_train_union_interface']
            and c['mode'] in ['raw', 'delta'] and c['ranking'] == 'variance'
            and c['pca'] == 0 and not c['cp_trm']
            and c['retain'] in ['all', 25, 100, 250, 1000]]
    thresholds = cfg['grids']['C388']
    seed = core.seed('C388-binary-fixed-five', 20260923)
    av = pd.read_csv(args.v4 / 'audit/matrix_availability.csv')
    av = av[av.endpoint == 'C388']
    args.output.mkdir(parents=True)
    protocol = dict(mode=args.mode, classes=CLASSES,
                    label_rule='delta_pp < -a: decrease; delta_pp >= -a: unchanged (includes increase)',
                    candidates_a_pp=thresholds, primary_metric=cfg['primary_metric']['C388'],
                    decision_rule='native classifier boundary; probability 0.5 or margin 0',
                    seed=seed, train=d.iloc[tr].construct_id.tolist(),
                    test=d.iloc[te].construct_id.tolist(), test_outcomes_masked=True,
                    batch_counts=d.batch.value_counts().to_dict(),
                    availability=av.groupby(['batch', 'PP_available', 'PR_available', 'RR_available']).size().reset_index(name='n').to_dict('records'),
                    representations=reps, screening_model=cfg['screening_classifier'],
                    refinement_models=cfg['refinement_models'],
                    refinement='top two screened configurations per feature family; max six total',
                    source_sha256={rel: digest(args.v4 / rel) for rel in provenance['source_sha256']},
                    input_sha256={str(q.relative_to(args.v4)): digest(q) for q in
                                  sorted((args.v4 / 'cp_features').glob('C388_*.npy'))},
                    runner_sha256=digest(__file__),
                    versions=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__,
                                  scipy=scipy.__version__, sklearn=sklearn.__version__, joblib=joblib.__version__))
    save(args.output, 'protocol.json', protocol)
    if args.mode == 'audit':
        print(json.dumps({k: protocol[k] for k in ['batch_counts', 'availability']}))
        print('screen configurations', len(reps) * len(thresholds))
        return
    if args.mode == 'pilot':
        # Small, stratified training subset only; no fixed-test predictions.
        yy = binary_labels(d, 'C388', 12.5)
        tr = np.sort(np.concatenate([tr[yy[tr] == k][:12] for k in [0, 1]]))
        thresholds = [12.5]
        reps = [next(c for c in reps if family(c) == f and c['mode'] == 'raw')
                for f in ['all_PR', 'high_variance_PR', 'RNA_protein_interface']]
    start = time.monotonic()
    records, audits, feasible = [], [], {}
    for a in thresholds:
        y = binary_labels(d, 'C388', a)
        support = np.bincount(y[tr], minlength=2)
        if support.min() < cfg['min_outer_class']['C388']:
            audits.append(dict(a=a, status='insufficient class support', support=support.tolist()))
            continue
        inner, kind = core.inner_splits(d, tr, y, 'LOCO', seed)
        if not inner:
            audits.append(dict(a=a, status=kind))
            continue
        for aa, bb, _ in inner:
            assert not set(aa) & set(te) and not set(bb) & set(te)
            assert not set(d.iloc[aa].sequence_identity) & set(d.iloc[bb].sequence_identity)
        feasible[a] = inner
        audits.append(dict(a=a, status='ok', support=support.tolist(), kind=kind,
                           folds=[dict(train=d.iloc[aa].construct_id.tolist(),
                                       validation=d.iloc[bb].construct_id.tolist()) for aa, bb, _ in inner]))
        for c in reps:
            rec = dict(a=a, b=0, representation=c, family=family(c),
                       model=cfg['screening_classifier'], stage='screen', status='invalid')
            try:
                metrics, counts = train.evaluate(ds, c, rec['model'], y, inner, seed)
                rec.update(status='ok', metrics=metrics, actual_counts=counts)
            except (ValueError, np.linalg.LinAlgError) as e:
                rec['reason'] = str(e)
            records.append(rec)
        train.clear_fitted(ds)
        print('screened', a, len(records), flush=True)
    rank = lambda r: train.rank(r, 'C388')
    shortlist = []
    for f in sorted({family(c) for c in reps}):
        shortlist += sorted([r for r in records if r['status'] == 'ok' and r['family'] == f], key=rank)[:2]
    models = cfg['refinement_models'] if args.mode != 'pilot' else [cfg['refinement_models'][3], cfg['refinement_models'][8]]
    for base in shortlist:
        y = binary_labels(d, 'C388', base['a'])
        for mc in models:
            rec = dict(a=base['a'], b=0, representation=base['representation'], family=base['family'],
                       model=mc, stage='refine', status='invalid')
            try:
                metrics, counts = train.evaluate(ds, rec['representation'], mc, y, feasible[rec['a']], seed)
                rec.update(status='ok', metrics=metrics, actual_counts=counts)
            except (ValueError, np.linalg.LinAlgError) as e:
                rec['reason'] = str(e)
            records.append(rec)
        train.clear_fitted(ds)
    save(args.output, 'inner_splits.json', audits)
    save(args.output, 'search.json.gz', records)
    successful = [r for r in records if r['status'] == 'ok']
    if not successful:
        raise RuntimeError('No valid candidates; see search records')
    best = min(successful, key=rank)
    save(args.output, 'selected.json', best)
    pd.DataFrame([dict(family=r['family'], threshold_pp=-r['a'], stage=r['stage'],
                       representation=r['representation']['id'], model=json.dumps(r['model']),
                       **r['metrics']) for r in successful]).to_csv(args.output / 'leaderboard.csv', index=False)
    if args.mode == 'pilot':
        save(args.output, 'pilot.json', dict(n_train=len(tr), elapsed_seconds=time.monotonic()-start,
                                           scientific_result=False, fixed_test_evaluated=False))
        print('Pilot completed; this is a software check, not a model benchmark', flush=True)
        return
    y = binary_labels(d, 'C388', best['a'])
    model, state, z, embedding = core.model_fit(ds, best['representation'], best['model'], tr, y, seed)
    if state['bad'][te].any():
        raise ValueError('Fixed test has missing selected CP coordinates')
    pred, scores, kind = core.model_scores(model, z[te])
    sealed = [dict(construct_id=d.iloc[i].construct_id, predicted_class=CLASSES[int(q)],
                   scores=sc.tolist(), score_kind=kind) for i, q, sc in zip(te, pred, scores)]
    save(args.output, 'predictions_sealed.json', sealed)
    raw_z = ds.original_z(state, best['representation'])
    basis = raw_z[tr].T @ state['projection'] if embedding else None
    artifact = dict(classifier=model, representation=best['representation'],
                    state={k: v for k, v in state.items() if k not in ['embedding', 'projection', 'bad']},
                    basis=basis, classes=CLASSES, threshold_pp=-best['a'], embedding=embedding,
                    reference=ds.base(best['representation'])[5][state['pick']],
                    feature_ids=core.PAIR.iloc[state['selected']].feature_id.tolist(),
                    training_ids=d.iloc[tr].construct_id.tolist())
    joblib.dump(artifact, args.output / 'model.joblib', compress=3)
    core.PAIR.iloc[state['selected']].to_csv(args.output / 'selected_CP_pairs.csv', index=False)
    actual = binary_labels(observed, 'C388', best['a'])
    majority = int(np.bincount(y[tr], minlength=2).argmax())
    metrics = core.scores_metric(actual, pred, scores, 'C388')
    save(args.output, 'evaluation.json', dict(metrics=metrics, n_train=len(tr), n_test=len(te),
         majority_correct=int(sum(actual == majority)), selected_threshold_pp=-best['a'],
         sealed_sha256=digest(args.output / 'predictions_sealed.json'),
         evaluated_after_seal=True, elapsed_seconds=time.monotonic()-start,
         limitation='Historical fixed-panel reconstruction; selection metrics are internal, not independent validation.'))
    rows = [dict(**s, true_class=CLASSES[int(a)], delta_pp=float(delta), correct=bool(q == a))
            for s, a, delta, q in zip(sealed, actual, observed.C388_delta_pp, pred)]
    pd.DataFrame(rows).to_csv(args.output / 'predictions.csv', index=False)
    print('Completed', json.dumps(metrics), flush=True)


if __name__ == '__main__':
    main()
