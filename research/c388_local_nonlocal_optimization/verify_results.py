"""Verify the curated C388 snapshot from saved predictions; no training or writes."""
from __future__ import annotations
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONFIGS = {
    'interface_rf_fixed': ('Interface_CP', 'RF_depth2', 'fixed_0.5'),
    'total_lr01_fixed': ('Total_CP', 'LR_C0.1', 'fixed_0.5'),
    'nested_all': ('nested_best_available', 'nested_model_selection', 'training_only_threshold'),
}

def read_csv(name: str) -> list[dict[str, str]]:
    with (ROOT / 'results' / name).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))

def main() -> None:
    cohort = {r['construct_key']: r for r in read_csv('frozen_cohort.csv')}
    if len(cohort) != 22:
        raise ValueError('Expected 22 frozen constructs')
    metrics = {(r['feature_set'], r['model'], r['score_policy'], r['validation']): r
               for r in read_csv('comparison_metrics.csv')}
    groups = defaultdict(list)
    for row in read_csv('key_predictions.csv'):
        groups[(row['config_id'], row['validation'])].append(row)
    if len(groups) != 6:
        raise ValueError('Expected three configurations and two validations')
    max_error = 0.0
    summaries = []
    for (config_id, validation), rows in groups.items():
        if len(rows) != 22 or {r['construct_key'] for r in rows} != set(cohort):
            raise ValueError('Missing or duplicated construct')
        y = [int(r['actual_class']) for r in rows]
        pred = [int(r['predicted_class']) for r in rows]
        scores = [float(r['predicted_score']) for r in rows]
        for row, label, prediction in zip(rows, y, pred):
            original = cohort[row['construct_key']]
            if label != int(float(original['delta_editing_pp']) >= -15):
                raise ValueError('Endpoint mismatch')
            if row['mutation_group'] != original['mutation_group']:
                raise ValueError('Group mismatch')
            if prediction != int(float(row['predicted_score']) >= float(row['score_threshold'])):
                raise ValueError('Score threshold mismatch')
        tp = sum(a == 1 and b == 1 for a, b in zip(y, pred))
        tn = sum(a == 0 and b == 0 for a, b in zip(y, pred))
        fp = sum(a == 0 and b == 1 for a, b in zip(y, pred))
        fn = sum(a == 1 and b == 0 for a, b in zip(y, pred))
        pos = [s for a, s in zip(y, scores) if a == 1]
        neg = [s for a, s in zip(y, scores) if a == 0]
        denominator = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        values = {
            'accuracy': (tp+tn)/len(y),
            'balanced_accuracy': (tp/(tp+fn)+tn/(tn+fp))/2,
            'MCC': (tp*tn-fp*fn)/denominator if denominator else 0.0,
            'AUC': sum((p > n)+0.5*(p == n) for p in pos for n in neg)/(len(pos)*len(neg)),
            'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn,
        }
        reference = metrics[(*CONFIGS[config_id], validation)]
        for name, value in values.items():
            error = abs(value-float(reference[name]))
            max_error = max(max_error, error)
            if error > 1e-12:
                raise ValueError(f'Metric mismatch: {config_id}/{validation}/{name}')
        summaries.append({'config_id': config_id, 'validation': validation, **values})
    print(json.dumps({'verified': True, 'prediction_rows': 132,
                      'max_absolute_error': max_error, 'results': summaries,
                      'scope': 'Saved key predictions only; no new model fitting or PP analysis'}, indent=2))

if __name__ == '__main__':
    main()
