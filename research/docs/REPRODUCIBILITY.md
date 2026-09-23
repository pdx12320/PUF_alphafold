# Reproducibility

Run the commands below from the repository root. Published tables are retained at their original precision. The reorganization changes navigation and documentation without retraining models.

## Environment and snapshot checks

```bash
python -m pip install -r research/requirements.txt
python research/tools/validate_snapshot.py
```

The snapshot validator checks file hashes, 24 unique PUF12 constructs, four work labels, prediction coverage and equivalent architecture/source-order holdouts. The current manifest uses repository-relative paths. Original input provenance records retain their historical source paths and hashes.

## Current scaffold model

```bash
python research/results_20260922/scaffold_RF/retrain.py
```

Use the [recorded environment](../results_20260922/scaffold_RF/audit/environment.json) to reproduce saved model results. This command overwrites the corresponding scaffold outputs. It does not retrain the archived TRM models.

## Historical architecture analysis

```bash
python research/tools/reproduce_current.py
```

This runs `stats.py`, `models.py`, `plot_static.py`, `finish.py` and `write_report.py` within the archived architecture workflow. It reads retained feature and mapping caches. Historical report generators preserve their original language. Use an isolated checkout for reruns because scripts overwrite their outputs.

## Wiki figures

```bash
python wiki/scripts/make_figures.py
python wiki/scripts/add_wiki_figures.py
```

These scripts read recorded values without model fitting. The first regenerates the endpoint confusion and non-TRM support figures, along with archived standalone panels. The second generates the combined scaffold validation and S12 distribution figures. Only the four figures listed in [current captions](../../wiki/figures/captions.md) appear in the Wiki.

## Additional analyses

- [C388 threshold analysis](../c388_threshold/REPRODUCE.md)
- [C295/C871 frozen-cohort analysis](../trm22_offtarget/README.md)
- [Local-contact comparison snapshot](../c388_local_nonlocal_optimization/README.md)
- [ProteinMPNN and LigandMPNN workflow](../aice_mpnn_20260922/README.md)

Large raw AF3 archives are retained separately and identified by the source manifests. Verify archive integrity before extraction. Model/seed predictions are technical replicates; reported generalization uses construct-level held-out predictions. Classifier scores are uncalibrated, and all-data training scores do not replace held-out evaluation.
