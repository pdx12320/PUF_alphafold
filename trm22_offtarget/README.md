# TRM22 C871 / C295 dynamic classification

Read REPORT.md for the full Chinese report.

Reproduce from this standalone folder:

```bash
python -m pip install -r requirements.txt
python run.py
python report.py
```

Default: four CPU processes; override with TRM_JOBS. Python 3.12 used. No AlphaFold rerun required. Training overwrites results.

Key outputs:
- results/label_scan.csv: experimental label scan and class counts
- results/group_selected_candidate_comparison.csv: same group-selected candidate across all validations
- results/strict_selected_candidate_comparison.csv: same strict-position-selected candidate across all validations
- results/performance.csv: all 5,115 metric rows
- results/outer_predictions.csv: all 119,350 fixed-label validation predictions
- results/joint_label_performance.csv: fully training-selected labels/features/models/score thresholds; endpoints vary by fold
- results/joint_label_selection.csv and joint_label_predictions.csv: adaptive choices and test outcomes
- results/split_audit.json and validation_checks.json: split membership and independent checks
- SOURCE_MANIFEST.json: immutable input provenance

Best observed configurations remain exploratory. No deployment model is claimed.
