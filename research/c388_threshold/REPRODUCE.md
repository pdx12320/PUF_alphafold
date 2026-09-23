# Reproduce joint threshold analyses

Run from the repository root:

```bash
python -m pip install -r research/c388_threshold/requirements.txt
python research/c388_threshold/thresholds.py
python research/c388_threshold/joint.py
python research/c388_threshold/validate_joint.py
```

`thresholds.py` selects classifier cutoffs within training folds for the fixed 50% experimental label. `joint.py` searches biological and classifier thresholds. `validate_joint.py` performs training-only score-threshold selection for the 30%, 35%, 40%, 50% and 55% label candidates. Commands overwrite their corresponding outputs.

The retained inputs support 28 distinct label partitions, four fixed models and two validation schemes. Global exploratory selection of labels and models remains distinct from training-fold cutoff selection. Mutation-combination holdouts can share component positions. No independent prospective validation or deployment weights are supplied for the exploratory 40% candidate.

Input migration checks are in `inputs/provenance.json`; original source inventory is in `inputs/raw_input_manifest.json`. Earlier fixed-label workflows are documented in the [engineering history](../docs/DRY_LAB_DBTL.md).
