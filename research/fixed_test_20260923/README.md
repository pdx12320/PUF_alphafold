# Fixed test panels for Models 1 and 2

**C388 status:** the three-class reconstruction below is retained as historical evidence. The [current C388 binary model and nested validation](../c388_binary/README.md) supersede its C388 presentation; C295/C871 records here remain unchanged.

[Project home](../../README.md) · [Wiki](../../wiki/PUF_Model_EN.md)

## Experimental chronology and computational reconstruction

The project owner clarified on 23 September 2026 that the four scaffold designs and five TRM candidates were intended as test-only constructs and that their wet-lab assays followed the original model-development and prediction stage. This report follows that experimental sequence.

The numerical results in this directory are a **23 September reconstruction with fixed test membership**, using the now-available measurements. This rerun does not establish the date of the original model freeze or recover a timestamped pre-assay prediction file. Test labels were masked during the reconstructed Model 2 search and fitting; predictions were written before scoring. Original cross-validation results remain available as development history, separately from the fixed-test metrics.

## Model 1: twenty training constructs and four test designs

The test designs are 1, 3, 7 and 8. All four are excluded together from preprocessing and fitting. The existing nine-feature random-forest recipe and 0.5 cutoff are retained without tuning on this panel. Agreement remains **3/4**: one true positive, two true negatives and one false positive (Design 7).

See [split and recipe](model1/split.json), [training inputs](model1/training_data.csv) and [test predictions](model1/predictions.csv). Run `python research/fixed_test_20260923/retrain_model1.py` from the repository root.

## Model 2: seventy-six training constructs and five test candidates

Test membership is fixed to P9-GNS, P9-NPS, P9-NTQ, P8-GVE and P4-R5-SNE+P7-R5-SNE. All matching sequence identities are excluded simultaneously. The combined protein–RNA cohort contains 76 training constructs and five test constructs. The input matrices are already construct-level aggregates, so individual seeds cannot cross the split.

The v4 search recipe is reused. Each endpoint selects its feature representation, classifier and biological class thresholds using only grouped internal validation among the training constructs. Test outcomes are masked throughout this process. Internal fold membership is retained in the endpoint `_inner_splits.json` files. Final models are fitted to all 76 training constructs; sealed predictions are then compared with the five measured outcomes.

| Endpoint | Test agreement | Training-majority reference | Selected biological boundary |
|---|---:|---:|---|
| C295 | 5/5 | 5/5 | Decrease ≤ −15 pp |
| C388 | 3/5 | 0/5 | Decrease ≤ −12.5 pp; increase ≥ +2.5 pp |
| C871 | 3/5 | 5/5 | Decrease ≤ −7.5 pp |

The aggregate **11/15** counts correlated endpoint calls. C295 has only no-decrease test labels at its chosen threshold; C871 has only decrease test labels. Report the counts and endpoint-specific errors rather than interpreting the aggregate as a broad generalization estimate. The previous 13/15 result used individual-construct folds and is not a result from this fixed-five split.

## Files and reproduction

The [combined prediction table](model2/predictions.csv) contains all fifteen calls. Endpoint `_selected.json` files contain the selected model and inner-validation metrics; `_search.json.gz` files retain every candidate search record. Protocol files preserve the complete source configuration, source-code hashes and train/test identities. `_predictions_sealed.json` files contain predictions saved before outcome comparison. The `.joblib` files preserve the fitted models and their transformation states.

Recover the verified `PUF_classification_v4_dynamic_fullCP_20260921_complete.zip` archive identified in [source_archives.json](../results_20260922/source_archives.json). Keep the extracted repository layout because the v4 source also reads its sibling v3 feature files. For each endpoint, run:

```bash
python research/fixed_test_20260923/retrain_model2.py \
  --v4 /path/to/PUF_alphafold/retrain_classification_v4_dynamic_fullCP_20260921 \
  --endpoint C295
```

Repeat with `C388` and `C871`, then concatenate the three endpoint prediction CSVs. No earlier fitted model or prediction table is read by the training script. The underlying representations and search recipe have historical development exposure; this computational reconstruction cannot remove that history. Later training updates that incorporate these test outcomes require a new test panel.
