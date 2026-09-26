# Methods and interpretation notes

## Data scope

The delivered analysis contains 168 experimental rows and 84 named groups. The main modeling cohort contains 81 mutant constructs with verified distinct PUF modeling-core sequences, plus three control groups.

Main-analysis valid labels:
- C295: 66
- C388: 79
- C871: 81
- Complete three-site specificity score S: 66

Records marked as not detected without independent coverage evidence were treated as missing in the main analysis. A separate zero-value sensitivity analysis was retained in the original delivery.

## Evaluation

The principal comparison uses construct-level out-of-fold predictions. Model selection uses the mean MAE across C295, C388 and C871. Ranking uses:

`S = C388 - (C295 + C871) / 2`

The selected quantitative model is M3a_global_Ridge. Ranking performance is reported separately because the model minimizing three-site MAE is not necessarily the model maximizing rank correlation.

## Baselines

- B0: no-feature mean baseline.
- B1_Ridge: sequence/mutation descriptor baseline used in the delivered analysis.
- M0_Ridge / M0_RF: CP-based baselines.
- M1/M3a families: pretrained protein-language-model representations and combinations described in the original analysis package.

## Work classifier

Work is defined as C388 > 0.50. Work classification is an independent logistic-regression comparison and should not be interpreted as a low-bystander or overall-good label.

## Structure subset

Only 20 constructs had matched coordinates usable for the structure/SaProt comparison. These constructs come from the first batch. Structure-subset performance therefore has a different cohort from the full 81-construct sequence analysis.

## Interpretation boundary

All results are retrospective small-sample evaluations. Fixed-OOF bootstrap intervals condition on the saved OOF predictions and do not include uncertainty from repeating model-family selection or retraining. Cross-batch performance is affected by the fact that batch and mutated repeat position are partially confounded.

The existing GitHub Model 2 fixed-five endpoint classification and this quantitative PLM analysis use different targets and evaluation protocols, so their numerical metrics should not be treated as a direct head-to-head benchmark.
