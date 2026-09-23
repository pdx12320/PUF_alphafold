# Joint C388 label and classifier threshold analysis

The analysis uses 31 independent constructs, with C388 labels defined from the mean of two experimental records per construct. Structural models and seeds are aggregated before evaluation. Changing the editing threshold changes the prediction task.

The search covers C388 thresholds from 5% to 95% in 5-percentage-point increments and adjacent observed-value midpoints. Equivalent label partitions are evaluated once; 28 partitions have at least two constructs in both classes. CP3 and CP+Structure9 features are compared using fixed logistic-regression and shallow-random-forest recipes under LOCO and mutation-group-out validation. Grouping by mutation-position combinations can retain shared positions between groups.

The exploratory group-out candidate uses C388 ≥40% and CP+Structure9 logistic regression with a score cutoff near 0.46. Its labels contain 23 work and eight non-work examples. After choosing the classifier cutoff inside training folds, group-out accuracy is 77.42% and balanced accuracy is 80.71%. The corresponding LOCO values are 64.52% and 67.93%, indicating sensitivity to threshold selection. The search does not establish a prospectively validated deployment threshold.

Logistic regression uses C=0.1, balanced class weights and liblinear. Random forests use 100 trees, depth 2, minimum leaf size 3, max_features=1.0, balanced class weights and seed 2026. See [reproduction instructions](REPRODUCE.md).

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/c388_threshold/README.md) preserves the historical prose and analysis context.

### Table 1. Exploratory threshold comparisons

| Selection objective | Validation | C388 threshold | Model | Score cutoff | Accuracy | Balanced accuracy |
|---|---|---:|---|---:|---:|---:|
| Highest accuracy | LOCO | 30% | CP+Structure9 LR | 0.376386 | 90.32% | 70.00% |
| Highest balanced accuracy | LOCO | 35% | CP3 LR | 0.453027 | 83.87% | 79.46% |
| Highest balanced accuracy | mutation-group-out | 40% | CP+Structure9 LR | 0.458533 | 80.65% | 82.88% |
