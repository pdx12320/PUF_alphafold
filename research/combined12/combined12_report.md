# Combined PUF12 scaffold analysis

The combined collection contains 24 independent 12-repeat constructs: four work and 20 non-work, represented by 96 seed tasks and 480 structure models. Each seed is averaged before equal-weight aggregation across seeds. The input RNA remains AUGGAGGACGUGC (13 nucleotides). Two PUF11 examples and three 15-repeat designs were excluded from this PUF12 comparison.

Fixed feature sets comprise three nonlocal contact densities, 14 contact summaries, and 17 contact-plus-structure features. All preprocessing is fitted within each training fold. LOCO withholds every structure and seed belonging to one construct. A common 0.5 score cutoff is used; final all-data fits are separate from held-out predictions.

The three-density random forest achieves AUC 0.9125 and 21/24 correct classifications (balanced accuracy 0.825). The 17-feature forest achieves AUC 0.9250 and 20/24 correct, with lower recall. Their ranking and classification trade-offs differ. Predicting every construct as non-work would already yield 20/24 accuracy. The tables retain all recorded predictions and comparisons.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/combined12/combined12_report.md) preserves the historical prose and analysis context.

### Table 1. Combined LOCO evaluation

| model | AUC | AP | accuracy | balanced_accuracy | precision | recall | F1 | TP | FP | FN | TN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3_LR | 0.9125 | 0.5845 | 0.6250 | 0.7750 | 0.3077 | 1.0000 | 0.4706 | 4 | 9 | 0 | 11 |
| CP_density3_RF | 0.9125 | 0.5845 | 0.8750 | 0.8250 | 0.6000 | 0.7500 | 0.6667 | 3 | 2 | 1 | 18 |
| CP_density3_XGB | 0.9000 | 0.6012 | 0.7917 | 0.6750 | 0.4000 | 0.5000 | 0.4444 | 2 | 3 | 2 | 17 |
| CP_density3_Stump | 0.7438 | 0.4167 | 0.8333 | 0.8000 | 0.5000 | 0.7500 | 0.6000 | 3 | 3 | 1 | 17 |
| CP_summary14_LR | 0.8250 | 0.4132 | 0.6250 | 0.7750 | 0.3077 | 1.0000 | 0.4706 | 4 | 9 | 0 | 11 |
| CP_summary14_RF | 0.8500 | 0.4444 | 0.8333 | 0.7000 | 0.5000 | 0.5000 | 0.5000 | 2 | 2 | 2 | 18 |
| CP_summary14_XGB | 0.8875 | 0.5179 | 0.7917 | 0.6750 | 0.4000 | 0.5000 | 0.4444 | 2 | 3 | 2 | 17 |
| CP_structure17_LR | 0.8250 | 0.5769 | 0.6250 | 0.6750 | 0.2727 | 0.7500 | 0.4000 | 3 | 8 | 1 | 12 |
| CP_structure17_RF | 0.9250 | 0.7679 | 0.8333 | 0.7000 | 0.5000 | 0.5000 | 0.5000 | 2 | 2 | 2 | 18 |
| CP_structure17_XGB | 0.8375 | 0.4433 | 0.7500 | 0.5500 | 0.2500 | 0.2500 | 0.2500 | 1 | 3 | 3 | 17 |

### Table 2. Held-out predictions for additional designs

| design_id | success | density_RF_score | density_RF_prediction | structure_RF_score | structure_RF_prediction |
| --- | --- | --- | --- | --- | --- |
| 1.0000 | 1 | 0.8122 | 1 | 0.7172 | 1 |
| 2.0000 | 0 | 0.0000 | 0 | 0.1008 | 0 |
| 3.0000 | 0 | 0.0612 | 0 | 0.1265 | 0 |
| 4.0000 | 0 | 0.7563 | 1 | 0.6090 | 1 |
| 5.0000 | 0 | 0.1440 | 0 | 0.1382 | 0 |
| 7.0000 | 0 | 0.9122 | 1 | 0.6234 | 1 |
| 8.0000 | 0 | 0.4574 | 0 | 0.3559 | 0 |
| 12.0000 | 0 | 0.1919 | 0 | 0.1870 | 0 |
| 13.0000 | 0 | 0.0000 | 0 | 0.0015 | 0 |
| 14.0000 | 0 | 0.0000 | 0 | 0.0127 | 0 |
| 15.0000 | 0 | 0.0000 | 0 | 0.0015 | 0 |
| 16.0000 | 0 | 0.0000 | 0 | 0.0048 | 0 |

### Table 3. Earlier-cohort training and later-cohort testing

| model | AUC | accuracy | balanced_accuracy | precision | recall | TP | FP | FN | TN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3_LR | 1.0000 | 0.3333 | 0.6364 | 0.1111 | 1.0000 | 1 | 8 | 0 | 3 |
| CP_density3_RF | 0.9091 | 0.6667 | 0.8182 | 0.2000 | 1.0000 | 1 | 4 | 0 | 7 |
| CP_density3_XGB | 0.7727 | 0.7500 | 0.4091 | 0.0000 | 0.0000 | 0 | 2 | 1 | 9 |
| CP_density3_Stump | 0.8182 | 0.6667 | 0.8182 | 0.2000 | 1.0000 | 1 | 4 | 0 | 7 |
| CP_summary14_LR | 0.9091 | 0.4167 | 0.6818 | 0.1250 | 1.0000 | 1 | 7 | 0 | 4 |
| CP_summary14_RF | 1.0000 | 0.6667 | 0.8182 | 0.2000 | 1.0000 | 1 | 4 | 0 | 7 |
| CP_summary14_XGB | 0.7727 | 0.7500 | 0.4091 | 0.0000 | 0.0000 | 0 | 2 | 1 | 9 |
| CP_structure17_LR | 1.0000 | 0.4167 | 0.6818 | 0.1250 | 1.0000 | 1 | 7 | 0 | 4 |
| CP_structure17_RF | 1.0000 | 0.6667 | 0.8182 | 0.2000 | 1.0000 | 1 | 4 | 0 | 7 |
| CP_structure17_XGB | 0.7727 | 0.7500 | 0.4091 | 0.0000 | 0.0000 | 0 | 2 | 1 | 9 |

### Table 4. Construct labels

| construct | batch | design_id | success |
| --- | --- | --- | --- |
| puf12_r123_r5678_r56788 | previous | nan | 0 |
| puf12_r123_r5678_r5loop6788 | previous | nan | 0 |
| puf12_r123_r567_r567_r678 | previous | nan | 1 |
| puf12_r123_r567_r5loop67_loopr678 | previous | nan | 1 |
| puf12_r123_r567_r5loop67_r678 | previous | nan | 0 |
| puf12_r123_r567_r6loopr7_r5678 | previous | nan | 0 |
| puf12_r123_r567_r6loopr7_r5loop678 | previous | nan | 0 |
| puf12_r123_r567_r6r7_r5678 | previous | nan | 1 |
| puf_12_loopr4loopr4 | previous | nan | 0 |
| puf_12_loopr4r4 | previous | nan | 0 |
| puf_12_r4loopr4 | previous | nan | 0 |
| puf_12_r4r4 | previous | nan | 0 |
| 1_130_r123_r567_r567_r67_no_loop_8 | new | 1.0000 | 1 |
| 2_130_r123_r567_r5loop67_r678 | new | 2.0000 | 0 |
| 3_130_r123_r567_r56loop7_r678 | new | 3.0000 | 0 |
| 4_130_r123_r567_r567_rloop678 | new | 4.0000 | 0 |
| 5_130_r123_r567_r567_r6loop78 | new | 5.0000 | 0 |
| 7_130_r123_r567_r567_rloop67loop8 | new | 7.0000 | 0 |
| 8_130_r123_r567_r567_r6loop7loop8 | new | 8.0000 | 0 |
| 12_130_r123456654loop328 | new | 12.0000 | 0 |
| 13_130_r12345665432loop8 | new | 13.0000 | 0 |
| 14_130_r123456734loop568 | new | 14.0000 | 0 |
| 15_130_r12345673456loop8 | new | 15.0000 | 0 |
| 16_130_r12365673156loop8 | new | 16.0000 | 0 |
