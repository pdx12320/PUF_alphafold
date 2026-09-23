# C871 and C295 classification in the frozen 22-variant cohort

All 22 variants follow the frozen C388 cohort membership; no variant is reselected using the current endpoint outcomes. Experimental rates are means of two records. Changes are mutant minus site-matched WT in percentage points; WT rates are 67.5015% at C871 and 32.958% at C295.

A positive label indicates a decrease beyond the specified negative threshold. The grid runs from −5 to −60 percentage points in steps of five. Eight feature families and nine classifier configurations are evaluated, including contact summaries, local contact changes, sequence changes and training-fold ContactSeek regions. Preprocessing and feature discovery are fitted within training folds.

Evaluation includes construct leave-one-out, mutation-combination holdout and strict-position holdout. Double mutants may appear in multiple strict-position tests and are inversely weighted by prediction count. Classifier thresholds maximize inner-fold balanced accuracy, with ties favouring cutoffs near 0.5. Adaptive analyses select biological labels, features, models and score thresholds entirely within each outer training set.

The tables retain fixed-candidate results, model comparisons and adaptive evaluation. Exploratory best candidates were selected from many alternatives; changes in biological labels change the task. Adaptive pooled metrics can combine different endpoint boundaries across folds. No independent prospective cohort or full-pipeline permutation test is supplied.

Run `python run.py` and `python report.py` within this directory to reproduce the numerical analysis; the original report generator retains its historical language. The input caches, source licenses and provenance records accompany the results.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/trm22_offtarget/REPORT.md) preserves the historical prose and analysis context.

### Table 1. Biological label partitions

| site | threshold_pp | absolute_editing_pct | positive | negative | eligible_headline |
| --- | --- | --- | --- | --- | --- |
| C871 | -5 | 62.5015 | 19 | 3 | False |
| C871 | -10 | 57.5015 | 15 | 7 | True |
| C871 | -15 | 52.5015 | 13 | 9 | True |
| C871 | -20 | 47.5015 | 10 | 12 | True |
| C871 | -25 | 42.5015 | 9 | 13 | True |
| C871 | -30 | 37.5015 | 7 | 15 | True |
| C871 | -35 | 32.5015 | 4 | 18 | False |
| C295 | -5 | 27.9580 | 18 | 4 | False |
| C295 | -10 | 22.9580 | 14 | 8 | True |
| C295 | -15 | 17.9580 | 8 | 14 | True |
| C295 | -20 | 12.9580 | 3 | 19 | False |

### Table 2. Exploratory best fixed candidates

| site | validation | label_threshold_pp | family | model | BA | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | LOOCV | -15 | Repository | RF_depth2 | 0.8846 | 0.8376 | 1.0000 | 0.7692 | 10.0000 | 9.0000 | 0.0000 | 3.0000 |
| C871 | Position_group_out | -20 | Repository | RF_depth1 | 0.8750 | 0.7917 | 0.7692 | 1.0000 | 10.0000 | 9.0000 | 3.0000 | 0.0000 |
| C871 | Strict_position_out | -30 | CP_summary | RF_depth1 | 0.8571 | 0.8714 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C295 | LOOCV | -15 | Sequence_Local | RF_depth1 | 0.7768 | 0.7500 | 0.8333 | 0.6250 | 5.0000 | 13.0000 | 1.0000 | 3.0000 |
| C295 | Position_group_out | -10 | Local_deltaCP | LR_0.01 | 0.7321 | 0.5536 | 0.8333 | 0.7143 | 10.0000 | 6.0000 | 2.0000 | 4.0000 |
| C295 | Strict_position_out | -10 | ContactSeek_CCR | LR_1 | 0.8571 | 0.7299 | 0.8710 | 0.9643 | 13.5000 | 6.0000 | 2.0000 | 0.5000 |

### Table 3. Group-selected candidate across validation schemes

| site | validation | policy | label_threshold_pp | family | model | BA | AUC | precision | recall |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | LOOCV | fixed_0.5 | -20 | Repository | RF_depth1 | 0.7833 | 0.8167 | 0.6923 | 0.9000 |
| C871 | LOOCV | training_only_threshold | -20 | Repository | RF_depth1 | 0.6917 | 0.8167 | 0.6154 | 0.8000 |
| C871 | Position_group_out | fixed_0.5 | -20 | Repository | RF_depth1 | 0.7750 | 0.7917 | 0.7273 | 0.8000 |
| C871 | Position_group_out | training_only_threshold | -20 | Repository | RF_depth1 | 0.8750 | 0.7917 | 0.7692 | 1.0000 |
| C871 | Strict_position_out | fixed_0.5 | -20 | Repository | RF_depth1 | 0.8083 | 0.7708 | 0.7037 | 0.9500 |
| C871 | Strict_position_out | training_only_threshold | -20 | Repository | RF_depth1 | 0.8083 | 0.7708 | 0.7037 | 0.9500 |
| C295 | LOOCV | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5982 | 0.5357 | 0.7273 | 0.5714 |
| C295 | LOOCV | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.5625 | 0.5357 | 0.7000 | 0.5000 |
| C295 | Position_group_out | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5536 | 0.5536 | 0.7143 | 0.3571 |
| C295 | Position_group_out | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.7321 | 0.5536 | 0.8333 | 0.7143 |
| C295 | Strict_position_out | fixed_0.5 | -10 | Local_deltaCP | LR_0.01 | 0.5536 | 0.5982 | 0.7143 | 0.3571 |
| C295 | Strict_position_out | training_only_threshold | -10 | Local_deltaCP | LR_0.01 | 0.4375 | 0.5982 | 0.5385 | 0.2500 |

### Table 4. Training-selected classifier cutoffs

| site | validation | score_threshold_median | score_threshold_min | score_threshold_max |
| --- | --- | --- | --- | --- |
| C871 | LOOCV | 0.3906 | 0.1881 | 0.6186 |
| C871 | Position_group_out | 0.5000 | 0.4260 | 0.5641 |
| C871 | Strict_position_out | 0.5000 | 0.4838 | 0.5573 |
| C295 | LOOCV | 0.4889 | 0.3893 | 0.5912 |
| C295 | Position_group_out | 0.5639 | 0.3943 | 0.6500 |
| C295 | Strict_position_out | 0.4666 | 0.2699 | 0.5608 |

### Table 5. Fully adaptive training-only selection

| site | validation | BA | accuracy | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C295 | LOOCV | 0.5729 | 0.6818 | 0.8021 | 0.4000 | 0.3333 | 2.0000 | 13.0000 | 3.0000 | 4.0000 |
| C295 | Position_group_out | 0.5104 | 0.5909 | 0.4271 | 0.7333 | 0.6875 | 11.0000 | 2.0000 | 4.0000 | 5.0000 |
| C295 | Strict_position_out | 0.4643 | 0.4773 | 0.5631 | 0.6522 | 0.5000 | 7.5000 | 3.0000 | 4.0000 | 7.5000 |
| C871 | LOOCV | 0.5128 | 0.5455 | 0.4444 | 0.4286 | 0.3333 | 3.0000 | 9.0000 | 4.0000 | 6.0000 |
| C871 | Position_group_out | 0.3482 | 0.4091 | 0.4375 | 0.5333 | 0.5714 | 8.0000 | 1.0000 | 7.0000 | 6.0000 |
| C871 | Strict_position_out | 0.6417 | 0.6364 | 0.6385 | 0.7000 | 0.5833 | 7.0000 | 7.0000 | 3.0000 | 5.0000 |

### Table 6. Exploratory classifier comparison

| site | model | label_threshold_pp | family | BA | AUC | precision | recall |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | RF_depth1 | -20 | Repository | 0.8750 | 0.7917 | 0.7692 | 1.0000 |
| C871 | RF_depth2 | -20 | Repository | 0.8750 | 0.7833 | 0.7692 | 1.0000 |
| C871 | ExtraTrees_depth2 | -20 | ContactSeek_CCR | 0.8500 | 0.8000 | 1.0000 | 0.7000 |
| C871 | LDA_shrinkage | -30 | CP_summary | 0.7524 | 0.5714 | 0.8000 | 0.5714 |
| C871 | LR_1 | -30 | Local_deltaCP | 0.7333 | 0.6190 | 0.4667 | 1.0000 |
| C295 | LR_0.01 | -10 | Local_deltaCP | 0.7321 | 0.5536 | 0.8333 | 0.7143 |
| C295 | LDA_shrinkage | -15 | Sequence | 0.7232 | 0.7143 | 0.5385 | 0.8750 |
| C295 | RF_depth2 | -10 | Interface_RMS | 0.7232 | 0.6875 | 0.8889 | 0.5714 |
| C295 | RF_depth1 | -10 | Interface_RMS | 0.7232 | 0.6027 | 0.8889 | 0.5714 |
| C871 | LR_0.01 | -25 | Total_CP | 0.7137 | 0.6752 | 0.5714 | 0.8889 |
| C871 | LR_0.1 | -25 | Total_CP | 0.6966 | 0.6667 | 0.5833 | 0.7778 |
| C295 | ExtraTrees_depth2 | -10 | Total_CP | 0.6964 | 0.3393 | 0.8182 | 0.6429 |
| C295 | SVM_RBF | -15 | Total_CP | 0.6875 | 0.5000 | 0.5000 | 0.8750 |
| C871 | GaussianNB | -10 | Sequence | 0.6857 | 0.6429 | 0.8000 | 0.8000 |
| C295 | GaussianNB | -10 | Sequence | 0.6786 | 0.6786 | 0.7500 | 0.8571 |
| C295 | LR_1 | -10 | ContactSeek_CCR | 0.6518 | 0.6607 | 0.8571 | 0.4286 |
| C295 | LR_0.1 | -15 | Sequence | 0.6518 | 0.6071 | 0.4667 | 0.8750 |
| C871 | SVM_RBF | -25 | Total_CP | 0.6368 | 0.6923 | 0.5000 | 0.8889 |

### Table 7. Strict-position-selected candidates

| site | label_threshold_pp | family | model | validation | BA | AUC | precision | recall | TP | TN | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C871 | -30 | CP_summary | RF_depth1 | LOOCV | 0.8286 | 0.9333 | 0.6667 | 0.8571 | 6.0000 | 12.0000 | 3.0000 | 1.0000 |
| C871 | -30 | CP_summary | RF_depth1 | Position_group_out | 0.8571 | 0.8190 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C871 | -30 | CP_summary | RF_depth1 | Strict_position_out | 0.8571 | 0.8714 | 1.0000 | 0.7143 | 5.0000 | 15.0000 | 0.0000 | 2.0000 |
| C295 | -10 | ContactSeek_CCR | LR_1 | LOOCV | 0.6429 | 0.6518 | 0.7333 | 0.7857 | 11.0000 | 4.0000 | 4.0000 | 3.0000 |
| C295 | -10 | ContactSeek_CCR | LR_1 | Position_group_out | 0.6518 | 0.6607 | 0.8571 | 0.4286 | 6.0000 | 7.0000 | 1.0000 | 8.0000 |
