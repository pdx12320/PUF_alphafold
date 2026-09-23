# Original local-contact comparison: English guide

The original analysis reproduced baseline comparisons and evaluated the available RNA-local feature families. Full protein–protein matrix inputs were unavailable, so the planned mutation-centred global and nonlocal communication analyses were incomplete.

C388 maps to nucleotide 15 of ACAUGGAGGACGUGC. Target-centred windows are clipped at this RNA boundary, with no artificial downstream padding. The 512-residue protein core corresponds to positions 8–519 of the full sequence. Features are computed after equal-weight seed aggregation. Interface_RMS retains the historical 519×15 denominator.

All preprocessing and ContactSeek feature discovery occur inside the appropriate training folds. Mutation-combination holdouts may share component positions. Fixed 0.5 cutoffs, training-selected cutoffs and adaptive model selection use distinct protocols. Feature importances describe the fitted models; they do not establish causal contact mechanisms.

For the full source package, run `python run_baselines.py --stage reproduce`, `python run_baselines.py --stage fixed --include-local`, and `python report_available.py` inside that package. The curated repository snapshot may omit those source entry points; recover the verified source archive before attempting a full rerun.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/c388_local_nonlocal_optimization/SOURCE_README.md) preserves the historical prose and analysis context.

### Table 1. Recorded baseline and feature results

| validation | balanced_accuracy | MCC | AUC |
|---|---|---|---|
| LOOCV | 0.5729 | 0.1550 | 0.7708 |
| mutation_group_out | 0.8583 | 0.7258 | 0.6833 |

### Table 2. Mutation-group errors

| model | score_policy | balanced_accuracy_group_out | balanced_accuracy_loocv |
|---|---|---|---|
| LR_C0.1 | fixed_0.5 | 0.6875 | 0.7188 |
| LR_C0.1 | training_only_threshold | 0.5000 | 0.7188 |
| nested_model_selection | training_only_threshold | 0.5000 | 0.7188 |

### Recorded table 3

| feature_set | model | balanced_accuracy_group_out | MCC_group_out | balanced_accuracy_loocv |
|---|---|---|---|---|
| Interface_CP | RF_depth2 | 0.9375 | 0.8101 | 0.8542 |
| CP_summary | LR_C1 | 0.8750 | 0.6708 | 0.7917 |
| ContactSeek_Top3 | LR_C1 | 0.7500 | 0.4629 | 0.7500 |
| Total_CP | LR_C1 | 0.7500 | 0.4629 | 0.7500 |
| ContactSeek_CCR | LR_C1 | 0.7500 | 0.4629 | 0.6979 |
| WT_Top5_pm0 | LR_C0.1 | 0.7188 | 0.4183 | 0.7500 |
| WT_Top5_pm1 | LR_C1 | 0.7188 | 0.4183 | 0.6042 |
| protein_RNA_CP | LR_C1 | 0.6875 | 0.3750 | 0.8125 |
| C388_pm0 | LR_C0.1 | 0.6875 | 0.3750 | 0.6875 |
| C388_pm1 | LR_C0.1 | 0.6875 | 0.3750 | 0.6875 |

### Recorded table 4

| feature_set | balanced_accuracy_group_out | MCC_group_out | balanced_accuracy_loocv |
|---|---|---|---|
| Interface_CP | 0.5000 | 0.0000 | 0.7708 |
| Total_CP | 0.5000 | 0.0000 | 0.7188 |
| nested_best_available | 0.5000 | 0.0000 | 0.6562 |
| nested_best_local | 0.4062 | -0.2433 | 0.6979 |
| nested_best_WT_contact | 0.4062 | -0.2433 | 0.3021 |

### Recorded table 5

| mutation_group | n | accuracy | balanced_accuracy | FP | FN |
|---|---|---|---|---|---|
| P4 | 2 | 1.0000 | nan | 0 | 0 |
| P4+P7 | 4 | 1.0000 | 1.0000 | 0 | 0 |
| P5 | 9 | 0.7778 | 0.7500 | 0 | 2 |
| P6 | 5 | 1.0000 | nan | 0 | 0 |
| P7 | 1 | 1.0000 | nan | 0 | 0 |
| P8 | 1 | 1.0000 | nan | 0 | 0 |
