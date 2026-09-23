# Scaffold architecture validation

This analysis evaluates whether scaffold contact features transfer across repeat arrangements. The principal set contains 24 PUF12 constructs, including four work and 20 non-work labels. Structures are aggregated within seeds and then across seeds; constructs are the independent evaluation units. Architecture-out and source-order-out define identical partitions in this dataset and provide one grouping assessment.

Nonlocal contact density contributes useful ranking information, with substantial architecture dependence. The three-density random forest has architecture-out AUC 0.86875, compared with LOCO AUC 0.9125. The nine-feature contact-plus-structure random forest has architecture-out AUC 0.9375, but detects two of the four work constructs at the fixed 0.5 cutoff. Local features and contact regions do not consistently improve held-out performance. With only four positive constructs, ranking and thresholded classification must be read together.

The tables below retain the recorded numerical results, including sensitivity analyses. Model scores are uncalibrated. Use the current Wiki for the later fixed-four comparison; its evaluation scope differs from this historical architecture analysis.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/architecture_validation/PUF_architecture_validation_report.md) preserves the historical prose and analysis context.

### Table 1. Dataset and architecture groups

| architecture | source_order | n_construct | n_success | n_failure |
| --- | --- | --- | --- | --- |
| G01 | 1,2,3,4,5,6,6,5,4,3,2,8 | 2 | 0 | 2 |
| G02 | 1,2,3,4,5,6,7,3,4,5,6,8 | 2 | 0 | 2 |
| G03 | 1,2,3,5,6,7,5,6,7,6,7,8 | 10 | 3 | 7 |
| G04 | 1,2,3,5,6,7,6,7,5,6,7,8 | 3 | 1 | 2 |
| G05 | 1,2,3,5,6,7,8,5,6,7,8,8 | 2 | 0 | 2 |
| G06 | 1,2,3,5,6,7,8,6,7,8,4,4 | 4 | 0 | 4 |
| G08 | 1,2,3,6,5,6,6,5,1,3,2,8 | 1 | 0 | 1 |

### Table 2. Feature definitions and importance

| feature | definition | minimum | maximum | LOCO_RF_MDI | architecture_out_RF_MDI |
| --- | --- | --- | --- | --- | --- |
| pr_cp_sum | sum of all protein-RNA CP_ij | 42.56 | 54.39 | 0.04391 | 0.04711 |
| pr_cp_per_nt | protein-RNA CP sum / RNA nucleotide count | 3.274 | 4.184 | 0.03716 | 0.04014 |
| pr_top20_mean | mean of 20 largest protein-RNA CP values | 0.9279 | 0.9959 | 0.02139 | 0.02203 |
| pr_rna_max_mean | mean across RNA nucleotides of maximum CP over protein residues | 0.7021 | 0.8444 | 0.05289 | 0.05424 |
| pr_rna_max_min | minimum across RNA nucleotides of maximum CP over protein residues | 0.0825 | 0.21 | 0.008 | 0.01086 |
| pr_rna_coverage_05 | fraction of RNA nucleotides with maximum protein CP>0.5 | 0.6923 | 0.8077 | 0.0333 | 0.03727 |
| pr_protein_coverage_05 | fraction of protein residues with maximum RNA CP>0.5 | 0.05408 | 0.07749 | 0.03005 | 0.02374 |
| pr_high_contacts_per_nt | count protein-RNA CP>0.5 / RNA nucleotide count | 2.654 | 4.058 | 0.0192 | 0.0185 |
| pp_nonlocal4_per_res | sum CP_ij / L, protein i<j and j-i>=4 | 2.345 | 2.44 | 0.05102 | 0.05648 |
| pp_nonlocal4_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=4 | 2.213 | 2.387 | 0.1967 | 0.1943 |
| pp_nonlocal12_per_res | sum CP_ij / L, protein i<j and j-i>=12 | 1.286 | 1.36 | 0.06751 | 0.0749 |
| pp_nonlocal12_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=12 | 1.189 | 1.333 | 0.1704 | 0.1636 |
| pp_nonlocal24_per_res | sum CP_ij / L, protein i<j and j-i>=24 | 1.104 | 1.162 | 0.05158 | 0.05092 |
| pp_nonlocal24_high_per_res | count(CP_ij>0.5) / L, protein i<j and j-i>=24 | 1.027 | 1.128 | 0.2169 | 0.206 |

### Table 3. Exact permutation and architecture effects

| feature | success_mean | failure_mean | delta | hedges_g | pearson_r | spearman_r | p_exact | q_family | architecture_R2 | adjusted_success_beta | p_within_architecture_exact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pp_nonlocal12_high_per_res | 1.323 | 1.275 | 0.04799 | 1.286 | 0.4604 | 0.5815 | 0.0191 | 0.08915 | 0.7446 | 0.01928 | 0.06389 |
| P5_P6_CP | 0.03895 | 0.03758 | 0.001363 | 0.7273 | 0.2814 | 0.323 | 0.1619 | 0.9071 | 0.9872 | 0.000194 | 0.2528 |
| K204_slot_density | 3.189 | 2.282 | 0.907 | 0.8433 | 0.3219 | 0.5653 | 0.1023 | 0.1534 | 0.9908 | 0.1769 | 0.04167 |
| A235_slot_density | 3.344 | 2.71 | 0.6333 | 0.862 | 0.3282 | 0.4522 | 0.08743 | 0.1526 | 0.9902 | 0.08949 | 0.1222 |
| CCR01_fixed | 3.067 | 3.107 | -0.03946 | -0.08938 | -0.03601 | -0.09691 | 0.9145 | 0.9145 | 0.9975 | -0.03066 | 0.04167 |

### Table 4. Architecture-controlled comparisons

| feature | architecture | n_success | n_failure | delta |
| --- | --- | --- | --- | --- |
| pp_nonlocal12_high_per_res | G03 | 3 | 7 | 0.01178 |
| pp_nonlocal12_high_per_res | G04 | 1 | 2 | 0.04288 |
| P5_P6_CP | G03 | 3 | 7 | 0.0001096 |
| P5_P6_CP | G04 | 1 | 2 | 0.0004599 |
| K204_slot_density | G03 | 3 | 7 | 0.09762 |
| K204_slot_density | G04 | 1 | 2 | 0.4267 |
| A235_slot_density | G03 | 3 | 7 | 0.02888 |
| A235_slot_density | G04 | 1 | 2 | 0.2804 |
| CCR01_fixed | G03 | 3 | 7 | -0.03153 |
| CCR01_fixed | G04 | 1 | 2 | -0.02792 |

### Table 5. Repeat interfaces and contact regions

| feature | reference_residues | repeat_slot | min_internal_Pearson | min_internal_Spearman | pearson_r | p_exact | q_family | p_within_architecture_exact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCR01_fixed | 178,179,180 | 5 | -0.8922 | -0.1001 | -0.03601 | 0.9145 | 0.9145 | 0.04167 |
| CCR_new_01 | 162,163,164 | 4 | 0.7405 | 0.7548 | -0.05276 | 0.659 | 0.659 | 0.002778 |
| CCR_new_02 | 195,196,197,198 | 5 | 0.7622 | 0.6139 | -0.1915 | 0.5305 | 0.659 | 0.2222 |
| CCR_new_03 | 411,412,413 | 11 | 0.6435 | 0.7722 | 0.3762 | 0.07086 | 0.1856 | 0.9889 |
| CCR_new_04 | 414,415,416 | 11 | 0.7828 | 0.6144 | 0.3657 | 0.07425 | 0.1856 | 0.7 |
| CCR_new_05 | 447,448,449 | 12 | 0.6513 | 0.709 | 0.2521 | 0.2372 | 0.3954 | 0.075 |

### Table 6. Random-forest validation

| feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3 | RF | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| CP_interface | RF | LOCO | 0.925 | 0.622 | 0.5193 | 0.8 | 0.5 | 0.75 | 17 | 3 | 1 | 3 |
| CP_structural | RF | LOCO | 0.9125 | 0.7875 | 0.7688 | 0.85 | 0.75 | 0.75 | 19 | 1 | 1 | 3 |
| CP_summary14 | RF | LOCO | 0.85 | 0.4444 | 0.3608 | 0.7 | 0.5 | 0.5 | 18 | 2 | 2 | 2 |
| S12_only | RF | LOCO | 0.875 | 0.4833 | 0.3994 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | RF | LOCO | 0.85 | 0.6104 | 0.5701 | 0.725 | 0.6667 | 0.5 | 19 | 1 | 2 | 2 |
| structural_only | RF | LOCO | 0.75 | 0.7917 | 0.7871 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| CP_density3 | RF | architecture_out | 0.8688 | 0.475 | 0.6542 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_interface | RF | architecture_out | 0.9062 | 0.6625 | 0.7479 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_structural | RF | architecture_out | 0.9375 | 0.8611 | 0.8524 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_summary14 | RF | architecture_out | 0.7625 | 0.4196 | 0.3182 | 0.65 | 0.3333 | 0.5 | 16 | 4 | 2 | 2 |
| S12_only | RF | architecture_out | 0.8 | 0.3409 | 0.583 | 0.725 | 0.3333 | 0.75 | 14 | 6 | 1 | 3 |
| full | RF | architecture_out | 0.85 | 0.7019 | 0.6822 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |
| structural_only | RF | architecture_out | 0.7625 | 0.7935 | 0.7888 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |

### Table 7. All classifier configurations

| feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CP_density3 | LR | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_density3 | RF | LOCO | 0.9125 | 0.5845 | 0.4798 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| CP_density3 | XGB | LOCO | 0.9 | 0.6012 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| CP_interface | LR | LOCO | 0.925 | 0.6917 | 0.65 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_interface | RF | LOCO | 0.925 | 0.622 | 0.5193 | 0.8 | 0.5 | 0.75 | 17 | 3 | 1 | 3 |
| CP_interface | XGB | LOCO | 0.9 | 0.6012 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| CP_structural | LR | LOCO | 0.8375 | 0.8088 | 0.8028 | 0.625 | 0.2308 | 0.75 | 10 | 10 | 1 | 3 |
| CP_structural | RF | LOCO | 0.9125 | 0.7875 | 0.7688 | 0.85 | 0.75 | 0.75 | 19 | 1 | 1 | 3 |
| CP_structural | XGB | LOCO | 0.65 | 0.3738 | 0.2761 | 0.55 | 0.25 | 0.25 | 17 | 3 | 3 | 1 |
| CP_summary14 | LR | LOCO | 0.825 | 0.4132 | 0.3309 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_summary14 | RF | LOCO | 0.85 | 0.4444 | 0.3608 | 0.7 | 0.5 | 0.5 | 18 | 2 | 2 | 2 |
| CP_summary14 | XGB | LOCO | 0.8875 | 0.5179 | 0.4131 | 0.675 | 0.4 | 0.5 | 17 | 3 | 2 | 2 |
| S12_only | LR | LOCO | 0.9 | 0.65 | 0.6077 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| S12_only | RF | LOCO | 0.875 | 0.4833 | 0.3994 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| S12_only | XGB | LOCO | 0.675 | 0.4417 | 0.2996 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | LR | LOCO | 0.725 | 0.3125 | 0.2394 | 0.675 | 0.2727 | 0.75 | 12 | 8 | 1 | 3 |
| full | RF | LOCO | 0.85 | 0.6104 | 0.5701 | 0.725 | 0.6667 | 0.5 | 19 | 1 | 2 | 2 |
| full | XGB | LOCO | 0.3375 | 0.1655 | 0.1261 | 0.375 | 0 | 0 | 15 | 5 | 4 | 0 |
| structural_only | LR | LOCO | 0.7875 | 0.6568 | 0.6201 | 0.6 | 0.2143 | 0.75 | 9 | 11 | 1 | 3 |
| structural_only | RF | LOCO | 0.75 | 0.7917 | 0.7871 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| structural_only | XGB | LOCO | 0.75 | 0.7917 | 0.7871 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_density3 | LR | architecture_out | 0.9 | 0.7042 | 0.6646 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_density3 | RF | architecture_out | 0.8688 | 0.475 | 0.6542 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_density3 | XGB | architecture_out | 0.8938 | 0.5179 | 0.5298 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_interface | LR | architecture_out | 0.8875 | 0.7409 | 0.7205 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| CP_interface | RF | architecture_out | 0.9062 | 0.6625 | 0.7479 | 0.75 | 0.375 | 0.75 | 15 | 5 | 1 | 3 |
| CP_interface | XGB | architecture_out | 0.9062 | 0.5679 | 0.5548 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_structural | LR | architecture_out | 0.925 | 0.85 | 0.8417 | 0.75 | 0.2857 | 1 | 10 | 10 | 0 | 4 |
| CP_structural | RF | architecture_out | 0.9375 | 0.8611 | 0.8524 | 0.75 | 1 | 0.5 | 20 | 0 | 2 | 2 |
| CP_structural | XGB | architecture_out | 0.8812 | 0.5 | 0.5119 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| CP_summary14 | LR | architecture_out | 0.8 | 0.4623 | 0.3606 | 0.7 | 0.25 | 1 | 8 | 12 | 0 | 4 |
| CP_summary14 | RF | architecture_out | 0.7625 | 0.4196 | 0.3182 | 0.65 | 0.3333 | 0.5 | 16 | 4 | 2 | 2 |
| CP_summary14 | XGB | architecture_out | 0.8875 | 0.5179 | 0.4131 | 0.575 | 0.3333 | 0.25 | 18 | 2 | 3 | 1 |
| S12_only | LR | architecture_out | 0.8875 | 0.7409 | 0.7205 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| S12_only | RF | architecture_out | 0.8 | 0.3409 | 0.583 | 0.725 | 0.3333 | 0.75 | 14 | 6 | 1 | 3 |
| S12_only | XGB | architecture_out | 0.8125 | 0.5088 | 0.6562 | 0.825 | 0.6 | 0.75 | 18 | 2 | 1 | 3 |
| full | LR | architecture_out | 0.875 | 0.5417 | 0.4375 | 0.775 | 0.3077 | 1 | 11 | 9 | 0 | 4 |
| full | RF | architecture_out | 0.85 | 0.7019 | 0.6822 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |
| full | XGB | architecture_out | 0.5438 | 0.2917 | 0.3165 | 0.55 | 0.25 | 0.25 | 17 | 3 | 3 | 1 |
| structural_only | LR | architecture_out | 0.8875 | 0.8269 | 0.8197 | 0.725 | 0.2667 | 1 | 9 | 11 | 0 | 4 |
| structural_only | RF | architecture_out | 0.7625 | 0.7935 | 0.7888 | 0.875 | 1 | 0.75 | 20 | 0 | 1 | 3 |
| structural_only | XGB | architecture_out | 0.3812 | 0.375 | 0.3393 | 0.625 | 1 | 0.25 | 20 | 0 | 3 | 1 |

### Table 8. PUF11 and TRM sensitivity

| subset | feature_set | model | validation | AUC | AP | PR_AUC | balanced_accuracy | precision | recall | TN | FP | FN | TP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| plus_PUF11 | CP_density3 | RF | LOCO | 0.9205 | 0.5845 | 0.4798 | 0.8295 | 0.6 | 0.75 | 20 | 2 | 1 | 3 |
| plus_PUF11 | S12_only | RF | LOCO | 0.8864 | 0.4833 | 0.3994 | 0.8295 | 0.6 | 0.75 | 20 | 2 | 1 | 3 |
| plus_PUF11 | full | RF | LOCO | 0.8977 | 0.6528 | 0.6108 | 0.7273 | 0.6667 | 0.5 | 21 | 1 | 2 | 2 |
| plus_PUF11 | CP_density3 | RF | architecture_out | 0.8807 | 0.475 | 0.6542 | 0.7614 | 0.375 | 0.75 | 17 | 5 | 1 | 3 |
| plus_PUF11 | S12_only | RF | architecture_out | 0.8182 | 0.3409 | 0.583 | 0.7386 | 0.3333 | 0.75 | 16 | 6 | 1 | 3 |
| plus_PUF11 | full | RF | architecture_out | 0.9318 | 0.85 | 0.8417 | 0.5 | 0 | 0 | 22 | 0 | 4 | 0 |
| plus_TRM | CP_density3 | RF | LOCO | 0.925 | 0.6829 | 0.6088 | 0.8667 | 0.7143 | 0.8333 | 18 | 2 | 1 | 5 |
| plus_TRM | S12_only | RF | LOCO | 0.8833 | 0.569 | 0.5058 | 0.8667 | 0.7143 | 0.8333 | 18 | 2 | 1 | 5 |
| plus_TRM | full | RF | LOCO | 0.925 | 0.8357 | 0.8238 | 0.8083 | 0.8 | 0.6667 | 19 | 1 | 2 | 4 |
| plus_TRM | CP_density3 | RF | architecture_out | 0.8875 | 0.6042 | 0.7566 | 0.7917 | 0.5 | 0.8333 | 15 | 5 | 1 | 5 |
| plus_TRM | S12_only | RF | architecture_out | 0.8167 | 0.4557 | 0.6792 | 0.7667 | 0.4545 | 0.8333 | 14 | 6 | 1 | 5 |
| plus_TRM | full | RF | architecture_out | 0.975 | 0.9444 | 0.941 | 0.75 | 1 | 0.5 | 20 | 0 | 3 | 3 |
