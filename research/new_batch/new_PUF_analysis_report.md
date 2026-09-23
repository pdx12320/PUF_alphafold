# Additional scaffold designs: frozen scores and structural checks

This historical scoring snapshot covers 15 designs, 60 seed tasks and 300 structure models. Each design has four seeds and five models per seed; design 6 is absent. All RNA inputs are AUGGAGGACGUGC (13 nucleotides). Experimental labels were unavailable at the time of this snapshot, so no accuracy or AUC is calculated here.

Among 12-repeat candidates, designs 7 and 8 were highlighted for follow-up. Their construct-level XGBoost scores were 0.756 and 0.619; respectively, three and four of the four seed scores exceeded 0.5. Designs 1 and 4 showed disagreement between models. Designs 9–11 contain 15 repeats and require extrapolation beyond the training scaffold class. Scores are uncalibrated and were generated from construct-level aggregated features.

Sequence and repeat-origin checks accompany the scores. Preserve the recorded sequence mappings when interpreting design names. Later experimental classifications are presented in the current Wiki and should not be retroactively treated as labels available to this historical scoring run.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/new_batch/new_PUF_analysis_report.md) preserves the historical prose and analysis context.

### Table 1. Frozen-model scores

| design_id | repeat_count | pp_nonlocal12_high_per_res | XGB_score | RF_score | LR_score | XGB_seed_positive_count | iptm | plddt_protein_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 12.0000 | 1.3327 | 0.3807 | 0.7797 | 0.7733 | 3.0000 | 0.9230 | 91.6007 |
| 2 | 12.0000 | 1.2979 | 0.2443 | 0.2504 | 0.6617 | 0.0000 | 0.9255 | 89.9931 |
| 3 | 12.0000 | 1.3081 | 0.2443 | 0.3918 | 0.6433 | 0.0000 | 0.9260 | 90.7871 |
| 4 | 12.0000 | 1.3232 | 0.3807 | 0.7717 | 0.7642 | 0.0000 | 0.9150 | 90.2774 |
| 5 | 12.0000 | 1.3164 | 0.2443 | 0.5128 | 0.6812 | 0.0000 | 0.9215 | 90.6663 |
| 7 | 12.0000 | 1.3237 | 0.7557 | 0.7374 | 0.7847 | 3.0000 | 0.9230 | 90.5177 |
| 8 | 12.0000 | 1.3165 | 0.6193 | 0.7660 | 0.7358 | 4.0000 | 0.9250 | 90.9730 |
| 9 | 15.0000 | 1.3698 | 0.7557 | 0.7616 | 0.9557 | 4.0000 | 0.9200 | 91.9485 |
| 10 | 15.0000 | 1.3536 | 0.7557 | 0.7524 | 0.9247 | 4.0000 | 0.9210 | 90.6612 |
| 11 | 15.0000 | 1.3437 | 0.7557 | 0.7637 | 0.9092 | 4.0000 | 0.9230 | 90.8630 |
| 12 | 12.0000 | 1.2987 | 0.2670 | 0.3972 | 0.6301 | 0.0000 | 0.9215 | 90.1263 |
| 13 | 12.0000 | 1.2398 | 0.2443 | 0.0435 | 0.1627 | 0.0000 | 0.9250 | 87.1102 |
| 14 | 12.0000 | 1.2822 | 0.2443 | 0.0461 | 0.4144 | 0.0000 | 0.9220 | 89.8512 |
| 15 | 12.0000 | 1.2714 | 0.2443 | 0.0461 | 0.2991 | 0.0000 | 0.9255 | 90.1504 |
| 16 | 12.0000 | 1.2183 | 0.2443 | 0.0196 | 0.0957 | 0.0000 | 0.9215 | 87.0432 |

### Table 2. Sequence and repeat-origin checks

| Physical repeat | Name-implied source | Sequence-matched source | Original core positions |
|---|---|---|---|
| P7 | R7 | R6 | 250–285 |
| P8 | R3 | R5 | 286–321 |
| P10 | R5 | R3 | 358–393 |
| P11 | R6 | R2 | 394–436 (including loop) |

### Table 3. Local-contact checks

| design_id | CP_204_239 | CP_196_235 | CCR01_density | CCR02_density |
| --- | --- | --- | --- | --- |
| 1 | 0.3850 | 0.5850 | 3.2175 | 3.0850 |
| 2 | 0.2825 | 0.4675 | 3.2383 | 3.0350 |
| 3 | 0.2825 | 0.5200 | 3.2408 | 3.0725 |
| 4 | 0.2900 | 0.5125 | 3.2058 | 3.0617 |
| 5 | 0.3225 | 0.5425 | 3.2300 | 3.0767 |
| 7 | 0.4025 | 0.5950 | 3.2442 | 3.0517 |
| 8 | 0.3900 | 0.5875 | 3.2267 | 3.0917 |
| 9 | 0.6175 | 0.6475 | 3.1733 | 3.0442 |
| 10 | 0.4875 | 0.6000 | 3.2100 | 3.0992 |
| 11 | 0.5375 | 0.6125 | 3.1967 | 3.0992 |

### Table 4. Output files

| design_id | construct | domain |
| --- | --- | --- |
| 1 | 1_130_r123_r567_r567_r67_no_loop_8 | related 12-repeat |
| 2 | 2_130_r123_r567_r5loop67_r678 | related 12-repeat |
| 3 | 3_130_r123_r567_r56loop7_r678 | related 12-repeat |
| 4 | 4_130_r123_r567_r567_rloop678 | related 12-repeat |
| 5 | 5_130_r123_r567_r567_r6loop78 | related 12-repeat |
| 7 | 7_130_r123_r567_r567_rloop67loop8 | related 12-repeat |
| 8 | 8_130_r123_r567_r567_r6loop7loop8 | related 12-repeat |
| 9 | 9_130_r123_r567_r567_r567_r67loop8 | 15-repeat extrapolation |
| 10 | 10_130_r123_r567_r567_r56loop7_r67_no8 | 15-repeat extrapolation |
| 11 | 11_130_r123_r567_r567_r56loop7_r67loop8 | 15-repeat extrapolation |
| 12 | 12_130_r123456654loop328 | new module order |
| 13 | 13_130_r12345665432loop8 | new module order |
| 14 | 14_130_r123456734loop568 | new module order |
| 15 | 15_130_r12345673456loop8 | new module order |
