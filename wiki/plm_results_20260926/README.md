# Protein language model prediction results — 26 September 2026

This directory records the sequence/structure representation analysis added after the original contact-probability Model 2. The analysis predicts continuous editing fractions at **C295, C388 and C871**, derives a specificity score, and separately evaluates C388 work/non-work classification.

## Main result

The pre-specified model-selection criterion was mean out-of-fold MAE across the three editing sites. **M3a_global_Ridge** was selected with a mean MAE of **9.979 percentage points**, compared with **10.619 pp** for the CP-only Ridge baseline and **11.244 pp** for the no-feature mean baseline.

| Model | C295 MAE | C388 MAE | C871 MAE | Mean MAE |
|---|---:|---:|---:|---:|
| M3a_global_Ridge | 7.201 | 10.070 | 12.666 | **9.979** |
| M0_Ridge (CP baseline) | 8.354 | 10.847 | 12.655 | 10.619 |
| B0 (mean baseline) | 8.117 | 11.563 | 14.052 | 11.244 |

For the selected model versus B0, the fixed-OOF paired MAE difference was −0.917 pp at C295, −1.494 pp at C388 and −1.386 pp at C871. The fixed-OOF C388 bootstrap interval was −2.751 to −0.120 pp.

## Ranking objective

Specificity is defined as:

`S = C388 - (C295 + C871) / 2`

On the main 66-construct complete-label cohort:

| Model | S MAE (pp) | Spearman | Kendall |
|---|---:|---:|---:|
| B1_Ridge | 6.691 | **0.490** | 0.344 |
| M1_local_Ridge | 7.163 | 0.459 | 0.308 |
| M3a_local_Ridge | 6.645 | 0.399 | 0.278 |
| M3a_global_Ridge | 6.957 | 0.348 | 0.226 |

For **M3a_global_Ridge**, NDCG was **0.867 at K=5** and **0.893 at K=10**. Four of the five highlighted favourable TRMs fall inside its predicted top 13: P9-NPS (#4), P9-NTQ (#9), P9-GNS (#10) and P8-GVE (#13). P4-R5-SNE + P7-R5-SNE is a major ranking miss for this model.

## C388 work classification

Work is defined as `C388 > 0.50`. The strongest AP in the logistic-classifier comparison was **M3a_local_Logistic**:

- AP: **0.778**
- AUROC: **0.656**
- Precision: **0.750**
- Recall: **0.692**
- Balanced accuracy: **0.624**

The positive fraction is 0.658.

## Structure subset

Only 20 variants had matched coordinates in this delivery. On that restricted subset, the CP/sequence baseline M0_Ridge had mean MAE **13.444 pp** and the tested SaProt/structure families did not improve on it. These 20 variants come from the first batch and should not be compared directly with the full-cohort results.

## Relationship to the existing Model 2

The original Model 2 classifies control-relative endpoint changes using protein–RNA contact features. This directory adds a quantitative task: continuous editing prediction and ranking by a shared specificity score. The existing fixed-five Model 2 test remains the classification evidence reported on the main Wiki page.

## Files

- `three_site_mae_models.csv` — pooled three-site MAE comparison.
- `ranking_models.csv` — pooled S ranking metrics.
- `top13_M3a_global_Ridge.csv` — predicted top 13 with measured S and measured rank.
- `topk_M3a_global_Ridge.csv` — NDCG and measured top-K summaries.
- `work_classification.csv` — C388 work/not-work classifier metrics.
- `batch_transfer_M3a_global_Ridge.csv` — cross-batch transfer metrics.
- `model_selection.csv` — selected final model and selection criterion.
- `methods_notes.md` — data scope, evaluation and interpretation notes.

The original delivery also contains large embedding matrices, fitted model objects, structural caches and run logs. Those intermediate binaries are intentionally not duplicated inside the Wiki directory.
