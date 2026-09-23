# Quantitative figure contract

Backend: Python / Matplotlib, matching the repository workflow. All graphics are rendered from existing recorded data. No model is retrained by the figure script.

| Figure | Question and claim | Evidence role | Archetype | Export size |
|---|---|---|---|---|
| 2 | How well does the CP + structure random forest classify the fixed four constructs? Three of four labels agree; one of three non-work designs is a false positive. | Fixed-panel assessment | Quantitative single panel | 89 × 78 mm |
| 1 | How well does the same feature recipe classify the expanded scaffold collection under LOCO? Twenty-two of 24 labels agree. | Broader cross-validation | Quantitative single panel | 89 × 78 mm |
| 3 | Which endpoint predictions agree for the selected five constructs? Agreement is 5/5, 4/5 and 4/5 for C295, C388 and C871. | Endpoint-specific error decomposition | Quantitative grid | 183 × 84 mm |
| 4 | What edits were measured and how do they differ from the matched controls? The three sites have distinct absolute levels and changes. | Absolute measurements plus paired-control effects | Quantitative grid | 183 × 94 mm |
| 5 | Which non-TRM substitutions receive consistent sequence-model support? Seventeen substitutions share the most frequent residue across both models. | Computational candidate nomination | Quantitative single panel | 120 × 126 mm |

Figures 1 and 2 remain separate because the four-construct fixed assessment and 24-construct LOCO analysis use distinct evaluation designs. Neither is attributed to the historical decision stump.

Figure 3 uses 15 recorded predictions for five constructs at three endpoints. Each construct is individually withheld in LOCO; the selected five were not held out together. The label thresholds vary with the training fold and are retained in source data. The five C871 observations all have the decrease label, so this panel cannot estimate specificity. The aggregate 13/15 is a descriptive count across correlated endpoint predictions.

Figure 4 displays all 15 corresponding measurements and all 15 paired differences. The matched-control value is reconstructed as the recorded measurement minus the recorded difference. No observations are removed. These are construct-level summaries; replicate-level spread is unavailable in the selected source, so no invented error bars or inferential tests are added.

Figure 5 selects all 17 rows with `consensus == True`; all are non-TRM positions and nominate the same amino acid in both models. Frequencies are counts among 10,000 generated sequences per model. They describe model preference and do not estimate a mutation's experimental benefit. The full 493-position source remains unchanged.

Every figure has editable SVG/PDF text and a 600-dpi PNG. Minimum text is 6 pt. Multi-panel axes have measured final alignment within 1.5 pt. Final PDF text and collision audits are retained alongside the figures. Source-data CSVs and source hashes permit direct reproduction. Export dimensions serve this Wiki; no journal submission compliance is claimed.
