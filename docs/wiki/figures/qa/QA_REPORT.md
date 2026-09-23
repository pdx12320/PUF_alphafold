# Figure verification

All five final figures pass rendered PDF collision checking with **zero FAIL and zero WARN findings**. All PDF text is at least **6.5 pt**. The two multi-panel figures pass the final plot-area alignment audit at the default **1.5 pt** tolerance without exemptions. The three single-panel figures record alignment as not applicable.

| Figure / panel | Evidence inspected | Summary and uncertainty | Observation unit | Visual and geometry result |
|---|---|---|---|---|
| 2 | Work/non-work confusion for the fixed four constructs | Counts; no interval | Construct | Labels, one false positive, count sum 4 and final bounds verified |
| 1 | Work/non-work confusion for LOCO24 | Counts; no interval | Construct held out from its fold | Labels, both error cells, count sum 24 and final bounds verified |
| 3a | C295 class agreement | Counts; no interval | Construct–endpoint prediction | Five observations; both classes shown; row alignment passes |
| 3b | C388 class agreement | Counts; no interval | Construct–endpoint prediction | Three classes shown, including the empty measured decrease row; all five observations retained |
| 3c | C871 class agreement | Counts; no interval | Construct–endpoint prediction | Five decreases and zero measured no-decreases are visible; one error retained |
| 4a | Absolute editing | Recorded construct summaries; replicate spread unavailable | Construct–endpoint summary | All 15 values visible; sequential 0–100% scale; exact units |
| 4b | Matched-control differences | Recorded paired differences; replicate spread unavailable | Construct–endpoint difference | All 15 differences visible; symmetric −60 to +60 pp scale; zero at midpoint |
| 5 | Two-model sequence preferences | Empirical generated-sequence frequencies; no inferential interval | Residue preference in 10,000 sequences per model | All 17 consensus substitutions; distinct colors and marker shapes; 0–1 frequency axis |

Full-figure and panel-by-panel visual inspection was performed on the final Python-rendered previews. The first Figure 3 draft had a left-side axis label beyond the page. The left margin was enlarged and the complete rendering and audit sequence was repeated. The final output has no clipping or collisions. Heatmap cell labels are intentional contained overlays and remain legible.

Two source-preflight warnings were reviewed:

- PNG is exported without TIFF: the intended destination is a Wiki, and SVG/PDF plus a 600-dpi PNG are the selected deliverables.
- The static width detector reads the literal `89` in `89 * MM` as an inch value. Direct inspection of final PDF page geometry confirms the intended dimensions: 89 × 78 mm for Figures 1–2, 183 × 84 mm for Figure 3, 183 × 94 mm for Figure 4, and 120 × 126 mm for Figure 5.

Data checks in the plotting script independently reconstruct confusion matrices from recorded predictions, compare the Model 1 counts with the metrics table, require all 15 Model 2 records, verify 13 agreements, verify all 17 Model 3 consensus rows are outside the protected TRM set, and verify the two models nominate the same amino acid at each selected position. No observations were sampled or removed from the requested panels. The Model 3 display includes the complete 17-row consensus subset of the 493-position source, as specified by its biological question.

The plotting script remains usable without extra QA tools. Passing `--qa-tools-dir` enables the extended final-layout audit used here. JSON reports and diagnostic overlays accompany this report. Diagnostic overlays are QA artifacts and do not replace the scientific figures.

## Interpretation boundaries

- Figures 1–2 show the recorded random-forest analyses and should not be attributed to the historical decision stump.
- Figure 3 shows five selected individually withheld LOCO records at each endpoint. It cannot establish prospective chronology or evaluate a simultaneous five-construct holdout.
- Fold-dependent label thresholds remain in the source CSV. “Within interval” does not imply a universal ±5-percentage-point criterion.
- Figure 4 presents construct-level summaries. Its sources do not support new replicate-level significance tests.
- Figure 5 shows computational mutation nominations. Its frequencies do not measure editing improvements or successful wet-lab validation.
