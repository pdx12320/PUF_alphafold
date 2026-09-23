# Wiki figure captions

Figures are provided as SVG, PDF and 600-dpi PNG. Source tables are in [figure_data](../figure_data/). Run [make_figures.py](../scripts/make_figures.py) for the original endpoint and sampling-support figures, and [add_wiki_figures.py](../scripts/add_wiki_figures.py) for the combined validation and added figures. These scripts use recorded data without retraining. No inferential statistics or uncertainty estimates are shown.

## Figure 1

Scaffold classification in two recorded evaluations. (a) Leave-one-construct-out (LOCO) evaluation in the expanded collection: 22/24 correct. (b) Four constructs withheld jointly: 3/4 correct, comprising one true positive, two true negatives and one false positive. Rows are measured classes and columns are predictions of the nine-feature random forest at a score threshold of 0.5. Both panels use the same colour scale for counts. The historical decision stump is a separate analysis.

![scaffold_validation_combined](scaffold_validation_combined.svg)

## Figure 2

Observed S12 distribution across the expanded scaffold collection. Each point represents one construct (20 non-work; four work); horizontal lines indicate medians. S12 is the per-residue high-confidence nonlocal contact feature, pp_nonlocal12_high_per_res. Horizontal offsets separate observations and have no quantitative meaning. The overlapping distributions are descriptive and do not establish a universal decision threshold.

![scaffold_s12_distribution](scaffold_s12_distribution.svg)

## Figure 3

Construct-level prediction–measurement comparison at C295, C388 and C871. Each construct was withheld individually in the saved LOCO evaluation; the five were not withheld jointly. Colours encode classes defined by thresholds selected within each training fold: D, decrease; ND, no threshold-crossing decrease; I, within the selected interval; U, increase. Outlines identify the two disagreements, P9-GNS at C388 and P8-GVE at C871. Agreement is 13/15 across correlated endpoint predictions. These five constructs were selected retrospectively from the existing experimental results.

![five_construct_class_pairs](five_construct_class_pairs.svg)

## Figure 4

Endpoint-level confusion matrices for the same five constructs. Agreement is 5/5 at C295, 4/5 at C388 and 4/5 at C871. All five C871 measurements belong to the decrease class, so C871 specificity cannot be estimated from this panel. Each construct was withheld individually; class thresholds are fold-specific.

![fig3_model2_five_construct_confusion](fig3_model2_five_construct_confusion.svg)

## Figure 5

Sampling frequencies for the exact 17 consensus substitutions. Each value is the fraction of 10,000 generated sequences carrying the indicated residue in one model. Frequency describes structural sequence preference and does not estimate the probability of improved RNA editing.

![fig5_model3_nontrm_consensus](fig5_model3_nontrm_consensus.svg)

## Figure 6

Locations of the 17 consensus non-TRM substitutions on the input PUF–RNA structure. Two orthographic projections show the protein Cα trace (grey), RNA phosphorus trace (blue) and nominated positions (orange). Numbers in panel b map to the mutation key. Coordinates come directly from aice_mpnn_20260922/inputs/complex.pdb; residue numbers follow its 493-residue protein chain A. The figure locates computational candidates and does not establish their effects on binding, stability or editing.

![nontrm_structure_locations](nontrm_structure_locations.svg)
