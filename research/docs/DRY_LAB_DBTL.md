# Design–Build–Test–Learn

[Project home](../../README.md) · [Full model narrative](../../wiki/PUF_Model_EN.md) · [Research index](../README.md)

## Engineering objective

Retain C388 target editing while reducing C295 and C871 bystander editing. Scaffold arrangement, recognition motifs and non-TRM sequence changes are evaluated as distinct design decisions.

| Stage | Design and build | Test | Learn and next experiment |
|---|---|---|---|
| Scaffold screening | Use contact features to establish a decision-tree baseline, then combine contact and structural descriptors in a shallow random forest. | The fixed four-construct panel is excluded jointly from fitting. Recorded labels comprise one work and three non-work constructs; three predictions agree. | Contact features can prioritize scaffolds; test additional arrangements using a frozen recipe. |
| Recognition-motif engineering | Predict control-relative editing classes separately at C295, C388 and C871. | Five retrospectively selected constructs contribute 15 individually held-out endpoint calls; 13 agree with measurements. | Compare complementary editing profiles and carry useful backgrounds into subsequent designs. |
| Non-TRM sequence exploration | Sample 10,000 sequences each from ProteinMPNN and LigandMPNN and filter recognition-code positions after sampling. | Seventeen substitutions receive same-residue support from both models. | Compare individual substitutions and combinations with their matched TRM backgrounds in the wet lab. |

## Development history

The project progressed from contact-based scaffold classification to expanded scaffold evaluation, architecture-held-out testing, endpoint-specific TRM models and inverse-folding nominations. Each archived result belongs to its recorded cohort and split definition.

| Historical stage | Finding | Consequence |
|---|---|---|
| Initial scaffold baseline | Small contact-feature models separated the initial work/non-work records. | Expand the scaffold collection and evaluate repeat-arrangement transfer. |
| Local contacts and contact regions | Local hypotheses did not show consistent corrected evidence or transferable gains. | Retain global nonlocal contact density as a compact baseline. |
| Expanded scaffold evaluation | Three-density RF LOCO AUC was 0.9125; architecture-out AUC was 0.86875. | Report both construct and arrangement holdouts. |
| Structure-feature comparison | Contact-plus-structure RF architecture-out AUC was 0.9375, with two of four positives detected at 0.5. | Read ranking and thresholded classification together. |
| Endpoint-specific TRM modelling | C388, C295 and C871 require their own labels, controls and validation scopes. | Use endpoint-specific predictions and explicit threshold definitions. |
| Inverse-folding proposals | Structural sequence preferences nominate non-TRM candidates. | Measure editing effects experimentally before assigning functional benefit. |

The [complete earlier engineering record](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/docs/DRY_LAB_DBTL.md) preserves the original detailed history. [Analysis history](ANALYSIS_HISTORY.md) and the [research index](../README.md) locate the retained numerical snapshots.
