# Module 2: Design–Build–Test–Learn

[Project home](../../README.md) · [Full model narrative](../../wiki/PUF_Model_EN.md) · [Research index](../README.md)

Our engineering objective is to retain C388 target editing while reducing C295 and C871 bystander editing. We organize Module 2 into three linked Design–Build–Test–Learn (DBTL) cycles, moving from scaffold selection to recognition-motif engineering and non-TRM design. Each cycle connects a design decision to measured outcomes and a specific next step. This iterative use of experimental feedback follows the DBTL framework described by Radivojević et al. (2020).

### Cycle 1. Establishing a functional scaffold

**Design.** We first asked which repeat arrangements could support a functional PUF scaffold. Structural contact features provided an interpretable starting point for screening arrangements before experimental testing.

**Build.** We established a decision-tree baseline and subsequently combined contact and structural descriptors in a shallow random forest. The final classifier used twenty training constructs. Designs 1, 3, 7 and 8 were excluded together from preprocessing, fitting and model selection, forming a separate test panel.

**Test.** Predictions were compared with the subsequent wet-lab outcomes for the four designs. Design 1 worked, whereas Designs 3, 7 and 8 did not work. The classifier agreed with three outcomes and incorrectly predicted Design 7 as functional.

**Learn.** Contact-based screening supported prioritization of functional scaffolds, while the Design 7 error identified a limit of structural descriptors. Scaffold function also left a separate engineering question unresolved: how would recognition-motif changes affect target and bystander editing? This motivated the endpoint-specific analysis in Cycle 2.

### Cycle 2. Balancing target and bystander editing

**Design.** We examined recognition-motif substitutions to preserve target editing while reducing unwanted editing at nearby reporter sites. The programmability of Pumilio recognition residues provided the biological rationale for this design space (Cheong & Hall, 2006).

**Build.** We trained separate classifiers for C295, C388 and C871 using contact features and, where selected, sequence descriptors. P9-GNS, P9-NPS, P9-NTQ, P8-GVE and P4-SNE + P7-SNE formed the fixed test panel. All five candidates and identical sequences were excluded together before preprocessing and model selection. Training used the remaining 76 constructs, with endpoint thresholds selected within the training data.

**Test.** Subsequent wet-lab measurements agreed with 11 of 15 endpoint classifications: five at C295, three at C388 and three at C871. Experimental editing profiles also identified useful backgrounds. Relative to its matched control, P4-SNE + P7-SNE increased C388 editing by 7.77 percentage points and reduced C295 and C871 editing by 13.15 and 34.90 points, respectively.

**Learn.** Retaining all three endpoints exposed trade-offs that a single work/non-work label could not capture. Contact attribution highlighted P3/P5/P2 for C295, P3/P4/P7 for C388 and P1/P10/P11 for C871. These associations suggest regions for further investigation, while the measured editing profiles identify backgrounds for the next design round.

### Cycle 3. Extending the search beyond recognition motifs

**Design.** We next explored substitutions outside the recognition motifs while preserving the intended recognition code. ProteinMPNN and LigandMPNN provided complementary structure-conditioned sequence proposals (Dauparas et al., 2022, 2025).

**Build.** Each model generated 10,000 sequences from the PUF structural input. After excluding recognition-code positions, we retained 17 non-TRM substitutions with the same alternative residue supported at frequencies of at least 80% by both models. These substitutions constitute computational designs; their wet-lab construction and testing remain the next cycle.

**Test.** We compared sampled preferences with exact TRM substitutions already represented in the experimental dataset. Several favourable experimental variants had low or zero sampling counts, demonstrating limited correspondence between structural sequence preference and editing performance. None of the 17 non-TRM nominations had an exact experimentally measured substitution match in that dataset.

**Learn.** Sampling frequency can nominate structurally compatible candidates, but provides insufficient evidence for improved editing. The next test should compare individual non-TRM substitutions with their unchanged, experimentally characterized TRM backgrounds before evaluating combinations. Measuring C388, C295 and C871 together would determine whether each substitution improves the desired editing profile and whether combinations retain individual benefits.

### Feeding experimental evidence into the next design round

The three cycles progressively refine the question from scaffold function to editing selectivity and then sequence context. New measurements should be linked to exact substitutions, matched controls and the frozen predictions that preceded testing. After evaluation, those measurements can inform an updated training set, with new candidates reserved for the next independent test panel. This closes the learning loop while preserving a clear record of what each model predicted before experimental feedback.

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

## References

Cheong, C.-G., & Hall, T. M. T. (2006). Engineering RNA sequence specificity of Pumilio repeats. *Proceedings of the National Academy of Sciences, 103*(37), 13635–13639. https://doi.org/10.1073/pnas.0606294103

Dauparas, J., Anishchenko, I., Bennett, N., Bai, H., Ragotte, R. J., Milles, L. F., Wicky, B. I. M., Courbet, A., de Haas, R. J., Bethel, N., Leung, P. J. Y., Huddy, T. F., Pellock, S., Tischer, D., Chan, F., Koepnick, B., Nguyen, H., Kang, A., Sankaran, B., … Baker, D. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science, 378*(6615), 49–56. https://doi.org/10.1126/science.add2187

Dauparas, J., Lee, G. R., Pecoraro, R., An, L., Anishchenko, I., Glasscock, C., & Baker, D. (2025). Atomic context-conditioned protein sequence design using LigandMPNN. *Nature Methods, 22*(4), 717–723. https://doi.org/10.1038/s41592-025-02626-1

Radivojević, T., Costello, Z., Workman, K., & Garcia Martin, H. (2020). A machine learning Automated Recommendation Tool for synthetic biology. *Nature Communications, 11*, Article 4879. https://doi.org/10.1038/s41467-020-18008-4
