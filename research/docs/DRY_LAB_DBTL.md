# Module 2 | PUF modelling and design: Design–Build–Test–Learn

## Learning from model failures and revising the design workflow

Our PUF modelling workflow developed through repeated attempts to resolve specific engineering problems. We compared structural representations, investigated inconsistent predictions and revised the modelling objectives when scaffold scores failed to describe editing activity. Experimental feedback also changed how we used inverse folding to propose new designs.

The cycles below describe these decisions and the evidence that motivated them. Each Learn section identifies a change in the workflow. The module connects scaffold screening, endpoint-specific TRM models and non-TRM design. Iterative use of experimental feedback follows the Design–Build–Test–Learn (DBTL) framework described by Radivojević et al. (2020). Historical diagnostics and the current fixed-panel evaluations are identified within their respective cycles.

| Cycle | Problem encountered | Intervention tested | Decision after testing |
|---|---|---|---|
| 1 | The useful structural representation was unclear | Compare contact summaries, local structural features and model complexity | Retain a compact baseline and strengthen validation across repeat arrangements |
| 2 | A known functional scaffold received low scores on new inputs | Audit sequence boundaries and recompute features over the matched core | Standardise preprocessing and group records by protein identity |
| 3 | Scaffold scores did not predict target-editing performance | Test direct transfer and retrain against editing measurements | Separate scaffold function from editing activity and use matched-control changes |
| 4 | Editing selectivity required simultaneous consideration of three sites | Train endpoint models, evaluate a fixed five-candidate test panel and identify influential repeat regions | Combine three-site measurements with model attribution to choose backgrounds and sites for further engineering |
| 5 | Inverse folding failed to recover useful recognition motifs | Compare sequence preferences with existing wet-lab measurements | Preserve experimentally supported TRMs and redirect inverse folding towards non-TRM candidates |

## Cycle 1 | Choosing a structural representation that remains useful beyond the initial dataset

### Design

We initially lacked a clear basis for choosing structural features for PUF scaffold assessment. Protein-internal contacts and protein–RNA contacts were both plausible inputs. We therefore compared these representations before committing to a more complex model.

### Build

We extracted contact-probability (CP) summaries and RNA-contact profiles from AlphaFold 3 predictions for an initial scaffold collection (Abramson et al., 2024). Structures were aggregated by construct to prevent repeated predictions from becoming independent training examples. The protein-internal representation included S4, S12 and S24, which summarised high-probability contacts at different sequence separations.

We subsequently added structural confidence, geometry, residue-level contacts and correlated contact regions. As the collection expanded, we tested whether feature expansion and random-forest tuning improved the compact baseline. A subsequent contact-plus-structure classifier was fitted to twenty training constructs. Designs 1, 3, 7 and 8 were excluded together from preprocessing, fitting and model selection for the fixed test evaluation.

### Test

The initial protein-internal summaries distinguished functional and non-functional constructs more effectively than the RNA-contact profiles. However, adding more features did not consistently improve performance. Local contact candidates lacked support after multiple-testing correction, and an initially identified correlated region did not reproduce in the expanded collection.

Tuning the three-density random forest also failed to improve its classification counts in the expanded dataset. More elaborate feature and parameter selection gave variable results. These tests exposed the limited value of continuing to increase model complexity on the available data.

In the current fixed four-design evaluation, predictions agreed with three subsequent wet-lab outcomes. Design 1 worked, whereas Designs 3, 7 and 8 did not work. Design 7 was the false-positive prediction. The initial decision tree established the baseline; these fixed-panel predictions came from the later shallow random forest.

### Learn

We retained the compact structural baseline as a reference and evaluated additional descriptors against it. Validation was extended to hold out entire repeat arrangements, addressing whether performance depended on shared scaffold backgrounds.

**Change implemented:** retain a simple reference model, assess added features against it and evaluate transfer across repeat arrangements. The next problem was ensuring that new structures were represented consistently with the training inputs.

## Cycle 2 | Diagnosing unexpectedly low scores for a known functional scaffold

### Design

When we applied the scaffold workflow to new TRM complexes, even a known functional reference received low scores. This prompted an input audit. We asked whether differences in the analysed protein region could explain the disagreement with the existing experimental label.

### Build

The new inputs comprised two protein sequences, each paired with three RNA contexts. Their 519-residue sequences included an additional seven-residue N-terminal segment relative to the 512-residue analysis region. We recalculated the features over the matched core while keeping the trained classifier unchanged.

### Test

All six full-length complexes initially scored below the classification threshold of 0.5. After matching the analysis region, five of six scores exceeded this threshold. For the reference paired with C388 RNA, the score changed from 0.000 to 0.786.

No structures were refolded and the classifier was not retrained. The score changes therefore diagnosed sensitivity to preprocessing. The six complexes represented only two proteins, so they also exposed the need to distinguish structural records from independent constructs.

### Learn

We learned that length-normalised contact features could still be sensitive to inconsistent sequence boundaries. The preprocessing rule needed to be defined before scoring new inputs. RNA contexts and repeated structures also needed to remain grouped under their protein identity.

**Change implemented:** align protein analysis regions and preserve construct-level grouping. Recovering a scaffold score still left a separate question unresolved: could that score predict editing activity?

## Cycle 3 | Replacing direct scaffold-model transfer with an editing-specific task

### Design

A functional scaffold score did not directly describe how much editing a TRM variant retained. We tested whether the existing scaffold model could nevertheless classify C388 activity, then compared direct transfer with models trained for that endpoint.

### Build

The historical test included 31 TRM proteins with experimental C388 measurements. Mean editing of at least 50% defined the initial activity class. We compared the frozen scaffold model with retrained classifiers using CP and structural features, under construct-held-out and mutation-position-combination evaluations.

### Test

The frozen scaffold model incorrectly classified eight of nine below-threshold variants as meeting the activity criterion. Its balanced accuracy was 0.396. Retraining against the activity labels improved some comparisons, confirming that the original scaffold classifier could not simply be reused for this task.

The fixed absolute cutoff also represented a different question from retaining activity relative to a matched control. This distinction became important when combining experimental batches and interpreting reductions at the bystander sites.

### Learn

We separated scaffold functionality from reporter editing and defined editing changes relative to the appropriate control. The workflow then evaluated C388 retention and C295/C871 reduction as distinct endpoints. This revised the model's objective as well as its training labels.

**Change implemented:** train against the editing outcome of interest and express changes in percentage points relative to matched controls. The next evaluation asked whether the three endpoint outputs could jointly support candidate selection.

## Cycle 4 | Turning endpoint models into experimentally informed design choices

### Design

Scaffold functionality left the target–bystander trade-off unresolved. A recognition-motif variant could reduce C295 and C871 editing while also losing C388 activity. We therefore treated the three reporter sites as separate prediction targets and considered their measured changes together. Engineering Pumilio recognition residues provided the biological basis for exploring this sequence space (Cheong & Hall, 2006).

### Build

We trained endpoint-specific classifiers using protein–RNA contact features and, where selected, recognition-motif sequence descriptors. P9-GNS, P9-NPS, P9-NTQ, P8-GVE and P4-SNE + P7-SNE formed a fixed test panel. All five constructs and identical sequences were excluded together before preprocessing, feature selection, threshold selection and fitting. The remaining 76 constructs supplied the training data.

The original project workflow trained models and predicted candidates before obtaining their wet-lab measurements. We subsequently revised C388 as a binary computational reanalysis, retaining the original five-construct exclusion. This rerun does not establish new pre-assay predictions. C295 and C871 retain their original fixed-panel models here.

For C388, we combined protein–RNA CP from both batches and compared interface, high-variance and full-CP representations. Training-only selection used seven candidate boundaries from −2.5 to −20 pp. The final Ridge classifier selected −5 pp: lower changes were labelled decrease, while changes at or above the boundary were labelled unchanged, including increases. The model used 876 eligible CP entries and no sequence descriptors. Its nominal top-1,000 variance configuration retained every eligible entry, so it did not establish a feature-selection advantage.

For the retained C295 and C871 classifiers, contact importance was aggregated by residue and repeat using impurity importance or absolute standardized coefficients. C388 binary attribution has not yet been recomputed.

### Test

With binary C388 and the retained C295/C871 models, predictions agreed with **11 of 15** endpoint classes: **5/5 at C295, 3/5 at C388 and 3/5 at C871**. The five constructs were excluded together for each model. The C295 test panel contained only the no-decrease class, so its agreement does not establish sensitivity to decreases.

Binary C388 correctly classified P9-NPS and P8-GVE as decreases and P4-SNE + P7-SNE as unchanged. P9-GNS and P9-NTQ were incorrectly classified as decreases. Its 3/5 agreement exceeded the training-majority count of 2/5; balanced accuracy was 66.7% versus 50.0%. This small difference did not establish a robust advantage.

Nested leave-one-identity-out validation across the 76 development constructs gave **48/76 correct (63.2%)** and **60.4% balanced accuracy**. The matched majority baseline achieved **50/76 (65.8%)** and **61.8% balanced accuracy**. Thresholds, features and classifiers were reselected within each outer training fold, while the five fixed tests remained excluded. Selected thresholds ranged from −2.5 to −20 pp. Thus, the complete automatic selection procedure did not outperform the baseline, despite the higher count on the fixed panel.

The measured profiles nevertheless identified useful engineering backgrounds. P4-SNE + P7-SNE increased C388 editing by **7.77 percentage points** relative to its matched control, while reducing C295 and C871 editing by **13.15 and 34.90 points**, respectively. These continuous changes provide information beyond the thresholded prediction classes.

Contact attribution highlighted different repeat regions at each endpoint. The table uses the common 493-residue reference; the verified mapping from the earlier Model 2 reference subtracts 16 residues.

| Endpoint | Leading repeat regions | Representative high-importance residues |
|---|---|---|
| C295 | P3, P5, P2 and P4 | R109, Y181, Y73 and Y145 |
| C388 (binary) | Not re-estimated | Previous three-class residue priorities are not transferred |
| C871 | P1, P10 and P11, with additional signal outside repeat cores | R361, Q364 and Y397; E7 outside repeat cores |

Earlier retrospective ranking tests exposed a separate problem: the P4/P7-SNE combination ranked seventh computationally but first experimentally in a ten-construct panel. That historical comparison used a different split and evaluated overall ordering. It remains evidence that endpoint classification and candidate ranking require separate assessment.

### Learn

We used the three-site experimental profiles to identify useful TRM backgrounds and retained C295/C871 attribution to guide regional hypotheses. C388 binary attribution has not been re-estimated. Its nested validation did not establish a baseline advantage, so experimental profiles remain essential for candidate selection. Importance describes model associations rather than evidence that a substitution improves editing.

**Decision carried forward:** preserve favourable experimental backgrounds, inspect endpoint-specific errors and test residue-level hypotheses on matched backgrounds. The resulting handoff includes a named construct, its measured editing profile and candidate regions for additional design.

The [C388 binary results and validation records](../c388_binary/) preserve both evaluations and their matched baselines. The former three-class results remain available as historical development evidence.

## Cycle 5 | Redirecting inverse folding after disagreement with experimental TRM results

### Design

We next explored whether structure-conditioned sequence preferences could guide further PUF design. Before using those preferences to select recognition motifs, we checked whether the workflow recovered TRMs with favourable experimental editing profiles.

### Build

ProteinMPNN and LigandMPNN each sampled 10,000 sequences using the predicted PUF12 structure, with LigandMPNN incorporating the surrounding atomic context (Dauparas et al., 2022, 2025). A position-wise scan covered the 493-residue protein. The updated comparison used 81 experimentally measured variant identities, reconstructed through a verified residue mapping. Recognition positions were tracked separately from other residues. This replaced the older 84-record comparison with the same matched-control measurements used for the endpoint analysis.

### Test

Several experimentally useful recognition motifs were rarely sampled. P8-GVE, corresponding to S288G and Y289V, had zero joint occurrences in both sets of 10,000 sequences. Experimentally, its C388 change was −6.70 percentage points, accompanied by C295 and C871 reductions of 7.74 and 57.27 points.

The comparison below matches exact substitutions. Each sampling count requires all listed changes to occur within the same sequence, with other positions unconstrained. WT counts refer to the unchanged residues at those same positions.

| Experimental variant | Exact substitution set | ProteinMPNN mutant / WT count | LigandMPNN mutant / WT count |
|---|---|---:|---:|
| P1-SHE | Y37H | 22 / 5 | 2 / 0 |
| P4-SNE | Y145N | 569 / 839 | 0 / 1889 |
| P4-SNE + P7-SNE | Y145N; Y253N | 120 / 359 | 0 / 461 |
| P9-NTQ | Y325T | 19 / 1017 | 1 / 7464 |
| P8-GVE | S288G; Y289V | 0 / 305 | 0 / 3006 |

These are selected favourable measured examples, not an estimate of overall prediction accuracy. The complete 81-variant comparison is retained in the linked project evidence. P1-SHE provides limited agreement between sampling preference and the measured profile, although its mutant frequencies were only 0.22% and 0.02%.

The mismatch is consistent with the models' objectives: they learn structure-conditioned sequence preferences, without training here on editing outcomes. The PUF structural input also omits the complete PUF–APOBEC fusion and cellular context. These differences plausibly contribute to the disagreement, but have not been experimentally established as its causes. A zero count denotes absence from the sampled sequences.

### Learn

We revised how inverse-folding outputs were used for design. Experimentally supported TRMs were preserved in the proposed constructs, and sequence nominations were directed towards non-TRM residues. Recognition-code positions were filtered after sampling in the present workflow.

Both models supported the same alternative residue at **17 non-TRM positions**, each at a sampling frequency of at least 80%. P10 contained four nominations: T350V, E351L, V366I and D374E. Additional examples included R45T in P1, S100E in P3 and L274I in P7. These repeat assignments identify concrete sites for experiments, without assigning functional benefit from sampling frequency alone.

**Change implemented:** preserve favourable TRM backgrounds and use inverse folding to propose additional non-TRM substitutions. None of the 17 nominated substitutions has an exact measured match in the current experimental cohort. Their experimental testing forms the next cycle.

## Next cycle | Testing non-TRM additions on matched backgrounds

The next cycle should determine which individual substitutions improve editing and whether their effects persist in combinations. The 17 consensus candidates provide a defined starting set. Designs should retain the parental TRM code so that each comparison isolates the added non-TRM change.

| Stage | Proposed action | Question resolved |
|---|---|---|
| Design | Select individual consensus substitutions on experimentally characterized TRM backgrounds | Which structural nominations preserve or improve the parental editing profile? |
| Build | Prepare the individual variants and unchanged parental comparators before constructing combinations | Can effects be assigned to a defined substitution and background? |
| Test | Measure C388, C295 and C871 in matched comparisons, reporting individual measurements and uncertainty | Does target retention accompany reduced bystander editing? |
| Learn | Retain supported substitutions, remove disruptive changes and compare selected combinations with their individual components | Are the benefits reproducible and compatible across substitutions? |

New results should be linked to exact protein sequences, matched controls and the predictions recorded before testing. Once a test round is complete, its measurements can inform subsequent model updates. Newly selected candidates should form a fresh test panel for that next iteration.

## Project evidence

The early diagnostic results in Cycles 1–3 retain the evidence described in the supplied engineering record. The 519/512-residue preprocessing audit is distinct from the 493-residue inverse-folding reference. Current fixed-panel metrics and exact-mutation comparisons use the updated analyses below. The 23 September rerun reconstructs the training/test separation; its artifacts do not independently establish the timing of the original experiments.

- [Historical scaffold diagnostics and activity-model transfer](https://github.com/pdx12320/PUF_alphafold/blob/c38943d46e4e08754a01b3160f53a599505878e4/docs/DRY_LAB_DBTL.md).
- [Fixed-panel training and evaluation](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/fixed_test_20260923/README.md).
- [Four-design scaffold predictions](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/fixed_test_20260923/model1/predictions.csv).
- [Five-construct endpoint predictions](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/fixed_test_20260923/model2/predictions.csv).
- [Repeat and residue attribution methods](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/model_interpretation_20260923/README.md).
- [Complete exact-mutation comparison](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/model_interpretation_20260923/model3_exact_mutation_wetlab_comparison.csv).
- [Historical ten-construct ranking](https://github.com/pdx12320/PUF_alphafold/blob/b5a7aaa7ef8426d1c98243d1a09c0db9a5db1c1a/research/results_20260922/ten_construct_holdout/ranking.csv).

## References

Abramson, J., Adler, J., Dunger, J., Evans, R., Green, T., Pritzel, A., Ronneberger, O., Willmore, L., Ballard, A. J., Bambrick, J., Bodenstein, S. W., Evans, D. A., Hung, C.-C., O’Neill, M., Reiman, D., Tunyasuvunakool, K., Wu, Z., Žemgulytė, A., Arvaniti, E., … Jumper, J. M. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature, 630*(8016), 493–500. https://doi.org/10.1038/s41586-024-07487-w

Cheong, C.-G., & Hall, T. M. T. (2006). Engineering RNA sequence specificity of Pumilio repeats. *Proceedings of the National Academy of Sciences, 103*(37), 13635–13639. https://doi.org/10.1073/pnas.0606294103

Dauparas, J., Anishchenko, I., Bennett, N., Bai, H., Ragotte, R. J., Milles, L. F., Wicky, B. I. M., Courbet, A., de Haas, R. J., Bethel, N., Leung, P. J. Y., Huddy, T. F., Pellock, S., Tischer, D., Chan, F., Koepnick, B., Nguyen, H., Kang, A., Sankaran, B., … Baker, D. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science, 378*(6615), 49–56. https://doi.org/10.1126/science.add2187

Dauparas, J., Lee, G. R., Pecoraro, R., An, L., Anishchenko, I., Glasscock, C., & Baker, D. (2025). Atomic context-conditioned protein sequence design using LigandMPNN. *Nature Methods, 22*(4), 717–723. https://doi.org/10.1038/s41592-025-02626-1

Radivojević, T., Costello, Z., Workman, K., & Garcia Martin, H. (2020). A machine learning Automated Recommendation Tool for synthetic biology. *Nature Communications, 11*, Article 4879. https://doi.org/10.1038/s41467-020-18008-4
