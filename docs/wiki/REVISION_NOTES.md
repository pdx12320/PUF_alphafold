# Wiki revision record — 23 September 2026

## Scope and evidence

This revision reorganizes the supplied English draft and repository documentation into three engineering modules. Each explains its biological question, computational method, evaluation and wet-lab handoff. Original experiments, model outputs, sampled sequences and mutation FASTA files are unchanged. The new figures are derived views of recorded source tables.

| Result | Narrative role | Location |
|---|---|---|
| Initial decision tree, 14/14 LOCO | Discovery-stage evidence | Model 1 fitting and pinned historical report |
| Expanded random forest, 22/24 LOCO | Broader evaluation | Model 1, Figure 1 |
| Simultaneous four-design RF panel, 3/4 | Prediction–experiment case comparison | Model 1, Figure 2 and all four rows |
| Combined-PR 81-construct evaluation and matched baseline | Context that bounds generalization | Model 2 evaluation table |
| Selected five constructs, 13/15 endpoint agreement | Construct-level illustration | Model 2, Figure 3 |
| Measured editing and control-relative changes | Experimental design consequence | Model 2, Figure 4 and complete five-row table |
| 17 exact MPNN consensus substitutions | Computational nominations | Model 3, Figure 5 |
| G0–G10 | Explicit next-round design | Model 3 and unchanged source FASTA/CSV |

## Terminology ledger

| Term | Definition and editorial decision |
|---|---|
| Work / non-work | Original scaffold outcome labels; the archived label does not identify the failed experimental stage |
| Decision tree | Historical depth-one baseline; its results remain distinct from the random forest |
| Random forest | Ensemble of shallow decision trees used for the recorded four-design comparison |
| CP | Contact probability; density summaries and individual probabilities are distinct quantities |
| TRM | Tripartite recognition motif within a PUF repeat |
| C388 | Reporter target RNA cytidine; distinct from protein residue 388 |
| C295 / C871 | Reporter bystander RNA cytidines; the panel does not establish transcriptome-wide specificity |
| Percentage points (pp) | 100 × (variant fraction − matched-control fraction) |
| Held out | Excluded from the corresponding computational fit; does not assert new prospective wet-lab chronology |
| N12 | Residue 12 in the full 493-residue design reference |
| Sampling frequency | Fraction of generated sequences containing the exact alternative amino acid; not an editing-improvement probability |

## Replaced, retained and relocated material

- Replaced the two-module opening with a three-model engineering objective and input/output map.
- Retained the historical decision-tree discovery and named the random forest used in the four-design assessment.
- Replaced the prior ten-panel fixed classifier definitions with the actual dynamic classifier/threshold workflow behind the five displayed LOCO predictions.
- Retained both errors in the five-construct panel and contextualized that selected subset using complete-cohort results and matched baselines.
- Added exact non-TRM mutation identities, pretrained sampling details and G0–G10 package definitions.
- Corrected RNA conditioning, the stage of TRM protection, AiCE's full name, and the N12S/N12G export discrepancy.
- Removed unsupported causal claims about core stability, conserved hinges, catalytic geometry and improved editing from generated sequence frequencies.
- Moved detailed source interpretation and reanalysis commands to the Model 3 methods page; retained original historical DBTL text under its date.
- Added verified APA references and acknowledged the official iGEM pages used as structural examples.

## Evidence chain and repetition check

The overview introduces the three decisions. Each model section supplies its own decisive evidence. Captions define labels, units and the evaluation population. The integration section states the next design consequences without repeating the full performance tables. The complete tables and figure source CSVs retain numerical provenance.

The revision reports descriptive counts and source metrics. It adds no new significance tests, confidence intervals, wet-lab measurements or model training claims. Full-cohort baseline evidence remains visible in the main Wiki because it changes how the favourable five-construct panel should be interpreted.

## Reference verification

Bibliographic fields were checked against publisher, PubMed and official author/repository records. Crossref was attempted but unavailable through the search service. Metadata verification and claim support were assessed separately. See the Wiki bibliography for seven APA references and the page-structure sources for three official iGEM precedents.

## Validation

Figure generation reconstructs counts from prediction rows and cross-checks saved metrics. It verifies all 15 selected construct–endpoint rows and all 17 same-amino-acid consensus nominations. Model 3 consensus frequencies were also checked directly against the two archived 10,000-sequence ensembles. Figure layout checks and visual inspection are recorded in [the figure QA report](figures/qa/QA_REPORT.md).

The existing result tables and FASTA outputs remain byte-identical to the base commit. Relative documentation links were checked after English and Chinese synchronization. This edit contains no new model fits or wet-lab experiments.

## Length record

Base repository commit: `c38943d46e4e08754a01b3160f53a599505878e4`. Counts use whitespace-delimited tokens, including tables and captions.

| Model section | Before | After |
|---|---:|---:|
| Model 1 | 345 | 607 |
| Model 2 | 552 | 775 |
| Model 3 | 251 | 877 |
