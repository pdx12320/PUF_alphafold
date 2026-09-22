---
title: "PUF scaffold function and TRM-dependent editing"
section: "Dry Lab / Model"
revision_date: "2026-09-22"
language: "en"
---

# PUF scaffold function and TRM-dependent editing

## Overview

We developed two complementary modeling modules to support PUF–APOBEC engineering. The first relates protein-internal contact organization to PUF scaffold function. The second evaluates how tripartite recognition motif (TRM) substitutions affect editing at the reporter target C388 and the reporter bystander sites C295 and C871. Together, these modules connect structural features with construct-level measurements and guide successive rounds of PUF design.[^project]

## Model 1. PUF scaffold functionality

### Structural representation and model construction

The scaffold model used the association between functional PUF arrangements and their protein-internal contact networks. Each construct was represented by three contact-probability (CP) densities, S4, S12 and S24. Each density counts unique protein residue pairs with CP greater than 0.5 and a sequence separation of at least 4, 12 or 24 residues, normalized by protein length. Features were calculated for each predicted structure and averaged within and across seeds, yielding one feature vector per construct.

The primary random forest used these three densities, 300 trees, a maximum depth of two, a minimum leaf size of two and balanced class weights. A companion forest added six descriptors of structural confidence and geometry: core pLDDT mean and minimum, core mean PAE, contact-weighted PAE, length-normalized radius of gyration and anisotropy. Scores of at least 0.5 were assigned to the work class.[^scaffold]

### Model performance and experimental characterization

The three-density random forest correctly classified all 14 constructs in the initial dataset during leave-one-construct-out validation, identifying all three work and eleven non-work constructs. Restricting the dataset to the twelve PUF12 constructs retained 12/12 correct predictions. In the expanded collection of 24 PUF12 constructs, the three-density forest correctly classified 21 constructs, while the nine-feature forest correctly classified 22. The nine-feature model achieved a balanced accuracy of 0.850 and a ROC-AUC of 0.913.[^scaffold]

We then examined four designs with the companion model trained on the remaining twenty constructs.

| Candidate | Work score | Predicted class | Experimental class |
|---|---:|---|---|
| Design 1: R123/R567/R567/R67-no-loop-8 | 0.854 | Work | Work |
| Design 3: R123/R567/R56-loop-7/R678 | 0.172 | Non-work | Non-work |
| Design 7: R123/R567/R567/R-loop-67-loop-8 | 0.680 | Work | Non-work |
| Design 8: R123/R567/R567/R6-loop-7-loop-8 | 0.431 | Non-work | Non-work |

The predictions agreed with three experimental labels. Design 1 received the highest work score and was experimentally functional. Designs 3 and 8 were predicted non-work and showed the corresponding experimental outcomes. These results established a structural basis for prioritizing PUF scaffold arrangements for functional characterization.[^scaffold]

## Model 2. TRM-dependent reporter editing

### Endpoint definitions and model construction

The TRM module used an experimental collection of 81 variants. Each editing measurement was compared with the matched control from the same experimental batch. Editing changes were expressed in percentage points, calculated as 100 times the difference between mutant and control editing fractions. Experimental replicates were aggregated by construct and condition, and model inputs were derived from structural contacts rather than measured editing outcomes.[^trmdata]

We defined three endpoint-specific classification tasks. The C295 classifier identified reductions greater than 10 percentage points, using balanced logistic regression with regularization parameter C = 0.01 and 110 protein–RNA contact-change descriptors. The C388 classifier identified retention of target activity, defined as a decrease of no more than 15 percentage points. It used a depth-two balanced random forest and the root-mean-square change in the whole protein–RNA CP matrix. The C871 classifier identified reductions greater than 20 percentage points, using a depth-one balanced random forest and nineteen contact, confidence and geometry descriptors.[^trmprotocol]

### Candidate prediction and experimental characterization

The models were fitted to 71 variants and applied to a ten-construct panel containing five designated candidates and five comparison constructs. The C295, C388 and C871 predictions agreed with 7/10, 5/10 and 5/10 experimental classes, respectively. Target-activity retention was correctly identified for P9-NTQ and P8-GVE, and C871 reduction was correctly identified for P9-NPS.[^trmresults]

The experimental characterization panel comprised P9-GNS, P9-NPS, P9-NTQ, P8-GVE and P4-R5-SNE+P7-R5-SNE. The P9 variants compared different TRMs at the same repeat position. P8-GVE provided a distinct positional context, while the double-SNE construct combined substitutions at P4 and P7.

| Construct | C295 editing | C388 editing | C871 editing |
|---|---:|---:|---:|
| P9-GNS | 30.55% | 62.68% | 26.37% |
| P9-NPS | 25.52% | 60.46% | 20.39% |
| P9-NTQ | 23.03% | 62.98% | 18.46% |
| P8-GVE | 25.22% | 65.86% | 10.23% |
| P4-R5-SNE+P7-R5-SNE | 19.80% | 80.34% | 32.60% |

P9-NTQ combined 62.98% C388 editing with 18.46% C871 editing. P8-GVE produced the lowest C871 editing in the panel, 10.23%, while retaining 65.86% C388 editing. The double-SNE construct reached 80.34% C388 editing and reduced C295 and C871 by 13.15 and 34.90 percentage points, respectively. These profiles provided complementary starting points for the next engineering round.[^candidates]

## Model 3. AI-guided PUF design outside TRM positions

Models 1 and 2 were evaluated against wet-lab measurements of scaffold function and reporter editing, respectively. Model 3 extends this work to protein residues outside the TRM recognition positions and provides mutation guidance for the next wet-lab design round.

We used an AiCE-inspired inverse-folding workflow on the AlphaFold 3 PUF12–APOE4 mRNA complex. ProteinMPNN and LigandMPNN each sampled 10,000 sequences while retaining the RNA context. A 493-position scan nominated 17 substitutions supported by both models. The wet-lab dataset identified P8, P9 and P7 as the main TRM positions for selectivity design, including the P8-GVE, P9-NTQ, P9-NPS, P9-GNS and P7-SYVIRR backbones. Model 3 therefore keeps these experimentally validated TRMs and introduces additional mutations only outside the recognition code.[^aice]

| Structural design class | Priority residues | Design purpose |
|---|---|---|
| Core stability and packing | T350V, D374E, V294I, M458L, L274I, T278I, L166I, V366I, R93K, E351L, R462L | Improve scaffold packing on favourable TRM backbones |
| C388-proximal fusion-end geometry | R45T, S30A, L41R, N12S | Tune the geometry near the deaminase fusion |
| Protein–RNA interface | A438G, H392D | Refine RNA-interface packing |
| Conserved scaffold hinges | G35, G71, G107, G179, G251, G287, V290, V362, P348, L101 | Preserve the native scaffold framework |

Eleven third-generation constructs, G0–G10, combine P8-GVE, P9-NTQ or P7-SYVIRR with the AI-guided core, geometry and interface substitutions. This set will be tested by measuring C295, C388 and C871 editing, providing direct experimental validation of non-TRM residue guidance.

## Conclusions

The scaffold module connected nonlocal protein contact density with recorded PUF function and correctly identified the functional Design 1 in the four-design assessment. The TRM module organized editing predictions by site and identified experimentally useful profiles at P9, P8 and the P4/P7 double-SNE design. The expanded AiCE–MPNN workflow now directs the next PUF engineering cycle toward non-recognition residues that can be combined with the validated TRM backbones. Together, the two modules support a design strategy that preserves experimentally favourable RNA recognition while optimizing the surrounding PUF scaffold.

## Data and methods

The source package contains complete model settings, experimental measurements, predictions and evaluation tables. The AiCE–MPNN release contains the PUF12–RNA complex, sampling configuration, full 493-position scan, repeat-level experimental summary and the G0–G10 design set.[^record]

[^project]: Team source, *dry_lab_figures_v1_clean(1).pptx*, “PUF Function and TRM Modeling.” C388 is the reporter target; C295 and C871 are reporter bystander sites.

[^scaffold]: [Scaffold metrics](source_data/scaffold_metrics.csv), [candidate and cross-validation predictions](source_data/scaffold_predictions.csv), and [random-forest configuration](source_data/scaffold_protocol.json).

[^trmdata]: [Ten-construct panel](source_data/TRM_ten_construct_manifest.csv) and [construct-level experimental records](source_data/TRM_candidate_measurements_and_ranks.csv).

[^trmprotocol]: [TRM model definitions and training settings](source_data/TRM_protocol.json).

[^trmresults]: [Complete endpoint evaluation](source_data/TRM_comparison_metrics.csv) and [candidate-level predictions](source_data/TRM_candidate_predictions.csv).

[^candidates]: [Measured editing fractions and matched-control changes](source_data/TRM_candidate_measurements_and_ranks.csv).

[^aice]: [AiCE × PUF12 inverse-folding workflow, experimental comparison and third-generation designs](https://github.com/pdx12320/PUF_alphafold/tree/main/aice_mpnn_20260922).

[^record]: [Repository results index](https://github.com/pdx12320/PUF_alphafold/tree/main/results_20260922). This revision adds no new experimental measurement.
