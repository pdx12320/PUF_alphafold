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

### Important PUF repeats identified from experimental editing data

We next used the complete set of 84 construct-level wet-lab records to identify repeat positions that consistently informed TRM design. P8, P9 and P7 emerged as the principal selectivity-design positions. P8-GVE gave the lowest C871/C388 ratio (0.153) in the collection, whereas P9-NTQ, P9-NPS and P9-GNS retained high C388 editing. P7-SYVIRR also retained substantial target activity. These positions therefore provide the main sequence-design space for maintaining C388 editing while reducing bystander editing.

P3 acted as a useful C295 control point. TRM changes at this repeat reduced median C295 editing to zero and included the selective VFQ profile, although target retention varied more widely. In contrast, P2 and P5 were sensitive positions: their variants showed lower median C388 retention and contributed few favourable C871/C388 profiles. P1, P4, P6, P10 and P11 were generally tolerated but produced fewer selectivity gains. P12 has not yet been represented in the experimental dataset and will be added to the next design cycle.[^aice]

| Design role | Repeat positions | Experimental design use |
|---|---|---|
| Primary selectivity positions | P8, P9, P7 | Prioritize TRM exploration and combined designs |
| C295 modulation position | P3 | Test selective C295 reduction with C388 retention monitoring |
| Sensitive positions | P2, P5 | Retain the native recognition code unless a targeted rationale is available |
| Tolerated positions | P1, P4, P6, P10, P11 | Use as secondary design positions |
| Unexplored position | P12 | Add systematic TRM coverage in the next round |

### AI-guided PUF design beyond TRM substitutions

To extend Model 2 beyond the recognition code, we established an AiCE-inspired inverse-folding workflow on the AlphaFold 3 PUF12–APOE4 mRNA complex. The input complex contained the 493-residue PUF12 protein and a 17-nt APOE4 RNA context centred on C388. ProteinMPNN and LigandMPNN each generated 10,000 sequences while the RNA context was retained during design. For every protein position, we compared the sampled frequency of the wild-type residue with that of the most frequent substitution. A site was nominated when the sampled substitution was favoured over the wild-type residue. The two models jointly nominated 17 consensus substitutions from a 493-position scan.[^aice]

The workflow separates two design layers. The experimental TRM analysis defines recognition-code substitutions at repeat positions P8, P9, P7 and P3. The inverse-folding scan supplies orthogonal non-recognition substitutions that can be combined with these experimentally favourable TRMs. This division preserves the successful C388 and bystander-editing profiles while expanding the design space toward protein stability, RNA-interface packing and fusion-end geometry.

| Structural design class | Priority residues | Intended use |
|---|---|---|
| Core stability and packing | T350V, D374E, V294I, M458L, L274I, T278I, L166I, V366I, R93K, E351L, R462L | Test additive stabilization on favourable TRM backbones |
| C388-proximal fusion-end geometry | R45T; S30A, L41R and N12S | Modify the N-terminal geometry facing the deaminase fusion |
| Protein–RNA interface | A438G, H392D | Tune RNA-interface packing on P8-GVE and P9-NTQ backbones |
| Conserved scaffold hinges | G35, G71, G107, G179, G251, G287, V290, V362, P348 and L101 | Preserve the native residues during combinatorial design |

The experimental comparison defined how these predictions should be used. Recognition-triplet residues at repeat positions 12, 13 and 16 were not advanced from inverse-folding frequency alone. The experimentally favourable TRMs, including GVE, NTQ, NPS, GNS and SYVIRR, were selected from the editing dataset and PUF recognition code. The AI scan instead nominated non-recognition residues for combination with these TRM backbones. This strategy directly links position-resolved wet-lab selection with structure-guided protein optimization.[^aice]

### Next PUF design set

We assembled eleven third-generation constructs that combine the P8-GVE or P9-NTQ experimental backbones with the AI-guided substitutions. G0 retains P8-GVE as the reference. G1 adds R45T, G2 adds the core-stability package, G3 combines R45T with the core package, G4 adds A438G and H392D, and G5 adds the N-terminal geometry package. G6 combines all design layers. G7 and G8 introduce the core and interface packages into P9-NTQ. G9 combines P8-GVE, P9-NTQ and the core package, and G10 adds the core package to the high-C388 SYVIRR backbone.[^aice]

The next wet-lab cycle will quantify C295, C388 and C871 editing for G0–G10. These measurements will test whether the experimentally selected recognition layer and the AI-guided non-recognition layer act additively, while expanding coverage at P8 and P12.

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
