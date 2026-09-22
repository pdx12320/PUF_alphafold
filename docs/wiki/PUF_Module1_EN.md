# Module 1 — Predicting PUF scaffold functionality

**Engineering question:** can contact organization help us distinguish PUF repeat arrangements with different experimentally recorded work/non-work outcomes?

**Main result:** a shallow random forest correctly classified all **14 initial constructs in leave-one-construct-out evaluation**. Restricting that initial set to the **12 PUF12 constructs also gave 12/12 correct predictions**. Expansion to a second scaffold batch revealed false positives and a narrower applicability domain.

> **Evidence scope.** These are retrospective computational evaluations against existing wet-lab records. Each reported test construct was excluded from its corresponding fit. The original labels, structures and earlier results had already informed model development; no new wet-lab experiment was conducted in this update. “Work” retains the original scaffold label and is separate from the later TRM C388 activity threshold.

## 1. From experimental observations to an interpretable hypothesis

The wet lab supplied work/non-work labels for redesigned PUF scaffolds. The initial collection contained **14 constructs: 12 PUF12 and two PUF11**, with three work and eleven non-work labels. Predicted structure models supplied contact-probability (CP) features. Features were aggregated model → seed → construct, leaving one observation per independently labeled construct.

The historical decision-stump analysis had found that successful constructs had higher protein-internal, high-probability nonlocal-contact density. Its training-fold feature selection produced **14/14 correct leave-one-construct-out predictions**. That result is retained as the discovery-stage baseline; the present update adds a random forest rather than relabeling the stump as a forest.[S1]

Our structural hypothesis was that a functional scaffold is associated with a sufficiently organized protein-internal contact network. This is an association to test; a high CP density does not demonstrate thermodynamic stability or identify the experimental failure mechanism.

## 2. CP features and random-forest design

For sequence separations d = 4, 12 and 24 residues, we used:

\[
S_d=\frac{1}{L}\sum_{i<j,\;j-i\geq d}\mathbf{1}(CP_{ij}>0.5).
\]

L is the analyzed protein length. Each contact pair is counted once; “nonlocal” refers to separation along the sequence. S can exceed one because it is a count of qualifying pairs per residue, not a probability. The CP threshold is applied per structure model before the established aggregation.

The primary **RF-CP3** model uses S4, S12 and S24. Its fixed settings are **300 trees, maximum depth 2, minimum leaf size 2, sqrt feature subsampling, balanced class weights, seed 2026**. Median imputation and scaling are fitted on training rows only. A work score of at least 0.5 gives the work class. The score is uncalibrated.

A predeclared **RF-CP+structure** companion uses the same forest settings plus core pLDDT mean/minimum, core mean PAE, contact-weighted PAE, length-normalized radius of gyration and anisotropy. These nine inputs follow the repository's expanded scaffold analysis.[S2] No parameter or test subset was selected to make this rerun reach 100%.

## 3. Initial-batch validation: perfect internal classification

For every initial construct, we trained RF-CP3 on the remaining constructs and predicted the omitted one. All structures and seeds for an identity remained together.

| Initial dataset | Model | Correct predictions | Work detected | Non-work rejected | Balanced accuracy |
|---|---|---:|---:|---:|---:|
| 14 constructs, including two PUF11 | Historical depth-one decision tree | 14/14 | 3/3 | 11/11 | 1.000 |
| **Same initial 14-construct collection** | **RF-CP3, newly fitted** | **14/14** | **3/3** | **11/11** | **1.000** |
| **Initial 12 PUF12 only** | **RF-CP3, newly fitted** | **12/12** | **3/3** | **9/9** | **1.000** |

The RF predictions reproduce the initial separation without relying on inclusion of the two shorter PUF11 failures. A separately declared single-S12 forest sensitivity also gave 14/14. Because the feature family was informed by previous exploration, these results describe strong internal agreement, not a universal success rate for new PUF architectures.

[All new metrics](../../results_20260922/scaffold_RF/results/metrics.csv) · [Every held-out prediction](../../results_20260922/scaffold_RF/results/predictions_indexed.csv)

## 4. Learning from a second scaffold batch

The PUF12 expansion contains **24 constructs: 4 work and 20 non-work**, comprising 3/9 in the first batch and 1/11 in the second. We kept the model recipes and class-score cutoff fixed and evaluated progressively different questions.

| Evaluation | RF-CP3 | RF-CP+structure |
|---|---|---|
| Leave one construct out of all 24 | 21/24; BA 0.825; AUC 0.913 | **22/24; BA 0.850; AUC 0.913** |
| Train on first 12; test all second-batch 12 | 8/12; BA 0.818 | **9/12; BA 0.864** |
| Leave a repeat-arrangement group out | 18/24; BA 0.750; AUC 0.869 | **22/24; BA 0.750; AUC 0.938** |
| Train on 20; test the frozen four-construct case study | 2/4; BA 0.667 | **3/4; BA 0.833** |

The architecture-group score is a distinct, harder generalization question. At the fixed 0.5 cutoff, RF-CP+structure detects two of the four work constructs while rejecting all twenty non-work constructs. Its 0.938 AUC must not be described as 93.8% classification accuracy.

### A concrete prediction–measurement comparison

The four-construct case study reuses the previously frozen second-batch panel: its unique work construct plus three randomly selected non-work constructs. All four were excluded simultaneously from training; the panel was not redrawn after prediction.

| Second-batch construct | Existing wet-lab label | RF-CP3 score / prediction | RF-CP+structure score / prediction |
|---|---|---|---|
| **Design 1: R123/R567/R567/R67-no-loop-8** | **Work** | **0.903 / Work** | **0.854 / Work** |
| Design 3: R123/R567/R56-loop-7/R678 | Non-work | 0.158 / Non-work | 0.172 / Non-work |
| Design 7: R123/R567/R567/R-loop-67-loop-8 | Non-work | 0.909 / Work | 0.680 / Work |
| Design 8: R123/R567/R567/R6-loop-7-loop-8 | Non-work | 0.524 / Work | 0.431 / Non-work |

**Design 1 is a positive held-out agreement case:** the model predicts work and the excluded wet-lab record documents work. Design 7 is an important counterexample: strong contact-density scores alone do not guarantee the recorded functional outcome. Adding structural summaries corrects Design 8 but leaves Design 7 unresolved.

## 5. What this module contributes to the engineering cycle
\n
This module provides an interpretable scaffold prioritization signal and identifies cases requiring experimental investigation. The early RF achieves perfect initial-batch classification, while expanded evaluation defines its limits. All 24 main scaffolds retain a shared front R123 background; transfer beyond that design space remains untested. The specificity of reporter editing is assessed separately in the TRM module.

For the next wet-lab cycle, the relevant output is a candidate and a structural hypothesis, accompanied by uncertainty and failure cases. A prospective claim requires predictions to be frozen before new measurements are generated.

## Reproduce and inspect

```bash
python results_20260922/scaffold_RF/retrain.py
```

The rerun uses the small published feature tables; raw AF files are unnecessary for reproducing the classifier fits. [Protocol](../../results_20260922/scaffold_RF/audit/protocol.json), [initial CP inputs](../../results_20260922/scaffold_RF/data/initial14_CP3.csv), and [expanded inputs](../../results_20260922/scaffold_RF/data/scaffolds24.csv) specify the data and decisions. The full downloadable update also retains trained model objects and complete split manifests; these can be regenerated by the script.

[S1] [Historical CP report and decision-stump validation, pinned source](https://github.com/pdx12320/PUF_alphafold/blob/958de2cc3567671c8f9c452cf819aa2133b630d5/previous/results/PUF_CP_report.md).

[S2] [Expanded scaffold feature and validation definitions, pinned source](https://github.com/pdx12320/PUF_alphafold/blob/0d68259e15c3d8d26b308bd938ef8acde16083f5/architecture_validation/PUF_architecture_validation_report.md).
