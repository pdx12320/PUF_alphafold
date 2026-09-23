# Repository guide

## Recommended reading order

1. [Project overview](../README.md): the three biological questions and their experimental connections.
2. [Wiki — English](wiki/PUF_Model_EN.md) or [Wiki — 中文](wiki/PUF_Model_ZH.md): background, model rationale, model results, experimental comparisons and the next design round, with APA references.
3. [Wiki figures and captions](wiki/figures/captions.md): confusion matrices, editing measurements and non-TRM candidates.
4. [Design–Build–Test–Learn](DRY_LAB_DBTL.md): the current engineering cycle followed by the dated development history.
5. [Current results index](../results_20260922/README.md): numerical outputs, split definitions and reproducibility details.

## Model 1 — Scaffold function

The initial decision-tree analysis motivated a low-dimensional contact-based model. The current four-construct comparison uses a shallow random forest trained with all four panel members excluded. Its recorded experimental labels are one work and three non-work; the CP + structure model makes three correct calls. This panel and the broader construct-held-out analysis have separate rows in the result files.

| File | Contents |
|---|---|
| [Protocol](../results_20260922/scaffold_RF/audit/protocol.json) | Features, random-forest settings, threshold and frozen four-construct membership |
| [Metrics](../results_20260922/scaffold_RF/results/metrics.csv) | Performance by dataset, model and validation setting |
| [Predictions](../results_20260922/scaffold_RF/results/predictions_indexed.csv) | Construct-level scores, calls and training-set size |
| [Construct lookup](../results_20260922/scaffold_RF/results/construct_ids.csv) | Mapping from indexed results to construct identities |
| [Reproduction script](../results_20260922/scaffold_RF/retrain.py) | Re-runs the fixed random-forest recipes and saved validation settings |
| [Architecture validation](../architecture_validation/PUF_architecture_validation_report.md) | Supporting assessment with repeat arrangements held out |

The focused [English](wiki/PUF_Module1_EN.md) and [Chinese](wiki/PUF_Module1_ZH.md) pages reproduce the current scaffold section. They link to the earlier detailed evaluation record. Use the full Wiki pages for the three-model narrative.

## Model 2 — TRM-dependent reporter editing

The main illustrative panel comprises P8-GVE, P9-GNS, P9-NPS, P9-NTQ and P4-R5-SNE+P7-R5-SNE. Each construct is held out in its own fold, and experimental editing is compared with its matched control. The reported 13/15 agreement counts the three endpoints for each construct. This retrospective panel was selected for its experimentally favourable behaviour.

| File | Contents |
|---|---|
| [Five-construct panel](../results_20260922/five_construct_holdout/predictions.csv) | Measured rates, control-relative changes, selected class thresholds and predicted classes for all 15 endpoints |
| [All-construct predictions](../results_20260922/v4/all81_LOCO.csv) | Full construct-held-out result table from which the five-construct panel is drawn |
| [Main metrics](../results_20260922/v4/main_metrics.csv) | Broader performance and matched baseline comparisons |
| [Ten-construct comparison](../results_20260922/ten_construct_holdout/ranking.csv) | A separate, simultaneous ten-construct holdout ranking analysis |
| [Ten-construct membership](../results_20260922/ten_construct_holdout/manifest.csv) | Membership of that separate ranking panel |

The five-construct classification panel and the ten-construct ranking analysis use different splits and outcomes. Their results should be read with the corresponding table and evaluation definition.

## Model 3 — Non-TRM mutation nomination

Pretrained ProteinMPNN and LigandMPNN each sample 10,000 sequences from the PUF12 structural input. Empirical amino-acid frequencies identify alternative residues, after which recognition-code positions are filtered. Seventeen substitutions receive dual-model support. The G0–G10 sequences combine selected nominations with experimentally characterized recognition-motif backbones; their editing effects remain to be measured.

| File | Contents |
|---|---|
| [Workflow](../aice_mpnn_20260922/README.md) | Structure preparation, sampling, screening, experimental context and next-round tests |
| [Structure and configuration](../aice_mpnn_20260922/inputs/) | PUF12–RNA complex and the 493-residue reference mapping |
| [Position-level screen](../aice_mpnn_20260922/results/aice_single_ranked.csv) | Full sampling-frequency scan |
| [Non-TRM nominations](../aice_mpnn_20260922/results/recommended_mutations.csv) | Candidate substitutions, support and priority categories |
| [Experimental comparison](../aice_mpnn_20260922/results/wetlab_vs_prediction.csv) | Previously measured TRM variants used to contextualize model output |
| [G0–G10 design table](../aice_mpnn_20260922/results/gen3_constructs.csv) and [FASTA](../aice_mpnn_20260922/results/gen3_constructs.fasta) | Exact proposed mutations and corresponding sequences |

## Earlier analysis snapshots

| Directory | Focus |
|---|---|
| [`c388_threshold/`](../c388_threshold/README.md) | C388 dynamic-threshold classification |
| [`c388_local_nonlocal_optimization/`](../c388_local_nonlocal_optimization/README.md) | Local and non-local contact-feature comparisons |
| [`trm22_offtarget/`](../trm22_offtarget/README.md) | C295/C871 analysis in the frozen 22-variant cohort |
| [`combined12/`](../combined12/) | Earlier combined PUF12 scaffold analyses |
| [`new_batch/`](../new_batch/) | Additional experimental-batch analyses |
| [`data/architecture_inputs/`](../data/architecture_inputs/README.md) | Preserved scaffold inputs and provenance |

The dated section of [DRY_LAB_DBTL.md](DRY_LAB_DBTL.md) records how these analyses informed subsequent choices. Numerical values in an earlier snapshot belong to that snapshot's dataset, labels and validation scheme.
