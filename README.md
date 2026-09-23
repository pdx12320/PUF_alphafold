# PUF–APOBEC Structure-Guided Engineering

We use structural modelling and reporter measurements to engineer PUF–APOBEC RNA editing. The dry-lab workflow addresses three linked questions: which PUF repeat arrangements retain function, which recognition-motif changes improve editing selectivity, and which residues outside the recognition code merit experimental testing.

**Read the Wiki:** [English](docs/wiki/PUF_Model_EN.md) · [中文](docs/wiki/PUF_Model_ZH.md) · [Design–Build–Test–Learn](docs/DRY_LAB_DBTL.md)

## Three models, one design cycle

| Stage | Biological question | Model and evaluation | Experimental connection |
|---|---|---|---|
| **Model 1 — Scaffold function** | Can a redesigned PUF repeat arrangement retain function? | A decision-tree baseline identifies informative protein-contact features. A shallow random forest combines contact and structural features for the four-construct comparison. | The frozen panel contains one experimentally labelled work construct and three labelled non-work constructs; the model agrees with three of the four recorded outcomes. |
| **Model 2 — TRM-dependent editing** | Can C388 editing be retained while reducing C295/C871 bystander editing? | Site-specific classifiers use changes in protein–RNA contacts to predict editing classes relative to matched controls. | Five experimentally favourable TRM constructs provide 15 site-level comparisons; 13 predictions agree with the measurements. |
| **Model 3 — Non-TRM design** | Which changes outside the RNA-recognition code should be tested next? | Pretrained ProteinMPNN and LigandMPNN sample sequences from the PUF12 structure; their frequency profiles nominate non-TRM substitutions. | Seventeen dual-model consensus substitutions inform an 11-construct G0–G10 test set on experimentally characterized TRM backbones. These designs await experimental assessment. |

The Wiki follows **background → model rationale → model fitting or sampling → experimental comparison → next design**. Models 1 and 2 are evaluated against existing experimental records using held-out predictions. The five Model 2 constructs are each omitted from their own training fold; they are an outcome-selected retrospective panel. Model 3 generates hypotheses for the next wet-lab round.

## Results to display

| Result | Recorded outcome | Source |
|---|---|---|
| Model 1 four-construct panel | Trained on the remaining 20 constructs; TP = 1, TN = 2, FP = 1, FN = 0; **3/4 agreement** | [Predictions](results_20260922/scaffold_RF/results/predictions_indexed.csv) · [Fixed protocol](results_20260922/scaffold_RF/audit/protocol.json) |
| Model 1 broader construct-held-out assessment | CP + structure random forest: **22/24 correct** | [Metrics](results_20260922/scaffold_RF/results/metrics.csv) |
| Model 2 five-construct panel | **13/15 endpoint calls correct**: C295 5/5, C388 4/5, C871 4/5 | [Predictions and measurements](results_20260922/five_construct_holdout/predictions.csv) |
| Model 3 sequence-based nomination | **17 dual-model consensus substitutions**; **G0–G10** proposed constructs | [Candidate table](aice_mpnn_20260922/results/recommended_mutations.csv) · [Design sequences](aice_mpnn_20260922/results/gen3_constructs.fasta) |

See the [Wiki figure set and captions](docs/wiki/figures/captions.md) for confusion matrices, measured editing rates, matched-control changes and the non-TRM candidate summary. The separate [ten-construct ranking analysis](results_20260922/ten_construct_holdout/ranking.csv) is retained with its own evaluation scope in the [results index](results_20260922/README.md).

## From predictions to the next experiments

The Model 2 panel comprises **P8-GVE**, **P9-GNS**, **P9-NPS**, **P9-NTQ** and **P4-R5-SNE+P7-R5-SNE**. Their measured C388 editing spans 60.46–80.34%, while C871 editing decreases by 30.70–57.27 percentage points relative to matched controls. The dual P4/P7 construct also increases C388 editing and reduces both bystander sites.

Model 3 builds on these experimentally characterized recognition motifs and on the high-C388 **P7-SYVIRR** backbone. Proposed additions include **R45T**; the six-residue package **T350V, D374E, M458L, L274I, T278I and V366I**; **A438G/H392D**; and **N12S/S30A/L41R**. The [Model 3 workflow](aice_mpnn_20260922/README.md) documents the 493-residue reference numbering, sampling, filters, sequence checks and proposed tests. Structural roles assigned to these substitutions are design hypotheses.

## Repository map

| Directory | Purpose |
|---|---|
| [`docs/wiki/`](docs/wiki/) | Current English and Chinese Wiki narrative, APA references and display-ready figures |
| [`results_20260922/`](results_20260922/README.md) | Model 1 reproduction script, model outputs, experimental comparisons and split records |
| [`aice_mpnn_20260922/`](aice_mpnn_20260922/README.md) | Model 3 structure inputs, sampling workflow, non-TRM nominations and G0–G10 sequences |
| [`architecture_validation/`](architecture_validation/) | Supporting scaffold analysis with repeat arrangements held out |
| [`data/`](data/) and [`tools/`](tools/) | Input records, provenance and reproduction utilities |

Earlier analyses remain accessible through the [repository guide](docs/REPOSITORY_GUIDE.md) and the dated history in [DRY_LAB_DBTL.md](docs/DRY_LAB_DBTL.md).

## Reproduce the scaffold models

```bash
python -m pip install -r requirements.txt
python results_20260922/scaffold_RF/retrain.py
```

Use the [recorded package versions](results_20260922/scaffold_RF/audit/environment.json) when reproducing the saved results. The [results index](results_20260922/README.md) documents source archives and the distinct validation settings; the [Model 3 README](aice_mpnn_20260922/README.md) provides its sampling and design commands.
