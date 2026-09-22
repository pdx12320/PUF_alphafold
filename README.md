# PUF–APOBEC Structure-Guided Engineering

Computational design and validation framework for PUF–APOBEC RNA editing. The repository connects protein structure, TRM editing measurements and AI-guided non-TRM mutation design.

## Project entry points

| Module | Question | Main outcome | Documentation |
|---|---|---|---|
| **Model 1 — PUF scaffold function** | Which repeat arrangements retain scaffold function? | Contact-probability and structural features prioritize work/non-work PUF designs; predictions were assessed with wet-lab function measurements. | [Model 1](docs/wiki/PUF_Module1_EN.md) |
| **Model 2 — TRM-dependent editing** | Which TRM changes retain C388 editing and reduce C295/C871 bystanders? | Candidate TRMs were predicted and assessed using reporter editing measurements. | [Full Model page](docs/wiki/PUF_Model_EN.md#model-2-trm-dependent-reporter-editing) |
| **Model 3 — AI-guided non-TRM design** | Which residues outside the recognition code should be modified next? | AiCE-inspired ProteinMPNN/LigandMPNN scanning defines core, geometry and interface mutations for the next wet-lab design round. | [Full Model page](docs/wiki/PUF_Model_EN.md#model-3-ai-guided-puf-design-outside-trm-positions) · [AiCE–MPNN workflow](aice_mpnn_20260922/README.md) |

## Key experimental design backbones

- **P8-GVE**: low C871 editing with retained C388 editing.
- **P9-NTQ, P9-NPS and P9-GNS**: complementary P9 TRM backbones with retained target editing.
- **P4-R5-SNE+P7-R5-SNE**: increased C388 editing with reductions at C295 and C871.
- **P7-SYVIRR**: high-C388 backbone for non-TRM AI-guided combinations.

## Current results

| Analysis | Result | Files |
|---|---|---|
| Initial scaffold Random Forest | 14/14 correct construct-level predictions | [Results](results_20260922/scaffold_RF/results/metrics.csv) |
| Expanded scaffold assessment | CP + structure Random Forest: 22/24 correct | [Results](results_20260922/scaffold_RF/results/metrics.csv) |
| Five candidate TRM panel | 13/15 site-level predictions | [Predictions](results_20260922/five_construct_holdout/predictions.csv) |
| Ten-construct comparison panel | Prediction and experimental rankings | [Ranking](results_20260922/ten_construct_holdout/ranking.csv) |
| AiCE–MPNN non-TRM design | 17 dual-model consensus substitutions; G0–G10 design set | [Workflow and designs](aice_mpnn_20260922/README.md) |

## Repository map

```text
docs/
  wiki/                   Model pages and Wiki-ready narrative
  REPOSITORY_GUIDE.md     Detailed directory guide
results_20260922/         Current reproducible model outputs and split manifests
aice_mpnn_20260922/       Model 3: AF3 complex, MPNN workflow and G0–G10 designs
architecture_validation/  Scaffold validation analysis
data/                     Input tables and registries
tools/                    Reproducibility and validation utilities

c388_threshold/
c388_local_nonlocal_optimization/
trm22_offtarget/
combined12/
new_batch/                Earlier analysis snapshots retained for provenance
```

## Reproduce

```bash
python -m pip install -r requirements.txt
python results_20260922/scaffold_RF/retrain.py
```

See [the repository guide](docs/REPOSITORY_GUIDE.md) for the role of each analysis directory and direct links to the current model outputs.
