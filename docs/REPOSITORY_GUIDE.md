# Repository guide

## Recommended reading order

1. [Root README](../README.md) for the project overview and current results.
2. [Full Model page](wiki/PUF_Model_EN.md) for the three-model dry-lab narrative.
3. [Model 1 page](wiki/PUF_Module1_EN.md) for scaffold-function prediction.
4. [AiCE–MPNN workflow](../aice_mpnn_20260922/README.md) for non-TRM mutation design.
5. [Current results index](../results_20260922/README.md) for model outputs and reproducibility files.

## Active project directories

| Directory | Contents | Use |
|---|---|---|
| `docs/wiki/` | Wiki-ready Model pages | Primary project narrative |
| `results_20260922/` | Current Random Forest, TRM and holdout outputs | Current result tables and scripts |
| `aice_mpnn_20260922/` | AlphaFold 3 PUF12–RNA complex, MPNN sampling, repeat-level comparison and G0–G10 constructs | Model 3 and next-round designs |
| `architecture_validation/` | Scaffold contact-probability validation | Model 1 supporting analysis |
| `data/` | Input tables and construct records | Source data |
| `tools/` | Validation and reproduction utilities | Reproducibility |

## Earlier analysis snapshots

The following directories preserve earlier model-development runs and their original file paths:

| Directory | Focus |
|---|---|
| `c388_threshold/` | C388 dynamic-threshold classification |
| `c388_local_nonlocal_optimization/` | Local and non-local contact optimization |
| `trm22_offtarget/` | C295/C871 off-target analysis |
| `combined12/` | Combined PUF12 analysis |
| `new_batch/` | Additional experimental batch analysis |

## Result selection

Use `results_20260922/` for current reported outputs. Use the earlier directories when reproducing the corresponding historical analysis or tracing a result to its original workflow.
