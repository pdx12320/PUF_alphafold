# PUF results update — 22 September 2026

This directory publishes a **curated, outcome-complete result snapshot** from the latest uploaded archives, together with a newly executed scaffold random-forest rerun. It does not declare one universal “best” model across incompatible endpoints or validation schemes. Existing repository analyses are retained.

## Start here

- [Rewritten Module 1 — English](../docs/wiki/PUF_Module1_EN.md) / [中文](../docs/wiki/PUF_Module1_ZH.md)
- [Scaffold RF metrics](scaffold_RF/results/metrics.csv), [all prediction rows](scaffold_RF/results/predictions_indexed.csv), [construct-ID lookup](scaffold_RF/results/construct_ids.csv)
- [v4 all-81 construct-held-out predictions and measurements](v4/all81_LOCO.csv)
- [v4 main and matched-Dummy metrics](v4/main_metrics.csv)
- [Designated-five predictions](five_construct_holdout/predictions.csv)
- [Simultaneous ten-construct holdout ranking](ten_construct_holdout/ranking.csv)
- [Separate fixed-recipe first-batch versus two-batch comparison](wiki_two_batch/TRM_heldout_metrics.csv)

## Results, with their actual evaluation scopes

| Analysis | Result | Interpretation |
|---|---|---|
| Initial 14 scaffolds, newly fitted RF-CP3, construct-held-out | **14/14 correct; BA/AUC 1.000** | Internal retrospective evaluation; includes two PUF11 failures |
| Initial 12 PUF12 only, newly fitted RF-CP3 | **12/12 correct** | Confirms the result does not require the two PUF11 examples |
| Expanded 24, CP+structure RF, construct-held-out | 22/24 correct; BA 0.850; AUC 0.913 | Same fixed shallow-RF recipe, broader dataset |
| Expanded 24, CP+structure RF, architecture-held-out | 22/24 correct; BA 0.750; AUC 0.938 | Detects 2/4 work, rejects 20/20 non-work at cutoff 0.5 |
| Frozen one-work/three-non-work scaffold panel | CP3 2/4; CP+structure 3/4 | All four excluded simultaneously; one false positive remains |
| All 81 TRMs, v4 combined-PR LOCO | C295 57/81; C388 36/81; C871 49/81 correct | Fold-specific dynamic experimental thresholds |
| Five designated TRMs, separate construct holdouts | **13/15 endpoint calls correct** | Same-position training allowed; C295 5/5, C388 4/5, C871 4/5 |
| Ten TRMs held out simultaneously | Composite Spearman 0.103; designated Top-5 overlap 3/5 | Outcome-selected retrospective comparison; enrichment p=0.50 |

The five-construct experiment and ten-construct experiment use different training sets. The former must not be substituted for the latter. The five constructs were already known to the team; none of these analyses constitutes newly completed prospective wet-lab validation.

### High-scoring CP families need their own baseline

v4's independently adaptive nonlocal C871 branch has PP-cohort position-held-out BA **0.885**, but a training-majority rule under that branch's own selected thresholds reaches **0.962**. The analogous C295 values are **0.669 versus 0.869**. These training-only baselines were independently reconstructed from saved split/label records in this update; [all values](v4/nonlocal_own_threshold_baseline.csv) are published alongside the attractive scores.

Different folds may choose different biological boundaries. Neither a large pooled BA nor a large pooled AUC alone establishes improvement attributable to CP. The main v4 PP+PR+RR matched comparison is recorded in the source archive; its high BA must likewise be evaluated against the corresponding matched baseline.

## What is actually included

Small readable CSVs, the rewritten Wiki module, an executable fixed-recipe scaffold rerun, original input feature snapshots, a frozen four-scaffold panel, and source archive hashes are tracked. The all81 CSV preserves measured fractions, matched controls, delta values, selected thresholds, true/predicted classes and native scores for all three sites. Values remain at source precision.

v4 label coding: C295/C871 `0=decrease`, `1=not_decrease`; C388 `0=decrease`, `1=within_selected_interval`, `2=increase`. `a` is the positive decrease magnitude in pp and `b` the positive increase magnitude. Scores may be uncalibrated probabilities or native decision scores; do not combine them as a calibrated common-scale ranking.

The ten-construct ranking uses its separately defined common-event scoring workflow. Its `S_pred = score_C388 + (score_C295 + score_C871)/2`; observed `S_exp = delta_C388 - (delta_C295 + delta_C871)/2`. Dynamic class labels and ranking scores are separate outputs.

**Large raw AF/CP archives and all historical model pickles are not uploaded into Git.** They are identified by [source_archives.json](source_archives.json). Use [extract_source_archives.py](extract_source_archives.py) with the locally retained ZIPs to recover original reports, code, predictions, matrix inputs and models. The extraction checks ZIP hashes and prevents path traversal. Missing matrices are never fabricated. The full downloadable delivery accompanying this update retains the new RF weights and complete split logs.

## Reproduce the new RF fits

```bash
python results_20260922/scaffold_RF/retrain.py
```

Run with the recorded package versions in [environment.json](scaffold_RF/audit/environment.json). The script writes the full prediction CSV, metrics, split records and fitted holdout models inside `scaffold_RF/`. It does not change the archived TRM results or the user-supplied experimental labels.

## Source and validation boundaries

The repository base was commit `0d68259e15c3d8d26b308bd938ef8acde16083f5`. The original decision-stump/PUF11 evidence is pinned to `958de2cc3567671c8f9c452cf819aa2133b630d5`. The RF rerun uses the original model recipe, with fixed seed and no post-hoc test-panel replacement. Models and labels have historical development exposure. Shared-R123 scaffold background, incomplete old-TRM PP, unmatched new-batch WT structure, rare C388 increases and batch differences remain relevant limitations.
