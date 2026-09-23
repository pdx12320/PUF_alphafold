# C388 contact-feature comparisons

This snapshot compares contact-feature families for C388 activity retention in the frozen 22-variant cohort. The endpoint is a control-relative editing change of at least −15 percentage points. Fixed classifier cutoffs and cutoffs selected within training folds are reported separately.

At the fixed 0.5 cutoff, Interface_CP with a depth-2 random forest achieves mutation-group-out balanced accuracy 0.9375 and LOOCV balanced accuracy 0.854167. Its group-out balanced accuracy falls to 0.583333 when the score cutoff is selected within training folds. Joint training-only selection of features, models and cutoffs achieves group-out balanced accuracy 0.5. The high fixed-candidate result belongs to an exploratory comparison.

Interface_CP denotes the RMS difference of the protein–RNA contact-probability matrix relative to the matching WT, using the retained 519×15 normalization. The cached core has 512 residues. It summarizes protein–RNA perturbation and does not measure protein–protein nonlocal coupling. Missing full protein–protein matrices prevent several originally proposed nonlocal comparisons.

Use `python verify_results.py` from this directory to check the retained snapshot. Full reruns require the separately retained source archive and its verified inputs. The tables preserve the recorded comparisons.

This condensed English guide retains the numerical tables below. The [complete original report at the pre-reorganization commit](https://github.com/pdx12320/PUF_alphafold/blob/fe61ff5cd05e53f7e3b3dab169017f940aa46067/c388_local_nonlocal_optimization/README.md) preserves the historical prose and analysis context.

### Table 1. Fixed and training-selected cutoff comparisons

| Features and model | Score protocol | Group-out BA | LOO BA | Group-out AUC |
|---|---|---:|---:|---:|
| Interface_CP / RF_depth2 | Fixed 0.5 | 0.937500 | 0.854167 | 0.947917 |
| Total_CP / LR_C0.1 | Fixed 0.5 | 0.687500 | 0.718750 | 0.947917 |
| Total_CP / LR_C1 | Fixed 0.5 | 0.750000 | 0.750000 | 0.927083 |
| Total_CP / RF_depth2 | Fixed 0.5 | 0.708333 | 0.708333 | 0.796875 |
| CP_summary / LR_C1 | Fixed 0.5 | 0.875000 | 0.791667 | 0.947917 |
| C388_pm0 / LR_C0.1 | Fixed 0.5 | 0.687500 | 0.687500 | See CSV |
| Interface_CP / RF_depth2 | Training-selected score cutoff | 0.583333 | 0.854167 | 0.947917 |
| Interface_CP / nested_model_selection | Training-selected model and cutoff | 0.500000 | 0.770833 | 0.875000 |
| nested_best_available | Training-selected features, model and cutoff | 0.500000 | 0.656250 | 0.687500 |
