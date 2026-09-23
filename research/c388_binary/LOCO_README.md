# Nested leave-one-identity-out validation

`loco.py` evaluates the 76 development constructs from the approved fixed-five
protocol. All five fixed test constructs stay excluded from every outer and
inner fold. This is deliberately not an 81-construct LOOCV analysis.

Each outer fold masks the held-out outcomes and repeats all seven biological
thresholds, 46 representations, and training-only top-two-per-family refinement
with the same 11 model configurations as the fixed-five search. Feature selection,
scaling and identity grouping use the outer-training data only. Outer predictions
are sealed before their labels are evaluated. Outputs retain all inner splits,
candidate records, selected settings and per-sample predictions.

The primary result evaluates an adaptive-threshold procedure. Different outer
folds may use different biological boundaries; pooled metrics must not be
presented as performance at one common -5 pp threshold. Pooled AUC is intentionally
omitted because margins/probabilities and biological thresholds differ across
models. Matched feature-family comparators use each primary fold's selected
biological threshold and training-only selection within that family. They are
conditional comparisons, not independent fixed endpoints.

Run from the repository with the same isolated environment and verified v4 inputs:

```bash
python research/c388_binary/loco.py --v4 /path/to/v4 \
  --protocol /path/to/fixed_five/protocol.json --output /new/loco_output
```

Execution is serial with one BLAS thread. The output directory must not exist.
The script does not refit a final deployment model, change source inputs or
rewrite the existing fixed-test results.
