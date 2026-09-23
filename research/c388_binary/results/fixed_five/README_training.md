# C388 binary two-batch search

This is a new runner; it does not overwrite the archived three-class analysis,
source arrays, fitted models, labels or frozen splits. It requires the verified
v4 source archive identified in `../results_20260922/source_archives.json` and its
sibling v3 `data/features_C388_TRM.csv` (loaded by the v4 dataset class, although
TRM features are excluded from this search).

## Task and cohort

The response is the recorded C388 change relative to the batch-matched assay
control, in percentage points. For a selected magnitude `a`, the new labels are:

- `decrease = 0`: delta < −a.
- `unchanged = 1`: delta ≥ −a, **including increases**.

Equality belongs to unchanged. This differs deliberately from the archived
`≤ −a` decrease boundary. Missing outcomes are −1 and never silently labelled
unchanged. The existing magnitudes 2.5, 5, 7.5, 10, 12.5, 15 and 20 are used;
the increase threshold and third class are removed.

All 81 construct-level records have PR (protein–RNA) CP: 30 old_TRM and 51
new_TRM. Preserve the existing five-construct fixed test, including sequence
identity exclusion; both batches contribute to the 76 training records.
Individual seeds/models are already aggregated and are not independent samples.

Only 58 constructs have PP and RR matrices. The missing 23 include two fixed
test constructs. This runner's **all_PR** means all eligible PR entries, not
all PP+PR+RR entries. It does not impute missing blocks or silently restrict
the cohort. A full-complex subset experiment requires a separate approved
protocol and cannot reuse a five-construct performance claim.

## Selection

The 46 archived representations include raw and WT-relative delta CP:

- RNA–protein interface: reference-CP masks or training-union masks at archived
  support values, with all entries or top-variance retention.
- High variance PR: top 25, 100, 250 or 1000 entries.
- All PR: every finite, nonconstant training-eligible PR entry.

Variance ranking, interface union, centering, scaling and exact dual linear
embedding are fitted only on each inner training fold. No truncated PCA or
sequence/TRM inputs are used. Selected missing validation coordinates make a
candidate invalid. Class support follows the archived minimum of three.

For each of seven biological thresholds, screen all 46 representations with
the archived Ridge recipe using three identity-grouped internal folds. Refine
the two best configurations per family with the archived 11 recipes spanning
Ridge, logistic regression, linear/RBF SVM, RF, ExtraTrees and LDA. Archived
dimension limits still apply. Select by macro F1, then MCC, model simplicity
and actual feature count, following the C388 ranking rule in v4. The maximum
budget is 388 candidate evaluations, each with up to three fits, plus refit.

Classifier probability/margin cutoffs remain native (0.5/0); the **biological
editing threshold** is the automatically selected threshold. These are
different quantities. Different thresholds define different prediction tasks;
the search leaderboard is internal development evidence, not independent
validation or proof of biological equivalence.

Fixed-test assay columns are masked before selection. Predictions are saved
and hashed before scoring; the artifact contains the selected CP coordinates,
training means/scales, reference values, exact linear basis and classifier.
The five test outcomes have historical development exposure, so this is a
reconstruction rather than a new prospective validation. A small test panel
also cannot establish broad generalization.

## Commands

Use an isolated environment with the existing repository versions:
numpy 2.3.5, pandas 2.2.3, scipy 1.17.0, scikit-learn 1.8.0, and joblib.
The runner records actual versions and source/input hashes. It uses one BLAS
thread and one tree-estimator worker. Every output directory must be new.

```bash
python research/c388_binary/train.py --v4 /path/to/v4 --output /new/audit --mode audit
C388_V4=/path/to/v4 python -m unittest discover -s research/c388_binary -p 'test_*.py'
python research/c388_binary/train.py --v4 /path/to/v4 --output /new/pilot --mode pilot
# Only after approval for the complete-data calculation:
python research/c388_binary/train.py --v4 /path/to/v4 --output /new/results --mode fixed-five
```

The pilot uses 24 training constructs, three feature families, one threshold
and Ridge/LR/RF recipes. It does not evaluate the fixed test and is solely a
software check. The full run emits protocol, split audits, compressed search
records, internal leaderboard, selected model, selected CP pair table, sealed
predictions, test predictions, metrics and a matched training-majority count.
No Git commit, push or PR is performed by these commands.
