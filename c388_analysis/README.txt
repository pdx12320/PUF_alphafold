C388 >= 50% reanalysis, 2026-09-08

Start with C388_analysis_report.md and per_construct_comparison.csv.
metrics.csv contains frozen, fixed-model, nested-selection and sensitivity results.
No classifier was trained with seeds, models, or multiple RNA contexts as independent samples.

Reproduce from feature caches (run from parent of c388_analysis):
python c388_analysis/analyze.py
python c388_analysis/summarize.py
python c388_analysis/finish.py

Re-extract raw uploads: place five original ZIP files and experimental XLSX in upload/.
python c388_analysis/inventory.py
python c388_analysis/extract_features.py
Only C388 matrices are extracted. The mmCIF parser is strict for supplied AF3 atom_site rows;
it validates 519 protein CA residues and exact atom column counts. It is not a general mmCIF library.
Protein crop is the previously identified SGSETPG N-terminal seven-residue prefix.

Required package versions are recorded in environment_versions.json.
Original ZIP/MSA/template/structure bytes are not duplicated in this results archive.
The provided model and feature caches are sufficient to repeat statistical and ML analysis.

Scoring new construct-level features:
python c388_analysis/predict.py --features new_features.csv --model CP3 --out new_scores.csv
Scores are uncalibrated, exploratory, and specific to this C388 label and scaffold.
JSON/CSV names retain the uploaded construct spellings. _plus_ maps to +;
WT_PUF12-9 maps to experimental PUF12-9. No new architecture validation is claimed.
