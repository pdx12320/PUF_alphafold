# Scaffold screening: fixed test panel

[Full model page](PUF_Model_EN.md) · [GitHub](https://github.com/pdx12320/PUF_alphafold)

## Model 1. Contact-based scaffold screening

### Why repeat arrangement matters

Changing repeat order or inserting loops can alter the organization of the PUF scaffold. We asked whether protein-internal contacts could distinguish arrangements with different experimental work/non-work outcomes. Here, “work” follows the recorded work/non-work label; C388 editing activity is evaluated separately in Model 2.

### Learning a structural decision rule

We represented each structure by contact-probability (CP) summaries. For a sequence separation of at least $d$ residues, the density of high-probability contacts was

$$
S_d=\frac{1}{L}\sum_{i<j,\;j-i\geq d}\mathbf{1}(CP_{ij}>0.5),\qquad d\in\{4,12,24\}.
$$

Each residue pair contributes once, and $L$ is the analyzed protein length. A larger value indicates more qualifying contacts per residue. Structure-level features were averaged within seeds and then across seeds, giving one feature vector per construct.

A depth-one decision tree established the initial classification rule. Feature selection was repeated within each training fold and consistently selected S12. We then used shallow ensembles of decision trees to combine contact density with additional structural information.

The three-feature random forest used S4, S12 and S24. Its nine-feature companion added core pLDDT mean and minimum, mean PAE, contact-weighted PAE, normalized radius of gyration and anisotropy. Both used 300 trees, maximum depth two, minimum leaf size two, balanced class weights and random seed 2026. Scores of at least 0.5 produced a work prediction. These scores are uncalibrated classifier outputs.

### Training with a separate four-design test panel

The four test designs were set aside before fitting. We trained the nine-feature random forest on the remaining twenty scaffold constructs. Missing-value handling, feature scaling and model fitting used only those training records. The forest recipe and 0.5 classification cutoff were fixed; none was selected using the four test outcomes.

The earlier decision-tree analysis motivated contact-based screening. Its cross-validation results are retained in the research archive. The result displayed here is the fixed four-design test.

### Experimental comparison: one work and three non-work designs

The four designs advanced to wet-lab testing after prediction in the project workflow. The reproducible fixed-test analysis predicts all four with one model fitted to the twenty training constructs, then compares their predictions with the experimental work/non-work outcomes.

| Construct | Work score | Model prediction | Experimental outcome |
|---|---:|---|---|
| Design 1: R123/R567/R567/R67-no-loop-8 | 0.854 | Work | **Work** |
| Design 3: R123/R567/R56-loop-7/R678 | 0.172 | Non-work | Non-work |
| Design 7: R123/R567/R567/R-loop-67-loop-8 | 0.680 | Work | Non-work |
| Design 8: R123/R567/R567/R6-loop-7-loop-8 | 0.431 | Non-work | Non-work |

The experimental panel contained **one work construct and three non-work constructs**. Predictions agreed with three of the four outcomes. Design 1 received the highest score and received a work outcome in the wet-lab record. Design 7 remained a false positive, showing why scaffold prioritization still requires experimental testing.

![scaffold_validation_combined](figures/fixed4_scaffold_test.png)

*Figure 1. Four scaffold designs evaluated as one fixed test panel. (a) Work scores from the model fitted to twenty separate training constructs; the dashed line marks the fixed 0.5 cutoff. Colours show the subsequently measured class. (b) Comparison with experimental labels: one true positive, two true negatives and one false positive. All four designs were excluded together from preprocessing and fitting.*

![scaffold_s12_distribution](figures/training20_s12.png)

*Figure 2. S12 distribution in the twenty-construct training set: seventeen non-work and three work constructs. Each point represents one training construct; horizontal lines show medians. The four test designs are absent. S12 denotes high-confidence nonlocal contacts per residue; horizontal offsets only separate points.*

### Scaffold selection for experimental testing

The model provides a structural criterion for prioritizing repeat arrangements and identifies Design 1 as a candidate with a matching recorded work outcome. The false positive motivates checking additional determinants of construct function. Once a functional scaffold is available, the next question concerns its editing profile.

Current fixed-test records: [train/test split](../research/fixed_test_20260923/model1/split.json) and [four test predictions](../research/fixed_test_20260923/model1/predictions.csv).

Earlier development records: [scaffold protocol](../research/results_20260922/scaffold_RF/audit/protocol.json), [metrics](../research/results_20260922/scaffold_RF/results/metrics.csv), [predictions](../research/results_20260922/scaffold_RF/results/predictions_indexed.csv) and [historical decision-tree analysis](https://github.com/pdx12320/PUF_alphafold/blob/958de2cc3567671c8f9c452cf819aa2133b630d5/previous/results/PUF_CP_report.md).
