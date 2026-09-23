# PUF–APOBEC: structure-guided RNA editing

### From scaffold selection to selective editing and non-TRM design

We combine protein–RNA structure analysis, interpretable classifiers and inverse folding to guide PUF–APOBEC engineering. Our objective is to retain editing at **C388** while reducing bystander editing at **C295 and C871**.

**[Read the full Dry Lab Wiki →](wiki/PUF_Model_EN.md)** · [Figures and captions](wiki/figures/captions.md) · [Research files](research/README.md) · [GitHub](https://github.com/pdx12320/PUF_alphafold)

## Our modelling workflow

| Design question | Approach | Experimental connection |
|---|---|---|
| Which repeat arrangements retain function? | Contact-based scaffold screening | Compare work/non-work predictions with four recorded outcomes. |
| Which recognition motifs improve editing profiles? | Endpoint-specific TRM classifiers | Compare predictions with C295, C388 and C871 measurements for five constructs. |
| Which surrounding residues merit testing? | ProteinMPNN and LigandMPNN sequence sampling | Nominate non-TRM substitutions for matched-background experiments. |

## 1 · Contact-based scaffold screening

A decision-tree baseline identified informative contact features. A shallow random forest subsequently combined contact and structural descriptors. It agrees with **22/24** labels under leave-one-construct-out evaluation and **3/4** labels when the four comparison constructs are withheld jointly. That panel contains one work and three non-work constructs.

![Scaffold classification under LOCO and joint four-construct holdout](wiki/figures/scaffold_validation_combined.png)

[Model rationale and experimental comparison →](wiki/PUF_Model_EN.md#model-1-contact-based-scaffold-screening)

## 2 · Balancing target and bystander editing

Separate models predict control-relative editing classes at the three reporter sites. For five experimentally favourable constructs, **13 of 15 endpoint predictions** agree with measured classes: C295 **5/5**, C388 **4/5** and C871 **4/5**. These are retrospectively selected examples; each construct was withheld in its own training fold.

![Endpoint-specific classification for five TRM constructs](wiki/figures/fig3_model2_five_construct_confusion.png)

[Editing profiles and model interpretation →](wiki/PUF_Model_EN.md#model-2-balancing-target-and-bystander-editing)

## 3 · Exploring non-TRM sequence space

ProteinMPNN and LigandMPNN each generated **10,000 sequences** from the PUF structural input. After recognition-code filtering, both models nominated the same alternative amino acid at **17 non-TRM positions**. These candidates connect structural sequence preferences to the next round of wet-lab comparisons.

[Sampling workflow and exact substitutions →](wiki/PUF_Model_EN.md#model-3-exploring-non-trm-sequence-space-with-aice)

## Explore the evidence

| Page | What you will find |
|---|---|
| [Full model narrative](wiki/PUF_Model_EN.md) | Background, model explanations, evaluation, experimental connections and APA references |
| [Figure gallery](wiki/figures/captions.md) | Four current figures with captions and editable vector exports |
| [Engineering cycle](research/docs/DRY_LAB_DBTL.md) | Design–Build–Test–Learn across the three models |
| [Research index](research/README.md) | Source data, scripts, outputs, provenance and historical analyses |
| [Reproducibility](research/docs/REPRODUCIBILITY.md) | Environment, validation and execution commands |

The repository is organized into two folders: **`wiki/` for presentation** and **`research/` for supporting analyses**. All displayed measurements and predictions come from the recorded source tables.
