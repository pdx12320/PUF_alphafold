# Residue, repeat and exact-mutation interpretation

This analysis reads the frozen Model 2 artifacts and the original Model 3 samples. It does not refit classifiers, change thresholds or alter the fixed test panel.

## Model 2

Random-forest contact importance uses mean decrease in impurity. LinearSVM and Ridge coefficients are mapped from the dual embedding back to the standardized original feature coordinates, then summarized by mean absolute coefficient across classes. Importance is summed over protein–RNA contacts for each residue and each repeat P, and normalized within each endpoint's CP block. Mean importance per selected contact and contact counts are also supplied because summed importance depends on feature count.

These are fitted-model attributions, not causal mutation effects. Endpoints use different estimators, so normalized profiles are interpreted separately. C388 also includes sequence descriptors: CP accounts for 43.59% of total absolute standardized coefficient weight, with the remaining 56.41% in the sequence block. Its P chart displays only the CP block.

Model 2 reference residues 17–509 map exactly to the 493-residue Model 3 sequence after subtracting 16. The script verifies the complete WT sequence match before mapping. P=0 denotes residues outside the twelve canonical repeat cores, including terminal regions and the inter-repeat insertion. The crosswalk retains both numberings.

## Model 3

All 10,000 generated sequences per model are read; native FASTA records are excluded by header identifier. Repeat-level WT retention is the mean of position-specific WT frequencies. TRM retention uses only the three recognition positions in each repeat. High retention describes a strong structural sequence preference and does not measure functional essentiality.

The 17 consensus non-TRM substitutions are mapped to P using the same verified sequence crosswalk. Their exact substitutions are compared with all 81 experimentally measured variant identities reconstructed from the residue-alignment table. There are **zero exact non-TRM consensus matches** in the current experimental cohort.

A separate table compares the exact TRM substitutions observed experimentally with the MPNN samples. Counts require all changed residues of a construct to co-occur within one generated sequence; WT counts use the same changed positions. Other sampled positions are unconstrained, so this is an exact substitution-set comparison, not complete variant-sequence recovery. Zero counts mean not observed among 10,000 sequences.

The displayed positive assay examples are selected descriptively from the complete table using C388 change ≥ −10 pp, C295 change ≤ 0 pp and C871 change ≤ −10 pp. This selection is not a new test or a model-performance estimate. Raw mutant and control values, every comparison, exact substitution identities and the selection flag are retained in `model3_exact_mutation_wetlab_comparison.csv`.

The original `wetlab_vs_prediction.csv` mixes an earlier assay snapshot and an approximate triplet parser (including an abbreviated label for SYVIRR). Current comparisons instead use the verified v4 residue mapping and the same matched-control measurements used by Model 2, and recompute sample frequencies from all 10,000 sequences. The original table is preserved for provenance.

## Reproduce

```bash
python research/model_interpretation_20260923/analyze.py \
  --v4 /path/to/PUF_alphafold/retrain_classification_v4_dynamic_fullCP_20260921
```

The source archive and model freeze are documented in [fixed_test_20260923](../fixed_test_20260923/README.md). The script never reads test labels when computing Model 2 feature importance. Experimental outcomes enter only the separate descriptive Model 3 comparison.
