# Model 3 — AiCE-inspired design beyond the PUF recognition code

Model 3 extends PUF engineering to residues outside the tripartite recognition motifs (TRMs). It uses pretrained inverse-folding models to propose sparse changes that can be combined with experimentally characterized TRM backgrounds. The current output is a computational design set for the next wet-lab round.

## Background and design question

TRM substitutions alter RNA recognition, while the surrounding PUF scaffold provides the structural setting in which recognition occurs. After examining functional scaffolds in Model 1 and editing profiles in Model 2, we asked which additional scaffold positions could be explored without directly replacing the selected TRMs.

We adapted **AiCE — AI-informed constraints for protein engineering** (Fei et al., 2025). AiCE uses sequence samples from inverse-folding models to nominate mutations. Our implementation adds a PUF recognition-code filter and a comparison with the project's existing editing measurements. Neither ProteinMPNN nor LigandMPNN was fine-tuned on these measurements.

## Workflow

1. **Prepare the structural input.** Use the AlphaFold 3-predicted PUF12–RNA complex in [`inputs/complex.pdb`](inputs/complex.pdb): protein chain A contains 493 residues and RNA chain R contains the 17-nt context `ACAUGGAGGACGUGCGC`. In the project numbering, RNA C388 is nucleotide 15 of this context. Protein residue 388 and RNA C388 refer to different coordinate systems.
2. **Sample two sequence ensembles.** ProteinMPNN conditions on the protein backbone. LigandMPNN additionally represents the RNA atom environment. Each generates 10,000 protein sequences at temperature 0.5 with seed 111. Each archived FASTA also includes one native sequence, giving 10,001 records per file.
3. **Calculate residue frequencies.** At each protein position, calculate the empirical frequency of the native amino acid, `f_WT`, and the most frequent alternative, `f_alt`, separately for each model. A model nominates a position when `f_alt > f_WT` and `f_alt >= 0.8`. The script also defines a relaxed threshold of 0.5 for flexible positions; all exported flexibility flags are False in this release, so the reported nominations use the 0.8 threshold.
4. **Apply the PUF recognition-code filter.** The sampling run allows all protein positions to vary. During screening, changes at the 36 recognition positions must retain an allowed triplet for the configured target base. No recognition-position substitution passes the final screen. The final sparse designs retain the chosen experimental TRM background and add selected non-TRM substitutions.
5. **Prioritize sites and combinations.** Retain positions nominated by either model, highlight exact substitutions supported by both, and annotate RNA proximity and the N-terminal region. An exploratory combination module calculates weighted covariance and an amino-acid Cramér's V association score across the LigandMPNN samples.
6. **Prepare the wet-lab design set.** Combine selected non-TRM proposals with the experimentally characterized P8-GVE, P9-NTQ and P7-SYVIRR backgrounds. [`gen3_constructs.fasta`](results/gen3_constructs.fasta) and [`gen3_constructs.csv`](results/gen3_constructs.csv) contain 11 sequences, G0–G10.

The combination module is an adaptation of the AiCE idea: its Cramér's V calculation differs from the original nucleotide-level linkage-disequilibrium workflow. The G0–G10 packages are explicitly assembled in `make_gen3.py`; they are not an experimentally optimized output of the combination score.

## Sampling results

The 493-position scan nominated **81 non-TRM positions**, including **17 exact substitutions supported by both models**. The table reports the fraction of generated sequences containing the indicated substitution at its position. These frequencies describe a model's sequence preference under the supplied structure; they do not measure an editing rate, a probability of improvement, or a causal contribution to function.

| Consensus substitution | LigandMPNN frequency | ProteinMPNN frequency |
|---|---:|---:|
| T350V | 0.9968 | 0.9346 |
| D374E | 0.9915 | 0.8926 |
| V294I | 0.9915 | 0.8443 |
| M458L | 0.8826 | 0.9863 |
| L274I | 0.9845 | 0.8644 |
| T278I | 0.9770 | 0.8679 |
| L166I | 0.9662 | 0.9255 |
| V366I | 0.9636 | 0.8585 |
| R93K | 0.9203 | 0.9371 |
| R45T | 0.8948 | 0.9221 |
| E351L | 0.8282 | 0.9164 |
| R462L | 0.9109 | 0.8730 |
| V208E | 0.9101 | 0.8659 |
| A136E | 0.9073 | 0.8339 |
| V388E | 0.8983 | 0.8152 |
| A244E | 0.8400 | 0.8941 |
| S100E | 0.8907 | 0.8266 |

Source: [`results/aice_single_ranked.csv`](results/aice_single_ranked.csv), filtered to `consensus=True`. All 17 rows are outside the configured TRM positions and have the same preferred alternative amino acid in both models.

The stored position-ranking score is:

```text
rank_score = max(f_alt_LigandMPNN, f_alt_ProteinMPNN)
             + 0.2 × consensus
             + 0.1 × RNA_interface_flag
```

This is a heuristic prioritization score. `RNA_interface_flag` uses a protein Cα-to-RNA atom distance below 10 Å in the input structure. The `near_c388_zone` label is a configured N-terminal interval, residues 1–45; it is not a measured effect on C388 editing.

### Selected non-TRM proposals and their experimental questions

| Proposal or package | Residues | Question for wet-lab testing |
|---|---|---|
| N-terminal consensus site | R45T | Does an N-terminal scaffold change alter the editing profile on the selected TRM background? |
| Six-site consensus package (`core`) | T350V, D374E, M458L, L274I, T278I, V366I | Can these jointly sampled preferences be combined while preserving expression and editing activity? |
| RNA-proximal package (`iface`) | A438G, H392D | Do changes near RNA in the input structure improve the measured editing profile? |
| N-terminal package (`geo`) | N12S, S30A, L41R | How do selected N-terminal changes affect the editing profile? |

The package names `core`, `iface` and `geo` identify design groups. Stability, affinity and editing geometry remain experimental questions. In this structure, the nearest non-hydrogen atom distances to RNA are approximately 2.19 Å for N12, 6.73 Å for S30, 3.76 Å for L41, 5.26 Å for R45, 6.66 Å for H392 and 7.14 Å for A438. These distances refer to the native input structure. R45 has a minimum heavy-atom distance of 5.26 Å but a Cα distance of 11.24 Å, explaining its False interface flag under the Cα-based rule.

**Single-model export detail.** The 81-position summary sometimes labels a row with the LigandMPNN-preferred amino acid even when only the ProteinMPNN alternative passes the threshold. At position 12, LigandMPNN prefers S with frequency 0.4436, whereas ProteinMPNN nominates G with frequency 0.8035. Thus **N12S in G5/G6 is an explicit design choice; the threshold-passing ProteinMPNN nomination is N12G**. The same distinction applies to the printed M433A and V457L rows, whose threshold-passing ProteinMPNN alternatives are M433L and V457I. For any single-model nomination, consult the model-specific residue and frequency columns in `aice_single_ranked.csv`. This issue does not affect the 17 consensus substitutions above.

## Connection to experimental evidence

[`wetlab_vs_prediction.csv`](results/wetlab_vs_prediction.csv) compares model-derived frequencies with **84 existing variant records**. This is a retrospective comparison that contextualizes TRM background selection. The data show useful editing profiles among variants at P7, P8 and P9, with unequal numbers of variants per position; they do not establish a position-wide causal effect.

No new wet-lab results for the nominated non-TRM substitutions or the G1–G10 combined designs are included in this release. The absence of TRM nominations after filtering also does not establish that the native TRMs are optimal: sampling and recognition-code filtering answer different questions from an editing assay. The experimentally characterized TRM backgrounds supply the reference points for the proposed next round.

## G0–G10 design set for wet-lab evaluation

All coordinates below are 1-based positions in the 493-residue protein. Each sequence was generated with assertions that verify the starting amino acid before substitution.

- **GVE background:** S288G + Y289V.
- **NTQ background:** Y325T.
- **SYVIRR background:** E256R; residue 257 is already R in the reference sequence.
- **Core package:** T350V + D374E + M458L + L274I + T278I + V366I.
- **Interface package:** A438G + H392D.
- **N-terminal package:** N12S + S30A + L41R.

| Construct | Experimental TRM background | Added non-TRM changes | Total substitutions from reference |
|---|---|---|---:|
| G0_P8GVE_ctrl | GVE | None | 2 |
| G1_GVE+R45T | GVE | R45T | 3 |
| G2_GVE+core | GVE | Core package | 8 |
| G3_GVE+R45T+core | GVE | R45T + core package | 9 |
| G4_GVE+iface | GVE | Interface package | 4 |
| G5_GVE+geo | GVE | N-terminal package | 5 |
| G6_GVE+all | GVE | R45T + core + interface + N-terminal packages | 14 |
| G7_NTQ+core | NTQ | R45T + core package | 8 |
| G8_NTQ+iface | NTQ | R45T + interface package | 4 |
| G9_GVE+NTQ+core | GVE + NTQ | R45T + core package | 10 |
| G10_SYVIRR+core | SYVIRR | Core package | 7 |

G7 and G8 include R45T even though their short identifiers do not mention it. G9 combines two experimentally characterized TRM changes; their combined background has not been established by the separate measurements. Comparisons against the corresponding parental background are needed to assign the effect of each added package. G0 provides the P8-GVE reference; the NTQ, combined GVE+NTQ and SYVIRR comparisons require their matched parental controls in the experiment.

The next round should measure C388, C295 and C871 editing under matched conditions, alongside expression or abundance checks where available. Single-substitution comparisons can resolve the contribution of individual residues, while package comparisons test whether their effects combine. G6 tests the largest combination in this design set.

## Reproduction and files

### Re-analyse the archived sequence samples

Run from `aice_mpnn_20260922/code/` in an environment with the Python dependencies installed. This route preserves the published result tables and writes a separate reanalysis directory.

```bash
mkdir -p workdir/archive_samples workdir/reanalysis
cp ../inputs/complex.pdb workdir/complex.pdb
cp ../inputs/design_config.json workdir/design_config.json
gzip -dc ../results/mpnn_samples_proteinmpnn.fa.gz > workdir/archive_samples/proteinmpnn.fa
gzip -dc ../results/mpnn_samples_ligandmpnn.fa.gz > workdir/archive_samples/ligandmpnn.fa
python step3_aice_screen.py \
  --pdb workdir/complex.pdb \
  --config workdir/design_config.json \
  --pmpnn workdir/archive_samples/proteinmpnn.fa \
  --lmpnn workdir/archive_samples/ligandmpnn.fa \
  --outdir workdir/reanalysis
```

Secondary-structure annotation depends on the available DSSP implementation and should be checked when comparing a reanalysis with the archived flags. The screen regenerates the position and exploratory-combination tables; `recommended_mutations.csv`, the retrospective wet-lab comparison and the repeat-level summary are archived downstream outputs.

### Generate a new sequence ensemble

```bash
cd aice_mpnn_20260922/code
bash setup_env.sh
conda activate puf_aice
mkdir -p workdir
cp ../inputs/complex.pdb workdir/complex.pdb
python step0_prepare.py --outdir workdir
bash step2_sample_mpnn.sh 10000 0.5
python step3_aice_screen.py --outdir workdir/reanalysis
```

Setup requires access to external source repositories, Python packages and model weights. These dependencies are downloaded by the setup scripts and are not bundled in this directory. Compatibility depends on the installed LigandMPNN and dependency versions. `python make_gen3.py` rebuilds the archived G0–G10 FASTA/CSV from the explicit package definitions.

`step4_make_yamls.py`, `step4_validate.sh` and `step5_collect_results.py` provide an optional structure-prediction route. No completed `final_ranked.csv` or new wet-lab validation dataset is included here. A predicted confidence score alone would not establish improved editing.

| File | Contents |
|---|---|
| `inputs/complex.pdb` | Predicted PUF12–17-nt RNA structural input |
| `inputs/design_config.json`, `inputs/repeat_parse.csv` | Recognition positions and sequence numbering |
| `results/mpnn_samples_proteinmpnn.fa.gz` | 10,000 ProteinMPNN samples plus native reference |
| `results/mpnn_samples_ligandmpnn.fa.gz` | 10,000 LigandMPNN samples plus native reference |
| `results/aice_single_ranked.csv` | All 493 positions, model-specific alternatives, frequencies and flags |
| `results/recommended_mutations.csv` | Archived 81-position nomination summary; see export detail above |
| `results/aice_multi_combos.csv` | Exploratory covariance/association-based combinations |
| `results/wetlab_vs_prediction.csv` | Retrospective comparison with 84 experimental variant records |
| `results/repeat_importance.csv` | Descriptive experimental summaries by mutated repeat |
| `results/gen3_constructs.fasta`, `results/gen3_constructs.csv` | G0–G10 sequences and exact mutation lists |
| `code/` | Preparation, sampling, screening and construct-generation scripts |

## Reference

Fei, H., Li, Y., Liu, Y., Wei, J., Chen, A., & Gao, C. (2025). Advancing protein evolution with inverse folding models integrating structural and evolutionary constraints. *Cell, 188*(17), 4674–4692.e19. https://doi.org/10.1016/j.cell.2025.06.014
