# CMPIS: programmed cell death states in cutaneous melanoma

Code and results for the study:

> **Programmed cell death states are associated with immune heterogeneity and therapeutic vulnerability in cutaneous melanoma**
> Angelos Kyriakoulis, Ilias Georgakopoulos-Soares, Apostolos Zaravinos. *Submitted.*

The pipeline integrates TCGA-SKCM multi-omics data with public GEO, single-cell and
immunotherapy cohorts. It has three parts:

1. **Subtypes.** It identifies two programmed-cell-death (PCD)-associated melanoma subtypes: **CS1** (immune-hot) and **CS2** (immune-cold, melanocytic).
2. **Characterisation.** It characterises their regulatory, immune and pathway landscapes.
3. **Signature.** It builds the 8-gene **Cutaneous Melanoma Prognostic Immune Signature (CMPIS)** and tests it in independent survival cohorts, in single-cell data, in five immune-checkpoint-inhibitor cohorts, and against in silico drug-sensitivity predictions.

CMPIS genes: `CFLAR, CYBB, CYLD, GYS1, IRF1, NFKBIA, PIM2, VCAM1`. Weights are in
`skcm_analysis_results/final_model_info_corrected.rds` and
`results/Supplementary_Table_6_CMPIS_weights.csv`; the scoring function is `R/cmpis_weights.R`.
The CMPIS was called MPCDS in earlier versions of the code.

---

## Pipeline stages

| Stage | Script | Purpose |
|---|---|---|
| 1a | `TCGA_Preprocessing.R` | Download and preprocess TCGA-SKCM multi-omics (TCGAbiolinks) |
| 1b | `GEO.R` | Download GEO cohorts, probe-to-gene mapping, ComBat, scaling |
| 2 | `TCGA_Clustering.R` | MOVICS consensus clustering (k = 2), NTP/PAM validation, DEGs, regulons, immune deconvolution, pathways (Figures 1–6) |
| 2+ | `UMAP&Clinical.R`, `Mutation&DESeq2.R` | UMAP, clinical heatmap, mutation and DESeq2 comparisons |
| 4 | `SKCM_ML_Algorithms.R` | CMPIS construction: PCD ∩ DEG genes, univariate Cox, 286-configuration model tournament |
| 5 | `Validation_PCDRGS.R` | Original interactive validation script (superseded by Stages 10–12) |
| 6 | `R/06_PCD_Sensitivity.R` | Re-clustering with each PCD pathway; robustness of CS1/CS2 |
| 7 | `R/07_Classifier_Benchmark.R` | CS1/CS2 and CMPIS vs IFN-γ, GEP, ESTIMATE, Tsoi, MITF signatures |
| 8 | `R/08_ICB_Validation.R` | Early anti-PD-1 check on two GEO cohorts (superseded by Stage 15) |
| 9 | `R/09_Mechanistic_Network.R` | Integrated PCD → regulator → immune → drug model (Figure 8) |
| 10 | `R/10_CMPIS_Reweight.R` | Final CMPIS weights and re-validation in TCGA and four GEO cohorts |
| 11 | `R/11_CMPIS_Clinical.R` | Kaplan–Meier, time-ROC, Cox, nomogram, immune correlates, drug sensitivity (GDSC2) |
| 12 | `R/12_Drug_oncoPredict_Check.R` | Drug predictions repeated with the oncoPredict/pRRophetic ridge workflow |
| 13 | `R/13_SingleCell_Validation.R` | Single-cell analysis: Tirosh, Jerby-Arnon, Sade-Feldman (Figure 7) |
| 14 | `R/14_Manuscript_Figures.R` | CMPIS weight and drug-effect figures from saved tables |
| 15 | `R/15_ICB_cBioPortal.R` | CS1/CS2 and CMPIS in five checkpoint-inhibitor cohorts (cBioPortal / CRI iAtlas) |

## Running

Requirements: R 4.4 (developed on R 4.4.3) with [`renv`](https://rstudio.github.io/renv/).

```r
renv::restore()        # install pinned packages (first time)
source("validate.R")   # pre-flight: wiring and which caches exist (no compute)
source("run.R")        # Stages 1–5 through {targets} (downloads ~6 GB on first run)
```

Stages 6–15 are standalone scripts that read the Stage 1–5 caches:

```r
source("config.R"); source("R/15_ICB_cBioPortal.R")   # one stage
```

On Windows, `run_all_stages.bat` runs Stages 6–15 in dependency order and writes logs to
`logs/`. The other `run_*.bat` files run single stages. Downloads made by Stages 8, 13 and 15
(GEO, cBioPortal) are cached under `data/`.

Stage 1–5 caching is controlled by `cfg$reuse_cache` in `config.R`. `R/helpers.R` provides
`rerun("<stage>")`, `status()` and `graph()` for the `{targets}` part.

## Where each paper item comes from

| Paper | File(s) in `results/` | Stage |
|---|---|---|
| Figures 1–6 | `figures/` | 2 |
| Figure 7, S3, S12–S14 | `Figure_7_singlecell.pdf` (pages 1+4 = Fig 7; 7 = S3; 2+5 = S12; 6 = S13; 3 = S14) | 13 |
| Figure 8 | `Figure_Mechanistic_Network.pdf` | 9 |
| Figure S1, Tables S1–S2 | `Figure_S_PCD_sensitivity.pdf`, `Supplementary_Table_1_*`, `_2_*` | 6 |
| Figure S2, Table S3 | `Figure_S_benchmark.pdf`, `Supplementary_Table_3_*` | 7 |
| Figure S5, Tables S4, S7 | `Figure_MS_CMPIS_weights.pdf`, `Supplementary_Table_5_*`, `_6_*`, `_7_*` | 9, 10, 14 |
| Figures S6–S8, Tables S8–S10 | `Figure_S3_CMPIS_KM.pdf`, `Figure_S4_CMPIS_clinical.pdf`, `Figure_S5_CMPIS_immune.pdf`, `Supplementary_Table_8_*` to `_10_*` | 11 |
| Figure S9, Table S11A–C | `Figure_S9_ICB_response.pdf`, `Figure_S9b_ICB_survival.pdf`, `Supplementary_Table_11_ICB_*`, `_11b_*`, `_11c_*` | 15 |
| Figure S10, Table S12 | `Figure_S6_CMPIS_drugs.pdf`, `Figure_MS_drug_effects.pdf`, `Supplementary_Table_11_CMPIS_drugs.csv` | 11, 14 |
| Figure S11, Table S13 | `Figure_MS_drug_implementations.pdf`, `Supplementary_Table_12_drug_oncoPredict.csv` | 12, 14 |
| Table S5A–G, S6 | `Supplementary_Table_13*` | 13 |

File names follow the stage that produced them; the table above gives the paper numbering
at submission.

## Data

All data are public. Large files are not committed; see `data/README.md`.

| Source | Data | Access |
|---|---|---|
| TCGA-SKCM | mRNA, miRNA, lncRNA, methylation, mutations, clinical | `TCGAbiolinks` (automatic) |
| GEO | GSE15605, GSE22153, GSE46517, GSE54467, GSE19234, GSE98394, GSE65904, GSE53118 | `GEOquery` (automatic) |
| GEO single cell | GSE72056 (Tirosh), GSE115978 (Jerby-Arnon), GSE120575 (Sade-Feldman) | automatic, Stage 13 |
| cBioPortal / CRI iAtlas | Liu 2019, Riaz 2017, Hugo 2016, Gide 2019, Van Allen 2015 (checkpoint inhibitors) | public API, Stage 15 |
| DepMap 25Q3 / GDSC2 | cell-line expression and drug AUC | manual download, see `data/README.md` |

## Repository layout

```
├── TCGA_Preprocessing.R  GEO.R  TCGA_Clustering.R      Stages 1–2
├── UMAP&Clinical.R  Mutation&DESeq2.R
├── SKCM_ML_Algorithms.R  Validation_PCDRGS.R           Stages 4–5
├── R/                    Stages 6–15, cmpis_weights.R, helpers
├── _targets.R  config.R  run.R  validate.R              pipeline and configuration
├── renv.lock  renv/  .Rprofile                          R environment
├── PCD GENES.xlsx  pcdrg_genes.txt  *_cohort.txt  tcga_skcm_*.txt   curated inputs
├── skcm_analysis_results/final_model_info_corrected.rds             CMPIS weights
├── data/README.md        large / external data
└── results/              figures, tables and stage outputs
```

## Citation

If you use this code, please cite the paper above (citation details will be added on
publication) and this repository.

## License

Code: MIT License (see `LICENSE`). Data obtained from TCGA, GEO, cBioPortal, DepMap and GDSC
remain under the terms of their original providers.
