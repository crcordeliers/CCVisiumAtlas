# CCVisiumAtlas

Code for the spatial transcriptomics analyses of

> Cazelles A., Pilati C. *et al.* Profiling colon cancer architecture with spatial transcriptomics identifies
> clinically relevant stromal ecotypes.

48 stage III colon tumours from the PETACC-8 trial profiled with 10x Genomics Visium (one tumour per capture area),
4 of them re-profiled with Visium HD; prognostic validation of the spatially derived signatures in three bulk RNA-seq
cohorts (PETACC-8, IDEA-France, PRODIGE-13).

## Repository layout

| Folder | Content |
|---|---|
| `R/` | core Visium pipeline (`Step1`-`Step7`: Space Ranger outputs -> meta-spots -> meta-clusters -> stromal sub-clusters) and the scripts of Fig. 1b, 2a, 3a, 3c, 4a and Supplementary Fig. S4e, S8a, S8c, S8f; names and colours in `R/00_config.R` |
| `analyses/` | one folder per analysis (supplementary figures, Fig. 2c, 3b and 6); each script starts with a header giving the panels, inputs, outputs and usage |
| `data/` | gene signatures, Supplementary Table S2 (`reference_TableS2_signatures.csv`), meta-spot reference labels, minimal annotation of the 48 samples (sex, MMR status, relapse), Pelka et al. 2021 cell-type hierarchy, sample list |

## Data

| Data | Access |
|---|---|
| Visium raw data | EGA, dataset EGAD50000002091 |
| H&E images | Zenodo, doi:10.5281/zenodo.17711406 |
| Single-cell reference (Pelka et al. 2021) | GEO GSE178341; cell annotations from the Broad Single Cell Portal (SCP1162) |
| PETACC-8 and IDEA-France bulk RNA-seq | as described in the Data availability section of the paper |
| PRODIGE-13 bulk RNA-seq | FFCD data-sharing policy |

## Requirements

- Core pipeline (`R/`): R 4.3.1 with the package versions of `renv.lock` (`renv::restore()` from the repository
  root): Seurat 5.4.0, harmony 1.2.4 (harmony 1.x is required), GSVA 1.50.5.
- Other R analyses: R >= 4.3 with the packages loaded in each script (e.g. CellChat 2.2, spacexr / RCTD, BayesSpace,
  BASS, aricode, spdep, survival, survminer, enrichR); the R version is given in the script header when it differs.
- Python 3: cell2location 0.1.5 and scvi-tools 1.3.3 (`analyses/04_*`); scvi-tools 1.4.2, squidpy 1.8.2,
  cellcharter 0.3.7, GraphST 1.1.1 and harmonypy 2.0.0 (`analyses/03_*`); pandas, numpy, scipy and matplotlib
  (`analyses/10_*`, `analyses/13_*`). cell2location, scVI / CellCharter and GraphST are run on a GPU.

## Configuration

All scripts are run from the repository root. Paths are set with environment variables:

| Variable | Content |
|---|---|
| `CRCVISIUM_WORKDIR` | working / output directory (`Figures/`, `Tables/`, intermediate objects) |
| `CRCVISIUM_RAW` | Visium Space Ranger outputs (`<Sample>/outs`) |
| `CRCVISIUM_CNV` | per-spot fastCNV fractions |
| `CRCVISIUM_PELKA` | Pelka et al. 2021 single-cell data (GSE178341 `.h5`, SCP1162 tSNE / cluster table) |
| `CRCVISIUM_C2L` | cell2location working directory |
| `CRCVISIUM_YAPREV` | spot-level YAP-Rev tables (`merged_spots/`), invasive-front annotations (`fi_lines/`, `pixel_sizes.csv`) |
| `CRCVISIUM_FRONT_OUT` | distances of the spots to the invasive front (default `<CRCVISIUM_YAPREV>/spot_distances`) |
| `CRCVISIUM_CELLCHAT` | CellChat objects and tables |
| `CRCVISIUM_GRADIENT` | per-spot gene counts for the stromal gradients |
| `CRCVISIUM_HD` | Visium HD data (`data/<sample>/outs`, Space Ranger 4.0.1 with cell segmentation), `per_sample/`, `reference/` |
| `CRCVISIUM_BULK` | RDS of the three bulk cohorts (`analyses/15_bulk_survival_Fig6_S13/00_build_bulk_dataset.R`) |
| `CRCVISIUM_SAMPLES` | sample list (default `data/samples_48.txt`) |
| `CRCVISIUM_WORKERS` | number of parallel workers for Step 2 (default 15) |

## Core pipeline

```
Rscript R/Step1_Object_creation.R        # per-sample Seurat objects, CNV fractions
Rscript R/Step2_Metaspots.R              # meta-spots (per-sample clustering, resolution 10)
Rscript R/Step3_Meta_clusters.R          # integration (Harmony), meta-clusters, Supplementary Table S2
Rscript R/Step4_Stromal_Subclusters.R    # stromal sub-clusters
Rscript R/Step5_Add_meta.R               # labels transferred to every spot
Rscript R/Step6_Pelka_Signatures.R       # Pelka et al. signatures
Rscript R/Step7_Check_results.R          # comparison with the reference results in data/
```

## Figures and scripts

| Figure | Script |
|---|---|
| Fig. 1b | `R/Fig1b_metaspots.R` |
| Fig. 2a | `R/Fig2a_S4e_ecotypes.R` |
| Fig. 2c, 3b | `analyses/04_cell2location_Fig2c_Fig3b/01-06` |
| Fig. 3a, 3c | `R/Fig3ac_stromal_subclusters.R` |
| Fig. 4a | `R/Fig4a_S8acf_relapse_stroma.R` |
| Fig. 6a-g | `analyses/15_bulk_survival_Fig6_S13/00_build_bulk_dataset.R`, `01_Fig6_survival.R` |
| Supplementary Fig. S1 | `analyses/01_metaspots_S1/FigS1_metaspots_and_resolution.R` |
| Supplementary Fig. S2a,b | `analyses/02_spatial_organisation_S2ab/FigS2ab_projection_and_MoranI.R` |
| Supplementary Fig. S2c,d | `analyses/03_clustering_benchmark_S2cd/benchmark_clustering.R`, `benchmark_clustering.py` |
| Supplementary Fig. S4a-d | `analyses/05_ecotype_hallmarks_S4/FigS4abcd_violins.R` |
| Supplementary Fig. S4e | `R/Fig2a_S4e_ecotypes.R` |
| Supplementary Fig. S5a,b, S6, S7 | `analyses/06_ecotype_composition_S5ab_S6_S7/FigS5ab_S6_S7_composition.R` |
| Supplementary Fig. S5c | `analyses/07_CMS_pseudobulk_S5c/01-02` |
| Supplementary Fig. S8a, S8c, S8f | `R/Fig4a_S8acf_relapse_stroma.R` |
| Supplementary Fig. S8b, S9, Supplementary Table S4 | `analyses/08_relapse_genes_S8b_S9_TableS4/FigS8b_S9_TableS4.R` |
| Supplementary Fig. S8d | `analyses/09_TF_binding_S8de/FigS8d_EnrichR_TF.R` |
| Supplementary Fig. S8g, S10a, S10b, S10e | `analyses/11_cellchat_S10abe_S8g/01-03` |
| Supplementary Fig. S8h-j | `analyses/10_YAPRev_spots_S8g-j/01-06` |
| Supplementary Fig. S10c,d | `analyses/12_ANXA1_correlations_S10cd/FigS10cd_ANXA1_correlations.R` |
| Supplementary Fig. S10f,g | `analyses/13_stromal_gradients_S10fg/01-04` |
| Supplementary Fig. S12 | `analyses/14_visium_hd_S12/00-11` |
| Supplementary Fig. S13 | `analyses/15_bulk_survival_Fig6_S13/02_FigS13_forest_by_arm_MMR.R` |
| Supplementary Table S2 | `R/Step3_Meta_clusters.R`, `R/Step4_Stromal_Subclusters.R` |

Scripts numbered `01`, `02`, ... within a folder are run in that order. The spot-level YAP-Rev score
(`analyses/10_YAPRev_spots_S8g-j/01-02`) is also used by the CellChat (`11_*`) and stromal gradient (`13_*`) analyses.
