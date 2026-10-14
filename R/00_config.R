# Shared configuration: paths, sample list, cluster names and colours.
# Sourced at the top of every Step/Figure script.

## ---- Paths (edit to your environment) ----
# Scripts are run from the repository root (Rscript R/StepX.R). Paths can be set with environment variables.
work.dir  <- Sys.getenv("CRCVISIUM_WORKDIR", "~/data")                                 # outputs are written here
raw.dir   <- Sys.getenv("CRCVISIUM_RAW",   file.path(work.dir, "processed"))            # Space Ranger outputs: <raw.dir>/<Sample>/outs
cnv.file  <- Sys.getenv("CRCVISIUM_CNV",   file.path(work.dir, "CNV_fractions_71samples.RDS"))  # fastCNV per-spot CNV fractions (named list)
data.dir  <- Sys.getenv("CRCVISIUM_DATA",  "data")                                      # the repository data/ folder

allseurat.dir  <- file.path(work.dir, "All_Seurat")         # Step1 outputs
processed.dir  <- file.path(work.dir, "Processed_seurat")   # Step2 outputs (per-sample objects with Cluster_res_10)
results_res.dir <- file.path(work.dir, "Results_res")       # Step2 outputs (meta-spot tables)
harmony.dir    <- file.path(work.dir, "Results_Harmony", "res_10")  # Step3-5 outputs
tables.dir     <- file.path(work.dir, "Tables")             # exported tables (Supplementary Table S2, volcano)
figures.dir    <- file.path(work.dir, "Figures")
for (d in c(allseurat.dir, processed.dir, results_res.dir, harmony.dir, tables.dir, figures.dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

n_workers <- as.integer(Sys.getenv("CRCVISIUM_WORKERS", "15"))   # parallel workers (Step 2)

## ---- Samples (48 tumours, one patient per capture area) ----
Sample.names <- c("Ech_052_A1","Ech_052_B1","Ech_052_C1","Ech_052_D1",
                  "Ech_073_A1","Ech_073_B1","Ech_073_C1","Ech_073_D1",
                  "Ech_074_A1","Ech_074_B1","Ech_074_C1","Ech_074_D1",
                  "Ech_075_A1","Ech_075_B1","Ech_075_C1","Ech_075_D1",
                  "Ech_077_A1","Ech_077_B1","Ech_077_C1","Ech_077_D1",
                  "Ech_133_A1","Ech_133_B1","Ech_133_C1","Ech_133_D1",
                  "Ech_134_A1","Ech_134_B1","Ech_134_C1","Ech_134_D1",
                  "Ech_21-77_A1","Ech_21-77_B1","Ech_21-77_C1","Ech_21-77_D1",
                  "Ech_295_A1","Ech_295_B1","Ech_295_C1","Ech_295_D1",
                  "Ech_296_A1","Ech_296_B1","Ech_296_C1","Ech_296_D1",
                  "Ech_395_A1","Ech_395_B1","Ech_395_C1","Ech_395_D1",
                  "Ech_405_A1","Ech_405_B1","Ech_405_C1","Ech_405_D1")

## ---- Clustering parameters ----
res_metaspot <- 10    # per-sample clustering resolution defining meta-spots
res_meta     <- 0.5   # meta-cluster resolution (Harmony space, dims 1:15)
res_stromal  <- 0.5   # stromal sub-cluster resolution (meta-clusters 2 and 4)
stromal_metaclusters <- c("2", "4")

## ---- Names (codes = Seurat cluster numbers) ----
meta_names <- c(
  "0"  = "Tumor-cycling",
  "1"  = "Tumor-intermediate",
  "7"  = "Tumor-infiltrated",
  "8"  = "Tumor-fibroblasts",
  "2"  = "Stroma-A",
  "4"  = "Stroma-B",
  "9"  = "Necrosis",
  "3"  = "Muscularis propria",
  "10" = "Lymphoid aggregates",
  "5"  = "Nerves - arteries",
  "6"  = "Normal colon epithelium"
)
meta_order <- names(meta_names)   # display order: 0, 1, 7, 8, 2, 4, 9, 3, 10, 5, 6

stromal_names <- c(
  "0" = "Stroma-A1",
  "2" = "Stroma-A2",
  "1" = "Stroma-B1",
  "3" = "Stroma-B2",
  "4" = "Stroma-B3"
)
stromal_order <- names(stromal_names)  # display order: 0, 2, 1, 3, 4

## ---- Colours ----
meta_colors <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3",
  "4" = "#FF7F00", "5" = "#FFFF33", "6" = "#A65628", "7" = "#F781BF",
  "8" = "#999999", "9" = "black",   "10" = "blue"
)
meta_colors_named <- setNames(meta_colors[names(meta_names)], meta_names)

stromal_colors <- c(
  "0" = "#A7C7E7",
  "2" = "#B0EACD",
  "1" = "#F7C6A3",
  "3" = "#F4A9A8",
  "4" = "#F9E79F"
)
stromal_colors_named <- setNames(stromal_colors[names(stromal_names)], stromal_names)

relapse_colors <- c("0" = "#9197FF", "1" = "#FC7171")   # No Relapse / Relapse

## ---- Helpers ----
# Meta-spot identifier used everywhere: <Sample>_res_10_Clust_<cluster>
metaspot_id <- function(sample, cluster) paste0(sample, "_res_", res_metaspot, "_Clust_", cluster)
