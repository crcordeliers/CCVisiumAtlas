#!/usr/bin/env Rscript
# =====================================================================================================================
# Spatial CellChat v2, one run per Visium sample (Methods "Cell-cell communication analyses"; Supplementary Fig. S10a,b,e)
#
# Spot identities = meta-cluster for non-stromal spots and stromal sub-cluster for stromal spots; stromal spots with
# YAPRev_raw >= 0.40 are relabelled "YAP-Rev high" (internal label "Regenerative", renamed in the figures).
# computeCommunProb(raw.use = TRUE, distance.use = TRUE, interaction.range = 250 µm, scale.distance = 0.4,
#                   contact.dependent = TRUE, contact.range = 100 µm); filterCommunication(min.cells = 10).
# Distances are calibrated so that the minimum centre-to-centre spot distance = 100 µm (Visium spot pitch).
#
# Inputs : Processed_seurat/processed_seurat_<Sample>_withmeta.RData (R/Step5_Add_meta.R),
#          <YAPREV_DIR>/merged_spots/<Sample>_merged.csv (analyses/10_YAPRev_spots_S8g-j/02_yaprev_pool_export.R)
# Output : <CELLCHAT_DIR>/objects/CellChat_spatial_regen_<Sample>.rds
# Usage  : Rscript analyses/11_cellchat_S10abe_S8g/01_cellchat_spatial_per_sample.R [Sample]   (all 48 if no argument)
# Env    : R 4.3.1, Seurat 5.4.0, CellChat 2.2.0.9001
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(CellChat); library(Seurat) })
options(stringsAsFactors = FALSE)

yaprev_dir   <- Sys.getenv("CRCVISIUM_YAPREV", "~/data/yaprev_spots")
cellchat_dir <- Sys.getenv("CRCVISIUM_CELLCHAT", file.path(work.dir, "CellChat"))
out_dir <- file.path(cellchat_dir, "objects"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else Sample.names

# Group labels used in the CellChat objects (the figure scripts expect these names)
cc_label <- c("Stroma-A1" = "Stroma A 1", "Stroma-A2" = "Stroma A 2", "Stroma-B1" = "Stroma B 1",
              "Stroma-B2" = "Stroma B 2", "Stroma-B3" = "Stroma B 3", "Tumor-fibroblasts" = "Tumor-fibroblasts")
stroma_labels <- c("Stroma A 1", "Stroma A 2", "Stroma B 1", "Stroma B 2", "Stroma B 3")
YAP_THR <- 0.40
spot.radius <- 65 / 2          # µm
target_min_dist <- 100         # µm, Visium spot pitch

for (sample_name in samples) {
  output_file <- file.path(out_dir, paste0("CellChat_spatial_regen_", sample_name, ".rds"))
  if (file.exists(output_file)) next
  load(file.path(processed.dir, paste0("processed_seurat_", sample_name, "_withmeta.RData")))   # -> sobj
  seu <- sobj; rm(sobj)

  lab <- seu$meta_stromal_label
  lab <- ifelse(lab %in% names(cc_label), cc_label[lab], lab)
  yap <- read.csv(file.path(yaprev_dir, "merged_spots", paste0(sample_name, "_merged.csv")))
  seu$yap <- setNames(yap$YAPRev_raw, yap$barcode)[colnames(seu)]
  is_stroma <- lab %in% stroma_labels

  keep <- !is.na(lab) & (!is_stroma | !is.na(seu$yap))
  seu <- subset(seu, cells = colnames(seu)[keep]); lab <- lab[keep]; is_stroma <- is_stroma[keep]
  seu$CellChat_var <- factor(ifelse(is_stroma & seu$yap >= YAP_THR, "Regenerative", lab))
  Idents(seu) <- seu$CellChat_var
  message(sample_name, ": ", sum(seu$CellChat_var == "Regenerative"), " YAP-Rev high / ", sum(is_stroma), " stromal spots")

  data.input <- GetAssayData(seu, assay = "Spatial", layer = "data")
  coords <- GetTissueCoordinates(seu, filter = "tissue")[colnames(data.input), c(1, 2)]
  colnames(coords) <- c("x", "y")
  d.raw <- computeCellDistance(coordinates = coords, ratio = 1, tol = 0)
  spatial.factors <- data.frame(ratio = target_min_dist / min(d.raw[d.raw > 0], na.rm = TRUE), tol = spot.radius)

  cellchat <- createCellChat(object = data.input, meta = seu@meta.data, group.by = "CellChat_var",
                             datatype = "spatial", coordinates = coords, spatial.factors = spatial.factors)
  cellchat <- setIdent(cellchat, ident.use = "CellChat_var")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = TRUE)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, distance.use = TRUE, interaction.range = 250,
                                scale.distance = 0.4, contact.dependent = TRUE, contact.range = 100)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, output_file)
  rm(seu, cellchat); gc()
}
