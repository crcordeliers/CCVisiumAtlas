#!/usr/bin/env Rscript
# =====================================================================================================================
# Spot-level YAP-Rev score (YAPRev_raw), one Visium sample per run (Methods "Spatial proximity analyses").
#
# Raw Space Ranger counts -> NormalizeData -> ScaleData (whole gene panel, ~18,000 genes) -> ssGSEA (GSVA, normalize =
# TRUE) of the 3-gene YAP-Rev set (AHNAK, ANXA1, CAPN2) on the scaled matrix. A spot is "YAP-Rev high" if
# YAPRev_raw >= 0.40 (fixed threshold, ~top 5% of the 193,418 pooled spots). Used for S8g-j, S10 (CellChat, gradients).
# Spot coordinates are also converted to the frame of the QuPath project (vertical mirror + 90° rotation of the
# full-resolution image) so that the spots overlay the H&E images and the invasive-front annotations.
#
# Input : <CRCVISIUM_RAW>/<SAMPLE_ID>/outs (filtered_feature_bc_matrix, spatial/)
# Output: <YAPREV_DIR>/pass1/<SAMPLE_ID>.rds
# Usage : SAMPLE_ID=Ech_074_D1 Rscript analyses/10_YAPRev_spots_S8g-j/01_yaprev_score_one_sample.R  (SLURM array, 48 tasks)
# Env   : R 4.5, Seurat 5, GSVA 2.x
# =====================================================================================================================
suppressPackageStartupMessages({ library(Seurat); library(GSVA); library(jsonlite); library(png) })

sample_id  <- Sys.getenv("SAMPLE_ID"); stopifnot(nzchar(sample_id))
raw_dir    <- Sys.getenv("CRCVISIUM_RAW", "~/data/processed")
yaprev_dir <- Sys.getenv("CRCVISIUM_YAPREV", "~/data/yaprev_spots")
outs   <- file.path(raw_dir, sample_id, "outs")
outdir <- file.path(yaprev_dir, "pass1"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
yaprev_genes <- c("AHNAK", "ANXA1", "CAPN2")

cnt <- Read10X(data.dir = file.path(outs, "filtered_feature_bc_matrix"))
stopifnot(all(yaprev_genes %in% rownames(cnt)))
tp <- read.csv(file.path(outs, "spatial", "tissue_positions_list.csv"), header = FALSE,
               col.names = c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"))
tp <- tp[tp$in_tissue == 1, ]; rownames(tp) <- tp$barcode
tp <- tp[colnames(cnt), ]

sobj <- CreateSeuratObject(counts = cnt, assay = "RNA")
sobj <- NormalizeData(sobj, verbose = FALSE)
sobj <- ScaleData(sobj, verbose = FALSE)
expr_norm <- FetchData(sobj, vars = yaprev_genes)

gs <- gsva(ssgseaParam(as.matrix(sobj@assays$RNA$scale.data), list(YAPRev = yaprev_genes), normalize = TRUE))
yaprev_raw <- setNames(as.numeric(gs["YAPRev", ]), colnames(gs))[colnames(sobj)]

# coordinates in the QuPath image frame
sf <- fromJSON(file.path(outs, "spatial", "scalefactors_json.json"))
hires_dims   <- dim(readPNG(file.path(outs, "spatial", "tissue_hires_image.png")))
image_width  <- hires_dims[2] / sf$tissue_hires_scalef
image_height <- hires_dims[1] / sf$tissue_hires_scalef
center_x <- image_width / 2; center_y <- image_height / 2

res <- data.frame(sample_id = sample_id, barcode = colnames(sobj),
                  pxl_row_in_fullres = tp$pxl_row_in_fullres, pxl_col_in_fullres = tp$pxl_col_in_fullres,
                  YAPRev_raw = yaprev_raw, expr_norm[colnames(sobj), yaprev_genes])
x_px <- res$pxl_col_in_fullres
y_px <- image_height - res$pxl_row_in_fullres            # vertical mirror (raw positions vs NDPI image)
res$x_rot  <- center_x + (y_px - center_y)
res$y_rot  <- center_y + (x_px - center_x)
res$radius <- sf$spot_diameter_fullres / 2
saveRDS(res, file.path(outdir, paste0(sample_id, ".rds")))
message(sample_id, ": ", nrow(res), " spots")
