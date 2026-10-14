#!/usr/bin/env Rscript
# =====================================================================================================================
# Visium HD, step 3 - fine cell states (98 Pelka categories) WITHIN the 17-category family assigned at step 2
#
# For family F (index = SLURM_ARRAY_TASK_ID, 1..17):
#  - reference = fine states of F + the 16 other families as broad categories (max 1,000 cells per category,
#    seed 42), so that signal from neighbouring cells is absorbed by another family rather than forced into a state of F;
#  - cells = non-rejected cells whose step-2 first_type is F; RCTD doublet mode, UMI_min = 50.
# Final state (subtype98) = first fine state of F among first_type/second_type, otherwise "<F> (unresolved)".
# 98 -> 17 mapping: majority vote over reference cells between ClusterFull and the 17-category label.
# Inputs : objects of steps 1-2 ; <HD_DIR>/reference/{pelka_reference_98cats.rds, cell_type_17groups.tsv.gz}
# Output : <HD_DIR>/per_sample/<SAMPLE>/objects/rctd98_within_family/subtypes_<family>.rds
# Usage  : SAMPLE=<id> SLURM_ARRAY_TASK_ID=<1..17> Rscript analyses/14_visium_hd_S12/03_rctd98_within_family.R
# =====================================================================================================================
suppressPackageStartupMessages({
  library(spacexr)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

sample_id <- Sys.getenv("SAMPLE"); stopifnot(nzchar(sample_id))
hd_dir <- Sys.getenv("CRCVISIUM_HD", "~/data/visium_hd")
sroot  <- file.path(hd_dir, "per_sample", sample_id)
ref98_path <- file.path(hd_dir, "reference", "pelka_reference_98cats.rds")
lab17_path <- file.path(hd_dir, "reference", "cell_type_17groups.tsv.gz")
lvl1_path <- file.path(sroot, "objects", paste0("rctd17_doublet_all_", sample_id, ".rds"))
outdir <- file.path(sroot, "objects", "rctd98_within_family")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
max_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# ---- 98 -> 17 mapping ----
ref98 <- readRDS(ref98_path)
lab17 <- read.delim(lab17_path, stringsAsFactors = FALSE)
f17 <- lab17$cell_type[match(colnames(ref98), lab17$cellID)]
ct98 <- as.character(ref98$cell_type)
tab <- table(ct98[!is.na(f17)], f17[!is.na(f17)])
map <- data.frame(subtype = rownames(tab),
                  family = colnames(tab)[apply(tab, 1, which.max)],
                  purity = round(apply(tab, 1, max) / rowSums(tab), 3))
write.csv(map, file.path(outdir, "mapping_98_to_17.csv"), row.names = FALSE)
if (any(map$purity < 0.9)) {
  message("WARNING: fine states with purity < 0.9:")
  print(map[map$purity < 0.9, ])
}

families <- sort(unique(map$family))
stopifnot(task >= 1, task <= length(families))
fam <- families[task]
fine <- map$subtype[map$family == fam]
message("family ", task, "/", length(families), ": ", fam, " (", length(fine), " fine states)")

# ---- HD cells of the family ----
lvl1 <- readRDS(lvl1_path)
cd1 <- as.data.frame(colData(lvl1))
sel <- rownames(cd1)[cd1$spot_class != "reject" & cd1$first_type == fam]
safe <- gsub("[^A-Za-z0-9]+", "_", fam)
if (length(fine) < 2 || length(sel) < 50) {
  message("family skipped (", length(fine), " fine state(s), ", length(sel), " cells)")
  res <- data.frame(barcode = sel, family = rep(fam, length(sel)),
                    subtype98 = if (length(fine) == 1) rep(fine, length(sel)) else rep(paste0(fam, " (unresolved)"), length(sel)),
                    spot_class2 = NA, first_type2 = NA, second_type2 = NA)
  saveRDS(res, file.path(outdir, paste0("subtypes_", safe, ".rds")))
  quit(save = "no")
}

# ---- hybrid reference: fine states for F, broad categories for the other families ----
hyb <- ifelse(ct98 %in% fine, ct98, f17)
keep_ref <- !is.na(hyb)
set.seed(42)
idx <- unlist(lapply(split(which(keep_ref), hyb[keep_ref]), function(i) if (length(i) > 1000) sample(i, 1000) else i))
ref_h <- SummarizedExperiment(assays = list(counts = assay(ref98, "counts")[, idx]),
                              colData = DataFrame(cell_type = factor(hyb[idx])))
rm(ref98); gc()

cells <- readRDS(file.path(sroot, "objects", paste0("cells_meta_", sample_id, ".rds")))
cnt <- readRDS(file.path(sroot, "objects", paste0("counts_", sample_id, ".rds")))
cnt <- cnt[, sel]
# genes with zero counts in these cells removed: no effect on UMI per cell, but avoids the failure of the
# "90% of genes do not match" check of spacexr::getInitialGeneList for small families (CD4 T, CD8 T...)
cnt <- cnt[Matrix::rowSums(cnt) > 0, ]
cells <- cells[match(sel, cells$barcode), ]
coords <- as.matrix(cells[, c("x_um", "y_um")]); colnames(coords) <- c("x", "y"); rownames(coords) <- sel
spe <- SpatialExperiment(assay = list(counts = cnt), spatialCoords = coords)

rctd_data <- createRctd(spe, ref_h, cell_type_col = "cell_type", UMI_min = 50)
result <- runRctd(rctd_data, rctd_mode = "doublet", max_cores = max_cores)
saveRDS(result, file.path(outdir, paste0("rctd98_", safe, ".rds")))

cd <- as.data.frame(colData(result))
ft <- as.character(cd$first_type); st <- as.character(cd$second_type)
sub <- ifelse(ft %in% fine, ft,
              ifelse(cd$spot_class %in% c("doublet_certain", "doublet_uncertain") & st %in% fine, st,
                     paste0(fam, " (unresolved)")))
res <- data.frame(barcode = rownames(cd), family = fam, subtype98 = sub,
                  spot_class2 = cd$spot_class, first_type2 = ft, second_type2 = st)
missing <- setdiff(sel, res$barcode)
if (length(missing)) {
  res <- rbind(res, data.frame(barcode = missing, family = fam, subtype98 = paste0(fam, " (unresolved)"),
                               spot_class2 = "below_UMI_min", first_type2 = NA, second_type2 = NA))
}
saveRDS(res, file.path(outdir, paste0("subtypes_", safe, ".rds")))
print(sort(table(res$subtype98), decreasing = TRUE))

