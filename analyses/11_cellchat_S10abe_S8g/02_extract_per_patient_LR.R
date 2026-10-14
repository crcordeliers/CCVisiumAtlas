#!/usr/bin/env Rscript
# =====================================================================================================================
# Extracts from the 48 CellChat objects every non-zero communication probability (source, target, ligand-receptor
# pair, pathway, probability, permutation P) and the number of spots per group. No recomputation of CellChat.
# Input : <CELLCHAT_DIR>/objects/CellChat_spatial_regen_<Sample>.rds (01_cellchat_spatial_per_sample.R)
# Output: <CELLCHAT_DIR>/tables/LR_prob_per_sample_nonzero.csv.gz, n_spots_per_group_per_sample.csv
# Usage : Rscript analyses/11_cellchat_S10abe_S8g/02_extract_per_patient_LR.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(CellChat))

cellchat_dir <- Sys.getenv("CRCVISIUM_CELLCHAT", file.path(work.dir, "CellChat"))
in_dir  <- file.path(cellchat_dir, "objects")
out_dir <- file.path(cellchat_dir, "tables"); dir.create(out_dir, showWarnings = FALSE)
files <- list.files(in_dir, pattern = "^CellChat_spatial_regen_.*\\.rds$", full.names = TRUE)
stopifnot(length(files) == 48)

lr_all <- list(); n_all <- list()
for (f in files) {
  s  <- sub("^CellChat_spatial_regen_(.*)\\.rds$", "\\1", basename(f))
  cc <- readRDS(f)
  prob <- cc@net$prob; pval <- cc@net$pval
  idx  <- which(prob > 0, arr.ind = TRUE)
  if (nrow(idx) > 0)
    lr_all[[s]] <- data.frame(sample = s, source = dimnames(prob)[[1]][idx[, 1]], target = dimnames(prob)[[2]][idx[, 2]],
                              interaction_name = dimnames(prob)[[3]][idx[, 3]], prob = prob[idx], pval = pval[idx])
  tab <- table(cc@idents)
  n_all[[s]] <- data.frame(sample = s, group = names(tab), n_spots = as.integer(tab))
}
lr_df <- do.call(rbind, lr_all)
ann <- unique(CellChatDB.human$interaction[, c("interaction_name", "pathway_name", "ligand", "receptor", "annotation")])
lr_df <- merge(lr_df, ann, by = "interaction_name", all.x = TRUE)
write.csv(do.call(rbind, n_all), file.path(out_dir, "n_spots_per_group_per_sample.csv"), row.names = FALSE)
data.table::fwrite(lr_df, file.path(out_dir, "LR_prob_per_sample_nonzero.csv.gz"))
