#!/usr/bin/env Rscript
# =====================================================================================================================
# Pools the 48 per-sample YAP-Rev tables (01_yaprev_score_one_sample.R) and exports:
#   <YAPREV_DIR>/yaprev_all48_pooled.rds/.csv    all spots (YAPRev_raw is never modified; a pooled z-score is added for
#                                               display only)
#   <YAPREV_DIR>/merged_spots/<sample>_merged.csv barcode, x, y (frame of the front annotations), YAPRev_raw - input of S8h-j, S10
#   <YAPREV_DIR>/geojson/<sample>.geojson         spots as circles with YAP-Rev measurements, for import into QuPath
# Usage : Rscript analyses/10_YAPRev_spots_S8g-j/02_yaprev_pool_export.R
# =====================================================================================================================
suppressPackageStartupMessages(library(jsonlite))

yaprev_dir <- Sys.getenv("CRCVISIUM_YAPREV", "~/data/yaprev_spots")
files <- list.files(file.path(yaprev_dir, "pass1"), pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) == 48)
all_res <- do.call(rbind, lapply(files, readRDS))
stopifnot(!anyNA(all_res$YAPRev_raw))
all_res$YAPRev_zscore <- (all_res$YAPRev_raw - mean(all_res$YAPRev_raw)) / sd(all_res$YAPRev_raw)
message(nrow(all_res), " spots; YAP-Rev high (>= 0.40): ", round(100 * mean(all_res$YAPRev_raw >= 0.40), 2), "%")
saveRDS(all_res, file.path(yaprev_dir, "yaprev_all48_pooled.rds"))
write.csv(all_res, file.path(yaprev_dir, "yaprev_all48_pooled.csv"), row.names = FALSE)

dir.create(file.path(yaprev_dir, "merged_spots"), showWarnings = FALSE)
dir.create(file.path(yaprev_dir, "geojson"), showWarnings = FALSE)

circle_ring <- function(x, y, r, n_points = 32) {
  a <- seq(0, 2 * pi, length.out = n_points + 1)
  sprintf("[%.3f,%.3f]", x + r * cos(a), y + r * sin(a))
}
for (sid in unique(all_res$sample_id)) {
  df <- all_res[all_res$sample_id == sid, ]
  # frame of the invasive-front annotations (fi_lines): Space Ranger full-resolution pixels, x = column + spot diameter / 66,
  # y = row
  write.csv(data.frame(barcode = df$barcode, x = df$pxl_col_in_fullres + 2 * df$radius / 66, y = df$pxl_row_in_fullres,
                       YAPRev_raw = df$YAPRev_raw),
            file.path(yaprev_dir, "merged_spots", paste0(sid, "_merged.csv")), row.names = FALSE)
  feats <- vapply(seq_len(nrow(df)), function(i) sprintf(
    '{"type":"Feature","id":"%s","geometry":{"type":"Polygon","coordinates":[[%s]]},"properties":{"objectType":"detection","measurements":{"YAPRev_raw":%s,"YAPRev_zscore":%s,"ANXA1":%s,"AHNAK":%s,"CAPN2":%s}}}',
    df$barcode[i], paste(circle_ring(df$x_rot[i], df$y_rot[i], df$radius[i]), collapse = ","),
    format(df$YAPRev_raw[i], digits = 8), format(df$YAPRev_zscore[i], digits = 8),
    format(df$ANXA1[i], digits = 8), format(df$AHNAK[i], digits = 8), format(df$CAPN2[i], digits = 8)), character(1))
  writeLines(paste0('{"type":"FeatureCollection","features":[', paste(feats, collapse = ","), ']}'),
             file.path(yaprev_dir, "geojson", paste0(sid, ".geojson")))
}
