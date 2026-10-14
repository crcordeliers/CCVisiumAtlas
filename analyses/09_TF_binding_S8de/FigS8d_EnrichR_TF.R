#!/usr/bin/env Rscript
# =====================================================================================================================
# Supplementary Fig. S8d - transcription factors binding AHNAK, ANXA1, CAPN2 and LRP10 (EnrichR, ENCODE TF ChIP-seq 2015;
# other libraries exported for completeness). Status: Significant (adjusted P < 0.05), Exploratory (P < 0.05 only).
# S8e (predicted TEAD4 binding sites, matrix MA0809.1, in the AHNAK / ANXA1 / CAPN2 promoters) was obtained with the
# JASPAR web server (https://jaspar.elixir.no/), no script.
# Output: Tables/FigS8d_EnrichR_<library>.csv (terms with P < 0.05)
# Usage : Rscript analyses/09_TF_binding_S8de/FigS8d_EnrichR_TF.R      (requires internet access; results depend on
#         the version of the EnrichR libraries at the time of the query)
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(enrichR); library(dplyr) })

dbs <- c("ChEA_2016", "TF_Perturbations_Followed_by_Expression", "ENCODE_TF_ChIP-seq_2015", "TRANSFAC_and_JASPAR_PWMs",
         "ARCHS4_TFs_Coexp", "ENCODE_Histone_Modifications_2015", "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019")
setEnrichrSite("Enrichr")
enriched <- enrichr(c("AHNAK", "ANXA1", "CAPN2", "LRP10"), dbs)

for (db in dbs) {
  tab <- enriched[[db]] |>
    mutate(n_genes = as.numeric(sub("/.*", "", Overlap)),
           status = case_when(Adjusted.P.value < 0.05 ~ "Significant", P.value < 0.05 ~ "Exploratory", TRUE ~ "Not significant")) |>
    arrange(P.value) |>
    select(Term, n_genes, Genes, P.value, Adjusted.P.value, status) |>
    filter(status != "Not significant")
  write.csv(tab, file.path(tables.dir, paste0("FigS8d_EnrichR_", db, ".csv")), row.names = FALSE)
}
print(head(read.csv(file.path(tables.dir, "FigS8d_EnrichR_ENCODE_TF_ChIP-seq_2015.csv")), 15))   # = table of S8d
