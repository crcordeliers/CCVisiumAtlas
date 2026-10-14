#!/usr/bin/env Rscript
# =====================================================================================================================
# Assembles the bulk RNA-seq cohorts used in Fig. 6 and Supplementary Fig. S13 into one RDS file
#   PETACC8 (n = 1,733, stage III), IDEA-France (n = 1,248, stage III), PRODIGE13 (n = 529, stage II)
# Expression (one column per gene) + clinical variables: time to recurrence, grade (recgrade), WHO performance status
# (recOMS), baseline high-risk features (recRisk: pT4 and/or pN2; pT4 only for PRODIGE13), MMR status, bowel
# complication at diagnosis (occlusion or perforation), treatment arm.
# The cohort data are not public (PETACC8 / IDEA-France: as described in the paper; PRODIGE13: FFCD data-sharing
# policy); set the paths below (environment variables) to the local copies.
# Output: <WORKDIR>/bulk_cohorts.rds (list PETACC, IDEA, PRODIGE)
# Usage : Rscript analyses/15_bulk_survival_Fig6_S13/00_build_bulk_dataset.R
# =====================================================================================================================
source("R/00_config.R")
suppressPackageStartupMessages(library(dplyr))
IDEA_PETACC_RDATA <- Sys.getenv("BULK_IDEA_PETACC", "IDEA_PETACC.RData")            # IDEA.tot, IDEA.annot, PETACC.tot
PRODIGE_RDATA     <- Sys.getenv("BULK_PRODIGE", "PRODIGE.RData")                        # PRODIGE.tot
PRODIGE_ANNOT     <- Sys.getenv("BULK_PRODIGE_ANNOT", "PRODIGE13_annotation.rdata")   # prodige13.annot
PRODIGE_COMPLIC   <- Sys.getenv("BULK_PRODIGE_COMPLIC", "PRODIGE13_complications.csv")

# ---- IDEA-France and PETACC8 ----
load(IDEA_PETACC_RDATA)

patients_sans_cst <- c(
  "FR020012","FR228006","FR009028","FR014044","FR021038",
  "FR033039","FR110001","FR140003","FR143008","FR214006",
  "FR227026","FR228007","FR228008","FR003005","FR192004"
)

# IDEA-France: patients without consent for translational research removed
IDEA.tot <- IDEA.tot[!IDEA.annot$IDSubject %in% patients_sans_cst, ]

# ---- PRODIGE13 ----
load(PRODIGE_RDATA)


load(PRODIGE_ANNOT)

prodige13_complications <- read.csv(PRODIGE_COMPLIC)

PRODIGE.annot <- merge(
  prodige13_complications,
  prodige13.annot,
  by.x = "d_id_patient",
  by.y = "D_ID_PATIENT"
)

PRODIGE.tot <- merge(
  PRODIGE.annot,
  PRODIGE.tot,
  by.y = "Row.names",
  by.x = "d_id_patient"
)

# treatment variables
PRODIGE.tot$name_trt <- ifelse(PRODIGE.tot$F2_CT_POSTOP == 2, "chimio", "no_chimio")
PETACC.tot$name_trt  <- ifelse(PETACC.tot$trtn_m == 0, "Folfox", "Folfox+Cetux")
IDEA.tot$name_trt    <- ifelse(IDEA.tot$trt == 0, "6_mois", "3_mois")

PRODIGE.tot$rectrt <- as.factor(ifelse(PRODIGE.tot$F2_CT_POSTOP == 2, 1, 0))
PETACC.tot$rectrt  <- as.factor(PETACC.tot$trtn_m)
IDEA.tot$rectrt    <- as.factor(IDEA.tot$trt)

PETACC.tot$etat_TTR <- PETACC.tot$recurr_cens
PETACC.tot$TTR_time <- PETACC.tot$recurr
# =========================================================
# Derived clinical variables -------
# =========================================================

PETACC.tot$recRisk <- as.factor(
  ifelse(PETACC.tot$recrecT == "pT4" | PETACC.tot$recN == "pN2", "HighRisk", "LowRisk")
)
IDEA.tot$recRisk <- as.factor(
  ifelse(IDEA.tot$recrecT == "pT4" | IDEA.tot$recN == "pN2", "HighRisk", "LowRisk")
)
PRODIGE.tot$recRisk <- as.factor(
  ifelse(PRODIGE.tot$recrecT == "pT4", "HighRisk", "LowRisk")
)

# unique column names (needed by dplyr::mutate)
names(PETACC.tot) <- make.unique(names(PETACC.tot))
names(IDEA.tot)   <- make.unique(names(IDEA.tot))
names(PRODIGE.tot) <- make.unique(names(PRODIGE.tot))

# bowel complication at diagnosis (occlusion or perforation)
IDEA.tot <- IDEA.tot |>
  dplyr::mutate(
    bowel_complication = dplyr::case_when(
      bowel == "bowel obstruction and or perforation" ~ 1L,
      bowel == "No bowel obstruction and no perforation" ~ 0L,
      .default = NA_integer_
    )
  )
IDEA.tot$bowel_complication <- factor(
  IDEA.tot$bowel_complication,
  levels = c(0, 1),
  labels = c("NoComplication", "Complication")
)
PRODIGE.tot <- PRODIGE.tot |>
  dplyr::mutate(
    bowel_complication = dplyr::case_when(
      complication %in% c("Occlusion", "Perforation", "Perforation+Occlusion") ~ 1L,
      complication == "Pas de complication" ~ 0L,
      .default = NA_integer_
    )
  )
PRODIGE.tot$bowel_complication <- factor(
  PRODIGE.tot$bowel_complication,
  levels = c(0, 1),
  labels = c("NoComplication", "Complication")
)
PETACC.tot <- PETACC.tot |>
  dplyr::mutate(
    bowel_complication = dplyr::case_when(
      bowel == "bowel obstruction and or perforation" ~ 1L,
      bowel == "No bowel obstruction and no perforation" ~ 0L,
      .default = NA_integer_
    )
  )
PETACC.tot$bowel_complication <- factor(
  PETACC.tot$bowel_complication,
  levels = c(0, 1),
  labels = c("NoComplication", "Complication")
)

IDEA.tot$recgrade <- factor(IDEA.tot$recgrade, levels = c("G1-2","G3-4"))
IDEA.tot$recOMS   <- factor(IDEA.tot$recOMS, levels = c("OMS0","OMS1-2"))
IDEA.tot$recRisk  <- factor(IDEA.tot$recRisk, levels = c("LowRisk","HighRisk"))
IDEA.tot$mmr_status <- factor(IDEA.tot$mmr_status, levels = c("pMMR","dMMR"))

PRODIGE.tot$recgrade <- factor(PRODIGE.tot$recgrade, levels = c("G1-2","G3-4"))
PRODIGE.tot$recOMS   <- factor(PRODIGE.tot$recOMS, levels = c("OMS0","OMS1-2"))
PRODIGE.tot$recRisk  <- factor(PRODIGE.tot$recRisk, levels = c("LowRisk","HighRisk"))
PRODIGE.tot$mmr_status <- factor(PRODIGE.tot$mmr_status, levels = c("pMMR","dMMR"))

PETACC.tot$recgrade <- factor(PETACC.tot$recgrade, levels = c("G1-2","G3-4"))
PETACC.tot$recOMS   <- factor(PETACC.tot$recOMS, levels = c("OMS0","OMS1-2"))
PETACC.tot$recRisk  <- factor(PETACC.tot$recRisk, levels = c("LowRisk","HighRisk"))
PETACC.tot$mmr_status <- factor(PETACC.tot$mmr_status, levels = c("pMMR","dMMR"))

surv_data <- list(PETACC = PETACC.tot, IDEA = IDEA.tot, PRODIGE = PRODIGE.tot)
saveRDS(surv_data, file.path(work.dir, "bulk_cohorts.rds"))
