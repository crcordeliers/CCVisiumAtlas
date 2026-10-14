# Step 7 - Compare the outputs of Steps 3-4 with the reference results provided in data/:
#   data/reference_metaspot_labels.csv     meta-cluster, stromal sub-cluster and UMAP coordinates of the 3,231 meta-spots
#   data/reference_TableS2_signatures.csv  top-20 marker genes of Supplementary Table S2

source("R/00_config.R")

library(Seurat)
library(dplyr)

ari <- function(a, b) {
  keep <- !is.na(a) & !is.na(b); tab <- table(a[keep], b[keep]); c2 <- function(x) sum(x * (x - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / c2(sum(tab))
  (s - e) / ((sa + sb) / 2 - e)
}

ref <- read.csv(file.path(data.dir, "reference_metaspot_labels.csv"), stringsAsFactors = FALSE, colClasses = c(metacluster = "character", stromal_subcluster = "character"))
obj <- readRDS(file.path(harmony.dir, "seurat_object_res_10_with5subsets.rds"))
m <- match(ref$metaspot, colnames(obj))
stopifnot("meta-spots differ from the reference" = !anyNA(m) && ncol(obj) == nrow(ref))

mc <- as.character(obj$merged_meta_clusters)[m]
sc <- as.character(obj$stromal_subcluster)[m]
U  <- Embeddings(obj, "umap")[m, ]

cat("Meta-clusters     : identical =", identical(mc, ref$metacluster), "| ARI =", round(ari(ref$metacluster, mc), 4), "\n")
cat("Stromal subclust. : identical =", identical(sc, ref$stromal_subcluster), "| ARI =", round(ari(ref$stromal_subcluster, sc), 4), "\n")
cat("UMAP              : max |diff| =", signif(max(abs(U - as.matrix(ref[, c("umap_1", "umap_2")]))), 3), "\n")
if (!identical(mc, ref$metacluster)) print(table(reference = ref$metacluster, new = mc))
if (!identical(sc, ref$stromal_subcluster)) print(table(reference = ref$stromal_subcluster, new = sc, useNA = "ifany"))

s2 <- read.csv(file.path(data.dir, "reference_TableS2_signatures.csv"), stringsAsFactors = FALSE)
new <- bind_rows(
  read.csv(file.path(tables.dir, "TableS2_metaclusters_top20.csv"))        %>% mutate(level = "metacluster"),
  read.csv(file.path(tables.dir, "TableS2_stromal_subclusters_top20.csv")) %>% mutate(level = "stromal_subcluster"))
chk <- s2 %>% group_by(level, signature) %>% summarise(
  identical = identical(gene, new$gene[new$level == first(level) & new$cluster_name == first(signature)]),
  n_common  = length(intersect(gene, new$gene[new$level == first(level) & new$cluster_name == first(signature)])),
  .groups = "drop")
print(as.data.frame(chk))
cat("Supplementary Table S2: ", sum(chk$identical), "/", nrow(chk), "signatures identical (genes and order)\n")
