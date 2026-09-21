############################################################
# SCRIPT 13
# SLINGSHOT PSEUDOTIME ANALYSIS
# FLG / DSG1 / DSG3 EXPRESSION ALONG PSEUDOTIME
# GSE72056 MELANOMA
############################################################

############################################################
# 1. LOAD PACKAGES
############################################################

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(patchwork)
library(dplyr)

############################################################
# 2. PROJECT DIRECTORY
############################################################

setwd("~/Desktop/Melanoma_Project2")

dir.create(
  "Results/SingleCell",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/SingleCell",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "RData",
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3. LOAD ANNOTATED SEURAT OBJECT
############################################################

melanoma <- readRDS(
  "RData/GSE72056_Annotated_Seurat.rds"
)

cat("\n====================================================\n")
cat("SCRIPT 15: SLINGSHOT PSEUDOTIME ANALYSIS\n")
cat("====================================================\n")

############################################################
# 4. CHECK CELL TYPES
############################################################

cat("\nCell types:\n")

print(
  table(melanoma$cell_type)
)

############################################################
# 5. CHECK UMAP
############################################################

if (!"umap" %in% names(melanoma@reductions)) {
  stop("UMAP reduction is missing from the Seurat object.")
}

############################################################
# 6. CHECK REQUIRED GENES
############################################################

genes <- c(
  "DSG3",
  "FLG",
  "DSG1"
)

cat("\nGenes available:\n")

print(
  genes %in% rownames(melanoma)
)

############################################################
# 7. CONVERT SEURAT TO SINGLECELLEXPERIMENT
############################################################

sce <- as.SingleCellExperiment(
  melanoma
)

############################################################
# 8. TRANSFER UMAP COORDINATES
############################################################

reducedDim(
  sce,
  "UMAP"
) <- Embeddings(
  melanoma,
  "umap"
)

############################################################
# 9. ADD CELL-TYPE LABELS
############################################################

sce$cell_type <- melanoma$cell_type

cat("\nCell types in SCE:\n")

print(
  table(sce$cell_type)
)

############################################################
# 10. RUN SLINGSHOT
#
# Starting cluster = Keratinocytes
############################################################

sce <- slingshot(
  sce,
  clusterLabels = "cell_type",
  reducedDim = "UMAP",
  start.clus = "Keratinocytes"
)

############################################################
# 11. EXTRACT PSEUDOTIME
############################################################

pseudo_matrix <- slingPseudotime(
  sce
)

cat("\n====================================================\n")
cat("SLINGSHOT LINEAGES\n")
cat("====================================================\n")

cat(
  "Number of pseudotime lineages:",
  ncol(pseudo_matrix),
  "\n"
)

print(
  colnames(pseudo_matrix)
)

############################################################
# 12. USE FIRST SLINGSHOT LINEAGE
# Same approach as original analysis
############################################################

melanoma$pseudotime <-
  pseudo_matrix[, 1]

############################################################
# 13. PSEUDOTIME SUMMARY
############################################################

cat("\n====================================================\n")
cat("PSEUDOTIME SUMMARY\n")
cat("====================================================\n")

print(
  summary(
    melanoma$pseudotime
  )
)

cat(
  "\nCells with pseudotime:",
  sum(!is.na(melanoma$pseudotime)),
  "\n"
)

cat(
  "Cells without pseudotime:",
  sum(is.na(melanoma$pseudotime)),
  "\n"
)

############################################################
# 14. SAVE PSEUDOTIME VALUES
############################################################

pseudo_table <- data.frame(
  Cell = colnames(melanoma),
  Cell_Type = melanoma$cell_type,
  Pseudotime = melanoma$pseudotime
)

write.csv(
  pseudo_table,
  "Results/SingleCell/GSE72056_Slingshot_Pseudotime.csv",
  row.names = FALSE
)

############################################################
# 15. PANEL A — PSEUDOTIME UMAP
############################################################

p_pseudotime <- FeaturePlot(
  melanoma,
  features = "pseudotime",
  reduction = "umap",
  pt.size = 0.6,
  cols = c(
    "lightgrey",
    "blue"
  )
) +
  ggtitle(
    "Pseudotime (Slingshot)"
  ) +
  theme_classic()

print(
  p_pseudotime
)

############################################################
# 16. SAVE PSEUDOTIME UMAP
############################################################

ggsave(
  "Figures/SingleCell/Figure_Pseudotime_Slingshot.png",
  plot = p_pseudotime,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  "Figures/SingleCell/Figure_Pseudotime_Slingshot.pdf",
  plot = p_pseudotime,
  width = 7,
  height = 6
)

############################################################
# 17. FETCH GENE EXPRESSION + PSEUDOTIME
############################################################

pseudo_expr <- FetchData(
  melanoma,
  vars = c(
    "pseudotime",
    "DSG3",
    "FLG",
    "DSG1"
  )
)

############################################################
# 18. REMOVE CELLS WITHOUT PSEUDOTIME
############################################################

pseudo_expr <- pseudo_expr[
  !is.na(
    pseudo_expr$pseudotime
  ),
]

cat(
  "\nCells used for expression trajectory:",
  nrow(pseudo_expr),
  "\n"
)

############################################################
# 19. DSG3 EXPRESSION ALONG PSEUDOTIME
############################################################

p_DSG3 <- ggplot(
  pseudo_expr,
  aes(
    x = pseudotime,
    y = DSG3
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.35
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    title = "DSG3 expression along pseudotime",
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic()

############################################################
# 20. FLG EXPRESSION ALONG PSEUDOTIME
############################################################

p_FLG <- ggplot(
  pseudo_expr,
  aes(
    x = pseudotime,
    y = FLG
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.35
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    title = "FLG expression along pseudotime",
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic()

############################################################
# 21. DSG1 EXPRESSION ALONG PSEUDOTIME
############################################################

p_DSG1 <- ggplot(
  pseudo_expr,
  aes(
    x = pseudotime,
    y = DSG1
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.35
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    title = "DSG1 expression along pseudotime",
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_classic()

############################################################
# 22. COMBINE THREE TRAJECTORY PLOTS
#
# Order matches reference figure:
# DSG3 | FLG | DSG1
############################################################

p_gene_pseudotime <-
  p_DSG3 +
  p_FLG +
  p_DSG1 +
  plot_layout(
    ncol = 3
  )

print(
  p_gene_pseudotime
)

############################################################
# 23. SAVE COMBINED TRAJECTORY FIGURE
############################################################

ggsave(
  "Figures/SingleCell/Figure_FLG_DSG1_DSG3_Pseudotime.png",
  plot = p_gene_pseudotime,
  width = 13,
  height = 5,
  dpi = 300
)

ggsave(
  "Figures/SingleCell/Figure_FLG_DSG1_DSG3_Pseudotime.pdf",
  plot = p_gene_pseudotime,
  width = 13,
  height = 5
)

############################################################
# 24. SAVE INDIVIDUAL DSG3 PLOT
############################################################

ggsave(
  "Figures/SingleCell/DSG3_Expression_Along_Pseudotime.png",
  plot = p_DSG3,
  width = 5,
  height = 5,
  dpi = 300
)

############################################################
# 25. SAVE INDIVIDUAL FLG PLOT
############################################################

ggsave(
  "Figures/SingleCell/FLG_Expression_Along_Pseudotime.png",
  plot = p_FLG,
  width = 5,
  height = 5,
  dpi = 300
)

############################################################
# 26. SAVE INDIVIDUAL DSG1 PLOT
############################################################

ggsave(
  "Figures/SingleCell/DSG1_Expression_Along_Pseudotime.png",
  plot = p_DSG1,
  width = 5,
  height = 5,
  dpi = 300
)

############################################################
# 27. SAVE UPDATED SEURAT OBJECT
############################################################

saveRDS(
  melanoma,
  "RData/GSE72056_Annotated_With_Pseudotime.rds"
)

############################################################
# 28. SAVE SCRIPT OBJECTS
############################################################

save(
  melanoma,
  sce,
  pseudo_matrix,
  pseudo_table,
  pseudo_expr,
  p_pseudotime,
  p_DSG3,
  p_FLG,
  p_DSG1,
  file =
    "RData/Script15_GSE72056_Slingshot_Pseudotime.RData"
)

############################################################
# 29. FINAL OUTPUT
############################################################

cat("\n====================================================\n")
cat("SCRIPT 15 COMPLETED\n")
cat("====================================================\n")

cat(
  "\nNumber of Slingshot lineages:",
  ncol(pseudo_matrix),
  "\n"
)

cat(
  "Cells assigned pseudotime:",
  sum(!is.na(melanoma$pseudotime)),
  "\n"
)

cat(
  "Pseudotime range:",
  range(
    melanoma$pseudotime,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nSaved figures:\n"
)

cat(
  "Figures/SingleCell/Figure_Pseudotime_Slingshot.png\n"
)

cat(
  "Figures/SingleCell/Figure_FLG_DSG1_DSG3_Pseudotime.png\n"
)

############################################################
# 30. SAVE SCRIPT HISTORY
############################################################

savehistory(
  "Scripts/13_GSE72056_Slingshot_Pseudotime.R"
)
