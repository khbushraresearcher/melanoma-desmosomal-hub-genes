############################################################
# SCRIPT 11
# GSE72056 SINGLE-CELL RNA-seq
# PREPROCESSING, CLUSTERING AND CELL-TYPE ANNOTATION
############################################################

############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(SingleCellExperiment)
library(slingshot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggpubr)

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
# 3. LOAD GSE72056 DATA
############################################################

expr_raw <- read.table(
  "Data/GSE72056_melanoma_single_cell_revised_v2.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

############################################################
# 4. EXTRACT EXPRESSION DATA
############################################################

metadata <- expr_raw[1:3, ]

expr <- expr_raw[-c(1:3), ]

############################################################
# 5. EXTRACT GENE NAMES
############################################################

gene_names <- make.unique(
  as.character(expr[, 1])
)

expr <- expr[, -1]

rownames(expr) <- gene_names

############################################################
# 6. CONVERT EXPRESSION MATRIX TO NUMERIC
############################################################

expr_matrix <- as.matrix(expr)

mode(expr_matrix) <- "numeric"

############################################################
# 7. BASIC DATA CHECK
############################################################

cat("\n====================================================\n")
cat("GSE72056 DATA CHECK\n")
cat("====================================================\n")

cat(
  "Genes:",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Cells:",
  ncol(expr_matrix),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(expr_matrix)),
  "\n"
)

############################################################
# 8. CREATE SEURAT OBJECT
############################################################

melanoma <- CreateSeuratObject(
  counts = expr_matrix,
  project = "Melanoma_scRNA"
)

############################################################
# 9. NORMALIZATION
############################################################

melanoma <- NormalizeData(
  melanoma,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

############################################################
# 10. IDENTIFY VARIABLE FEATURES
############################################################

melanoma <- FindVariableFeatures(
  melanoma,
  selection.method = "vst",
  nfeatures = 2000
)

############################################################
# 11. SCALE DATA
############################################################

melanoma <- ScaleData(
  melanoma
)

############################################################
# 12. PCA
############################################################

melanoma <- RunPCA(
  melanoma
)

############################################################
# 13. ELBOW PLOT
############################################################

png(
  "Figures/SingleCell/GSE72056_ElbowPlot.png",
  width = 2000,
  height = 1600,
  res = 300
)

print(
  ElbowPlot(melanoma)
)

dev.off()

############################################################
# 14. NEIGHBOR GRAPH
############################################################

melanoma <- FindNeighbors(
  melanoma,
  dims = 1:15
)

############################################################
# 15. CLUSTERING
############################################################

melanoma <- FindClusters(
  melanoma,
  resolution = 0.5
)

############################################################
# 16. UMAP
############################################################

melanoma <- RunUMAP(
  melanoma,
  dims = 1:15
)

############################################################
# 17. UMAP BY CLUSTER
############################################################

png(
  "Figures/SingleCell/GSE72056_UMAP_Clusters.png",
  width = 2400,
  height = 1800,
  res = 300
)

print(
  DimPlot(
    melanoma,
    reduction = "umap",
    label = TRUE
  )
)

dev.off()

############################################################
# 18. CELL-TYPE ANNOTATION
############################################################

cluster_to_celltype <- c(
  "0"  = "Melanoma cells",
  "1"  = "Keratinocytes",
  "2"  = "B cells",
  "3"  = "Macrophages",
  "4"  = "Endothelial",
  "5"  = "Fibroblasts",
  "6"  = "Pericytes",
  "7"  = "Fibroblasts",
  "8"  = "Other",
  "9"  = "NK cells",
  "10" = "Other",
  "11" = "Endothelial",
  "12" = "Other",
  "13" = "Other",
  "14" = "Pericytes",
  "15" = "Other"
)

############################################################
# 19. RENAME IDENTITIES
############################################################

melanoma <- RenameIdents(
  melanoma,
  cluster_to_celltype
)

############################################################
# 20. STORE CELL-TYPE ANNOTATION
############################################################

melanoma$cell_type <- Idents(
  melanoma
)

############################################################
# 21. CELL-TYPE COUNTS
############################################################

cell_type_counts <- as.data.frame(
  table(
    melanoma$cell_type
  )
)

colnames(
  cell_type_counts
) <- c(
  "CellType",
  "CellCount"
)

print(
  cell_type_counts
)

write.csv(
  cell_type_counts,
  "Results/SingleCell/GSE72056_CellType_Counts.csv",
  row.names = FALSE
)

############################################################
# 22. ANNOTATED UMAP
############################################################

png(
  "Figures/SingleCell/GSE72056_UMAP_Annotated_CellTypes.png",
  width = 2400,
  height = 1800,
  res = 300
)

print(
  DimPlot(
    melanoma,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(
      "GSE72056 melanoma single-cell annotation"
    )
)

dev.off()

############################################################
# 23. CHECK DESMOSOMAL GENES
############################################################

genes <- c(
  "FLG",
  "DSG1",
  "DSG3"
)

genes_present <- genes[
  genes %in% rownames(melanoma)
]

cat(
  "\nGenes present in the dataset:\n"
)

print(
  genes_present
)

############################################################
# 24. DESMOSOMAL GENE UMAP
############################################################

png(
  "Figures/SingleCell/GSE72056_FLG_DSG1_DSG3_UMAP.png",
  width = 2400,
  height = 1800,
  res = 300
)

print(
  FeaturePlot(
    melanoma,
    features = genes_present,
    ncol = 3
  )
)

dev.off()

############################################################
# 25. SAVE SEURAT OBJECT
############################################################

saveRDS(
  melanoma,
  "RData/GSE72056_Annotated_Seurat.rds"
)

############################################################
# 26. SAVE WORKSPACE
############################################################

save(
  melanoma,
  cluster_to_celltype,
  cell_type_counts,
  genes_present,
  file =
    "RData/Script13_GSE72056_Preprocessing_Annotation.RData"
)

############################################################
# 27. FINAL SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "SCRIPT 13 COMPLETED\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "Number of genes:",
  nrow(melanoma),
  "\n"
)

cat(
  "Number of cells:",
  ncol(melanoma),
  "\n"
)

cat(
  "\nCell-type distribution:\n"
)

print(
  cell_type_counts
)

cat(
  "\nDesmosomal genes detected:\n"
)

print(
  genes_present
)

############################################################
# END OF SCRIPT 11
############################################################

savehistory(
  "Scripts/11_scRNA_Preprocessing_Annotation.R"
)
