############################################################
# SCRIPT 12
# MELANOMA-CELL SUBSET AND DSG3+ / DSG3- ANALYSIS
# GSE72056
############################################################

############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

############################################################
# 2. SET PROJECT DIRECTORY
############################################################

setwd("~/Desktop/Melanoma_Project2")

############################################################
# 3. CREATE OUTPUT DIRECTORIES
############################################################

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
# 4. LOAD ANNOTATED SEURAT OBJECT
############################################################

melanoma <- readRDS(
  "RData/GSE72056_Annotated_Seurat.rds"
)

cat("\n====================================================\n")
cat("SCRIPT 14: DSG3+ / DSG3- MELANOMA ANALYSIS\n")
cat("====================================================\n")

print(melanoma)

############################################################
# 5. CHECK CELL TYPES
############################################################

cat("\nCell-type distribution:\n")

print(
  table(melanoma$cell_type)
)

############################################################
# 6. SUBSET MELANOMA CELLS
############################################################

mel <- subset(
  melanoma,
  subset = cell_type == "Melanoma cells"
)

cat("\n====================================================\n")
cat("MELANOMA CELL SUBSET\n")
cat("====================================================\n")

cat(
  "Number of melanoma cells:",
  ncol(mel),
  "\n"
)

############################################################
# 7. CHECK DSG3 EXPRESSION
############################################################

cat("\nDSG3 expression summary:\n")

print(
  summary(
    FetchData(
      mel,
      "DSG3"
    )[, 1]
  )
)

############################################################
# 8. DEFINE DSG3+ AND DSG3- CELLS
############################################################

mel$DSG3_group <- ifelse(
  FetchData(
    mel,
    "DSG3"
  )[, 1] > 0,
  "DSG3_pos",
  "DSG3_neg"
)

############################################################
# 9. CHECK GROUP SIZES
############################################################

cat("\n====================================================\n")
cat("DSG3 GROUP DISTRIBUTION\n")
cat("====================================================\n")

print(
  table(
    mel$DSG3_group
  )
)

############################################################
# 10. SAVE GROUP INFORMATION
############################################################

group_table <- as.data.frame(
  table(
    mel$DSG3_group
  )
)

colnames(group_table) <- c(
  "DSG3_Group",
  "Cell_Count"
)

write.csv(
  group_table,
  "Results/SingleCell/GSE72056_DSG3_Group_Cell_Counts.csv",
  row.names = FALSE
)

############################################################
# 11. SET DSG3 GROUP AS ACTIVE IDENTITY
############################################################

Idents(mel) <- "DSG3_group"

############################################################
# 12. DIFFERENTIAL EXPRESSION
# DSG3-positive vs DSG3-negative
############################################################

cat("\n====================================================\n")
cat("DSG3+ VS DSG3- DIFFERENTIAL EXPRESSION\n")
cat("====================================================\n")

deg_dsg3 <- FindMarkers(
  mel,
  ident.1 = "DSG3_pos",
  ident.2 = "DSG3_neg",
  logfc.threshold = 0,
  min.pct = 0.1
)

############################################################
# 13. ADD GENE NAMES
############################################################

deg_dsg3$gene <- rownames(
  deg_dsg3
)

############################################################
# 14. ADD -LOG10 ADJUSTED P VALUE
############################################################

deg_dsg3$negLogP <-
  -log10(
    deg_dsg3$p_val_adj + 1e-300
  )

############################################################
# 15. DEFINE SIGNIFICANT GENES
############################################################

deg_dsg3$sig <- ifelse(
  deg_dsg3$p_val_adj < 0.05 &
    abs(deg_dsg3$avg_log2FC) > 0.25,
  "Significant",
  "NS"
)

############################################################
# 16. ORDER RESULTS
############################################################

deg_dsg3 <- deg_dsg3 %>%
  arrange(
    p_val_adj
  )

############################################################
# 17. SAVE COMPLETE DEG TABLE
############################################################

write.csv(
  deg_dsg3,
  "Results/SingleCell/GSE72056_DSG3_Pos_vs_Neg_DEG.csv",
  row.names = FALSE
)

############################################################
# 18. SIGNIFICANT GENE TABLE
############################################################

sig_dsg3 <- deg_dsg3 %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25
  )

write.csv(
  sig_dsg3,
  "Results/SingleCell/GSE72056_DSG3_Significant_DEGs.csv",
  row.names = FALSE
)

############################################################
# 19. DISPLAY DEG SUMMARY
############################################################

cat(
  "\nTotal genes tested:",
  nrow(deg_dsg3),
  "\n"
)

cat(
  "Significant genes:",
  nrow(sig_dsg3),
  "\n"
)

cat(
  "Upregulated in DSG3+ cells:",
  sum(
    sig_dsg3$avg_log2FC > 0
  ),
  "\n"
)

cat(
  "Downregulated in DSG3+ cells:",
  sum(
    sig_dsg3$avg_log2FC < 0
  ),
  "\n"
)

############################################################
# 20. TOP DIFFERENTIALLY EXPRESSED GENES
############################################################

cat("\nTop 20 genes:\n")

print(
  head(
    deg_dsg3[
      ,
      c(
        "gene",
        "avg_log2FC",
        "p_val_adj"
      )
    ],
    20
  )
)

############################################################
# 21. SELECT TOP 40 GENES FOR HEATMAP
############################################################

top_genes <- deg_dsg3 %>%
  filter(
    p_val_adj < 1e-5,
    abs(avg_log2FC) > 0.25
  ) %>%
  arrange(
    desc(
      abs(avg_log2FC)
    )
  ) %>%
  head(40) %>%
  pull(gene)

cat(
  "\nNumber of genes selected for heatmap:",
  length(top_genes),
  "\n"
)

############################################################
# 22. SAVE TOP 40 GENE LIST
############################################################

write.csv(
  data.frame(
    Gene = top_genes
  ),
  "Results/SingleCell/GSE72056_DSG3_Top40_Genes.csv",
  row.names = FALSE
)

############################################################
# 23. SCALE TOP GENES
############################################################

mel <- ScaleData(
  mel,
  features = top_genes,
  verbose = FALSE
)

############################################################
# 24. DSG3-ASSOCIATED HEATMAP
############################################################

if (
  length(top_genes) > 0
) {
  
  pdf(
    "Figures/SingleCell/Figure_DSG3_Associated_Heatmap.pdf",
    width = 9,
    height = 10
  )
  
  print(
    DoHeatmap(
      mel,
      features = top_genes,
      group.by = "DSG3_group"
    ) +
      ggtitle(
        "DSG3-associated transcriptional program"
      )
  )
  
  dev.off()
  
}

############################################################
# 25. SAVE PNG VERSION
############################################################

png(
  "Figures/SingleCell/Figure_DSG3_Associated_Heatmap.png",
  width = 2400,
  height = 2800,
  res = 300
)

print(
  DoHeatmap(
    mel,
    features = top_genes,
    group.by = "DSG3_group"
  ) +
    ggtitle(
      "DSG3-associated transcriptional program"
    )
)

dev.off()

############################################################
# 26. GO BIOLOGICAL PROCESS — OVER-REPRESENTATION
############################################################

genes_up <- deg_dsg3 %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  ) %>%
  pull(gene)

genes_up <- unique(
  genes_up[
    genes_up %in%
      keys(
        org.Hs.eg.db,
        keytype = "SYMBOL"
      )
  ]
)

cat(
  "\nGenes used for GO over-representation analysis:",
  length(genes_up),
  "\n"
)

############################################################
# 27. GO ENRICHMENT
############################################################

if (
  length(genes_up) > 0
) {
  
  ego <- enrichGO(
    gene = genes_up,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
} else {
  
  ego <- NULL
  
}

############################################################
# 28. SAVE GO RESULTS
############################################################

if (
  !is.null(ego)
) {
  
  write.csv(
    as.data.frame(ego),
    "Results/SingleCell/GSE72056_DSG3_GO_BP_ORA.csv",
    row.names = FALSE
  )
  
}

############################################################
# 29. RANKED GENE LIST FOR GSEA
############################################################

gene_ranks <- deg_dsg3$avg_log2FC

names(gene_ranks) <-
  deg_dsg3$gene

gene_ranks <- gene_ranks[
  !is.na(gene_ranks)
]

gene_ranks <- gene_ranks[
  !duplicated(
    names(gene_ranks)
  )
]

gene_ranks <- sort(
  gene_ranks,
  decreasing = TRUE
)

############################################################
# 30. GO BIOLOGICAL PROCESS GSEA
############################################################

gsea_go <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 300,
  pvalueCutoff = 0.25,
  verbose = FALSE
)

############################################################
# 31. SAVE GO GSEA
############################################################

write.csv(
  as.data.frame(gsea_go),
  "Results/SingleCell/GSE72056_DSG3_GSEA_GO_BP.csv",
  row.names = FALSE
)

############################################################
# 32. HALLMARK GENE SETS
############################################################

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_t2g <- hallmark[
  ,
  c(
    "gs_name",
    "gene_symbol"
  )
]

############################################################
# 33. HALLMARK GSEA
############################################################

gsea_hallmark <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 0.25
)

############################################################
# 34. SAVE HALLMARK GSEA
############################################################

write.csv(
  as.data.frame(
    gsea_hallmark
  ),
  "Results/SingleCell/GSE72056_DSG3_GSEA_Hallmark.csv",
  row.names = FALSE
)

############################################################
# 35. DISPLAY TOP GO RESULTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "TOP GO BIOLOGICAL PROCESS RESULTS\n"
)

cat(
  "====================================================\n"
)

if (
  !is.null(ego) &&
  nrow(as.data.frame(ego)) > 0
) {
  
  print(
    head(
      as.data.frame(ego),
      10
    )
  )
  
} else {
  
  cat(
    "No significant GO terms detected.\n"
  )
  
}

############################################################
# 36. DISPLAY TOP HALLMARK RESULTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "TOP HALLMARK GSEA RESULTS\n"
)

cat(
  "====================================================\n"
)

print(
  head(
    as.data.frame(
      gsea_hallmark
    ),
    10
  )
)

############################################################
# 37. SAVE ANALYSIS OBJECTS
############################################################

saveRDS(
  mel,
  "RData/GSE72056_Melanoma_DSG3_Analysis.rds"
)

save(
  melanoma,
  mel,
  deg_dsg3,
  sig_dsg3,
  top_genes,
  genes_up,
  ego,
  gene_ranks,
  gsea_go,
  hallmark,
  hallmark_t2g,
  gsea_hallmark,
  file =
    "RData/Script14_GSE72056_DSG3_Analysis.RData"
)

############################################################
# 38. FINAL MESSAGE
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "SCRIPT 14 COMPLETED\n"
)

cat(
  "====================================================\n"
)

cat(
  "\nDSG3+ cells:",
  sum(
    mel$DSG3_group == "DSG3_pos"
  ),
  "\n"
)

cat(
  "DSG3- cells:",
  sum(
    mel$DSG3_group == "DSG3_neg"
  ),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(sig_dsg3),
  "\n"
)

cat(
  "Top heatmap genes:",
  length(top_genes),
  "\n"
)

cat(
  "\nResults saved to:\n"
)

cat(
  "Results/SingleCell/\n"
)

cat(
  "\nFigures saved to:\n"
)

cat(
  "Figures/SingleCell/\n"
)

############################################################
# 39. SAVE COMMAND HISTORY
############################################################

savehistory(
  "Scripts/12_GSE72056_DSG3_Analysis.R"
)
