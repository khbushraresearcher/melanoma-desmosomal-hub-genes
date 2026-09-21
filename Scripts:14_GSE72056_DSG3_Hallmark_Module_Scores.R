############################################################
# SCRIPT 14
# DSG3-ASSOCIATED HALLMARK MODULE SCORES
# GSE72056 MELANOMA
############################################################

############################################################
# 1. LOAD PACKAGES
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
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
# 3. LOAD MELANOMA OBJECT
############################################################

mel <- readRDS(
  "RData/GSE72056_Melanoma_DSG3_Analysis.rds"
)

############################################################
# 4. CHECK DSG3 GROUP
############################################################

if (
  !"DSG3_group" %in%
  colnames(mel@meta.data)
) {
  
  mel$DSG3_group <- ifelse(
    FetchData(
      mel,
      "DSG3"
    )[, 1] > 0,
    "DSG3_pos",
    "DSG3_neg"
  )
  
}

Idents(mel) <- "DSG3_group"

############################################################
# 5. CHECK GROUPS
############################################################

cat("\n====================================================\n")
cat("DSG3 GROUPS\n")
cat("====================================================\n")

print(
  table(
    mel$DSG3_group
  )
)

############################################################
# 6. DOWNLOAD HALLMARK GENE SETS
############################################################

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

############################################################
# 7. SELECT PATHWAYS
#
# These correspond to the pathway/module-score analysis
# in the original single-cell workflow.
############################################################

hallmark_paths <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

############################################################
# 8. CREATE PATHWAY GENE LISTS
############################################################

pathway_genes <- list()

for (
  pathway in hallmark_paths
) {
  
  genes <- hallmark$gene_symbol[
    hallmark$gs_name == pathway
  ]
  
  genes <- unique(
    genes[
      genes %in% rownames(mel)
    ]
  )
  
  pathway_genes[[pathway]] <- genes
  
  cat(
    "\n",
    pathway,
    ": ",
    length(genes),
    " genes present",
    sep = ""
  )
}

############################################################
# 9. ADD MODULE SCORES
############################################################

for (
  pathway in names(pathway_genes)
) {
  
  genes <- pathway_genes[[pathway]]
  
  if (
    length(genes) >= 5
  ) {
    
    score_name <- paste0(
      "DSG3_",
      gsub(
        "HALLMARK_",
        "",
        pathway
      )
    )
    
    mel <- AddModuleScore(
      object = mel,
      features = list(genes),
      name = score_name,
      assay = DefaultAssay(mel)
    )
  }
}

############################################################
# 10. IDENTIFY SCORE COLUMNS
############################################################

score_columns <- grep(
  "^DSG3_.*1$",
  colnames(
    mel@meta.data
  ),
  value = TRUE
)

cat(
  "\n\nModule-score columns:\n"
)

print(
  score_columns
)

############################################################
# 11. SAVE MODULE SCORE DATA
############################################################

score_data <- FetchData(
  mel,
  vars = c(
    "DSG3_group",
    score_columns
  )
)

write.csv(
  score_data,
  "Results/SingleCell/GSE72056_DSG3_Hallmark_Module_Scores.csv",
  row.names = TRUE
)

############################################################
# 12. CREATE INDIVIDUAL VIOLIN PLOTS
############################################################

module_plots <- list()

for (
  score in score_columns
) {
  
  pathway_name <- gsub(
    "_1$",
    "",
    gsub(
      "^DSG3_",
      "",
      score
    )
  )
  
  pathway_name <- gsub(
    "_",
    " ",
    pathway_name
  )
  
  pathway_name <- tools::toTitleCase(
    pathway_name
  )
  
  p <- VlnPlot(
    mel,
    features = score,
    group.by = "DSG3_group",
    pt.size = 0
  ) +
    
    stat_compare_means(
      method = "wilcox.test"
    ) +
    
    labs(
      title = pathway_name,
      x = NULL,
      y = "Module score"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        size = 9,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  module_plots[[score]] <- p
}

############################################################
# 13. COMBINE MODULE-SCORE PLOTS
############################################################

p_modules <- wrap_plots(
  module_plots,
  ncol = 2
)

print(
  p_modules
)

############################################################
# 14. SAVE MODULE-SCORE FIGURE
############################################################

ggsave(
  "Figures/SingleCell/Figure_DSG3_Hallmark_Module_Scores.png",
  plot = p_modules,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  "Figures/SingleCell/Figure_DSG3_Hallmark_Module_Scores.pdf",
  plot = p_modules,
  width = 10,
  height = 10
)

############################################################
# 15. STATISTICAL COMPARISON
############################################################

wilcox_results <- data.frame()

for (
  score in score_columns
) {
  
  pos_values <- score_data[
    score_data$DSG3_group == "DSG3_pos",
    score
  ]
  
  neg_values <- score_data[
    score_data$DSG3_group == "DSG3_neg",
    score
  ]
  
  test <- wilcox.test(
    pos_values,
    neg_values
  )
  
  wilcox_results <- rbind(
    wilcox_results,
    data.frame(
      Module = score,
      DSG3_Pos_Median =
        median(
          pos_values,
          na.rm = TRUE
        ),
      DSG3_Neg_Median =
        median(
          neg_values,
          na.rm = TRUE
        ),
      P_value =
        test$p.value
    )
  )
}

############################################################
# 16. MULTIPLE-TEST CORRECTION
############################################################

wilcox_results$FDR <-
  p.adjust(
    wilcox_results$P_value,
    method = "BH"
  )

############################################################
# 17. SAVE STATISTICAL RESULTS
############################################################

write.csv(
  wilcox_results,
  "Results/SingleCell/GSE72056_DSG3_Hallmark_Wilcoxon_Results.csv",
  row.names = FALSE
)

############################################################
# 18. PRINT RESULTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "HALLMARK MODULE-SCORE COMPARISON\n"
)

cat(
  "====================================================\n\n"
)

print(
  wilcox_results
)

############################################################
# 19. SAVE OBJECT
############################################################

saveRDS(
  mel,
  "RData/GSE72056_DSG3_Hallmark_Scored.rds"
)

save(
  mel,
  hallmark,
  pathway_genes,
  score_columns,
  score_data,
  wilcox_results,
  file =
    "RData/Script16_DSG3_Hallmark_Module_Scores.RData"
)

############################################################
# 20. COMPLETION MESSAGE
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "SCRIPT 16 COMPLETED\n"
)

cat(
  "====================================================\n"
)

cat(
  "\nModule-score figure saved:\n"
)

cat(
  "Figures/SingleCell/Figure_DSG3_Hallmark_Module_Scores.png\n"
)

cat(
  "\nStatistical results saved:\n"
)

cat(
  "Results/SingleCell/GSE72056_DSG3_Hallmark_Wilcoxon_Results.csv\n"
)

############################################################
# 21. SAVE COMMAND HISTORY
############################################################

savehistory(
  "Scripts/14_GSE72056_DSG3_Hallmark_Module_Scores.R"
)