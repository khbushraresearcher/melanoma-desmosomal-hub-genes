############################################################
# PART 7
# INDIVIDUAL GENE COX REGRESSION
# TCGA-SKCM
############################################################

############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(survival)
library(ggplot2)

############################################################
# 2. CHECK DATA
############################################################

cat("\n====================================================\n")
cat("PART 7: INDIVIDUAL GENE COX REGRESSION\n")
cat("====================================================\n\n")

cat("Number of patients:", nrow(analysis_data_surv), "\n")

############################################################
# 3. CANDIDATE DESMOSOMAL GENES
############################################################

candidate_genes <- c(
  "FLG",
  "DSG3",
  "DSG1",
  "DSP",
  "DSC1"
)

print(candidate_genes)

############################################################
# 4. CHECK GENE AVAILABILITY
############################################################

missing_genes <- setdiff(
  candidate_genes,
  colnames(analysis_data_surv)
)

cat("\nMissing genes:\n")
print(missing_genes)

if (length(missing_genes) > 0) {
  stop(
    paste(
      "The following candidate genes are missing:",
      paste(missing_genes, collapse = ", ")
    )
  )
}

############################################################
# 5. INDIVIDUAL COX REGRESSION
############################################################

individual_cox_results <- data.frame()

for (gene in candidate_genes) {
  
  cat("\nAnalyzing:", gene, "\n")
  
  gene_data <- analysis_data_surv[
    complete.cases(
      analysis_data_surv[, c(
        "OS_time",
        "OS_event",
        gene
      )]
    ),
  ]
  
  formula_gene <- as.formula(
    paste(
      "Surv(OS_time, OS_event) ~",
      gene
    )
  )
  
  cox_gene <- coxph(
    formula_gene,
    data = gene_data
  )
  
  cox_summary <- summary(cox_gene)
  
  individual_cox_results <- rbind(
    individual_cox_results,
    data.frame(
      Gene = gene,
      N = nrow(gene_data),
      Events = sum(gene_data$OS_event),
      Coefficient = cox_summary$coefficients[1, "coef"],
      HR = cox_summary$coefficients[1, "exp(coef)"],
      SE = cox_summary$coefficients[1, "se(coef)"],
      Z = cox_summary$coefficients[1, "z"],
      P_value = cox_summary$coefficients[1, "Pr(>|z|)"],
      Lower_95CI = cox_summary$conf.int[1, "lower .95"],
      Upper_95CI = cox_summary$conf.int[1, "upper .95"],
      stringsAsFactors = FALSE
    )
  )
}

############################################################
# 6. MULTIPLE-TESTING CORRECTION
############################################################

individual_cox_results$FDR <- p.adjust(
  individual_cox_results$P_value,
  method = "BH"
)

############################################################
# 7. ORDER RESULTS
############################################################

individual_cox_results <- individual_cox_results[
  order(individual_cox_results$P_value),
]

rownames(individual_cox_results) <- NULL

############################################################
# 8. PRINT RESULTS
############################################################

cat("\n====================================================\n")
cat("INDIVIDUAL GENE COX RESULTS\n")
cat("====================================================\n\n")

print(individual_cox_results)

############################################################
# 9. SAVE RESULTS
############################################################

write.csv(
  individual_cox_results,
  "Results/Prognostic/TCGA_SKCM_Individual_Gene_Cox.csv",
  row.names = FALSE
)

############################################################
# 10. PREPARE FOREST PLOT DATA
############################################################

forest_data <- individual_cox_results

forest_data$Gene <- factor(
  forest_data$Gene,
  levels = rev(forest_data$Gene)
)

############################################################
# 11. FOREST PLOT
############################################################

forest_plot <- ggplot(
  forest_data,
  aes(
    x = Gene,
    y = HR
  )
) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  
  geom_errorbar(
    aes(
      ymin = Lower_95CI,
      ymax = Upper_95CI
    ),
    width = 0.15
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_y_log10() +
  
  coord_flip() +
  
  labs(
    title = "Individual Gene Cox Regression",
    x = "Gene",
    y = "Hazard Ratio (95% CI)"
  ) +
  
  theme_classic()

print(forest_plot)

############################################################
# 12. SAVE FOREST PLOT
############################################################

ggsave(
  "Figures/TCGA_SKCM_Individual_Gene_Cox_Forest.png",
  forest_plot,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
# 13. SAVE RDATA
############################################################

save(
  individual_cox_results,
  forest_data,
  forest_plot,
  file = "RData/Script10_TCGA_Individual_Gene_Cox.RData"
)

############################################################
# 14. FINAL SUMMARY
############################################################

cat("\n====================================================\n")
cat("PART 7 COMPLETED\n")
cat("====================================================\n")

cat("\nCandidate genes analyzed:\n")
print(candidate_genes)

cat("\nResults saved to:\n")
cat("Results/Prognostic/TCGA_SKCM_Individual_Gene_Cox.csv\n")

cat("\nForest plot saved to:\n")
cat("Figures/TCGA_SKCM_Individual_Gene_Cox_Forest.png\n")

############################################################
# END OF PART 7
############################################################

savehistory(
  "Scripts/10_TCGA_Prognostic_Analysis_Part7.R"
)