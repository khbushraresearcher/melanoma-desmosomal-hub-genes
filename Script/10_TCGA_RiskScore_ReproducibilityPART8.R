############################################################
# SCRIPT 10 PART 8 
# TCGA-SKCM RISKSCORE REPRODUCIBILITY AND MODEL EXPORT
############################################################

############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(glmnet)

############################################################
# 2. PROJECT DIRECTORIES
############################################################

setwd("~/Desktop/Melanoma_Project2")

dir.create(
  "Results/Prognostic",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "RData",
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3. CHECK REQUIRED OBJECTS
############################################################

cat("\n====================================================\n")
cat("SCRIPT 11: RISKSCORE REPRODUCIBILITY EXPORT\n")
cat("====================================================\n\n")

required_objects <- c(
  "analysis_data_surv",
  "selected_genes"
)

missing_objects <- required_objects[
  !sapply(required_objects, exists)
]

cat("Missing required objects:\n")
print(missing_objects)

if (length(missing_objects) > 0) {
  stop(
    paste(
      "The following required objects are missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}

############################################################
# 4. IDENTIFY THE FITTED LASSO MODEL
############################################################

# IMPORTANT:
# cv.glmnet is the glmnet FUNCTION.
# cv_lasso is the fitted LASSO model.

lasso_model <- cv_lasso

cat(
  "\nUsing fitted LASSO model object: cv_lasso\n"
)

############################################################
# 5. IDENTIFY LAMBDA.MIN
############################################################

if (exists("lambda_min")) {
  
  lambda_used <- lambda_min
  
} else if (!is.null(lasso_model$lambda.min)) {
  
  lambda_used <- lasso_model$lambda.min
  
} else {
  
  stop(
    "lambda.min could not be identified."
  )
  
}

cat(
  "\nLambda.min:",
  lambda_used,
  "\n"
)

############################################################
# 6. EXTRACT LASSO COEFFICIENTS
############################################################

lasso_coef <- coef(
  lasso_model,
  s = lambda_used
)

lasso_coef_matrix <- as.matrix(
  lasso_coef
)

lasso_coefficient_table <- data.frame(
  Gene = rownames(lasso_coef_matrix),
  Coefficient = as.numeric(
    lasso_coef_matrix[, 1]
  ),
  stringsAsFactors = FALSE
)

############################################################
# 7. REMOVE ZERO COEFFICIENTS
############################################################

lasso_selected_coefficients <-
  lasso_coefficient_table[
    lasso_coefficient_table$Coefficient != 0,
  ]

rownames(
  lasso_selected_coefficients
) <- NULL

############################################################
# 8. PRINT LASSO COEFFICIENTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "LASSO COEFFICIENTS\n"
)

cat(
  "====================================================\n\n"
)

print(
  lasso_selected_coefficients
)

############################################################
# 9. SAVE ALL COEFFICIENTS
############################################################

write.csv(
  lasso_coefficient_table,
  "Results/Prognostic/TCGA_SKCM_LASSO_All_Coefficients.csv",
  row.names = FALSE
)

############################################################
# 10. SAVE SELECTED COEFFICIENTS
############################################################

write.csv(
  lasso_selected_coefficients,
  "Results/Prognostic/TCGA_SKCM_LASSO_Selected_Coefficients.csv",
  row.names = FALSE
)

############################################################
# 11. CREATE RISKSCORE EQUATION
############################################################

equation_terms <- paste(
  lasso_selected_coefficients$Coefficient,
  "*",
  lasso_selected_coefficients$Gene,
  collapse = " + "
)

risk_score_equation <- paste(
  "RiskScore =",
  equation_terms
)

cat(
  "\n====================================================\n"
)

cat(
  "RISKSCORE EQUATION\n"
)

cat(
  "====================================================\n\n"
)

cat(
  risk_score_equation,
  "\n"
)

############################################################
# 12. SAVE RISKSCORE EQUATION
############################################################

writeLines(
  risk_score_equation,
  "Results/Prognostic/TCGA_SKCM_RiskScore_Equation.txt"
)

############################################################
# 13. SAVE LAMBDA INFORMATION
############################################################

lambda_information <- data.frame(
  Lambda_Min = lambda_used,
  Number_Selected_Genes =
    nrow(lasso_selected_coefficients)
)

write.csv(
  lambda_information,
  "Results/Prognostic/TCGA_SKCM_LASSO_Lambda_Information.csv",
  row.names = FALSE
)

############################################################
# 14. SAVE MODEL INFORMATION
############################################################

model_information <- data.frame(
  Dataset = "TCGA-SKCM",
  Model = "LASSO Cox proportional hazards",
  Alpha = 1,
  Lambda = lambda_used,
  Number_Candidate_Genes = length(
    selected_genes
  ),
  Number_Selected_Genes =
    nrow(lasso_selected_coefficients)
)

write.csv(
  model_information,
  "Results/Prognostic/TCGA_SKCM_LASSO_Model_Information.csv",
  row.names = FALSE
)

############################################################
# 15. SAVE LASSO OBJECTS
############################################################

save(
  lasso_model,
  lasso_coefficient_table,
  lasso_selected_coefficients,
  lambda_used,
  risk_score_equation,
  model_information,
  file =
    "RData/Script11_TCGA_RiskScore_Reproducibility.RData"
)

############################################################
# 16. FINAL SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "SCRIPT 11 COMPLETED\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "LASSO model object: cv_lasso\n"
)

cat(
  "Lambda.min:",
  lambda_used,
  "\n\n"
)

cat(
  "Selected genes:\n"
)

print(
  lasso_selected_coefficients$Gene
)

cat(
  "\nRiskScore equation:\n"
)

cat(
  risk_score_equation,
  "\n"
)

############################################################
# END OF SCRIPT 10
############################################################

savehistory(
  "Scripts/10_TCGA_RiskScore_ReproducibilityPART8.R"
)
