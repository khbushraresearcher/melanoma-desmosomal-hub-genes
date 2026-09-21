

############################################################
# SCRIPT 10
# TCGA-SKCM PROGNOSTIC ANALYSIS
#
# Candidate genes:
# FLG, DSG3, DSG1, DSP, DSC1
#
# Analyses:
# 1. TCGA-SKCM expression preprocessing
# 2. Clinical data processing
# 3. Overall survival construction
# 4. LASSO-Cox regression
# 5. RiskScore construction
# 6. Kaplan-Meier survival analysis
# 7. Univariate Cox regression
# 8. Multivariate Cox regression
# 9. Proportional hazards test
# 10. Time-varying stage sensitivity analysis
# 11. Time-dependent ROC
# 12. Bootstrap validation
# 13. Full bootstrap LASSO validation
# 14. Gene-selection frequency
# 15. Calibration analysis
# 16. Individual-gene Cox analysis
############################################################


############################################################
# 1. PROJECT DIRECTORY
############################################################

project_dir <- path.expand(
  "~/Desktop/Melanoma_Project2"
)

setwd(project_dir)

cat(
  "\nWorking directory:\n",
  getwd(),
  "\n"
)


############################################################
# 2. CREATE REQUIRED DIRECTORIES
############################################################

dir.create(
  "Data",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Prognostic",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "RData",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Scripts",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. LOAD REQUIRED PACKAGES
############################################################

library(survival)
library(glmnet)
library(ggplot2)
library(pROC)
library(survminer)
library(timeROC)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rms)


############################################################
# 4. RANDOM SEED
############################################################

set.seed(1234)


############################################################
# 5. INPUT FILES
############################################################

expression_file <- paste0(
  "Data/",
  "TCGA-SKCM.star_fpkm.tsv"
)

clinical_file <- paste0(
  "Data/",
  "TCGA-SKCM.clinical 3.tsv"
)


############################################################
# 6. CHECK INPUT FILES
############################################################

if (!file.exists(expression_file)) {
  stop(
    paste(
      "Expression file not found:",
      expression_file
    )
  )
}

if (!file.exists(clinical_file)) {
  stop(
    paste(
      "Clinical file not found:",
      clinical_file
    )
  )
}


############################################################
# 7. LOAD TCGA EXPRESSION DATA
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "LOADING TCGA-SKCM EXPRESSION DATA\n"
)

cat(
  "====================================================\n"
)


TCGA.exp <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Expression dimensions:\n"
)

print(
  dim(TCGA.exp)
)


############################################################
# 8. LOAD TCGA CLINICAL DATA
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "LOADING TCGA-SKCM CLINICAL DATA\n"
)

cat(
  "====================================================\n"
)


TCGA.clin <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Clinical dimensions:\n"
)

print(
  dim(TCGA.clin)
)


############################################################
# 9. INSPECT DATA
############################################################

cat(
  "\nExpression columns:\n"
)

print(
  head(
    colnames(TCGA.exp),
    20
  )
)


cat(
  "\nClinical columns:\n"
)

print(
  head(
    colnames(TCGA.clin),
    30
  )
)


############################################################
# 10. IDENTIFY GENE COLUMN
############################################################

gene_column <- colnames(
  TCGA.exp
)[1]


cat(
  "\nGene column:",
  gene_column,
  "\n"
)


############################################################
# 11. REMOVE ENSEMBL VERSION NUMBERS
############################################################

TCGA.exp[[gene_column]] <- sub(
  "\\..*$",
  "",
  TCGA.exp[[gene_column]]
)


############################################################
# 12. MAP ENSEMBL IDS TO GENE SYMBOLS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "MAPPING ENSEMBL IDS TO GENE SYMBOLS\n"
)

cat(
  "====================================================\n"
)


ensembl_ids <- TCGA.exp[[gene_column]]


gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)


TCGA.exp$GeneSymbol <- unname(
  gene_symbols
)


############################################################
# 13. REMOVE GENES WITHOUT SYMBOLS
############################################################

TCGA.exp <- TCGA.exp[
  !is.na(
    TCGA.exp$GeneSymbol
  ) &
    TCGA.exp$GeneSymbol != "",
]


############################################################
# 14. REMOVE DUPLICATE GENE SYMBOLS
############################################################

TCGA.exp <- TCGA.exp[
  !duplicated(
    TCGA.exp$GeneSymbol
  ),
]


############################################################
# 15. CREATE EXPRESSION MATRIX
############################################################

expression_columns <- setdiff(
  colnames(TCGA.exp),
  c(
    gene_column,
    "GeneSymbol"
  )
)


TCGA_expression <- as.matrix(
  TCGA.exp[
    ,
    expression_columns,
    drop = FALSE
  ]
)


rownames(
  TCGA_expression
) <- TCGA.exp$GeneSymbol


mode(
  TCGA_expression
) <- "numeric"


############################################################
# 16. CHECK EXPRESSION MATRIX
############################################################

cat(
  "\nExpression matrix dimensions:\n"
)

print(
  dim(TCGA_expression)
)


############################################################
# 17. CANDIDATE PROGNOSTIC GENES
############################################################

prognostic_genes <- c(
  "FLG",
  "DSG3",
  "DSG1",
  "DSP",
  "DSC1"
)


cat(
  "\nCandidate prognostic genes:\n"
)

print(
  prognostic_genes
)


############################################################
# 18. CHECK CANDIDATE GENES
############################################################

gene_presence <- prognostic_genes %in%
  rownames(TCGA_expression)


cat(
  "\nCandidate gene availability:\n"
)

print(
  data.frame(
    Gene = prognostic_genes,
    Present = gene_presence
  )
)


############################################################
# 19. EXTRACT CANDIDATE GENE EXPRESSION
############################################################

missing_genes <- prognostic_genes[
  !gene_presence
]


if (length(missing_genes) > 0) {
  
  stop(
    paste(
      "The following candidate genes are missing:",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  )
  
}


candidate_expression <- t(
  TCGA_expression[
    prognostic_genes,
    ,
    drop = FALSE
  ]
)


candidate_expression <- as.data.frame(
  candidate_expression
)


candidate_expression$sample_id <- rownames(
  candidate_expression
)


############################################################
# 20. CLINICAL SAMPLE IDENTIFIERS
############################################################

cat(
  "\nClinical column names:\n"
)

print(
  colnames(TCGA.clin)
)


############################################################
# 21. FIND SAMPLE / PATIENT IDENTIFIER
############################################################

possible_id_columns <- c(
  "submitter_id",
  "case_submitter_id",
  "bcr_patient_barcode",
  "patient_id",
  "sample_id"
)


id_column <- possible_id_columns[
  possible_id_columns %in%
    colnames(TCGA.clin)
][1]


if (is.na(id_column)) {
  
  stop(
    "No appropriate patient identifier column was found."
  )
  
}


cat(
  "\nClinical ID column:",
  id_column,
  "\n"
)


############################################################
# 22. CREATE PATIENT IDS
############################################################

candidate_expression$patient_id <- substr(
  candidate_expression$sample_id,
  1,
  12
)


TCGA.clin$patient_id <- substr(
  TCGA.clin[[id_column]],
  1,
  12
)


############################################################
# 23. REMOVE DUPLICATE CLINICAL PATIENTS
############################################################

TCGA.clin <- TCGA.clin[
  !duplicated(
    TCGA.clin$patient_id
  ),
]


############################################################
# 24. MERGE EXPRESSION AND CLINICAL DATA
############################################################

analysis_data <- merge(
  candidate_expression,
  TCGA.clin,
  by = "patient_id",
  all = FALSE
)


############################################################
# 25. CHECK MERGED DATA
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "MERGED TCGA-SKCM DATA\n"
)

cat(
  "====================================================\n"
)

cat(
  "Number of patients:",
  nrow(analysis_data),
  "\n"
)


############################################################
# 26. IDENTIFY VITAL STATUS
############################################################

vital_status_columns <- grep(
  "vital_status",
  colnames(analysis_data),
  ignore.case = TRUE,
  value = TRUE
)


print(
  vital_status_columns
)


############################################################
# 27. CREATE OS EVENT
############################################################

analysis_data$OS_event <- ifelse(
  analysis_data$vital_status.demographic ==
    "Dead",
  1,
  0
)


############################################################
# 28. CREATE OVERALL SURVIVAL TIME
############################################################

analysis_data$OS_time <- ifelse(
  analysis_data$OS_event == 1,
  
  as.numeric(
    analysis_data$days_to_death.demographic
  ),
  
  as.numeric(
    analysis_data$days_to_last_follow_up.diagnoses
  )
)


############################################################
# 29. CHECK SURVIVAL VARIABLES
############################################################

cat(
  "\nOS event distribution:\n"
)

print(
  table(
    analysis_data$OS_event,
    useNA = "ifany"
  )
)


cat(
  "\nOS time summary:\n"
)

print(
  summary(
    analysis_data$OS_time
  )
)


############################################################
# 30. SAVE MERGED DATA
############################################################

write.csv(
  analysis_data,
  "Results/Prognostic/TCGA_SKCM_Merged_Analysis_Data.csv",
  row.names = FALSE
)


############################################################
# 31. CHECK CANDIDATE GENE EXPRESSION
############################################################

cat(
  "\nCandidate gene expression summary:\n"
)


for (
  gene in prognostic_genes
) {
  
  cat(
    "\n",
    gene,
    ":\n"
  )
  
  print(
    summary(
      analysis_data[[gene]]
    )
  )
  
}


############################################################
# 32. CREATE FINAL VALID SURVIVAL COHORT
############################################################

analysis_data_surv <- analysis_data[
  !is.na(
    analysis_data$OS_time
  ) &
    analysis_data$OS_time > 0,
  ,
  drop = FALSE
]


cat(
  "\n====================================================\n"
)

cat(
  "FINAL SURVIVAL COHORT\n"
)

cat(
  "====================================================\n"
)

cat(
  "Total merged cohort:",
  nrow(analysis_data),
  "\n"
)

cat(
  "Valid OS cohort:",
  nrow(analysis_data_surv),
  "\n"
)

cat(
  "Patients removed:",
  nrow(analysis_data) -
    nrow(analysis_data_surv),
  "\n"
)


############################################################
# 33. EVENT DISTRIBUTION
############################################################

cat(
  "\nOS event distribution:\n"
)

print(
  table(
    analysis_data_surv$OS_event,
    useNA = "ifany"
  )
)


############################################################
# 34. OS TIME SUMMARY
############################################################

cat(
  "\nOS time summary:\n"
)

print(
  summary(
    analysis_data_surv$OS_time
  )
)


############################################################
# 35. CHECK INVALID VALUES
############################################################

cat(
  "\nMissing OS time:\n"
)

print(
  sum(
    is.na(
      analysis_data_surv$OS_time
    )
  )
)


cat(
  "\nNon-positive OS time:\n"
)

print(
  sum(
    analysis_data_surv$OS_time <= 0,
    na.rm = TRUE
  )
)


############################################################
# 36. FIND SURVIVAL-RELATED COLUMNS
############################################################

survival_cols <- grep(
  "follow|death|vital|survival|last",
  colnames(analysis_data),
  ignore.case = TRUE,
  value = TRUE
)


print(
  survival_cols
)


############################################################
# 37. CHECK SURVIVAL VALUES
############################################################

for (
  x in survival_cols
) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    x,
    "\n"
  )
  
  cat(
    "========================================\n"
  )
  
  print(
    head(
      analysis_data[[x]],
      10
    )
  )
  
}


############################################################
# 38. FINAL SURVIVAL COHORT SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "FINAL SURVIVAL COHORT SUMMARY\n"
)

cat(
  "====================================================\n"
)

cat(
  "Total merged cohort:",
  nrow(analysis_data),
  "\n"
)

cat(
  "Valid OS cohort:",
  nrow(analysis_data_surv),
  "\n"
)

cat(
  "Deaths:",
  sum(
    analysis_data_surv$OS_event == 1
  ),
  "\n"
)

cat(
  "Censored:",
  sum(
    analysis_data_surv$OS_event == 0
  ),
  "\n"
)

cat(
  "Median OS time:",
  median(
    analysis_data_surv$OS_time
  ),
  "days\n"
)


############################################################
# 39. SAVE PRE-LASSO WORKSPACE
############################################################

save(
  analysis_data,
  analysis_data_surv,
  candidate_expression,
  candidate_expression,
  TCGA_expression,
  prognostic_genes,
  file =
    "RData/Script10_TCGA_PreLASSO_Workspace.RData"
)


############################################################
# 40. LASSO-COX ANALYSIS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "LASSO-COX REGRESSION\n"
)

cat(
  "====================================================\n"
)


lasso_data <- analysis_data_surv[
  ,
  c(
    prognostic_genes,
    "OS_time",
    "OS_event"
  ),
  drop = FALSE
]


lasso_data <- lasso_data[
  complete.cases(
    lasso_data
  ),
  ,
  drop = FALSE
]


x_lasso <- as.matrix(
  lasso_data[
    ,
    prognostic_genes,
    drop = FALSE
  ]
)


y_lasso <- Surv(
  lasso_data$OS_time,
  lasso_data$OS_event
)


set.seed(1234)


cv_lasso <- cv.glmnet(
  x = x_lasso,
  y = y_lasso,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance",
  standardize = TRUE
)


############################################################
# 41. SAVE LASSO CROSS-VALIDATION PLOT
############################################################

pdf(
  "Figures/TCGA_SKCM_LASSO_CV.pdf",
  width = 8,
  height = 7
)

plot(
  cv_lasso
)

dev.off()


############################################################
# 42. LAMBDA VALUES
############################################################

lambda_min <- cv_lasso$lambda.min

lambda_1se <- cv_lasso$lambda.1se


cat(
  "\nLambda.min:",
  lambda_min,
  "\n"
)

cat(
  "Lambda.1se:",
  lambda_1se,
  "\n"
)


############################################################
# 43. LASSO COEFFICIENTS
############################################################

lasso_coef <- coef(
  cv_lasso,
  s = "lambda.min"
)


lasso_coef_matrix <- as.matrix(
  lasso_coef
)


lasso_coef_df <- data.frame(
  Gene = rownames(
    lasso_coef_matrix
  ),
  Coefficient =
    as.numeric(
      lasso_coef_matrix[, 1]
    ),
  row.names = NULL
)


write.csv(
  lasso_coef_df,
  "Results/Prognostic/TCGA_SKCM_LASSO_Coefficients.csv",
  row.names = FALSE
)


############################################################
# 44. SELECT NONZERO GENES
############################################################

selected_genes <- lasso_coef_df$Gene[
  lasso_coef_df$Coefficient != 0
]


selected_genes


############################################################
# 45. SELECTED COEFFICIENTS
############################################################

selected_coefficients <- lasso_coef_df[
  lasso_coef_df$Coefficient != 0,
  ,
  drop = FALSE
]


print(
  selected_coefficients
)


############################################################
# 46. SAVE SELECTED GENES
############################################################

write.csv(
  selected_coefficients,
  "Results/Prognostic/TCGA_SKCM_Selected_LASSO_Genes.csv",
  row.names = FALSE
)


############################################################
# 47. CREATE RISKSCORE
############################################################

risk_expression <- as.matrix(
  lasso_data[
    ,
    selected_genes,
    drop = FALSE
  ]
)


risk_coefficients <- selected_coefficients$Coefficient


RiskScore <- as.numeric(
  risk_expression %*%
    risk_coefficients
)


analysis_data_surv$RiskScore <- RiskScore


############################################################
# 48. COMPLETE RISKSCORE EQUATION
############################################################

risk_terms <- paste(
  selected_coefficients$Coefficient,
  "*",
  selected_coefficients$Gene,
  collapse = " + "
)


risk_equation <- paste(
  "RiskScore =",
  risk_terms
)


cat(
  "\nRiskScore equation:\n"
)

cat(
  risk_equation,
  "\n"
)


writeLines(
  risk_equation,
  "Results/Prognostic/TCGA_SKCM_RiskScore_Equation.txt"
)


############################################################
# 49. MEDIAN RISKSCORE CUTOFF
############################################################

risk_cutoff <- median(
  analysis_data_surv$RiskScore,
  na.rm = TRUE
)


analysis_data_surv$RiskGroup <- ifelse(
  analysis_data_surv$RiskScore >=
    risk_cutoff,
  "High-risk",
  "Low-risk"
)


analysis_data_surv$RiskGroup <- factor(
  analysis_data_surv$RiskGroup,
  levels = c(
    "Low-risk",
    "High-risk"
  )
)


cat(
  "\nRiskScore cutoff:",
  risk_cutoff,
  "\n"
)


print(
  table(
    analysis_data_surv$RiskGroup
  )
)


############################################################
# 50. SAVE RISKSCORE DATA
############################################################

write.csv(
  analysis_data_surv,
  "Results/Prognostic/TCGA_SKCM_RiskScore_Data.csv",
  row.names = FALSE
)

############################################################
# SAVE PART 1 SCRIPT HISTORY
############################################################

savehistory(
  "Scripts/10_TCGA_Prognostic_Analysis.R"
)



