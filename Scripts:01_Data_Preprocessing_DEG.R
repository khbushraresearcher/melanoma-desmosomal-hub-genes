############################################################
# SCRIPT 01
# DATA PREPROCESSING AND DIFFERENTIAL EXPRESSION SETUP
# Melanoma Transcriptomic Analysis
############################################################

############################################################
# 1. SET PROJECT DIRECTORY
############################################################

setwd("~/Desktop/Melanoma_Project2")

############################################################
# 2. CREATE REQUIRED DIRECTORIES
############################################################

dir.create("Data", recursive = TRUE, showWarnings = FALSE)
dir.create("Scripts", recursive = TRUE, showWarnings = FALSE)
dir.create("Results", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Prognostic", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("RData", recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. LOAD REQUIRED PACKAGES
############################################################

library(GEOquery)
library(limma)
library(affy)
library(Biobase)
library(dplyr)

############################################################
# 4. DEFINE GEO DATASETS
############################################################

geo_ids <- c(
  "GSE3189",
  "GSE7553",
  "GSE46517"
)

############################################################
# 5. DOWNLOAD / LOAD GEO DATA
############################################################

gse3189 <- getGEO(
  "GSE3189",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

gse7553 <- getGEO(
  "GSE7553",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

gse46517 <- getGEO(
  "GSE46517",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

############################################################
# 6. EXTRACT EXPRESSION MATRICES
############################################################

expr3189 <- exprs(gse3189[[1]])
expr7553 <- exprs(gse7553[[1]])
expr46517 <- exprs(gse46517[[1]])

############################################################
# 7. EXTRACT PHENOTYPE DATA
############################################################

pheno3189 <- pData(gse3189[[1]])
pheno7553 <- pData(gse7553[[1]])
pheno46517 <- pData(gse46517[[1]])

############################################################
# 8. CHECK DATA DIMENSIONS
############################################################

cat("\nGSE3189 dimensions:\n")
print(dim(expr3189))

cat("\nGSE7553 dimensions:\n")
print(dim(expr7553))

cat("\nGSE46517 dimensions:\n")
print(dim(expr46517))

############################################################
# 9. SAVE RAW GEO OBJECTS
############################################################

save(
  gse3189,
  gse7553,
  gse46517,
  file = "RData/Script01_GEO_Raw_Objects.RData"
)

############################################################
# 10. SAVE EXPRESSION MATRICES
############################################################

write.csv(
  expr3189,
  "Data/GSE3189_Expression.csv"
)

write.csv(
  expr7553,
  "Data/GSE7553_Expression.csv"
)

write.csv(
  expr46517,
  "Data/GSE46517_Expression.csv"
)

############################################################
# 11. SAVE PHENOTYPE DATA
############################################################

write.csv(
  pheno3189,
  "Data/GSE3189_Phenotype.csv",
  row.names = TRUE
)

write.csv(
  pheno7553,
  "Data/GSE7553_Phenotype.csv",
  row.names = TRUE
)

write.csv(
  pheno46517,
  "Data/GSE46517_Phenotype.csv",
  row.names = TRUE
)

############################################################
# 12. CHECK SAMPLE INFORMATION
############################################################

cat("\nGSE3189 phenotype columns:\n")
print(colnames(pheno3189))

cat("\nGSE7553 phenotype columns:\n")
print(colnames(pheno7553))

cat("\nGSE46517 phenotype columns:\n")
print(colnames(pheno46517))

############################################################
# 13. CHECK PLATFORM INFORMATION
############################################################

cat("\nGSE3189 platform:\n")
print(unique(pheno3189$platform_id))

cat("\nGSE7553 platform:\n")
print(unique(pheno7553$platform_id))

cat("\nGSE46517 platform:\n")
print(unique(pheno46517$platform_id))

############################################################
# 14. CHECK FOR MISSING VALUES
############################################################

cat("\nMissing values in GSE3189 expression:",
    sum(is.na(expr3189)), "\n")

cat("Missing values in GSE7553 expression:",
    sum(is.na(expr7553)), "\n")

cat("Missing values in GSE46517 expression:",
    sum(is.na(expr46517)), "\n")

############################################################
# 15. SAVE PREPROCESSING WORKSPACE
############################################################

save(
  expr3189,
  expr7553,
  expr46517,
  pheno3189,
  pheno7553,
  pheno46517,
  file = "RData/Script01_Preprocessing_Workspace.RData"
)

############################################################
# END OF SCRIPT 01
############################################################