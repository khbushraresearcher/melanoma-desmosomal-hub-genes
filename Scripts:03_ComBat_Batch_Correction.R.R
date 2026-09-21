


###############################################################
## SCRIPT 03: ComBat Batch Effect Correction
## Project: Melanoma_Project2
###############################################################

setwd("~/Desktop/Melanoma_Project2")

###############################################################
## Load packages
###############################################################

library(sva)
library(ggplot2)
library(limma)

###############################################################
## Create output folders
###############################################################

dir.create("Results", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("RData", recursive = TRUE, showWarnings = FALSE)

###############################################################
## Load Script 2 workspace
###############################################################

load("RData/Script2_DEG_Workspace.RData")

###############################################################
## Load expression matrices and groups
###############################################################

load("RData/Script1_DEG_Workspace.RData")

###############################################################
## Create common genes across all datasets
###############################################################

commonGenes <- Reduce(
  intersect,
  list(
    rownames(expr3189),
    rownames(expr46517),
    rownames(expr7553)
  )
)

length(commonGenes)

###############################################################
## Subset each dataset to common genes
###############################################################

expr3189.common <- expr3189[commonGenes, ]

expr46517.common <- expr46517[commonGenes, ]

expr7553.common <- expr7553[commonGenes, ]

###############################################################
## Merge expression matrices
###############################################################

merged.expr <- cbind(
  expr3189.common,
  expr46517.common,
  expr7553.common
)

dim(merged.expr)

###############################################################
## Create batch vector
###############################################################

batch <- c(
  rep("GSE3189", ncol(expr3189.common)),
  rep("GSE46517", ncol(expr46517.common)),
  rep("GSE7553", ncol(expr7553.common))
)

table(batch)

###############################################################
## Create phenotype vector
###############################################################

group3189 <- ifelse(
  grepl("^Normal", colnames(expr3189)),
  "Normal",
  ifelse(
    grepl("^Melanoma", colnames(expr3189)),
    "Melanoma",
    "Other"
  )
)

group46517 <- ifelse(
  grepl(
    "^Primary Melanoma|^Metastatic Melanoma",
    colnames(expr46517)
  ),
  "Melanoma",
  ifelse(
    grepl(
      "^Normal Skin|^Normal Epithelial Melanocytes",
      colnames(expr46517)
    ),
    "Normal",
    "Other"
  )
)

group7553 <- ifelse(
  grepl(
    "^Primary Melanoma|^Metastatic Melanoma",
    colnames(expr7553)
  ),
  "Melanoma",
  ifelse(
    grepl(
      "^Normal Skin|^normal human epidermal melanocytes",
      colnames(expr7553),
      ignore.case = TRUE
    ),
    "Normal",
    "Other"
  )
)

group.all <- c(
  group3189,
  group46517,
  group7553
)

table(group.all)

###############################################################
## Remove Other samples
###############################################################

keep <- group.all != "Other"

merged.expr2 <- merged.expr[, keep]

batch2 <- batch[keep]

group2 <- factor(
  group.all[keep],
  levels = c("Normal", "Melanoma")
)

table(group2)

dim(merged.expr2)

###############################################################
## Run ComBat on filtered dataset
###############################################################

mod <- model.matrix(
  ~group2
)

combat.expr <- ComBat(
  dat = merged.expr2,
  batch = batch2,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

###############################################################
## PCA after ComBat
###############################################################

pca.after <- prcomp(
  t(combat.expr),
  scale. = TRUE
)

plot.df <- data.frame(
  PC1 = pca.after$x[,1],
  PC2 = pca.after$x[,2],
  Batch = batch2,
  Group = group2
)

p <- ggplot(
  plot.df,
  aes(
    PC1,
    PC2,
    color = Batch,
    shape = Group
  )
) +
  geom_point(size = 3) +
  theme_classic()

print(p)

ggsave(
  "Results/Figures/PCA_After_ComBat_Filtered.png",
  plot = p,
  width = 7,
  height = 6
)

###############################################################
## Save corrected matrix
###############################################################

save(
  combat.expr,
  group2,
  batch2,
  file = "RData/ComBat_Filtered_Merged.RData"
)

write.csv(
  combat.expr,
  "Results/Tables/ComBat_Filtered_Merged.csv"
)

###############################################################
## DEG on ComBat-corrected merged dataset
###############################################################

design <- model.matrix(
  ~0 + group2
)

colnames(design) <- levels(group2)

contrast <- makeContrasts(
  Melanoma - Normal,
  levels = design
)

fit <- lmFit(
  combat.expr,
  design
)

fit <- contrasts.fit(
  fit,
  contrast
)

fit <- eBayes(
  fit
)

deg.combat <- topTable(
  fit,
  number = Inf,
  adjust.method = "BH"
)

###############################################################
## Save ComBat DEG results
###############################################################

write.csv(
  deg.combat,
  "Results/Tables/ComBat_Merged_DEGs.csv"
)

save(
  deg.combat,
  file = "RData/ComBat_Merged_DEGs.RData"
)

###############################################################
## Significant DEGs
###############################################################

DEG.combat <- subset(
  deg.combat,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

nrow(DEG.combat)

write.csv(
  DEG.combat,
  "Results/Tables/ComBat_Significant_DEGs.csv"
)

###############################################################
## Save complete Script 3 workspace
###############################################################

save(
  commonGenes,
  merged.expr,
  batch,
  group.all,
  merged.expr2,
  batch2,
  group2,
  combat.expr,
  pca.after,
  deg.combat,
  DEG.combat,
  file = "RData/Script3_ComBat_Workspace.RData"
)

###############################################################
## Final summary
###############################################################

cat("\n")
cat("============================================\n")
cat("SCRIPT 3 COMPLETED\n")
cat("============================================\n")

cat(
  "Common genes:",
  length(commonGenes),
  "\n"
)

cat(
  "Samples after removing Other:",
  ncol(merged.expr2),
  "\n"
)

cat(
  "Significant ComBat DEGs:",
  nrow(DEG.combat),
  "\n"
)

cat("============================================\n")

setwd("~/Desktop/Melanoma_Project2")
savehistory("Scripts/03_ComBat_Batch_Correction.R")















