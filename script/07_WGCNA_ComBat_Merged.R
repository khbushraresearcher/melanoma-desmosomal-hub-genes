

###############################################################
## SCRIPT 07: WGCNA ON MERGED COMBAT-CORRECTED DATASET
## Based directly on the original WGCNApipe.R pipeline
###############################################################

###############################################################
# Set project directory
###############################################################

setwd("~/Desktop/Melanoma_Project2")

###############################################################
# Load packages
###############################################################

library(WGCNA)
library(ggplot2)

options(stringsAsFactors = FALSE)

###############################################################
# Create output folders
###############################################################

dir.create("Results", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/RData", recursive = TRUE, showWarnings = FALSE)

###############################################################
# Load original ComBat-corrected merged dataset
###############################################################

load("RData/ComBat_Filtered_Merged.RData")

###############################################################
# Check loaded objects
###############################################################

print(ls())

###############################################################
# Confirm ComBat expression matrix
###############################################################

print(dim(combat.expr))

###############################################################
# Convert expression matrix for WGCNA
# Original pipeline:
# samples = rows
# genes = columns
###############################################################

datExpr <- as.data.frame(t(combat.expr))

print(dim(datExpr))

###############################################################
# STEP 1: DATA QUALITY CONTROL
###############################################################

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

print(gsg$allOK)

###############################################################
# Remove bad samples/genes if required
###############################################################

if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) {
    cat(
      "Removing",
      sum(!gsg$goodGenes),
      "genes\n"
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    cat(
      "Removing",
      sum(!gsg$goodSamples),
      "samples\n"
    )
  }
  
  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes
  ]
}

###############################################################
# Check final expression matrix
###############################################################

print(dim(datExpr))

###############################################################
# Save cleaned ComBat expression matrix
###############################################################

write.csv(
  datExpr,
  "Results/Tables/Cleaned_ComBat_Expression.csv",
  row.names = TRUE
)

###############################################################
# Sample clustering
###############################################################

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

###############################################################
# Save sample clustering figure
###############################################################

png(
  "Results/Figures/Sample_Clustering.png",
  width = 3000,
  height = 2000,
  res = 300
)

plot(
  sampleTree,
  main = "Sample Clustering",
  xlab = "",
  sub = "",
  cex = 0.6
)

dev.off()

###############################################################
# Save sample tree
###############################################################

save(
  sampleTree,
  file = "Results/RData/SampleTree.RData"
)

###############################################################
# STEP 2: SOFT-THRESHOLD POWER SELECTION
###############################################################

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

###############################################################
# Save soft-threshold results
###############################################################

write.csv(
  sft$fitIndices,
  "Results/Tables/SoftThreshold_Results.csv",
  row.names = FALSE
)

###############################################################
# Plot scale-free topology fit
###############################################################

png(
  "Results/Figures/Figure2_SoftThreshold.png",
  width = 3200,
  height = 2600,
  res = 300
)

par(mfrow = c(1, 2))

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) *
    sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale Independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) *
    sft$fitIndices[, 2],
  labels = powers,
  col = "red",
  cex = 1
)

abline(
  h = 0.80,
  col = "blue",
  lty = 2
)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  col = "red",
  cex = 1
)

dev.off()

###############################################################
# Original selected soft power
###############################################################

softPower <- 11

cat(
  "\nSelected Soft Power:",
  softPower,
  "\n"
)

###############################################################
# Save soft-threshold objects
###############################################################

save(
  sft,
  softPower,
  file = "Results/RData/Step2_SoftThreshold.RData"
)

###############################################################
# STEP 3: NETWORK CONSTRUCTION & MODULE DETECTION
###############################################################

net <- blockwiseModules(
  datExpr,
  
  power = softPower,
  
  networkType = "signed",
  
  TOMType = "signed",
  
  maxBlockSize = 20000,
  
  deepSplit = 2,
  
  minModuleSize = 30,
  
  mergeCutHeight = 0.25,
  
  reassignThreshold = 0,
  
  pamRespectsDendro = FALSE,
  
  numericLabels = TRUE,
  
  saveTOMs = FALSE,
  
  verbose = 3
)

###############################################################
# Convert numeric module labels to colors
###############################################################

moduleColors <- labels2colors(
  net$colors
)

print(table(moduleColors))

print(
  sort(
    table(moduleColors),
    decreasing = TRUE
  )
)

###############################################################
# Save module assignment table
###############################################################

Module_Assignment <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

write.csv(
  Module_Assignment,
  "Results/Tables/Module_Assignment.csv",
  row.names = FALSE
)

###############################################################
# Number of blocks
###############################################################

print(
  length(net$dendrograms)
)

print(
  sapply(
    net$blockGenes,
    length
  )
)

###############################################################
# Plot dendrogram(s)
###############################################################

for (i in seq_along(net$dendrograms)) {
  
  png(
    paste0(
      "Results/Figures/Figure3_Dendrogram_Block_",
      i,
      ".png"
    ),
    width = 4200,
    height = 2400,
    res = 300
  )
  
  plotDendroAndColors(
    net$dendrograms[[i]],
    moduleColors[
      net$blockGenes[[i]]
    ],
    "Module Colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  
  dev.off()
}

###############################################################
# Module eigengenes
###############################################################

MEs <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs)

###############################################################
# Save module eigengenes
###############################################################

write.csv(
  MEs,
  "Results/Tables/Module_Eigengenes.csv"
)

###############################################################
# Save network objects
###############################################################

save(
  net,
  moduleColors,
  MEs,
  Module_Assignment,
  file = "Results/RData/Step3_NetworkConstruction.RData"
)

###############################################################
# STEP 4: MODULE–TRAIT RELATIONSHIP ANALYSIS
###############################################################

###############################################################
# Check group2
###############################################################

print(table(group2))

###############################################################
# Create trait data
###############################################################

traitData <- data.frame(
  Melanoma = group2,
  row.names = rownames(datExpr)
)

###############################################################
# Convert phenotype to numeric
# Normal = 0
# Melanoma = 1
###############################################################

traitData$Melanoma <- ifelse(
  traitData$Melanoma == "Melanoma",
  1,
  0
)

###############################################################
# Check trait
###############################################################

str(traitData)

print(
  table(traitData$Melanoma)
)

###############################################################
# Calculate module eigengenes
###############################################################

MEs <- moduleEigengenes(
  datExpr,
  moduleColors
)$eigengenes

MEs <- orderMEs(MEs)

###############################################################
# Module-trait correlation
###############################################################

moduleTraitCor <- cor(
  MEs,
  traitData,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nrow(datExpr)
)

###############################################################
# Save correlation tables
###############################################################

write.csv(
  moduleTraitCor,
  "Results/Tables/ModuleTrait_Correlation.csv",
  row.names = TRUE
)

write.csv(
  moduleTraitPvalue,
  "Results/Tables/ModuleTrait_Pvalue.csv",
  row.names = TRUE
)

###############################################################
# Create heatmap text
###############################################################

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 2),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(
  moduleTraitCor
)

###############################################################
# Module-trait heatmap
###############################################################

png(
  "Results/Figures/Figure4_ModuleTraitHeatmap.png",
  width = 2600,
  height = 2400,
  res = 300
)

par(
  mar = c(
    8,
    10,
    3,
    3
  )
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)

dev.off()

###############################################################
# Module summary
###############################################################

moduleSummary <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation =
    moduleTraitCor[
      , "Melanoma"
    ],
  Pvalue =
    moduleTraitPvalue[
      , "Melanoma"
    ]
)

moduleSummary <- moduleSummary[
  order(
    abs(moduleSummary$Correlation),
    decreasing = TRUE
  ),
]

###############################################################
# Save module significance table
###############################################################

write.csv(
  moduleSummary,
  "Results/Tables/Module_Significance.csv",
  row.names = FALSE
)

###############################################################
# Save module-trait objects
###############################################################

save(
  traitData,
  MEs,
  moduleTraitCor,
  moduleTraitPvalue,
  moduleSummary,
  file = "Results/RData/Step4_ModuleTrait.RData"
)

###############################################################
# Display top module
###############################################################

cat(
  "\n=====================================\n"
)

cat(
  "TOP CORRELATED MODULE\n"
)

cat(
  "=====================================\n"
)

print(
  moduleSummary[1, ]
)

###############################################################
# STEP 5: GENE SIGNIFICANCE & MODULE MEMBERSHIP
###############################################################

###############################################################
# Select the original module of interest
###############################################################

module <- "brown"

###############################################################
# Gene significance
###############################################################

geneTraitSignificance <- as.data.frame(
  cor(
    datExpr,
    traitData$Melanoma,
    use = "p"
  )
)

GSPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(
      geneTraitSignificance
    ),
    nrow(datExpr)
  )
)

colnames(
  geneTraitSignificance
) <- "GS"

colnames(
  GSPvalue
) <- "GS.Pvalue"

###############################################################
# Module membership
###############################################################

MEs <- orderMEs(
  moduleEigengenes(
    datExpr,
    moduleColors
  )$eigengenes
)

moduleMembership <- as.data.frame(
  cor(
    datExpr,
    MEs,
    use = "p"
  )
)

MMPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(
      moduleMembership
    ),
    nrow(datExpr)
  )
)

colnames(
  moduleMembership
) <- paste0(
  "MM_",
  substring(
    names(MEs),
    3
  )
)

colnames(
  MMPvalue
) <- paste0(
  "MM.P_",
  substring(
    names(MEs),
    3
  )
)

###############################################################
# Select genes in brown module
###############################################################

moduleGenes <- (
  moduleColors == module
)

###############################################################
# GS vs MM scatter plot
###############################################################

png(
  paste0(
    "Results/Figures/Figure5_GS_vs_MM_",
    module,
    ".png"
  ),
  width = 2600,
  height = 2400,
  res = 300
)

verboseScatterplot(
  abs(
    moduleMembership[
      moduleGenes,
      paste0(
        "MM_",
        module
      )
    ]
  ),
  abs(
    geneTraitSignificance[
      moduleGenes,
      1
    ]
  ),
  xlab = paste(
    "Module Membership in",
    module,
    "module"
  ),
  ylab = "Gene Significance for Melanoma",
  main = paste(
    module,
    "Module"
  ),
  col = module
)

dev.off()

###############################################################
# Hub gene table
###############################################################

hubTable <- data.frame(
  
  Gene =
    colnames(datExpr)[
      moduleGenes
    ],
  
  Module = module,
  
  GS =
    geneTraitSignificance[
      moduleGenes,
      1
    ],
  
  GS.Pvalue =
    GSPvalue[
      moduleGenes,
      1
    ],
  
  MM =
    moduleMembership[
      moduleGenes,
      paste0(
        "MM_",
        module
      )
    ],
  
  MM.Pvalue =
    MMPvalue[
      moduleGenes,
      paste0(
        "MM.P_",
        module
      )
    ]
)

###############################################################
# Sort by module membership
###############################################################

hubTable <- hubTable[
  order(
    abs(hubTable$MM),
    decreasing = TRUE
  ),
]

###############################################################
# Save hub genes
###############################################################

write.csv(
  hubTable,
  paste0(
    "Results/Tables/",
    module,
    "_HubGenes.csv"
  ),
  row.names = FALSE
)

###############################################################
# Candidate hub genes
###############################################################

candidateHubGenes <- subset(
  hubTable,
  abs(MM) > 0.80 &
    abs(GS) > 0.20
)

###############################################################
# Save candidate hub genes
# Keep EXACT filename from original pipeline
###############################################################

write.csv(
  candidateHubGenes,
  paste0(
    "Results/Tables/",
    module,
    "merged_CandidateHubGenes.csv"
  ),
  row.names = FALSE
)

###############################################################
# STRING input
###############################################################

write.table(
  candidateHubGenes$Gene,
  file =
    "Results/Tables/CandidateHubGenes_STRING.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

###############################################################
# Save Step 5 objects
###############################################################

save(
  geneTraitSignificance,
  GSPvalue,
  moduleMembership,
  MMPvalue,
  hubTable,
  candidateHubGenes,
  file =
    "Results/RData/Step5_GS_MM.RData"
)

###############################################################
# FINAL SUMMARY
###############################################################

cat(
  "\n=========================================\n"
)

cat(
  "SCRIPT 07 COMPLETED\n"
)

cat(
  "=========================================\n"
)

cat(
  "ComBat expression dimensions: ",
  dim(combat.expr)[1],
  " x ",
  dim(combat.expr)[2],
  "\n"
)

cat(
  "WGCNA expression dimensions: ",
  nrow(datExpr),
  " samples x ",
  ncol(datExpr),
  " genes\n"
)

cat(
  "Selected soft power: ",
  softPower,
  "\n"
)

cat(
  "Selected module: ",
  module,
  "\n"
)

cat(
  "Genes in selected module: ",
  sum(moduleGenes),
  "\n"
)

cat(
  "Candidate hub genes: ",
  nrow(candidateHubGenes),
  "\n"
)

cat(
  "=========================================\n"
)

###############################################################
# Save complete Script 07 workspace
###############################################################

save(
  combat.expr,
  group2,
  batch2,
  datExpr,
  sampleTree,
  sft,
  softPower,
  net,
  moduleColors,
  MEs,
  traitData,
  moduleTraitCor,
  moduleTraitPvalue,
  moduleSummary,
  geneTraitSignificance,
  GSPvalue,
  moduleMembership,
  MMPvalue,
  hubTable,
  candidateHubGenes,
  file =
    "RData/Script7_WGCNA_ComBat_Merged_Workspace.RData"
)

###############################################################
# Save console history as requested
###############################################################

savehistory(
  "Scripts/07_WGCNA_ComBat_Merged.R"
)

###############################################################
# END OF SCRIPT 07
###############################################################






















