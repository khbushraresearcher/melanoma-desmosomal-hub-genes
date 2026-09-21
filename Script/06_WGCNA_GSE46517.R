


###############################################################
## SCRIPT 06
## WGCNA Analysis – GSE46517 Melanoma Dataset
###############################################################

setwd("~/Desktop/Melanoma_Project2")

###############################################################
## Load libraries
###############################################################

library(WGCNA)

options(stringsAsFactors = FALSE)

allowWGCNAThreads()

###############################################################
## Create output folders
###############################################################

dir.create("Results", showWarnings = FALSE)
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("RData", showWarnings = FALSE)

###############################################################
## Load workspace from Script 02
###############################################################

load("RData/Script2_DEG_Workspace.RData")

###############################################################
## Check GSE46517 data
###############################################################

print(dim(expr46517.DEG))
print(table(group46517.DEG))

###############################################################
## Prepare expression matrix
###############################################################

datExpr <- t(expr46517.DEG)

###############################################################
## Check samples and genes
###############################################################

gsg <- goodSamplesGenes(datExpr)

print(gsg)

datExpr <- datExpr[
  gsg$goodSamples,
  gsg$goodGenes
]

print(dim(datExpr))

###############################################################
## Sample clustering
###############################################################

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

png(
  "Results/Figures/GSE46517_Figure1_SampleClustering.png",
  width = 2200,
  height = 1800,
  res = 300
)

plot(
  sampleTree,
  main = "GSE46517 Sample Clustering",
  sub = "",
  xlab = ""
)

dev.off()

###############################################################
## Trait data
###############################################################

traitData <- data.frame(
  Melanoma = ifelse(
    group46517.DEG == "Melanoma",
    1,
    0
  )
)

rownames(traitData) <- rownames(datExpr)

###############################################################
## Soft-threshold selection
###############################################################

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

###############################################################
## Select soft power
###############################################################

softPower <- 5

###############################################################
## Soft-threshold plot
###############################################################

png(
  "Results/Figures/GSE46517_Figure2_SoftThreshold.png",
  width = 3200,
  height = 1600,
  res = 300
)

par(mfrow = c(1,2))

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) *
    sft$fitIndices[,2],
  type = "n",
  xlab = "Soft Threshold (Power)",
  ylab = "Scale-Free Topology Fit (Signed R²)",
  main = "Scale Independence"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) *
    sft$fitIndices[,2],
  labels = powers,
  col = "red",
  cex = 1
)

abline(
  h = 0.90,
  col = "blue",
  lwd = 2
)

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  type = "n",
  xlab = "Soft Threshold (Power)",
  ylab = "Mean Connectivity",
  main = "Mean Connectivity"
)

text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  col = "red",
  cex = 1
)

dev.off()

###############################################################
## Construct signed co-expression network
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
## Module colors
###############################################################

moduleColors <- labels2colors(
  net$colors
)

print(table(moduleColors))

###############################################################
## Module eigengenes
###############################################################

MEs <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs)

###############################################################
## Gene dendrogram
###############################################################

for(i in seq_along(net$dendrograms)){
  
  png(
    paste0(
      "Results/Figures/GSE46517_Dendrogram_Block",
      i,
      ".png"
    ),
    width = 4500,
    height = 2500,
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
## Module-Trait Relationships
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

print(moduleTraitCor)
print(moduleTraitPvalue)

###############################################################
## Module-Trait Heatmap
###############################################################

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(moduleTraitCor)

png(
  "Results/Figures/GSE46517_ModuleTraitHeatmap.png",
  width = 2200,
  height = 3200,
  res = 300
)

par(
  mar = c(8,10,4,3)
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Melanoma",
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.9,
  zlim = c(-1,1),
  main = "Module-Trait Relationships"
)

dev.off()

###############################################################
## Gene Significance
###############################################################

GS <- as.data.frame(
  cor(
    datExpr,
    traitData$Melanoma,
    use = "p"
  )
)

colnames(GS) <- "GS"

###############################################################
## Module Membership
###############################################################

MM <- as.data.frame(
  cor(
    datExpr,
    MEs,
    use = "p"
  )
)

colnames(MM) <- paste0(
  "MM.",
  colnames(MM)
)

###############################################################
## Complete Hub Gene Table
###############################################################

hubTable46517 <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  GS = GS$GS,
  MM,
  stringsAsFactors = FALSE
)

###############################################################
## CSV 1 – Complete Hub Gene Table
###############################################################

write.csv(
  hubTable46517,
  "Results/Tables/HubTable_GSE46517.csv",
  row.names = FALSE
)

###############################################################
## Candidate hub genes
###############################################################

candidateGenes <- c(
  "FLG",
  "DSG3",
  "DSG1"
)

candidateHubTable46517 <- subset(
  hubTable46517,
  Gene %in% candidateGenes
)

print(candidateHubTable46517)

###############################################################
## CSV 2 – Candidate Hub Genes
###############################################################

write.csv(
  candidateHubTable46517,
  "Results/Tables/GSE46517_Candidate_Hub_Genes.csv",
  row.names = FALSE
)

###############################################################
## GS-MM Correlation
###############################################################

modules <- setdiff(
  unique(moduleColors),
  "grey"
)

GSMM46517 <- data.frame()

for(module in modules){
  
  moduleGenes <- moduleColors == module
  
  MMcolumn <- paste0(
    "MM.ME",
    module
  )
  
  if(!MMcolumn %in% colnames(MM)){
    next
  }
  
  ct <- cor.test(
    MM[
      moduleGenes,
      MMcolumn
    ],
    GS[
      moduleGenes,
      "GS"
    ]
  )
  
  GSMM46517 <- rbind(
    GSMM46517,
    data.frame(
      Module = module,
      Correlation = unname(
        ct$estimate
      ),
      Pvalue = ct$p.value,
      Genes = sum(moduleGenes)
    )
  )
}

GSMM46517 <- GSMM46517[
  order(
    -abs(
      GSMM46517$Correlation
    )
  ),
]

print(GSMM46517)

###############################################################
## CSV 3 – GS-MM Correlation
###############################################################

write.csv(
  GSMM46517,
  "Results/Tables/GSMM_Correlation_GSE46517.csv",
  row.names = FALSE
)

###############################################################
## GS vs MM plots
###############################################################

for(module in modules){
  
  moduleGenes <- moduleColors == module
  
  MMcolumn <- paste0(
    "MM.ME",
    module
  )
  
  if(!MMcolumn %in% colnames(MM)){
    next
  }
  
  r <- cor(
    MM[
      moduleGenes,
      MMcolumn
    ],
    GS[
      moduleGenes,
      "GS"
    ],
    use = "p"
  )
  
  p <- cor.test(
    MM[
      moduleGenes,
      MMcolumn
    ],
    GS[
      moduleGenes,
      "GS"
    ]
  )$p.value
  
  png(
    paste0(
      "Results/Figures/GSE46517_GSvsMM_",
      module,
      ".png"
    ),
    width = 2200,
    height = 2200,
    res = 300
  )
  
  verboseScatterplot(
    MM[
      moduleGenes,
      MMcolumn
    ],
    GS[
      moduleGenes,
      "GS"
    ],
    xlab = paste(
      "Module Membership (",
      module,
      ")",
      sep = ""
    ),
    ylab = "Gene Significance",
    main = paste0(
      tools::toTitleCase(module),
      " Module\nr = ",
      round(r, 3),
      ", P = ",
      signif(p, 3)
    ),
    col = module,
    pch = 19,
    cex = 0.8
  )
  
  dev.off()
}

###############################################################
## Significant DEGs
###############################################################

DEG46517 <- subset(
  deg46517,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

###############################################################
## Significant modules
###############################################################

moduleOrder <- names(
  sort(
    abs(moduleTraitCor[, "Melanoma"]),
    decreasing = TRUE
  )
)

###############################################################
## Module summary
###############################################################

ModuleSummary46517 <- data.frame()

for(moduleME in moduleOrder){
  
  module <- sub(
    "^ME",
    "",
    moduleME
  )
  
  moduleGenes <- hubTable46517$Gene[
    hubTable46517$Module == module
  ]
  
  moduleDEGs <- DEG46517[
    rownames(DEG46517) %in% moduleGenes,
  ]
  
  ModuleSummary46517 <- rbind(
    ModuleSummary46517,
    data.frame(
      Dataset = "GSE46517",
      Module = module,
      Trait = "Melanoma",
      Correlation = round(
        moduleTraitCor[
          moduleME,
          "Melanoma"
        ],
        3
      ),
      P_value = signif(
        moduleTraitPvalue[
          moduleME,
          "Melanoma"
        ],
        3
      ),
      Genes_in_Module = length(
        moduleGenes
      ),
      DEGs_in_Module = nrow(
        moduleDEGs
      ),
      Upregulated = sum(
        moduleDEGs$logFC > 1
      ),
      Downregulated = sum(
        moduleDEGs$logFC < -1
      )
    )
  )
}

print(ModuleSummary46517)

###############################################################
## CSV 4 – Significant Module Summary
###############################################################

write.csv(
  ModuleSummary46517,
  "Results/Tables/GSE46517_Significant_Modules.csv",
  row.names = FALSE
)

###############################################################
## Save WGCNA workspace
###############################################################

save(
  net,
  moduleColors,
  MEs,
  GS,
  MM,
  hubTable46517,
  candidateHubTable46517,
  GSMM46517,
  moduleTraitCor,
  moduleTraitPvalue,
  datExpr,
  traitData,
  ModuleSummary46517,
  file = "RData/GSE46517_WGCNA.RData"
)

###############################################################
## Save Script 06 as .R
###############################################################

savehistory(
  "Scripts/06_WGCNA_GSE46517.R"
)

###############################################################
## END OF SCRIPT 06
###############################################################


