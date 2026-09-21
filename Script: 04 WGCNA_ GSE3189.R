
###############################################################
## SCRIPT 04
## WGCNA Analysis – GSE3189 Melanoma Dataset
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
## Prepare expression matrix
###############################################################

datExpr <- t(expr3189.DEG)

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
  "Results/Figures/GSE3189_Figure1_SampleClustering.png",
  width = 2200,
  height = 1800,
  res = 300
)

plot(
  sampleTree,
  main = "Sample Clustering",
  sub = "",
  xlab = ""
)

dev.off()

###############################################################
## Trait data
###############################################################

traitData <- data.frame(
  Melanoma = ifelse(
    group3189.DEG == "Melanoma",
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
## Original GSE3189 analysis
## Soft power = 5
###############################################################

softPower <- 5

###############################################################
## Soft-threshold plot
###############################################################

png(
  "Results/Figures/GSE3189_Figure2_SoftThreshold.png",
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

print(
  sort(
    table(moduleColors),
    decreasing = TRUE
  )
)

###############################################################
## Module eigengenes
###############################################################

MEs <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs)

###############################################################
## Gene dendrograms
###############################################################

for(i in seq_along(net$dendrograms)){
  
  png(
    paste0(
      "Results/Figures/GSE3189_Dendrogram_Block",
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
  "Results/Figures/GSE3189_ModuleTraitHeatmap.png",
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
## Hub Gene Table
###############################################################

hubTable3189 <- data.frame(
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
  hubTable3189,
  "Results/Tables/HubTable_GSE3189.csv",
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

candidateHubTable3189 <- subset(
  hubTable3189,
  Gene %in% candidateGenes
)

print(candidateHubTable3189)

###############################################################
## CSV 2 – Candidate Hub Genes
###############################################################

write.csv(
  candidateHubTable3189,
  "Results/Tables/GSE3189_Candidate_Hub_Genes.csv",
  row.names = FALSE
)

###############################################################
## GS-MM Correlation
###############################################################

modules <- setdiff(
  unique(moduleColors),
  "grey"
)

GSMM3189 <- data.frame()

for(module in modules){
  
  moduleGenes <- moduleColors == module
  
  MMcolumn <- paste0(
    "MM.ME",
    module
  )
  
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
  
  GSMM3189 <- rbind(
    GSMM3189,
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

GSMM3189 <- GSMM3189[
  order(
    -abs(
      GSMM3189$Correlation
    )
  ),
]

print(GSMM3189)

###############################################################
## CSV 3 – GS-MM Correlation
###############################################################

write.csv(
  GSMM3189,
  "Results/Tables/GSMM_Correlation_GSE3189.csv",
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
      "Results/Figures/GSE3189_GSvsMM_",
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

DEG3189 <- subset(
  deg3189,
  adj.P.Val < 0.05 &
    abs(logFC) >= 1
)

###############################################################
## Significant modules from original analysis
###############################################################

sigModules <- c(
  "turquoise",
  "brown"
)

ModuleSummary3189 <- data.frame()

for(module in sigModules){
  
  moduleGenes <- hubTable3189$Gene[
    hubTable3189$Module == module
  ]
  
  moduleDEGs <- DEG3189[
    rownames(DEG3189) %in% moduleGenes,
  ]
  
  ModuleSummary3189 <- rbind(
    ModuleSummary3189,
    data.frame(
      Dataset = "GSE3189",
      Module = module,
      Trait = "Melanoma",
      Correlation = round(
        moduleTraitCor[
          paste0("ME", module),
          "Melanoma"
        ],
        3
      ),
      P_value = signif(
        moduleTraitPvalue[
          paste0("ME", module),
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

print(ModuleSummary3189)

###############################################################
## CSV 4 – Significant Module Summary
###############################################################

write.csv(
  ModuleSummary3189,
  "Results/Tables/GSE3189_Significant_Modules.csv",
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
  hubTable3189,
  candidateHubTable3189,
  GSMM3189,
  moduleTraitCor,
  moduleTraitPvalue,
  datExpr,
  traitData,
  ModuleSummary3189,
  file = "RData/GSE3189_WGCNA.RData"
)

###############################################################
## Save Script 04 as .R
###############################################################

savehistory(
  "Scripts/04_WGCNA_GSE3189.R"
)

###############################################################
## END OF SCRIPT 04
###############################################################



