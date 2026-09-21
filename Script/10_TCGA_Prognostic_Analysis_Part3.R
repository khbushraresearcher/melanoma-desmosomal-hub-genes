


















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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###############################################################
  ## SCRIPT 15: GSEA ANALYSIS
  ##
  ## GSEA performed for:
  ## 1. GSE3189
  ## 2. GSE7553
  ## 3. GSE46517
  ## 4. ComBat-corrected merged dataset
  ##
  ## Gene Ontology Biological Process
  ## Hallmark gene sets
  ###############################################################
  
  ###############################################################
  # Set project directory
  ###############################################################
  
  setwd("~/Desktop/Melanoma_Project2")
  
  ###############################################################
  # Load packages
  ###############################################################
  
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
  
  ###############################################################
  # Create output folders
  ###############################################################
  
  dir.create(
    "Results",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    "Results/Tables",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    "Results/Figures",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    "Results/RData",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ###############################################################
  # Load DEG workspace
  ###############################################################
  
  load(
    "RData/Script2_DEG_Workspace.RData"
  )
  
  ###############################################################
  # Load ComBat dataset
  ###############################################################
  
  load(
    "RData/ComBat_Filtered_Merged.RData"
  )
  
  ###############################################################
  # Check objects
  ###############################################################
  
  print(
    ls()
  )
  
  ###############################################################
  # Load Hallmark gene sets
  ###############################################################
  
  hallmark <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  )
  
  ###############################################################
  # Create Hallmark TERM2GENE table
  ###############################################################
  
  hallmark_t2g <- hallmark[
    ,
    c(
      "gs_name",
      "gene_symbol"
    )
  ]
  
  ###############################################################
  # Function for ranked GSEA
  ###############################################################
  
  run_GSEA <- function(
    deg_table,
    dataset_name
  ) {
    
    #############################################################
    # Check required columns
    #############################################################
    
    print(
      colnames(deg_table)
    )
    
    #############################################################
    # Remove missing statistics
    #############################################################
    
    deg_table <- deg_table[
      !is.na(deg_table$t),
    ]
    
    #############################################################
    # Remove duplicated genes
    #############################################################
    
    deg_table <- deg_table[
      !duplicated(
        rownames(deg_table)
      ),
    ]
    
    #############################################################
    # Create ranked gene list
    #############################################################
    
    gene_ranks <- deg_table$t
    
    names(gene_ranks) <- rownames(
      deg_table
    )
    
    #############################################################
    # Sort decreasing
    #############################################################
    
    gene_ranks <- sort(
      gene_ranks,
      decreasing = TRUE
    )
    
    #############################################################
    # Remove duplicated names
    #############################################################
    
    gene_ranks <- gene_ranks[
      !duplicated(
        names(gene_ranks)
      )
    ]
    
    #############################################################
    # GO Biological Process GSEA
    #############################################################
    
    gsea_GO <- gseGO(
      
      geneList = gene_ranks,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "SYMBOL",
      
      ont = "BP",
      
      minGSSize = 10,
      
      maxGSSize = 500,
      
      pvalueCutoff = 0.25,
      
      verbose = FALSE
    )
    
    #############################################################
    # Hallmark GSEA
    #############################################################
    
    gsea_Hallmark <- GSEA(
      
      geneList = gene_ranks,
      
      TERM2GENE = hallmark_t2g,
      
      minGSSize = 10,
      
      maxGSSize = 500,
      
      pvalueCutoff = 0.25,
      
      verbose = FALSE
    )
    
    #############################################################
    # Save GO results
    #############################################################
    
    if (!is.null(gsea_GO)) {
      
      write.csv(
        
        as.data.frame(gsea_GO),
        
        paste0(
          "Results/Tables/",
          dataset_name,
          "_GSEA_GO_BP.csv"
        ),
        
        row.names = FALSE
      )
    }
    
    #############################################################
    # Save Hallmark results
    #############################################################
    
    if (!is.null(gsea_Hallmark)) {
      
      write.csv(
        
        as.data.frame(gsea_Hallmark),
        
        paste0(
          "Results/Tables/",
          dataset_name,
          "_GSEA_Hallmark.csv"
        ),
        
        row.names = FALSE
      )
    }
    
    #############################################################
    # GO dotplot
    #############################################################
    
    if (
      !is.null(gsea_GO) &&
      nrow(as.data.frame(gsea_GO)) > 0
    ) {
      
      p_GO <- dotplot(
        gsea_GO,
        showCategory = 15
      ) +
        ggtitle(
          paste(
            dataset_name,
            "GSEA: GO Biological Process"
          )
        )
      
      ggsave(
        
        paste0(
          "Results/Figures/",
          dataset_name,
          "_GSEA_GO_BP.png"
        ),
        
        plot = p_GO,
        
        width = 10,
        
        height = 8,
        
        dpi = 300
      )
    }
    
    #############################################################
    # Hallmark dotplot
    #############################################################
    
    if (
      !is.null(gsea_Hallmark) &&
      nrow(as.data.frame(gsea_Hallmark)) > 0
    ) {
      
      p_Hallmark <- dotplot(
        gsea_Hallmark,
        showCategory = 15
      ) +
        ggtitle(
          paste(
            dataset_name,
            "GSEA: Hallmark Pathways"
          )
        )
      
      ggsave(
        
        paste0(
          "Results/Figures/",
          dataset_name,
          "_GSEA_Hallmark.png"
        ),
        
        plot = p_Hallmark,
        
        width = 10,
        
        height = 8,
        
        dpi = 300
      )
    }
    
    #############################################################
    # Return results
    #############################################################
    
    return(
      list(
        ranks = gene_ranks,
        GO = gsea_GO,
        Hallmark = gsea_Hallmark
      )
    )
  }
  
  ###############################################################
  # GSEA 1: GSE3189
  ###############################################################
  
  cat(
    "\n=========================================\n"
  )
  
  cat(
    "GSEA: GSE3189\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  ###############################################################
  # Reconstruct limma model
  ###############################################################
  
  design3189 <- model.matrix(
    ~0 + group3189.DEG
  )
  
  colnames(design3189) <- levels(
    group3189.DEG
  )
  
  contrast3189 <- makeContrasts(
    Melanoma - Normal,
    levels = design3189
  )
  
  fit3189 <- lmFit(
    expr3189.DEG,
    design3189
  )
  
  fit3189 <- contrasts.fit(
    fit3189,
    contrast3189
  )
  
  fit3189 <- eBayes(
    fit3189
  )
  
  deg_GSEA_3189 <- topTable(
    fit3189,
    number = Inf,
    sort.by = "none"
  )
  
  ###############################################################
  # Run GSEA
  ###############################################################
  
  GSEA_3189 <- run_GSEA(
    deg_GSEA_3189,
    "GSE3189"
  )
  
  ###############################################################
  # GSEA 2: GSE46517
  ###############################################################
  
  cat(
    "\n=========================================\n"
  )
  
  cat(
    "GSEA: GSE46517\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  ###############################################################
  # Reconstruct limma model
  ###############################################################
  
  design46517 <- model.matrix(
    ~0 + group46517.DEG
  )
  
  colnames(design46517) <- levels(
    group46517.DEG
  )
  
  contrast46517 <- makeContrasts(
    Melanoma - Normal,
    levels = design46517
  )
  
  fit46517 <- lmFit(
    expr46517.DEG,
    design46517
  )
  
  fit46517 <- contrasts.fit(
    fit46517,
    contrast46517
  )
  
  fit46517 <- eBayes(
    fit46517
  )
  
  deg_GSEA_46517 <- topTable(
    fit46517,
    number = Inf,
    sort.by = "none"
  )
  
  ###############################################################
  # Run GSEA
  ###############################################################
  
  GSEA_46517 <- run_GSEA(
    deg_GSEA_46517,
    "GSE46517"
  )
  
  ###############################################################
  # GSEA 3: GSE7553
  ###############################################################
  
  cat(
    "\n=========================================\n"
  )
  
  cat(
    "GSEA: GSE7553\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  ###############################################################
  # Reconstruct limma model
  ###############################################################
  
  design7553 <- model.matrix(
    ~0 + group7553.DEG
  )
  
  colnames(design7553) <- levels(
    group7553.DEG
  )
  
  contrast7553 <- makeContrasts(
    Melanoma - Normal,
    levels = design7553
  )
  
  fit7553 <- lmFit(
    expr7553.DEG,
    design7553
  )
  
  fit7553 <- contrasts.fit(
    fit7553,
    contrast7553
  )
  
  fit7553 <- eBayes(
    fit7553
  )
  
  deg_GSEA_7553 <- topTable(
    fit7553,
    number = Inf,
    sort.by = "none"
  )
  
  ###############################################################
  # Run GSEA
  ###############################################################
  
  GSEA_7553 <- run_GSEA(
    deg_GSEA_7553,
    "GSE7553"
  )
  
  ###############################################################
  # GSEA 4: COMBAT-MERGED DATASET
  ###############################################################
  
  cat(
    "\n=========================================\n"
  )
  
  cat(
    "GSEA: COMBAT-MERGED DATASET\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  ###############################################################
  # ComBat matrix
  ###############################################################
  
  print(
    dim(combat.expr)
  )
  
  ###############################################################
  # Design matrix
  ###############################################################
  
  design_combat <- model.matrix(
    ~0 + group2
  )
  
  colnames(design_combat) <- levels(
    factor(group2)
  )
  
  ###############################################################
  # Melanoma vs Normal contrast
  ###############################################################
  
  contrast_combat <- makeContrasts(
    Melanoma - Normal,
    levels = design_combat
  )
  
  ###############################################################
  # Limma model
  ###############################################################
  
  fit_combat <- lmFit(
    combat.expr,
    design_combat
  )
  
  fit_combat <- contrasts.fit(
    fit_combat,
    contrast_combat
  )
  
  fit_combat <- eBayes(
    fit_combat
  )
  
  ###############################################################
  # Complete ranked DEG table
  ###############################################################
  
  deg_GSEA_combat <- topTable(
    fit_combat,
    number = Inf,
    sort.by = "none"
  )
  
  ###############################################################
  # Run GSEA
  ###############################################################
  
  GSEA_ComBat <- run_GSEA(
    deg_GSEA_combat,
    "ComBat_Merged"
  )
  
  ###############################################################
  # Save limma DEG table used for GSEA
  ###############################################################
  
  write.csv(
    
    deg_GSEA_combat,
    
    "Results/Tables/ComBat_Merged_GSEA_Ranking.csv",
    
    row.names = TRUE
  )
  
  ###############################################################
  # Save individual GSEA ranking tables
  ###############################################################
  
  write.csv(
    
    deg_GSEA_3189,
    
    "Results/Tables/GSE3189_GSEA_Ranking.csv",
    
    row.names = TRUE
  )
  
  write.csv(
    
    deg_GSEA_46517,
    
    "Results/Tables/GSE46517_GSEA_Ranking.csv",
    
    row.names = TRUE
  )
  
  write.csv(
    
    deg_GSEA_7553,
    
    "Results/Tables/GSE7553_GSEA_Ranking.csv",
    
    row.names = TRUE
  )
  
  ###############################################################
  # Save all GSEA objects
  ###############################################################
  
  save(
    
    GSEA_3189,
    
    GSEA_46517,
    
    GSEA_7553,
    
    GSEA_ComBat,
    
    deg_GSEA_3189,
    
    deg_GSEA_46517,
    
    deg_GSEA_7553,
    
    deg_GSEA_combat,
    
    file =
      "RData/Script15_GSEA_Workspace.RData"
  )
  
  ###############################################################
  # Final summary
  ###############################################################
  
  cat(
    "\n=========================================\n"
  )
  
  cat(
    "SCRIPT 15 GSEA COMPLETED\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  cat(
    "GSE3189: GSEA completed\n"
  )
  
  cat(
    "GSE46517: GSEA completed\n"
  )
  
  cat(
    "GSE7553: GSEA completed\n"
  )
  
  cat(
    "ComBat merged: GSEA completed\n"
  )
  
  cat(
    "=========================================\n"
  )
  
  ###############################################################
  # Save console history
  ###############################################################
  
  savehistory(
    "Scripts/15_GSEA_Analysis.R"
  )
  
  ###############################################################
  # END OF SCRIPT 15
  ###############################################################
  
  
  
  
  
  
  
  
  
  
  list.files(
    "Results/Tables",
    pattern = "GSEA",
    full.names = TRUE
  )
  
  combat_hallmark <- read.csv(
    "Results/Tables/ComBat_Merged_GSEA_Hallmark.csv"
  )
  
  View(combat_hallmark)
  
  
  g3189_H <- read.csv("Results/Tables/GSE3189_GSEA_Hallmark.csv")
  
  colnames(g3189_H)
  
  head(
    g3189_H[
      order(g3189_H$p.adjust),
      c("Description", "NES", "pvalue", "p.adjust", "qvalue")
    ],
    20
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############################################################
  # SCRIPT 10 – PART 3
  # Time-Dependent ROC Analysis
  # TCGA-SKCM Prognostic Analysis
  ############################################################
  
  setwd("~/Desktop/Melanoma_Project2")
  
  ############################################################
  # 1. LOAD REQUIRED PACKAGES
  ############################################################
  
  library(survival)
  library(timeROC)
  library(ggplot2)
  
  ############################################################
  # 2. CHECK REQUIRED OBJECTS
  ############################################################
  
  # analysis_data_surv should contain:
  # OS_time
  # OS_event
  # RiskScore
  
  head(analysis_data_surv[, c("OS_time", "OS_event", "RiskScore")])
  
  summary(analysis_data_surv$OS_time)
  table(analysis_data_surv$OS_event)
  summary(analysis_data_surv$RiskScore)
  
  ############################################################
  # 3. DEFINE TIME POINTS
  ############################################################
  
  times <- c(365, 1095, 1825)
  
  time_labels <- c(
    "1-year",
    "3-year",
    "5-year"
  )
  
  ############################################################
  # 4. TIME-DEPENDENT ROC ANALYSIS
  ############################################################
  
  roc_td <- timeROC(
    T = analysis_data_surv$OS_time,
    delta = analysis_data_surv$OS_event,
    marker = analysis_data_surv$RiskScore,
    cause = 1,
    weighting = "marginal",
    times = times,
    iid = TRUE
  )
  
  ############################################################
  # 5. AUC VALUES
  ############################################################
  
  auc_values <- data.frame(
    Time = time_labels,
    Days = times,
    AUC = as.numeric(roc_td$AUC)
  )
  
  print(auc_values)
  
  ############################################################
  # 6. AUC STANDARD ERRORS AND 95% CI
  ############################################################
  
  auc_se <- as.numeric(
    roc_td$inference$vect_sd_1
  )
  
  auc_ci <- data.frame(
    Time = time_labels,
    Days = times,
    AUC = as.numeric(roc_td$AUC),
    SE = auc_se
  )
  
  auc_ci$Lower_95CI <- auc_ci$AUC -
    1.96 * auc_ci$SE
  
  auc_ci$Upper_95CI <- auc_ci$AUC +
    1.96 * auc_ci$SE
  
  # Restrict CI to possible AUC range
  auc_ci$Lower_95CI <- pmax(
    auc_ci$Lower_95CI,
    0
  )
  
  auc_ci$Upper_95CI <- pmin(
    auc_ci$Upper_95CI,
    1
  )
  
  print(auc_ci)
  
  ############################################################
  # 7. SAVE AUC RESULTS
  ############################################################
  
  write.csv(
    auc_values,
    "Results/Prognostic/TCGA_SKCM_TimeDependent_ROC_AUC.csv",
    row.names = FALSE
  )
  
  write.csv(
    auc_ci,
    "Results/Prognostic/TCGA_SKCM_TimeDependent_ROC_AUC_95CI.csv",
    row.names = FALSE
  )
  
  ############################################################
  # 8. ROC CURVES – 1, 3 AND 5 YEARS
  ############################################################
  
  png(
    "Figures/TCGA_SKCM_TimeDependent_ROC.png",
    width = 2400,
    height = 2000,
    res = 300
  )
  
  plot(
    roc_td,
    time = 365,
    col = 1,
    lwd = 2,
    title = FALSE
  )
  
  plot(
    roc_td,
    time = 1095,
    add = TRUE,
    col = 2,
    lwd = 2
  )
  
  plot(
    roc_td,
    time = 1825,
    add = TRUE,
    col = 3,
    lwd = 2
  )
  
  abline(
    a = 0,
    b = 1,
    lty = 2
  )
  
  legend(
    "bottomright",
    legend = c(
      paste0(
        "1-year AUC = ",
        round(auc_values$AUC[1], 3)
      ),
      paste0(
        "3-year AUC = ",
        round(auc_values$AUC[2], 3)
      ),
      paste0(
        "5-year AUC = ",
        round(auc_values$AUC[3], 3)
      )
    ),
    col = c(1, 2, 3),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
  
  ############################################################
  # 9. SAVE ROC OBJECTS
  ############################################################
  
  save(
    roc_td,
    auc_values,
    auc_ci,
    times,
    file = "RData/Script10_TCGA_TimeDependent_ROC.RData"
  )
  
  ############################################################
  # 10. SAVE SCRIPT HISTORY
  ############################################################
  
  savehistory(
    "Scripts/10_TCGA_Prognostic_Analysis_Part3.R"
  )
  
  ############################################################
  # END OF PART 3
  ############################################################
  
