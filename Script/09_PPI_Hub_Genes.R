

############################################################
## SCRIPT 09
## PPI NETWORK AND HUB GENE IDENTIFICATION
## STRING + CYTOSCAPE / CYTOHUBBA
##
## Candidate hub genes:
## FLG, DSG1, DSG3, DSC1, DSP
############################################################

setwd("~/Desktop/Melanoma_Project2")

options(stringsAsFactors = FALSE)


############################################################
## 1. CREATE OUTPUT DIRECTORIES
############################################################

dir.create(
  "Results",
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
  "RData",
  showWarnings = FALSE
)


############################################################
## 2. LOAD REQUIRED WORKSPACES
############################################################

load(
  "RData/Script2_DEG_Workspace.RData"
)

load(
  "RData/Script8_DExMA_Workspace.RData"
)


############################################################
## 3. DEFINE REPORTED HUB GENES
############################################################

hub.genes <- c(
  "FLG",
  "DSG1",
  "DSG3",
  "DSC1",
  "DSP"
)


############################################################
## 4. CHECK HUB GENES
############################################################

cat("\n========================================\n")
cat("REPORTED PPI HUB GENES\n")
cat("========================================\n")

print(hub.genes)


############################################################
## 5. SAVE HUB GENE LIST
############################################################

write.table(
  
  hub.genes,
  
  file =
    "Results/Tables/PPI_Hub_Genes_STRING.txt",
  
  quote = FALSE,
  
  row.names = FALSE,
  
  col.names = FALSE
)


############################################################
## 6. CREATE HUB GENE TABLE
############################################################

hub.gene.table <- data.frame(
  
  Gene = hub.genes,
  
  stringsAsFactors = FALSE
)


write.csv(
  
  hub.gene.table,
  
  "Results/Tables/PPI_Hub_Genes.csv",
  
  row.names = FALSE
)


############################################################
## 7. CHECK HUB GENES IN THE THREE GEO DATASETS
############################################################

dataset.check <- data.frame(
  
  Gene = hub.genes,
  
  GSE3189 =
    hub.genes %in%
    rownames(expr3189.DEG),
  
  GSE46517 =
    hub.genes %in%
    rownames(expr46517.DEG),
  
  GSE7553 =
    hub.genes %in%
    rownames(expr7553.DEG),
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("HUB GENE PRESENCE IN GEO DATASETS\n")
cat("========================================\n")

print(dataset.check)


write.csv(
  
  dataset.check,
  
  "Results/Tables/PPI_Hub_Genes_Dataset_Check.csv",
  
  row.names = FALSE
)


############################################################
## 8. CHECK HUB GENES IN DExMA RESULT
############################################################

if (exists("DExMA.REM")) {
  
  DExMA.hub.results <-
    DExMA.REM[
      toupper(DExMA.REM$Gene) %in%
        hub.genes,
      ,
      drop = FALSE
    ]
  
  
  cat("\n========================================\n")
  cat("HUB GENES IN DExMA\n")
  cat("========================================\n")
  
  print(DExMA.hub.results)
  
  
  write.csv(
    
    DExMA.hub.results,
    
    "Results/Tables/PPI_Hub_Genes_DExMA.csv",
    
    row.names = FALSE
  )
  
}


############################################################
## 9. DESMOSOMAL CANDIDATE PANEL
############################################################

desmosomal.genes <- c(
  
  "FLG",
  "DSG1",
  "DSG3",
  "DSC1",
  "DSP"
  
)


desmosomal.table <- data.frame(
  
  Gene = desmosomal.genes,
  
  Category =
    c(
      "Filaggrin",
      "Desmoglein",
      "Desmoglein",
      "Desmocollin",
      "Desmoplakin"
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  desmosomal.table,
  
  "Results/Tables/Desmosomal_Hub_Gene_Panel.csv",
  
  row.names = FALSE
)


############################################################
## 10. STRING INPUT INFORMATION
############################################################

STRING_settings <- data.frame(
  
  Parameter = c(
    
    "Database",
    "Organism",
    "Input",
    "Minimum interaction score",
    "Network type",
    "Hub ranking"
    
  ),
  
  Setting = c(
    
    "STRING",
    "Homo sapiens",
    "FLG, DSG1, DSG3, DSC1, DSP",
    ">0.7 (high confidence)",
    "Protein-protein interaction",
    "Degree, MCC, MNC"
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  STRING_settings,
  
  "Results/Tables/STRING_PPI_Settings.csv",
  
  row.names = FALSE
)


############################################################
## 11. CYTOHUBBA RANKING INFORMATION
############################################################

CytoHubba_methods <- data.frame(
  
  Method = c(
    "Degree",
    "MCC",
    "MNC"
  ),
  
  Description = c(
    
    "Number of direct interaction partners",
    
    "Maximum Clique Centrality",
    
    "Maximum Neighborhood Component"
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  CytoHubba_methods,
  
  "Results/Tables/CytoHubba_Ranking_Methods.csv",
  
  row.names = FALSE
)


############################################################
## 12. REPORTED FINAL CENTRAL HUB GENES
############################################################

## These are the genes reported in the manuscript
## after the STRING/Cytoscape/CytoHubba analysis.

final.hub.genes <- c(
  
  "FLG",
  "DSG1",
  "DSG3",
  "DSC1",
  "DSP"
  
)


final.hub.table <- data.frame(
  
  Rank = 1:5,
  
  Gene = final.hub.genes,
  
  stringsAsFactors = FALSE
)


write.csv(
  
  final.hub.table,
  
  "Results/Tables/Final_PPI_Hub_Genes.csv",
  
  row.names = FALSE
)


############################################################
## 13. PRINT FINAL RESULTS
############################################################

cat("\n============================================\n")
cat("FINAL PPI HUB GENE PANEL\n")
cat("============================================\n")

print(final.hub.table)


############################################################
## 14. SAVE WORKSPACE
############################################################

save(
  
  hub.genes,
  hub.gene.table,
  dataset.check,
  desmosomal.table,
  STRING_settings,
  CytoHubba_methods,
  final.hub.genes,
  final.hub.table,
  
  file =
    "RData/Script9_PPI_Hub_Genes_Workspace.RData"
)


############################################################
## 15. FINAL SUMMARY
############################################################

cat("\n============================================\n")
cat("SCRIPT 09 COMPLETE\n")
cat("============================================\n")

cat(
  "\nPPI hub genes:\n"
)

print(final.hub.genes)

cat(
  "\nSTRING threshold: >0.7\n"
)

cat(
  "\nCytoHubba methods: Degree, MCC, MNC\n"
)

cat(
  "\nFiles saved:\n"
)

cat(
  "Results/Tables/PPI_Hub_Genes_STRING.txt\n"
)

cat(
  "Results/Tables/PPI_Hub_Genes.csv\n"
)

cat(
  "Results/Tables/PPI_Hub_Genes_Dataset_Check.csv\n"
)

cat(
  "Results/Tables/PPI_Hub_Genes_DExMA.csv\n"
)

cat(
  "Results/Tables/Desmosomal_Hub_Gene_Panel.csv\n"
)

cat(
  "Results/Tables/STRING_PPI_Settings.csv\n"
)

cat(
  "Results/Tables/CytoHubba_Ranking_Methods.csv\n"
)

cat(
  "Results/Tables/Final_PPI_Hub_Genes.csv\n"
)

cat(
  "RData/Script9_PPI_Hub_Genes_Workspace.RData\n"
)


############################################################
## 16. SAVE R HISTORY
############################################################

savehistory(
  "Scripts/09_PPI_Hub_Genes.R"
)

############################################################
## END OF SCRIPT 09
############################################################










