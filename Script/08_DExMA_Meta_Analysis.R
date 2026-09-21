
############################################################
## SCRIPT 08
## DExMA 1.20.0
## RANDOM-EFFECTS META-ANALYSIS
## GSE3189 + GSE46517 + GSE7553
## FLG + DSG1 + DSG3
############################################################

setwd("~/Desktop/Melanoma_Project2")

options(stringsAsFactors = FALSE)


############################################################
## 1. LOAD LIBRARIES
############################################################

library(DExMA)
library(ggplot2)
library(ggrepel)


############################################################
## 2. CREATE OUTPUT FOLDERS
############################################################

dir.create("Results", showWarnings = FALSE)
dir.create("Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)
dir.create("RData", showWarnings = FALSE)


############################################################
## 3. LOAD PREVIOUS WORKSPACE
############################################################

load("RData/Script2_DEG_Workspace.RData")


############################################################
## 4. CHECK INPUT OBJECTS
############################################################

required.objects <- c(
  "expr3189.DEG",
  "expr46517.DEG",
  "expr7553.DEG",
  "group3189.DEG",
  "group46517.DEG",
  "group7553.DEG"
)

missing.objects <- required.objects[
  !sapply(required.objects, exists)
]

if (length(missing.objects) > 0) {
  
  stop(
    "Missing objects: ",
    paste(missing.objects, collapse = ", ")
  )
}


############################################################
## 5. CREATE EXPRESSION LIST
############################################################

expression.list <- list(
  GSE3189  = expr3189.DEG,
  GSE46517 = expr46517.DEG,
  GSE7553  = expr7553.DEG
)


############################################################
## 6. CHECK ACTUAL GROUP NAMES
############################################################

cat("\nGSE3189 groups:\n")
print(table(group3189.DEG))

cat("\nGSE46517 groups:\n")
print(table(group46517.DEG))

cat("\nGSE7553 groups:\n")
print(table(group7553.DEG))


############################################################
## 7. CREATE PHENOTYPE LIST
############################################################

phenotype.list <- lapply(
  
  list(
    group3189.DEG,
    group46517.DEG,
    group7553.DEG
  ),
  
  function(g) {
    
    data.frame(
      
      Group = factor(
        as.character(g),
        levels = c(
          "Melanoma",
          "Normal"
        )
      ),
      
      row.names = NULL
    )
  }
)


############################################################
## 8. ASSIGN CORRECT SAMPLE NAMES
############################################################

rownames(phenotype.list[[1]]) <-
  colnames(expr3189.DEG)

rownames(phenotype.list[[2]]) <-
  colnames(expr46517.DEG)

rownames(phenotype.list[[3]]) <-
  colnames(expr7553.DEG)


############################################################
## 9. CHECK PHENOTYPE / EXPRESSION MATCHING
############################################################

for (i in seq_along(expression.list)) {
  
  if (!identical(
    rownames(phenotype.list[[i]]),
    colnames(expression.list[[i]])
  )) {
    
    stop(
      "Sample names do not match for dataset ",
      names(expression.list)[i]
    )
  }
}


cat("\n========================================\n")
cat("EXPRESSION / PHENOTYPE MATCHING PASSED\n")
cat("========================================\n")


############################################################
## 10. CREATE DExMA OBJECT
############################################################

DExMA.object <- createObjectMA(
  
  listEX = expression.list,
  
  listPheno = phenotype.list,
  
  namePheno = c(
    "Group",
    "Group",
    "Group"
  ),
  
  expGroups = c(
    1,
    1,
    1
  ),       # Melanoma
  
  refGroups = c(
    2,
    2,
    2
  )        # Normal
)


cat("\n========================================\n")
cat("DExMA OBJECT CREATED\n")
cat("========================================\n")


############################################################
## 11. CALCULATE EFFECT SIZES
############################################################

DExMA.ES <- calculateES(
  objectMA = DExMA.object
)


cat("\n========================================\n")
cat("DExMA EFFECT SIZES CALCULATED\n")
cat("========================================\n")

print(names(DExMA.ES))


############################################################
## 12. RANDOM-EFFECTS META-ANALYSIS
############################################################

## DExMA 1.20.0 internal REM function
## This is the method used in the original analysis.

DExMA.REM <- DExMA:::.metaES(
  
  calESResults = DExMA.ES,
  
  metaMethod = "REM",
  
  proportionData = 0.5
)


cat("\n========================================\n")
cat("DExMA REM META-ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "Genes analyzed:",
  nrow(DExMA.REM),
  "\n"
)

print(colnames(DExMA.REM))


############################################################
## 13. SAVE COMPLETE DExMA RESULT
############################################################

DExMA.REM$Gene <-
  rownames(DExMA.REM)

DExMA.REM <-
  DExMA.REM[
    ,
    c(
      "Gene",
      setdiff(
        colnames(DExMA.REM),
        "Gene"
      )
    )
  ]


write.csv(
  
  DExMA.REM,
  
  "Results/Tables/DExMA_REM_AllGenes.csv",
  
  row.names = FALSE
)


############################################################
## 14. EXTRACT SIGNIFICANT GENES
############################################################

DExMA.sig <- subset(
  
  DExMA.REM,
  
  !is.na(FDR) &
    FDR < 0.05 &
    abs(Com.ES) >= 1
)


cat("\n========================================\n")
cat("DExMA SIGNIFICANT GENES\n")
cat("========================================\n")

cat(
  "Significant genes:",
  nrow(DExMA.sig),
  "\n"
)


write.csv(
  
  DExMA.sig,
  
  "Results/Tables/DExMA_REM_SignificantGenes.csv",
  
  row.names = FALSE
)


############################################################
## 15. EXTRACT FLG, DSG1 AND DSG3
############################################################

candidate.genes <- c(
  "FLG",
  "DSG1",
  "DSG3"
)


candidate.results <- DExMA.REM[
  
  toupper(
    DExMA.REM$Gene
  ) %in% candidate.genes,
  
  ,
  
  drop = FALSE
]


cat("\n========================================\n")
cat("FLG / DSG1 / DSG3 DExMA RESULTS\n")
cat("========================================\n")


print(
  
  candidate.results[
    ,
    c(
      "Gene",
      "Com.ES",
      "ES.var",
      "Qval",
      "Qpval",
      "tau2",
      "Zval",
      "Pval",
      "FDR",
      "AveFC",
      "propDataset"
    )
  ]
)


write.csv(
  
  candidate.results,
  
  "Results/Tables/DExMA_REM_FLG_DSG1_DSG3.csv",
  
  row.names = FALSE
)


############################################################
## 16. PREPARE VOLCANO-PLOT DATA
############################################################

volcano.df <- DExMA.REM


volcano.df$Meta_ES <-
  as.numeric(
    volcano.df$Com.ES
  )


volcano.df$FDR <-
  as.numeric(
    volcano.df$FDR
  )


## Protect against zero FDR

volcano.df$negLog10FDR <-
  -log10(
    pmax(
      volcano.df$FDR,
      .Machine$double.xmin
    )
  )


## Significance

volcano.df$Significance <-
  ifelse(
    
    volcano.df$FDR < 0.05 &
      abs(volcano.df$Meta_ES) >= 1,
    
    "Significant",
    
    "Not significant"
  )


volcano.df$Significance <-
  factor(
    
    volcano.df$Significance,
    
    levels = c(
      "Not significant",
      "Significant"
    )
  )


############################################################
## 17. LABEL THE THREE CANDIDATE GENES
############################################################

volcano.df$Label <-
  ifelse(
    
    toupper(
      volcano.df$Gene
    ) %in% candidate.genes,
    
    volcano.df$Gene,
    
    ""
  )


############################################################
## 18. CAP Y-AXIS FOR DISPLAY ONLY
############################################################

## The original FDR values remain unchanged.
## Only the displayed y-axis is capped.

volcano.df$plot_y <-
  pmin(
    volcano.df$negLog10FDR,
    35
  )


label.df <-
  subset(
    volcano.df,
    Label != ""
  )


label.df$label_y <-
  pmin(
    label.df$plot_y,
    30
  )


############################################################
## 19. PRINT CANDIDATE VALUES
############################################################

cat("\n========================================\n")
cat("VALUES USED IN VOLCANO PLOT\n")
cat("========================================\n")


print(
  
  label.df[
    ,
    c(
      "Gene",
      "Meta_ES",
      "FDR",
      "negLog10FDR",
      "Significance"
    )
  ]
)


############################################################
## 20. DExMA VOLCANO PLOT
############################################################

p.volcano <- ggplot(
  
  volcano.df,
  
  aes(
    x = Meta_ES,
    y = plot_y
  )
  
) +
  
  
  geom_point(
    
    aes(
      color = Significance
    ),
    
    alpha = 0.60,
    
    size = 1.8
  ) +
  
  
  geom_vline(
    
    xintercept = c(
      -1,
      1
    ),
    
    linetype = "dashed",
    
    linewidth = 0.7
  ) +
  
  
  geom_hline(
    
    yintercept = -log10(0.05),
    
    linetype = "dashed",
    
    linewidth = 0.7
  ) +
  
  
  geom_label_repel(
    
    data = label.df,
    
    aes(
      x = Meta_ES,
      y = label_y,
      label = Label
    ),
    
    size = 5,
    
    fontface = "bold",
    
    color = "black",
    
    fill = "white",
    
    box.padding = 0.8,
    
    point.padding = 0.5,
    
    force = 3,
    
    min.segment.length = 0,
    
    max.overlaps = Inf,
    
    seed = 123
  ) +
  
  
  scale_y_continuous(
    
    limits = c(
      0,
      35
    ),
    
    breaks = seq(
      0,
      35,
      5
    ),
    
    expand = expansion(
      mult = c(
        0.02,
        0.03
      )
    )
  ) +
  
  
  labs(
    
    title =
      "Random-effects meta-analysis of GSE3189, GSE46517 and GSE7553",
    
    x =
      "Meta-analysis effect size (Com.ES)",
    
    y =
      expression(
        -log[10](FDR)
      ),
    
    color =
      "Significance"
  ) +
  
  
  theme_classic(
    base_size = 16
  ) +
  
  
  theme(
    
    plot.title =
      element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
    
    axis.title =
      element_text(
        size = 17
      ),
    
    axis.text =
      element_text(
        size = 14
      ),
    
    legend.title =
      element_text(
        size = 15
      ),
    
    legend.text =
      element_text(
        size = 14
      )
  )


print(p.volcano)


############################################################
## 21. SAVE VOLCANO PLOT
############################################################

ggsave(
  
  "Results/Figures/DExMA_REM_Volcano_FLG_DSG1_DSG3.png",
  
  p.volcano,
  
  width = 12,
  
  height = 9,
  
  dpi = 600
)


############################################################
## 22. SAVE DExMA WORKSPACE
############################################################

save(
  
  expression.list,
  phenotype.list,
  DExMA.object,
  DExMA.ES,
  DExMA.REM,
  DExMA.sig,
  candidate.results,
  
  file =
    "RData/Script8_DExMA_Workspace.RData"
)


############################################################
## 23. FINAL SUMMARY
############################################################

cat("\n============================================\n")
cat("DExMA ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
  "Common/meta-analyzed genes:",
  nrow(DExMA.REM),
  "\n"
)

cat(
  "Significant genes:",
  nrow(DExMA.sig),
  "\n"
)

cat("\nCandidate genes:\n")

print(candidate.results)


cat("\nFiles saved:\n")

cat(
  "Results/Tables/DExMA_REM_AllGenes.csv\n"
)

cat(
  "Results/Tables/DExMA_REM_SignificantGenes.csv\n"
)

cat(
  "Results/Tables/DExMA_REM_FLG_DSG1_DSG3.csv\n"
)

cat(
  "Results/Figures/DExMA_REM_Volcano_FLG_DSG1_DSG3.png\n"
)

cat(
  "RData/Script8_DExMA_Workspace.RData\n"
)


############################################################
## 24. SAVE R HISTORY
############################################################

savehistory(
  "Scripts/08_DExMA_Meta_Analysis.R"
)

############################################################
## END OF SCRIPT 08
############################################################














