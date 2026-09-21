############################################################
# SCRIPT 10 – TCGA-SKCM PROGNOSTIC ANALYSIS
# PART 5 – BOOTSTRAP LASSO VALIDATION
# AND GENE-SELECTION FREQUENCY
############################################################

############################################################
# 1. SET PROJECT DIRECTORY
############################################################

setwd("~/Desktop/Melanoma_Project2")

############################################################
# 2. LOAD REQUIRED PACKAGES
############################################################

library(survival)
library(glmnet)

############################################################
# 3. LOAD TCGA PROGNOSTIC OBJECTS
############################################################

load(
  "RData/Script10_TCGA_Prognostic_Intermediate.RData"
)

############################################################
# 4. DEFINE CANDIDATE GENES
############################################################

candidate_genes <- c(
  "FLG",
  "DSG3",
  "DSG1",
  "DSP",
  "DSC1"
)

############################################################
# 5. CHECK CANDIDATE GENES
############################################################

missing_genes <- setdiff(
  candidate_genes,
  colnames(cox_multi_complete)
)

cat(
  "\nMissing candidate genes:\n"
)

print(
  missing_genes
)

if (
  length(missing_genes) > 0
) {
  
  stop(
    "One or more candidate genes are missing from cox_multi_complete."
  )
  
}

############################################################
# 6. PREPARE LASSO DATA
############################################################

lasso_data <- cox_multi_complete[
  ,
  c(
    "OS_time",
    "OS_event",
    candidate_genes
  ),
  drop = FALSE
]

############################################################
# 7. REMOVE INCOMPLETE CASES
############################################################

lasso_data <- lasso_data[
  complete.cases(
    lasso_data
  ),
  ,
  drop = FALSE
]

############################################################
# 8. CREATE MATRIX AND SURVIVAL OBJECT
############################################################

x_lasso <- as.matrix(
  lasso_data[
    ,
    candidate_genes,
    drop = FALSE
  ]
)

y_lasso <- Surv(
  lasso_data$OS_time,
  lasso_data$OS_event
)

############################################################
# 9. CHECK DATA
############################################################

cat(
  "\nLASSO validation samples:",
  nrow(lasso_data),
  "\n"
)

cat(
  "LASSO validation events:",
  sum(
    lasso_data$OS_event == 1
  ),
  "\n"
)

cat(
  "\nCandidate genes:\n"
)

print(
  candidate_genes
)

############################################################
# 10. FIT ORIGINAL LASSO
############################################################

set.seed(1234)

original_lasso <- cv.glmnet(
  x = x_lasso,
  y = y_lasso,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance",
  standardize = TRUE
)

############################################################
# 11. EXTRACT ORIGINAL LASSO COEFFICIENTS
############################################################

original_coef <- coef(
  original_lasso,
  s = "lambda.min"
)

original_coef_matrix <- as.matrix(
  original_coef
)

original_coef_table <- data.frame(
  Gene = rownames(
    original_coef_matrix
  ),
  Coefficient = as.numeric(
    original_coef_matrix[, 1]
  )
)

############################################################
# 12. ORIGINAL SELECTED GENES
############################################################

original_selected <- original_coef_table[
  original_coef_table$Coefficient != 0,
  ,
  drop = FALSE
]

cat(
  "\nOriginal LASSO-selected genes:\n"
)

print(
  original_selected
)

############################################################
# 13. BOOTSTRAP SETTINGS
############################################################

set.seed(1234)

n_boot <- 1000

############################################################
# 14. INITIALIZE STORAGE
############################################################

selection_matrix <- matrix(
  0,
  nrow = n_boot,
  ncol = length(
    candidate_genes
  )
)

colnames(
  selection_matrix
) <- candidate_genes

coefficient_matrix <- matrix(
  NA_real_,
  nrow = n_boot,
  ncol = length(
    candidate_genes
  )
)

colnames(
  coefficient_matrix
) <- candidate_genes

lambda_min_values <- rep(
  NA_real_,
  n_boot
)

############################################################
# 15. BOOTSTRAP LASSO
############################################################

for (
  i in seq_len(n_boot)
) {
  
  ##########################################################
  # Bootstrap sample
  ##########################################################
  
  boot_index <- sample(
    seq_len(
      nrow(lasso_data)
    ),
    size = nrow(
      lasso_data
    ),
    replace = TRUE
  )
  
  x_boot <- x_lasso[
    boot_index,
    ,
    drop = FALSE
  ]
  
  y_boot <- y_lasso[
    boot_index
  ]
  
  ##########################################################
  # Check event variation
  ##########################################################
  
  event_boot <- lasso_data$OS_event[
    boot_index
  ]
  
  if (
    length(
      unique(
        event_boot
      )
    ) < 2
  ) {
    
    next
    
  }
  
  ##########################################################
  # Bootstrap LASSO
  ##########################################################
  
  boot_lasso <- try(
    cv.glmnet(
      x = x_boot,
      y = y_boot,
      family = "cox",
      alpha = 1,
      nfolds = 10,
      type.measure = "deviance",
      standardize = TRUE
    ),
    silent = TRUE
  )
  
  ##########################################################
  # Skip failed iterations
  ##########################################################
  
  if (
    inherits(
      boot_lasso,
      "try-error"
    )
  ) {
    
    next
    
  }
  
  ##########################################################
  # Store lambda
  ##########################################################
  
  lambda_min_values[i] <-
    boot_lasso$lambda.min
  
  ##########################################################
  # Extract coefficients
  ##########################################################
  
  boot_coef <- as.matrix(
    coef(
      boot_lasso,
      s = "lambda.min"
    )
  )
  
  boot_coef_values <- as.numeric(
    boot_coef[, 1]
  )
  
  names(
    boot_coef_values
  ) <- rownames(
    boot_coef
  )
  
  ##########################################################
  # Store coefficients
  ##########################################################
  
  coefficient_matrix[
    i,
    candidate_genes
  ] <-
    boot_coef_values[
      candidate_genes
    ]
  
  ##########################################################
  # Store selection status
  ##########################################################
  
  selected_boot_genes <-
    candidate_genes[
      boot_coef_values[
        candidate_genes
      ] != 0
    ]
  
  if (
    length(
      selected_boot_genes
    ) > 0
  ) {
    
    selection_matrix[
      i,
      selected_boot_genes
    ] <- 1
    
  }
  
  ##########################################################
  # Progress
  ##########################################################
  
  if (
    i %% 100 == 0
  ) {
    
    cat(
      "Completed bootstrap:",
      i,
      "/",
      n_boot,
      "\n"
    )
    
  }
}

############################################################
# 16. DETERMINE VALID ITERATIONS
############################################################

valid_iterations <- rowSums(
  selection_matrix
) >= 0

valid_iterations[
  is.na(
    lambda_min_values
  )
] <- FALSE

selection_matrix_valid <-
  selection_matrix[
    valid_iterations,
    ,
    drop = FALSE
  ]

coefficient_matrix_valid <-
  coefficient_matrix[
    valid_iterations,
    ,
    drop = FALSE
  ]

lambda_min_valid <-
  lambda_min_values[
    valid_iterations
  ]

############################################################
# 17. NUMBER OF VALID BOOTSTRAPS
############################################################

n_valid_boot <- nrow(
  selection_matrix_valid
)

cat(
  "\nRequested bootstrap iterations:",
  n_boot,
  "\n"
)

cat(
  "Valid bootstrap iterations:",
  n_valid_boot,
  "\n"
)

############################################################
# 18. CALCULATE GENE-SELECTION FREQUENCY
############################################################

selection_frequency <- colSums(
  selection_matrix_valid
)

selection_percentage <-
  100 *
  selection_frequency /
  n_valid_boot

############################################################
# 19. CALCULATE COEFFICIENT SUMMARIES
############################################################

coefficient_mean <- apply(
  coefficient_matrix_valid,
  2,
  mean,
  na.rm = TRUE
)

coefficient_median <- apply(
  coefficient_matrix_valid,
  2,
  median,
  na.rm = TRUE
)

coefficient_lower <- apply(
  coefficient_matrix_valid,
  2,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

coefficient_upper <- apply(
  coefficient_matrix_valid,
  2,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

############################################################
# 20. CREATE GENE-SELECTION TABLE
############################################################

gene_selection_frequency <- data.frame(
  
  Gene = candidate_genes,
  
  Selection_Count =
    as.numeric(
      selection_frequency[
        candidate_genes
      ]
    ),
  
  Selection_Frequency =
    as.numeric(
      selection_percentage[
        candidate_genes
      ]
    ),
  
  Mean_Coefficient =
    as.numeric(
      coefficient_mean[
        candidate_genes
      ]
    ),
  
  Median_Coefficient =
    as.numeric(
      coefficient_median[
        candidate_genes
      ]
    ),
  
  Coefficient_2.5CI =
    as.numeric(
      coefficient_lower[
        candidate_genes
      ]
    ),
  
  Coefficient_97.5CI =
    as.numeric(
      coefficient_upper[
        candidate_genes
      ]
    )
)

############################################################
# 21. ADD ORIGINAL LASSO STATUS
############################################################

gene_selection_frequency$Original_LASSO_Selected <-
  gene_selection_frequency$Gene %in%
  original_selected$Gene

############################################################
# 22. ORDER BY SELECTION FREQUENCY
############################################################

gene_selection_frequency <-
  gene_selection_frequency[
    order(
      gene_selection_frequency$Selection_Frequency,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

############################################################
# 23. PRINT RESULTS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "BOOTSTRAP LASSO GENE-SELECTION FREQUENCY\n"
)

cat(
  "====================================================\n"
)

print(
  gene_selection_frequency
)

############################################################
# 24. SAVE GENE-SELECTION RESULTS
############################################################

write.csv(
  gene_selection_frequency,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_LASSO_Gene_Selection_Frequency.csv",
  row.names = FALSE
)

############################################################
# 25. SAVE BOOTSTRAP COEFFICIENTS
############################################################

bootstrap_coefficients_df <- data.frame(
  Iteration = seq_len(
    n_valid_boot
  ),
  coefficient_matrix_valid,
  check.names = FALSE
)

write.csv(
  bootstrap_coefficients_df,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_LASSO_Coefficients.csv",
  row.names = FALSE
)

############################################################
# 26. SAVE LAMBDA VALUES
############################################################

lambda_table <- data.frame(
  Iteration = seq_len(
    n_valid_boot
  ),
  Lambda_Min = lambda_min_valid
)

write.csv(
  lambda_table,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_LASSO_Lambda.csv",
  row.names = FALSE
)

############################################################
# 27. GENE-SELECTION FREQUENCY PLOT
############################################################

png(
  "Figures/TCGA_SKCM_Bootstrap_LASSO_Gene_Selection_Frequency.png",
  width = 2400,
  height = 1800,
  res = 300
)

barplot(
  gene_selection_frequency$Selection_Frequency,
  names.arg =
    gene_selection_frequency$Gene,
  ylim = c(
    0,
    100
  ),
  ylab =
    "Selection frequency (%)",
  xlab =
    "Candidate gene",
  main =
    "Bootstrap LASSO Gene-Selection Frequency"
)

abline(
  h = 50,
  lty = 2
)

dev.off()

############################################################
# 28. SAVE ORIGINAL LASSO RESULTS
############################################################

write.csv(
  original_coef_table,
  "Results/Prognostic/TCGA_SKCM_Original_LASSO_Coefficients.csv",
  row.names = FALSE
)

write.csv(
  original_selected,
  "Results/Prognostic/TCGA_SKCM_Original_LASSO_Selected_Genes.csv",
  row.names = FALSE
)

############################################################
# 29. SAVE WORKSPACE
############################################################

save(
  candidate_genes,
  lasso_data,
  x_lasso,
  y_lasso,
  original_lasso,
  original_coef_table,
  original_selected,
  selection_matrix_valid,
  coefficient_matrix_valid,
  lambda_min_valid,
  gene_selection_frequency,
  bootstrap_coefficients_df,
  lambda_table,
  n_boot,
  n_valid_boot,
  file =
    "RData/Script10_TCGA_Bootstrap_LASSO.RData"
)

############################################################
# 30. FINAL SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "FINAL BOOTSTRAP LASSO SUMMARY\n"
)

cat(
  "====================================================\n"
)

cat(
  "\nValid bootstrap iterations:",
  n_valid_boot,
  "\n"
)

cat(
  "\nGene-selection frequency:\n"
)

print(
  gene_selection_frequency[
    ,
    c(
      "Gene",
      "Selection_Count",
      "Selection_Frequency",
      "Original_LASSO_Selected"
    )
  ]
)

############################################################
# END OF PART 5
############################################################