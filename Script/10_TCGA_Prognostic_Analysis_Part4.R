############################################################
# SCRIPT 10 – TCGA-SKCM PROGNOSTIC ANALYSIS
# PART 4 – BOOTSTRAP INTERNAL VALIDATION
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
# LOAD TIME-DEPENDENT ROC OBJECTS IF AVAILABLE
############################################################

if (
  file.exists(
    "RData/Script10_TCGA_TimeDependent_ROC.RData"
  )
) {
  
  load(
    "RData/Script10_TCGA_TimeDependent_ROC.RData"
  )
  
}

############################################################
# 4. CHECK DATA
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "BOOTSTRAP INTERNAL VALIDATION\n"
)

cat(
  "====================================================\n"
)

cat(
  "\nNumber of samples:",
  nrow(analysis_data_surv),
  "\n"
)

cat(
  "Number of events:",
  sum(
    analysis_data_surv$OS_event == 1,
    na.rm = TRUE
  ),
  "\n"
)

############################################################
# 5. CHECK RISK SCORE
############################################################

summary(
  analysis_data_surv$RiskScore
)

############################################################
# 6. DEFINE BOOTSTRAP SETTINGS
############################################################

set.seed(1234)

n_boot <- 1000

############################################################
# 7. STORAGE OBJECTS
############################################################

bootstrap_coef <- numeric(
  n_boot
)

bootstrap_hr <- numeric(
  n_boot
)

bootstrap_cindex <- numeric(
  n_boot
)

############################################################
# 8. BOOTSTRAP COX VALIDATION
############################################################

for (
  i in seq_len(n_boot)
) {
  
  ##########################################################
  # Bootstrap sample
  ##########################################################
  
  boot_index <- sample(
    seq_len(
      nrow(analysis_data_surv)
    ),
    size = nrow(
      analysis_data_surv
    ),
    replace = TRUE
  )
  
  boot_data <- analysis_data_surv[
    boot_index,
    ,
    drop = FALSE
  ]
  
  ##########################################################
  # Check sufficient variation
  ##########################################################
  
  if (
    length(
      unique(
        boot_data$OS_event
      )
    ) < 2
  ) {
    
    next
    
  }
  
  ##########################################################
  # Bootstrap Cox model
  ##########################################################
  
  boot_model <- try(
    coxph(
      Surv(
        OS_time,
        OS_event
      ) ~ RiskScore,
      data = boot_data
    ),
    silent = TRUE
  )
  
  if (
    inherits(
      boot_model,
      "try-error"
    )
  ) {
    
    next
    
  }
  
  ##########################################################
  # Store coefficient and HR
  ##########################################################
  
  bootstrap_coef[i] <-
    coef(
      boot_model
    )[1]
  
  bootstrap_hr[i] <-
    exp(
      coef(
        boot_model
      )[1]
    )
  
  ##########################################################
  # Store concordance
  ##########################################################
  
  boot_summary <- summary(
    boot_model
  )
  
  bootstrap_cindex[i] <-
    as.numeric(
      boot_summary$concordance[1]
    )
}

############################################################
# 9. REMOVE FAILED ITERATIONS
############################################################

valid_boot <- !is.na(
  bootstrap_hr
)

bootstrap_hr_valid <-
  bootstrap_hr[
    valid_boot
  ]

bootstrap_coef_valid <-
  bootstrap_coef[
    valid_boot
  ]

bootstrap_cindex_valid <-
  bootstrap_cindex[
    valid_boot
  ]

############################################################
# 10. BOOTSTRAP HR SUMMARY
############################################################

bootstrap_hr_summary <- data.frame(
  Statistic = c(
    "Original HR",
    "Bootstrap mean HR",
    "Bootstrap median HR",
    "Bootstrap 2.5% CI",
    "Bootstrap 97.5% CI"
  ),
  Value = c(
    exp(
      coef(
        coxph(
          Surv(
            OS_time,
            OS_event
          ) ~ RiskScore,
          data = analysis_data_surv
        )
      )
    ),
    mean(
      bootstrap_hr_valid
    ),
    median(
      bootstrap_hr_valid
    ),
    quantile(
      bootstrap_hr_valid,
      0.025,
      na.rm = TRUE
    ),
    quantile(
      bootstrap_hr_valid,
      0.975,
      na.rm = TRUE
    )
  )
)

print(
  bootstrap_hr_summary
)

############################################################
# 11. BOOTSTRAP C-INDEX SUMMARY
############################################################

bootstrap_cindex_summary <- data.frame(
  Statistic = c(
    "Original C-index",
    "Bootstrap mean C-index",
    "Bootstrap median C-index",
    "Bootstrap 2.5% CI",
    "Bootstrap 97.5% CI"
  ),
  Value = c(
    summary(
      coxph(
        Surv(
          OS_time,
          OS_event
        ) ~ RiskScore,
        data = analysis_data_surv
      )
    )$concordance[1],
    mean(
      bootstrap_cindex_valid,
      na.rm = TRUE
    ),
    median(
      bootstrap_cindex_valid,
      na.rm = TRUE
    ),
    quantile(
      bootstrap_cindex_valid,
      0.025,
      na.rm = TRUE
    ),
    quantile(
      bootstrap_cindex_valid,
      0.975,
      na.rm = TRUE
    )
  )
)

print(
  bootstrap_cindex_summary
)

############################################################
# 12. SAVE BOOTSTRAP RESULTS
############################################################

write.csv(
  bootstrap_hr_summary,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_HR_Validation.csv",
  row.names = FALSE
)

write.csv(
  bootstrap_cindex_summary,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_CIndex_Validation.csv",
  row.names = FALSE
)

############################################################
# 13. SAVE ALL BOOTSTRAP VALUES
############################################################

bootstrap_results <- data.frame(
  Iteration = seq_len(
    n_boot
  ),
  Coefficient = bootstrap_coef,
  HR = bootstrap_hr,
  C_index = bootstrap_cindex
)

write.csv(
  bootstrap_results,
  "Results/Prognostic/TCGA_SKCM_Bootstrap_All_Iterations.csv",
  row.names = FALSE
)

############################################################
# 14. BOOTSTRAP HR DISTRIBUTION
############################################################

png(
  "Figures/TCGA_SKCM_Bootstrap_HR_Distribution.png",
  width = 2400,
  height = 1800,
  res = 300
)

hist(
  bootstrap_hr_valid,
  breaks = 40,
  main = "Bootstrap Distribution of RiskScore Hazard Ratio",
  xlab = "Hazard Ratio"
)

abline(
  v = median(
    bootstrap_hr_valid
  ),
  lty = 2,
  lwd = 2
)

dev.off()

############################################################
# 15. BOOTSTRAP C-INDEX DISTRIBUTION
############################################################

png(
  "Figures/TCGA_SKCM_Bootstrap_CIndex_Distribution.png",
  width = 2400,
  height = 1800,
  res = 300
)

hist(
  bootstrap_cindex_valid,
  breaks = 40,
  main = "Bootstrap Distribution of C-index",
  xlab = "C-index"
)

abline(
  v = median(
    bootstrap_cindex_valid
  ),
  lty = 2,
  lwd = 2
)

dev.off()

############################################################
# 16. SAVE BOOTSTRAP WORKSPACE
############################################################

save(
  bootstrap_results,
  bootstrap_hr_valid,
  bootstrap_coef_valid,
  bootstrap_cindex_valid,
  bootstrap_hr_summary,
  bootstrap_cindex_summary,
  n_boot,
  file =
    "RData/Script10_TCGA_Bootstrap_Validation.RData"
)

############################################################
# 17. FINAL SUMMARY
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "BOOTSTRAP VALIDATION SUMMARY\n"
)

cat(
  "====================================================\n"
)

cat(
  "\nRequested bootstrap iterations:",
  n_boot,
  "\n"
)

cat(
  "Valid bootstrap iterations:",
  length(
    bootstrap_hr_valid
  ),
  "\n"
)

cat(
  "\nBootstrap median HR:",
  median(
    bootstrap_hr_valid,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Bootstrap 95% CI:",
  quantile(
    bootstrap_hr_valid,
    0.025,
    na.rm = TRUE
  ),
  "–",
  quantile(
    bootstrap_hr_valid,
    0.975,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nBootstrap median C-index:",
  median(
    bootstrap_cindex_valid,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Bootstrap C-index 95% CI:",
  quantile(
    bootstrap_cindex_valid,
    0.025,
    na.rm = TRUE
  ),
  "–",
  quantile(
    bootstrap_cindex_valid,
    0.975,
    na.rm = TRUE
  ),
  "\n"
)

############################################################
# END OF PART 4
############################################################
