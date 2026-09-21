







############################################################
# LOAD REQUIRED PACKAGES
############################################################

library(survival)
library(survminer)
library(ggplot2)





############################################################
# SCRIPT 10 — TCGA-SKCM PROGNOSTIC ANALYSIS
# PART 2 — KM, COX REGRESSION, CLINICAL VARIABLES,
#          MULTIVARIABLE ANALYSIS, PH ASSUMPTION
############################################################

############################################################
# 10.1 SURVIVAL OBJECT
############################################################

surv_object <- Surv(
  analysis_data_surv$OS_time,
  analysis_data_surv$OS_event
)


############################################################
# 10.2 KAPLAN–MEIER SURVIVAL ANALYSIS
############################################################

fit_km <- survfit(
  surv_object ~ RiskGroup,
  data = analysis_data_surv
)

print(summary(fit_km))


############################################################
# Kaplan–Meier plot
############################################################

km_plot <- ggsurvplot(
  fit_km,
  data = analysis_data_surv,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Overall survival time (days)",
  ylab = "Overall survival probability",
  title = "TCGA-SKCM Overall Survival by Risk Group",
  legend.title = "Risk Group",
  legend.labs = c("Low Risk", "High Risk"),
  palette = c("#2E86AB", "#D1495B")
)

print(km_plot)


############################################################
# Save KM plot
############################################################

ggsave(
  filename = "Figures/TCGA_SKCM_Kaplan_Meier_RiskScore.png",
  plot = km_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)


############################################################
# Save KM risk table
############################################################

if (!is.null(km_plot$table)) {
  
  ggsave(
    filename = "Figures/TCGA_SKCM_Kaplan_Meier_RiskTable.png",
    plot = km_plot$table,
    width = 8,
    height = 3,
    dpi = 300
  )
}


############################################################
# 10.3 LOG-RANK TEST
############################################################

logrank_test <- survdiff(
  surv_object ~ RiskGroup,
  data = analysis_data_surv
)

print(logrank_test)

logrank_p <- 1 - pchisq(
  logrank_test$chisq,
  df = length(logrank_test$n) - 1
)

print(logrank_p)


############################################################
# Save log-rank result
############################################################

logrank_result <- data.frame(
  Chi_square = logrank_test$chisq,
  df = length(logrank_test$n) - 1,
  P_value = logrank_p
)

write.csv(
  logrank_result,
  "Results/Prognostic/TCGA_SKCM_LogRank_Test.csv",
  row.names = FALSE
)


############################################################
# 10.4 UNIVARIATE COX REGRESSION
# RiskScore as a continuous variable
############################################################

cox_uni_model <- coxph(
  Surv(OS_time, OS_event) ~ RiskScore,
  data = analysis_data_surv
)

print(summary(cox_uni_model))


############################################################
# Extract univariate Cox statistics
############################################################

cox_uni_summary <- summary(cox_uni_model)

cox_uni_result <- data.frame(
  Variable = "RiskScore",
  HR = cox_uni_summary$coefficients[,"exp(coef)"],
  HR_lower_95CI = cox_uni_summary$conf.int[,"lower .95"],
  HR_upper_95CI = cox_uni_summary$conf.int[,"upper .95"],
  P_value = cox_uni_summary$coefficients[,"Pr(>|z|)"]
)

write.csv(
  cox_uni_result,
  "Results/Prognostic/TCGA_SKCM_Univariate_Cox_RiskScore.csv",
  row.names = FALSE
)


############################################################
# 10.5 PREPARE CLINICOPATHOLOGICAL VARIABLES
############################################################

cox_multi <- analysis_data_surv


############################################################
# AGE
############################################################

if ("age_at_diagnosis" %in% colnames(cox_multi)) {
  
  cox_multi$Age <- as.numeric(
    cox_multi$age_at_diagnosis
  )
  
} else if ("age_at_index" %in% colnames(cox_multi)) {
  
  cox_multi$Age <- as.numeric(
    cox_multi$age_at_index
  )
}


############################################################
# GENDER
############################################################

if ("gender.demographic" %in% colnames(cox_multi)) {
  
  cox_multi$Gender <- as.factor(
    cox_multi$gender.demographic
  )
}


############################################################
# 10.6 STAGE GROUPING
############################################################

if ("ajcc_pathologic_stage.diagnoses" %in% colnames(cox_multi)) {
  
  cox_multi$Stage_original <-
    cox_multi$ajcc_pathologic_stage.diagnoses
  
} else if ("ajcc_pathologic_stage" %in% colnames(cox_multi)) {
  
  cox_multi$Stage_original <-
    cox_multi$ajcc_pathologic_stage
}


############################################################
# Clean stage terminology
############################################################

cox_multi$Stage_original <- gsub(
  "^Stage ",
  "",
  cox_multi$Stage_original
)

cox_multi$Stage_original <- trimws(
  cox_multi$Stage_original
)


############################################################
# Stage I–II versus Stage III–IV
############################################################

cox_multi$Stage_group <- NA

cox_multi$Stage_group[
  cox_multi$Stage_original %in%
    c("I", "IA", "IB", "II", "IIA", "IIB", "IIC")
] <- "Stage I-II"

cox_multi$Stage_group[
  cox_multi$Stage_original %in%
    c("III", "IIIA", "IIIB", "IIIC", "IV")
] <- "Stage III-IV"


cox_multi$Stage_group <- factor(
  cox_multi$Stage_group,
  levels = c(
    "Stage I-II",
    "Stage III-IV"
  )
)


############################################################
# 10.7 NUMERIC STAGE VARIABLE
# Used for the time-varying stage sensitivity analysis
############################################################

cox_multi$Stage_III_IV_numeric <- ifelse(
  cox_multi$Stage_group == "Stage III-IV",
  1,
  0
)









############################################################
# CREATE AGE VARIABLE
############################################################

cox_multi$Age <- as.numeric(
  cox_multi$age_at_index.demographic
)

############################################################
# CHECK AGE
############################################################

summary(cox_multi$Age)


setdiff(
  c(
    "OS_time",
    "OS_event",
    "RiskScore",
    "Age",
    "Gender",
    "Stage_group"
  ),
  colnames(cox_multi)
)





############################################################
# 10.8 REMOVE INCOMPLETE CASES FOR MULTIVARIABLE COX
############################################################

cox_multi_complete <- cox_multi[
  complete.cases(
    cox_multi[, c(
      "OS_time",
      "OS_event",
      "RiskScore",
      "Age",
      "Gender",
      "Stage_group"
    )]
  ),
]


############################################################
# Check final multivariable cohort
############################################################

cat(
  "\nNumber of samples in multivariable Cox model:",
  nrow(cox_multi_complete),
  "\n"
)

print(
  table(
    cox_multi_complete$Stage_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_multi_complete$Gender,
    useNA = "ifany"
  )
)


############################################################
# 10.9 MULTIVARIABLE COX REGRESSION
############################################################

cox_multi_model <- coxph(
  Surv(OS_time, OS_event) ~
    RiskScore +
    Age +
    Gender +
    Stage_group,
  data = cox_multi_complete
)

print(summary(cox_multi_model))


############################################################
# Extract multivariable Cox results
############################################################

cox_multi_summary <- summary(
  cox_multi_model
)

cox_multi_result <- data.frame(
  Variable = rownames(
    cox_multi_summary$coefficients
  ),
  HR = cox_multi_summary$coefficients[,"exp(coef)"],
  HR_lower_95CI =
    cox_multi_summary$conf.int[,"lower .95"],
  HR_upper_95CI =
    cox_multi_summary$conf.int[,"upper .95"],
  P_value =
    cox_multi_summary$coefficients[,"Pr(>|z|)"]
)

write.csv(
  cox_multi_result,
  "Results/Prognostic/TCGA_SKCM_Multivariable_Cox.csv",
  row.names = FALSE
)


############################################################
# 10.10 PROPORTIONAL HAZARDS ASSUMPTION
############################################################

ph_test <- cox.zph(
  cox_multi_model
)

print(ph_test)


############################################################
# Save PH test results
############################################################

ph_result <- data.frame(
  Variable = rownames(ph_test$table),
  Chisq = ph_test$table[,"chisq"],
  df = ph_test$table[,"df"],
  P_value = ph_test$table[,"p"]
)

write.csv(
  ph_result,
  "Results/Prognostic/TCGA_SKCM_PH_Assumption_Test.csv",
  row.names = FALSE
)


############################################################
# PH diagnostic plots
############################################################

png(
  "Figures/TCGA_SKCM_PH_Assumption_Test.png",
  width = 2400,
  height = 1800,
  res = 300
)

plot(ph_test)

dev.off()


############################################################
# 10.11 TIME-VARYING STAGE SENSITIVITY ANALYSIS
############################################################

cox_tv_stage <- coxph(
  Surv(OS_time, OS_event) ~
    RiskScore +
    Age +
    Gender +
    Stage_group +
    tt(Stage_III_IV_numeric),
  data = cox_multi_complete,
  tt = function(x, t, ...) {
    x * log(t + 1)
  }
)

print(
  summary(cox_tv_stage)
)


############################################################
# Extract time-varying stage model
############################################################

cox_tv_summary <- summary(
  cox_tv_stage
)

cox_tv_result <- data.frame(
  Variable = rownames(
    cox_tv_summary$coefficients
  ),
  HR = cox_tv_summary$coefficients[,"exp(coef)"],
  HR_lower_95CI =
    cox_tv_summary$conf.int[,"lower .95"],
  HR_upper_95CI =
    cox_tv_summary$conf.int[,"upper .95"],
  P_value =
    cox_tv_summary$coefficients[,"Pr(>|z|)"]
)

write.csv(
  cox_tv_result,
  "Results/Prognostic/TCGA_SKCM_TimeVarying_Stage_Cox.csv",
  row.names = FALSE
)


############################################################
# 10.12 SAVE INTERMEDIATE TCGA PROGNOSTIC OBJECTS
############################################################

save(
  analysis_data_surv,
  cox_multi_complete,
  fit_km,
  cox_uni_model,
  cox_multi_model,
  ph_test,
  cox_tv_stage,
  file = "RData/Script10_TCGA_Prognostic_Intermediate.RData"
)


############################################################
# END OF PART 2
############################################################


colnames(cox_multi)




exists("analysis_data_surv")
exists("RiskScore")
exists("analysis_data")
exists("lasso_data")

savehistory("Scripts/10_TCGA_Prognostic_Analysis_Part2.R")














