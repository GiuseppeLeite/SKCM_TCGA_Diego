# -------------------------------------------------------------------------------------
# Survival Analysis Pipeline
#
# This script performs survival analysis using Kaplan–Meier methods. It reads in clinical 
# data and relevant expression data, calculates survival estimates, and generates Kaplan–Meier 
# plots to compare survival outcomes between different groups.
#
# Author: Giuseppe Leite, PhD | ggfleite@gmail.com
# Date: 15/04/2025
# -------------------------------------------------------------------------------------

# ----------------------------- Load Required Packages ----------------------------------- 
library(ggpubr)
library(tidyverse)
library(rstatix)
library(parameters)
library(survival)
library(survminer)
library(ggcorrplot)
library(pheatmap)

set.seed(123)

# ----------------------------- Load and Merge Data -------------------------------------- 
sample_annotation_cleaned <- read_delim("sample_annotation_cleaned.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

clustered_samples <- read.csv("clustered_Samples_VSD.csv") %>%
  rename(SAMPLE_ID = sample_id)

sample_annotation <- merge(sample_annotation_cleaned, clustered_samples, by = "SAMPLE_ID", all = FALSE)

# --------------------------- Format Survival Status Variables --------------------------- 
sample_annotation <- sample_annotation %>%
  mutate(
    OS_STATUS = ifelse(OS_STATUS == "1:DECEASED", 1, 0),
    DSS_STATUS = ifelse(DSS_STATUS == "1:DEAD WITH TUMOR", 1, 0),
    cluster = factor(cluster)
  )

# ---------------------- Kaplan-Meier Survival Curve - OS ------------------------------- 
fit_os <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data = sample_annotation)
cluster_colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")

surv_plot_os <- ggsurvplot(
  fit_os,
  data = sample_annotation,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  ggtheme = theme_bw(),
  xlab = "Time in months",
  ylab = "Survival Probability",
  title = "Overall Survival per cluster",
  risk.table.col = "black",
  legend.labs = paste("cluster", 1:4),
  palette = cluster_colors
)
print(surv_plot_os)
dev.copy2pdf(file = "OS_STATUS.pdf", width = 5, height = 7)
dev.off()

# ---------------------- Kaplan-Meier Survival Curve - DSS ------------------------------ 
fit_dss <- survfit(Surv(DSS_MONTHS, DSS_STATUS) ~ cluster, data = sample_annotation)

surv_plot_dss <- ggsurvplot(
  fit_dss,
  data = sample_annotation,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  ggtheme = theme_bw(),
  xlab = "Time in months",
  ylab = "Survival Probability",
  title = "Disease-Specific Survival per cluster",
  risk.table.col = "black",
  legend.labs = paste("cluster", 1:4),
  palette = cluster_colors
)
print(surv_plot_dss)
dev.copy2pdf(file = "DSS_STATUS.pdf", width = 5, height = 7)
dev.off()


# -------------------------- Cox Proportional Hazards Model - OS ------------------------ 

# Unadjusted Model
cox_model_unadjusted_os <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ relevel(cluster, ref = "2"), data = sample_annotation)



cox_summary_unadjusted_os <- summary(cox_model_unadjusted_os)

cox_results_unadjusted_os <- data.frame(
  variavel = rownames(cox_summary_unadjusted_os$coefficients),
  hr = exp(cox_summary_unadjusted_os$coefficients[, "coef"]),
  lower_ci = exp(cox_summary_unadjusted_os$conf.int[, "lower .95"]),
  upper_ci = exp(cox_summary_unadjusted_os$conf.int[, "upper .95"]),
  p_value = cox_summary_unadjusted_os$coefficients[, "Pr(>|z|)"]
)
write.csv(cox_results_unadjusted_os, "resumo_modelo_cox_OS_unadjusted.csv", row.names = FALSE)

# Adjusted Model

cox_model_adjusted_os <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ relevel(cluster, ref = "2") + AJCC_PATHOLOGIC_TUMOR_STAGE, data = sample_annotation)
cox_summary_adjusted_os <- summary(cox_model_adjusted_os)
cox_results_adjusted_os <- data.frame(
  variavel = rownames(cox_summary_adjusted_os$coefficients),
  hr = exp(cox_summary_adjusted_os$coefficients[, "coef"]),
  lower_ci = exp(cox_summary_adjusted_os$conf.int[, "lower .95"]),
  upper_ci = exp(cox_summary_adjusted_os$conf.int[, "upper .95"]),
  p_value = cox_summary_adjusted_os$coefficients[, "Pr(>|z|)"]
)
write.csv(cox_results_adjusted_os, "resumo_modelo_cox_OS_adjusted.csv", row.names = FALSE)


# -------------------------- Cox Proportional Hazards Model - DSS ------------------------ 
# Unadjusted Model
cox_model_unadjusted_dss <- coxph(
  Surv(DSS_MONTHS, DSS_STATUS) ~ relevel(cluster, ref = "2"),
  data = sample_annotation
)

cox_summary_unadjusted_dss <- summary(cox_model_unadjusted_dss)

cox_results_unadjusted_dss <- data.frame(
  variavel = rownames(cox_summary_unadjusted_dss$coefficients),
  hr = exp(cox_summary_unadjusted_dss$coefficients[, "coef"]),
  lower_ci = exp(cox_summary_unadjusted_dss$conf.int[, "lower .95"]),
  upper_ci = exp(cox_summary_unadjusted_dss$conf.int[, "upper .95"]),
  p_value = cox_summary_unadjusted_dss$coefficients[, "Pr(>|z|)"]
)
write.csv(cox_results_unadjusted_dss, "resumo_modelo_cox_DSS_unadjusted.csv", row.names = FALSE)

# Adjusted Model

cox_model_adjusted_dss <- coxph(Surv(DSS_MONTHS, DSS_STATUS) ~ relevel(cluster, ref = "2") + AJCC_PATHOLOGIC_TUMOR_STAGE, data = sample_annotation)
cox_summary_adjusted_dss <- summary(cox_model_adjusted_dss)
cox_results_adjusted_dss <- data.frame(
  variavel = rownames(cox_summary_adjusted_dss$coefficients),
  hr = exp(cox_summary_adjusted_dss$coefficients[, "coef"]),
  lower_ci = exp(cox_summary_adjusted_dss$conf.int[, "lower .95"]),
  upper_ci = exp(cox_summary_adjusted_dss$conf.int[, "upper .95"]),
  p_value = cox_summary_adjusted_dss$coefficients[, "Pr(>|z|)"]
)
write.csv(cox_results_adjusted_dss, "resumo_modelo_cox_DSS_adjusted.csv", row.names = FALSE)
