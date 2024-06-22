library(survival)
library(ggplot2)
library(survminer)

# Get data
setwd("~/Tirocinio/Signature genes analysis/Calcolo signature nei database")
filename = "~/Tirocinio/Signature genes analysis/Calcolo signature nei database/RData/Signature2019.RData"
sig = load(file = filename)

#sample_data <- read.delim("~/Tirocinio/Signature genes analysis/Calcolo signature nei database/data_clinical_sample_2019.txt", header = T)
clinical_data <- read.delim("~/Tirocinio/Signature genes analysis/Calcolo signature nei database/data_clinical_patient_2019.txt", header = F)
clinical_data <- clinical_data[-c(1:4), ]
colnames(clinical_data) <- unlist(clinical_data[1, ])
clinical_data <- clinical_data[-1, ]

# Modify data to obtain the final df
sample_ids <- rownames(df_signature)
extract_numeric <- function(ids) {
  as.numeric(gsub("\\D", "", ids))
}

sample_numeric <- extract_numeric(sample_ids)
clinical_data$patient_numeric <- extract_numeric(clinical_data$PATIENT_ID)
shared_numeric <- intersect(sample_numeric, clinical_data$patient_numeric)
df_clinical <- clinical_data[clinical_data$patient_numeric %in% shared_numeric, ]
df_clinical <- df_clinical[, -ncol(df_clinical)]
rownames(df_clinical) <- sample_ids
rm(sample_ids, sample_numeric, shared_numeric)

df_clinical$SIGNATURE <- df_signature$signature

df <- df_clinical[, c("PATIENT_ID", "PFS_STATUS", 'OS_STATUS', 'PFS_MONTHS', 'OS_MONTHS', 'IO_THERAPY', 'SIGNATURE')]
df$PFS_MONTHS <- as.numeric(df$PFS_MONTHS)
df$OS_MONTHS <- as.numeric(df$OS_MONTHS)
df$SIGNATURE <- as.numeric(df$SIGNATURE)
df$OS_STATUS_LOGICAL <- ifelse(df$OS_STATUS == '1:DECEASED', T, F)
df$SIGNATURE_GROUP <- ifelse(df$SIGNATURE > 0, "Positive", "Negative")
#sixmonths_df <- df[df$PFS_MONTHS >= 6, ]

# Cox regression model, technique for assessing the association between signature and survival rate
cox_model <- coxph(Surv(OS_MONTHS, OS_STATUS_LOGICAL) ~ SIGNATURE, data = df)
summary(cox_model)

km_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_LOGICAL) ~ SIGNATURE_GROUP, data = df)
ggsurvplot(
  km_fit,
  data = df,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Signature Group",
  xlab = "Time",
  ylab = "Survival Probability",
  title = "Kaplan-Meier Curves by Signature Group"
)
