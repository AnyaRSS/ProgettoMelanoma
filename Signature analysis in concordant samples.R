setwd("~/Tirocinio")
library(tidyverse)
library(dplyr)
library(ggplot2)

#Take tpm and traits
tpm0 <- read.csv("~/Tirocinio/melanoma_gene_tpm.csv", header = T)
cellcolture <- "IRE.T1841"
tpm1 <- select(tpm0, -starts_with(cellcolture))
remove(cellcolture)
remove(tpm0)

samples <- names(tpm1)[grep("IRE", names(tpm1))]
aggregate_tpm <- aggregate(. ~ gene_name, data = tpm1[, c("gene_name", samples)], FUN = sum)
remove(tpm1)

rownames(aggregate_tpm) <- NULL
tpm <- column_to_rownames(aggregate_tpm, var = "gene_name")
remove(aggregate_tpm)

traits_0 <- read.csv("~/Tirocinio/traits.csv", header = T)
traits <- traits_0[order(traits_0$PATIENT_ID), ]
samples <- sort(samples)
rownames(traits) <- samples
remove(traits_0)
remove(samples)

#Take concordant samples
conc <- rownames(traits[traits$GRUPPO == "CONCORDANTE", ])
samples_conc <- select(tpm, all_of(conc))

#Take signature genes
signature1 <- read.csv("~/Tirocinio/signature_step1.csv", header = T)
signature2 <- read.csv("~/Tirocinio/signature_step2.csv", header = T)
sig_genes1 <- signature1$Gene
sig_genes2 <- signature2$gene_name

common_elements <- intersect(sig_genes1, sig_genes2)
sig_genes1 <- sig_genes1[!sig_genes1 %in% common_elements]
sig_genes <- append(sig_genes1, sig_genes2)
remove(signature1)
remove(signature2)
remove(sig_genes1)
remove(sig_genes2)
remove(common_elements)

#Select signature genes on concordant samples 
conc_sig_genes <- samples_conc[sig_genes, ]
remove(samples_conc)
write.csv(conc_sig_genes, file = "Concordant_tpm.csv")

#Statistics of signature genes in concordant samples
summary_stats <- apply(conc_sig_genes, 1, function(x) {
  c(mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    Q1 = quantile(x, probs = 0.25, na.rm = TRUE),
    Q3 = quantile(x, probs = 0.75, na.rm = TRUE))
})

summary_stats_df <- as.data.frame(t(as.data.frame(summary_stats)))
rownames(summary_stats_df) <- rownames(conc_sig_genes)
remove(summary_stats)

write.csv(summary_stats_df, file = "Concordant_tpm_statistics.csv")

#Plot a distribution
ggplot(summary_stats_df, aes(x = mean)) +
  geom_density(fill = "skyblue", color = "blue") +
  labs(title = "Distribution of Mean Expression",
       x = "Mean Expression",
       y = "Density") +
  theme_minimal()

#Write a table only with specific identified genes
selected_sig_genes <- c('PDPN', 'LEP', 'COL11A1', 'CXCL12', 'MXRA5', 'THBS2')

signature_genes <- conc_sig_genes[selected_sig_genes, ]
write.csv(signature_genes, file = "Signature_genes.csv")

summary_stats <- apply(signature_genes, 1, function(x) {
  c(mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    Q1 = quantile(x, probs = 0.25, na.rm = TRUE),
    Q3 = quantile(x, probs = 0.75, na.rm = TRUE))
})

selected_signature_summary_stats <- as.data.frame(t(as.data.frame(summary_stats)))
rownames(selected_signature_summary_stats) <- rownames(signature_genes)
remove(summary_stats)

write.csv(selected_signature_summary_stats, file = "Signature_statistics.csv")

ggplot(selected_signature_summary_stats, aes(x = mean)) +
  geom_density(fill = "skyblue", color = "blue") +
  labs(title = "Distribution of Mean Expression",
       x = "Mean Expression",
       y = "Density") +
  theme_minimal()