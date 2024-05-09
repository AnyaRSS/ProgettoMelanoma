setwd("~/Tirocinio")

library(dplyr)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

tpm0 <- read.csv("~/Tirocinio/melanoma_gene_tpm.csv", header = T)
cellcolture <- "IRE.T1841"
tpm1 <- select(tpm0, -starts_with(cellcolture))
remove(cellcolture)
remove(tpm0)

tpm1$base_gene <- gsub("-AS.*", "", tpm1$gene_name)
samples <- names(tpm1)[grep("IRE", names(tpm1))]
aggregate_tpm <- aggregate(. ~ base_gene, data = tpm1[, c("base_gene", samples)], FUN = sum)

names(aggregate_tpm)[1] <- "gene_name"
remove(tpm1)

sum_tpm <- rowSums(aggregate_tpm[, -1])
filtered_tpm_0 <- aggregate_tpm[sum_tpm != 0, ]
remove(sum_tpm)
remove(aggregate_tpm)

rownames(filtered_tpm_0) <- NULL
filtered_tpm_1 <- column_to_rownames(filtered_tpm_0, var = "gene_name")
remove(filtered_tpm_0)

iqr_values <- apply(filtered_tpm_1, 1, IQR)
threshold <- 0.1 # Adjust as needed
data <- filtered_tpm_1[iqr_values > threshold, ]
remove(filtered_tpm_1)
remove(iqr_values)
remove(threshold)

# Get traits
traits_0 <- read.csv("~/Tirocinio/traits.csv", header = T)

# Data manipulation for match
traits <- traits_0[order(traits_0$PATIENT_ID), ]
samples <- sort(samples)
rownames(traits) <- samples
remove(traits_0)
remove(samples)

check1 = all(rownames(traits) %in% colnames(data))
check2 = all(rownames(traits) == colnames(data))

if (!check2) 
{
  sorted_colnames <- names(data)[order(names(data))]
  data <- data[, sorted_colnames]
}
remove(check1)
remove(check2)
remove(sorted_colnames)

# Get outlayers and remove them from the analysis
gsg <- goodSamplesGenes(t(data))
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(rawcounts)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints outlier samples
  data <- data.frame(t(data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE])) # Removes the offending genes and samples from the data
  remove(rawcounts)
}

# Detect outlayer samples - hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")
plot(htree, xlab = "distance", cex = 0.8)

# Detect outlayer samples - PCA
pca <- prcomp(t(data))
pca.data <- as.data.frame(pca$x)
pca.var <- pca$sdev**2
pca.var.perc <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca.data, aes(PC1, PC2)) + 
  geom_point() +
  geom_text(label = rownames(pca.data)) +
  labs(x = paste0("PC1: ", pca.var.perc[1], '%'),
       y = paste0("PC2: ", pca.var.perc[2], '%'))

remove(pca.data)
remove(pca.var)
remove(pca.var.perc)

#Save needed data
data <- as.data.frame(t(data))

dirRes <- "C:/Users/HP/Documents/Tirocinio/Rdata/"
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "TPM.RData" , sep = "")
save(data, traits, file = filename)
