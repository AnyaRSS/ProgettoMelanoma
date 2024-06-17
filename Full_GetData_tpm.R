setwd("~/Tirocinio")

library(dplyr)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# GET DATA
tpm1 <- read.csv("~/Tirocinio/melanoma_gene_tpm_1.csv", header = T)
tpm2 <- read.csv("~/Tirocinio/melanoma_gene_tpm_2.csv", header = T)

tpm1 <- tpm1[, -c(1, 2)]
tpm2 <- tpm2[, -c(1, 2)]

cellcolture <- "IRE.T1841"
tpm1 <- select(tpm1, -starts_with(cellcolture)) #remove cell colture
remove(cellcolture)

tpm1$base_gene <- gsub("-AS.*", "", tpm1$gene_name)
tpm2$base_gene <- gsub("-AS.*", "", tpm2$gene_name)
samples1 <- names(tpm1)[grep("IRE", names(tpm1))]
samples2 <- names(tpm2)[grep("IRE", names(tpm2))]

#aggregate all genes with same name or splicing derivatives
aggregate_tpm1 <- aggregate(. ~ base_gene, data = tpm1[, c("base_gene", samples1)], FUN = sum)
aggregate_tpm2 <- aggregate(. ~ base_gene, data = tpm2[, c("base_gene", samples2)], FUN = sum) 
names(aggregate_tpm1)[1] <- "gene_name"
names(aggregate_tpm2)[1] <- "gene_name"

rownames(aggregate_tpm1) <- NULL
rownames(aggregate_tpm2) <- NULL
tpm1 <- column_to_rownames(aggregate_tpm1, var = "gene_name")
tpm2 <- column_to_rownames(aggregate_tpm2, var = "gene_name")

remove(aggregate_tpm1, aggregate_tpm2)

tpm <- cbind(tpm1, tpm2) #merge the two tpm files
rm(tpm1, tpm2)

sum_tpm <- rowSums(tpm)
filtered_tpm <- tpm[sum_tpm != 0, ]
remove(sum_tpm)

iqr_values <- apply(filtered_tpm, 1, IQR)
threshold <- 0.1 # Adjust as needed
data <- filtered_tpm[iqr_values > threshold, ]
remove(filtered_tpm, iqr_values, threshold, tpm)

# GET TRAITS
traits0 <- read.csv("~/Tirocinio/Traits_all.csv", header = T)

# DATA MANIPULATION FOR MATCH
rownames(traits0) <- NULL
traits <- column_to_rownames(traits0, var = "SAMPLE_ID")
rm(traits0)

check1 = all(rownames(traits) %in% colnames(data))
check2 = all(rownames(traits) == colnames(data))

if (!check2) 
{
  sorted_colnames <- names(data)[order(names(data))]
  data <- data[, sorted_colnames]
}
remove(check1, check2, sorted_colnames)

# BATCH EFFECT
batch <- ifelse(colnames(data) %in% samples1, "Batch1", "Batch2")

pca <- prcomp(t(data), scale = T)
pca.data <- data.frame(Sample = rownames(pca$x),
                       PC1 = pca$x[, 1],
                       PC2 = pca$x[, 2],
                       Batch = batch)
ggplot(pca.data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 4) +
  labs(title = "PCA of Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

remove(pca, pca.data)

traits$GRUPPO_COMBO <- ifelse(traits$GRUPPO_REAL == 'unknown', traits$GRUPPO_STIMA, traits$GRUPPO_REAL)
treatment <- traits$GRUPPO_COMBO
df <- data.frame(Sample = colnames(data), Batch = batch, Gruppo = treatment)

mod <- model.matrix(~ Gruppo, data = df)
mod0 <- model.matrix(~ 1, data = df)

exp_matrix <- as.matrix(data)

svobj <- sva(exp_matrix, mod, mod0)
svobj$n.sv # Number of significant surrogate variables (potential batch effects)
svobj$sv # Surrogate variables
mod_sv <- cbind(mod, svobj$sv)

sv_df <- as.data.frame(svobj$sv)
sv_df$Sample <- colnames(data)
sv_df_long <- reshape2::melt(sv_df, id.vars="Sample")

ggplot(sv_df_long, aes(x=Sample, y=value, color=variable)) +
  geom_point(size = 2) +
  labs(title="Surrogate Variables by Sample", x="Sample", y="Surrogate Variable Value")

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

#Save needed data
data <- as.data.frame(t(data))

dirRes <- "C:/Users/HP/Documents/Tirocinio/Rdata/"
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "TPM_all.RData" , sep = "")
save(data, traits, file = filename)
