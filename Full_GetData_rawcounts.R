setwd("~/Tirocinio")

library(dplyr)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(ggplot2)
library(factoextra)
library(sva)

#GET DATA
rc1 <- read.csv("~/Tirocinio/melanoma_gene_counts_1.csv", header = T)
rc2 <- read.delim("~/Tirocinio/melanoma_gene_counts_2.tsv", header = T)

rc1 <- rc1[, -c(1, 2)]
rc2 <- rc2[, -1]

cellcolture <- "IRE.T1841"
rc1 <- select(rc1, -starts_with(cellcolture)) #remove cell colture
remove(cellcolture)

rc1$base_gene <- gsub("-AS.*", "", rc1$gene_name)
rc2$base_gene <- gsub("-AS.*", "", rc2$gene_name)
samples1 <- names(rc1)[grep("IRE", names(rc1))]
samples2 <- names(rc2)[grep("IRE", names(rc2))]

#aggregate all genes with same name or splicing derivatives
aggregate_1 <- aggregate(. ~ base_gene, data = rc1[, c("base_gene", samples1)], FUN = sum)
aggregate_2 <- aggregate(. ~ base_gene, data = rc2[, c("base_gene", samples2)], FUN = sum) 
names(aggregate_1)[1] <- "gene_name"
names(aggregate_2)[1] <- "gene_name"

rownames(aggregate_1) <- NULL
rownames(aggregate_2) <- NULL
rc1 <- column_to_rownames(aggregate_1, var = "gene_name")
rc2 <- column_to_rownames(aggregate_2, var = "gene_name")

remove(aggregate_1, aggregate_2)

rawcounts <- cbind(rc1, rc2) #merge the two rc files
rm(rc1, rc2)

sum_rc <- rowSums(rawcounts)
filtered_rc <- rawcounts[sum_rc != 0, ]
remove(sum_rc)

iqr_values <- apply(filtered_rc, 1, IQR)
threshold <- 0.1 # Adjust as needed
data <- filtered_rc[iqr_values > threshold, ]
remove(filtered_rc, iqr_values, threshold, rawcounts)

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

# Normalization and DESeq2
data_int <- round(data)
dds <- DESeqDataSetFromMatrix(countData = data_int,
                              colData = traits,
                              design = ~ 1)

dds75 <- dds[rowSums(counts(dds) >= 10) >= 10, ] #keep only genes with more than 10 counts in more than 75% of the samples
dds_norm <- vst(dds75) #variance stabilization

matr_norm_counts <- assay(dds_norm) %>%
  t()

remove(data_int)
remove(dds)
remove(dds75)
remove(dds_norm)

data <- as.data.frame(t(matr_norm_counts))

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
  geom_text(aes(label = Sample), vjust = -1, hjust = 1, show.legend = F) +
  labs(title="Surrogate Variables by Sample", x="Sample", y="Surrogate Variable Value")

# Get outlayers and remove them from the analysis
gsg <- goodSamplesGenes(as.matrix(data))
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

data <- as.data.frame(t(data))

#Save needed data
dirRes <- "C:/Users/HP/Documents/Tirocinio/Rdata/"
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "NormRawcounts_all.RData" , sep = "")
save(data, traits, file = filename)
