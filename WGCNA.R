# Get data
setwd("~/Tirocinio")
filename_counts = "~/Tirocinio/RData/NormRawcounts_all.RData"
rawcounts = load(file = filename_counts)

#filename_tpm = "~/Tirocinio/RData/TPM_all.RData"
#tpm = load(file = filename_tpm)

# Soft thresholding power
power <- c(c(1:10), seq(from = 12, to = 20, by = 2))

sft <- pickSoftThreshold(data, dataIsExpr = T,
                         powerVector = power,
                         networkType = 'signed',
                         verbose = 5)
sft_data <- sft$fitIndices

Rplot <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale-free topology model fit, signed R^2') +
  theme_classic()

Meank_plot <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean connectivity') +
  theme_classic()

grid.arrange(Rplot, Meank_plot, nrow = 2)

remove(power, Rplot, Meank_plot)

soft_power <- 10
temp_corr <- cor
cor <- WGCNA::cor

# Modules
net <- blockwiseModules(data, 
                        power = soft_power,
                        TOMType = 'signed',
                        maxBlockSize = 15000,
                        minModuleSize = 40,
                        mergeCutHeight = 0.25,
                        numercLabels = F, 
                        pamRespectsDendro = F,
                        saveTOMs = T, 
                        saveTOMFileBase = "TOM",
                        randomSeed = 1234, 
                        verbose = 3)

cor <- temp_corr

moduleColors = net$colors
MEs = moduleEigengenes(data, moduleColors)$eigengenes
MEs = orderMEs(MEs)

ModuleColors_df <- data.frame(table(moduleColors))
ModuleColors_df <- ModuleColors_df[order(ModuleColors_df$Freq, decreasing = T),]

par(mar=c(10,5,3,1))
lab <- paste0(as.character(ModuleColors_df$moduleColors),"(",as.character(ModuleColors_df$Freq),")")
barplot(ModuleColors_df$Freq, col = as.character(ModuleColors_df$moduleColors),
        ylab = "# Genes", names.arg = lab, las = 2, cex.names = 1.1)

geneTree = net$dendrograms[[1]]
plotDendroAndColors(geneTree, moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Modules-traits association
status_map <- c("concordanti" = 0, "discordanti" = 1)
traits$GruppoNum <- status_map[traits$GRUPPO_COMBO]

nSamples <- nrow(data)
nGenes <- ncol(data)

module_trait_corr <- cor(MEs, traits$GruppoNum, use = 'p')
colnames(module_trait_corr) <- "Module-trait correlation"
module_trait_corr_pvalue <- corPvalueStudent(module_trait_corr, nSamples)
colnames(module_trait_corr_pvalue) <- "P-value"
module_trait_correlation <- data.frame(module_trait_corr, module_trait_corr_pvalue)

textMatrix = paste(signif(module_trait_corr, 2), "\n(", signif(module_trait_corr_pvalue, 1), ")", sep = "")
dim(textMatrix) = dim(module_trait_corr)

labeledHeatmap(Matrix = module_trait_corr,
               xLabels = 'GRUPPO',
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

heatmap.data <- cbind(MEs, Gruppo = traits$GruppoNum)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[31],
             y = names(heatmap.data)[1:30],
             col = c('blue', 'skyblue', 'white', 'pink', 'red'))

module_gene_mapping <- as.data.frame(moduleColors)
saddlebrown_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'saddlebrown']
lightcyan_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'lightcyan']
green_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'green']
purple_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'purple']
lightyellow_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'lightyellow']
midnightblue_genes <- row.names(module_gene_mapping)[module_gene_mapping$moduleColors == 'midnightblue']

# Representatives of each gene within its module
ModuleMembership <- cor(data, MEs, use = 'p')
ModuleMembership_pvalue <- corPvalueStudent(ModuleMembership, nSamples)

# Gene significance for traits
GeneSignificance <- cor(data, traits$GruppoNum, use = 'p')
GeneSignificance_pvalue <- corPvalueStudent(GeneSignificance, nSamples)

# Save data
dirRes <- "C:/Users/HP/Documents/Tirocinio/Rdata/"
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "WGCNA_results_all.RData" , sep = "")
save(net, MEs, ModuleColors_df, module_trait_correlation, module_gene_mapping, saddlebrown_genes, lightcyan_genes, green_genes, purple_genes, lightyellow_genes, midnightblue_genes, ModuleMembership, ModuleMembership_pvalue, GeneSignificance, GeneSignificance_pvalue, file = filename)
