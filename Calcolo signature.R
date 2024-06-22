setwd("~/Tirocinio/Signature genes analysis/Calcolo signature nei database")

tpm <- read.delim("~/Tirocinio/Signature genes analysis/Calcolo signature nei database/data_mrna_seq_tpm_2019.txt",header = T)
tpm <- tpm[-2]

df_unique <- tpm[!duplicated(tpm$Hugo_Symbol), ]
rownames(df_unique) <- NULL
tpm <- column_to_rownames(df_unique, var = "Hugo_Symbol")
rm(df_unique)

signature_genes <- tpm[c("LEP", "PDPN", "CXCL12", "COL11A1", "THBS2", "MXRA5", "NRAS", "BRAF"), ]
signature_genes <- signature_genes[, order(colnames(signature_genes))]
samples <- colnames(signature_genes)

sig_step1 <- c('LEP', 'PDPN')
sig_step2 <- c("CXCL12", "COL11A1", "THBS2", "MXRA5")
norm <- c('NRAS', 'BRAF')

df_signature <- data.frame()

for (col in seq_along(signature_genes)) {
  med1 <- median(signature_genes[sig_step1, col])
  med2 <- median(signature_genes[sig_step2, col])
  med_norm <- median(signature_genes[norm, col])
  sum <- sum(signature_genes[sig_step1, col], signature_genes[sig_step2, col], signature_genes[norm, col])
  signature <- (med1 + med2 - med_norm) / sum
  new_row <- data.frame(signature)
  df_signature <- rbind(df_signature, new_row)
}

rownames(df_signature) <- samples

#Save needed data
dirRes <- "C:/Users/HP/Documents/Tirocinio/Signature genes analysis/Calcolo signature nei database/RData/"
if(!file.exists(dirRes)){
  dir.create(dirRes)
}

filename = paste(dirRes, "Signature2019.RData" , sep = "")
save(df_signature, file = filename)
