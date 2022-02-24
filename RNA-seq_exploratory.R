# Import the necessary libraries
library(pheatmap)
library(stats)
library(ggplot2)
library(ggfortify)
library(corrplot)
library(DESeq2)
library(ggbeeswarm)
library(biomaRt)
library(dplyr)

# Import raw counts matrix for DEseq2
counts <- as.matrix(read.table('counts_from_SALMON.genes.tsv', 
                               header = T, sep = '\t'))
summary(counts)
head(counts, 3)

# Import TPM-normalised counts matrix for correlation plots and heatmaps
tpm <- read.table('TPM_counts_from_SALMON.genes.tsv', 
                  header = T, sep = '\t', row.names = 'Name')
tpm <- subset(tpm, select = -c(row.names))
tpm <- as.matrix(tpm)
summary(tpm)
head(tpm, 3)

# Provide the sample information sheet
colData <- read.table('colData.tsv', header = T, sep = '\t', 
                      stringsAsFactors = TRUE, row.names = 'sample_type')
print(colData)

# Ensure row names correspond to col names 
# if does not return TRUE, colData needs editing
all(rownames(colData) == colnames(counts))

# Clean up genes in tpm matrix with no reads
tpm[tpm == 0] <- NA
tpm_filtered <- tpm[!rowSums(!is.finite(tpm)),]

# Calculate the variance of each gene across samples and select the top 500
vars <- apply(tpm_filtered, 1, var)
top500 <- names(vars[order(vars, decreasing = T)][1:500])

# Plot the heatmap
pheatmap(tpm_filtered[top500,], scale = 'row', show_rownames = FALSE)

# Transpose the matrix and transform the counts to log2 scale, compute PCA
mtrans <- t(tpm_filtered[top500,])
mtrans <- log2(mtrans + 1)
pca <- prcomp(mtrans)

# Plot and print the PCA results
autoplot(pca, data = colData, colour = 'group', shape = 'IAA')
summary(pca)

# Calculate and plot the correlation matrix
correlation <- cor(tpm_filtered)
corrplot(correlation, type = 'lower', method = 'square', order = 'hclust',
         tl.col = 'black', cl.ratio = 0.2, addCoef.col = 'black', tl.srt = 45,
         col = COL2('RdBu', 10))

# Plot the heatmap of correlation matrix
pheatmap(correlation,  
         annotation_col = colData, 
         cutree_cols = 2)

# Construct DEseq dataset object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = colData, 
                              design = ~ group)

# Filter out genes with no reads
dds <- dds[rowSums(DESeq2::counts(dds)) >= 1,]

# Run the differential expression pipeline
dds <- DESeq(dds)

# Build the results table for the condition of interest
# Plot MA
res_def <- results(dds, contrast=c('group', 'no_IAA_control', 'WT_control'))
summary(res_def)
plotMA(res_def, ylim = c(-5,5))

# Build the results table for the condition of interest, lower alpha to 0.05
# Plot MA
res_p.05 <- results(dds, contrast=c('group', 'no_IAA_control', 'WT_control'), 
               alpha = 0.05)
summary(res_p.05)
plotMA(res_p.05, ylim = c(-5,5))

# Build the results table for the condition of interest, lower alpha to 0.05 and increase lfcthreshold to 1
# Plot MA
res_lfc1 <- results(dds, contrast=c('group', 'no_IAA_control', 'WT_control'), 
               alpha = 0.05, lfcThreshold = 1)
summary(res_lfc1)
plotMA(res_lfc1, ylim = c(-5,5))

# Build the results table for the next condition of interest with different parameters
# Plot MAs
res_def <- results(dds, contrast=c('group', 'CTCF_depletion', 'no_IAA_control'))
summary(res_def)
plotMA(res_def, ylim = c(-5,5))

res_p.05 <- results(dds, contrast=c('group', 'CTCF_depletion', 'no_IAA_control'), 
                    alpha = 0.05)
summary(res_p.05)
plotMA(res_p.05, ylim = c(-5,5))

res_lfc1 <- results(dds, contrast=c('group', 'CTCF_depletion', 'no_IAA_control'), 
                    alpha = 0.05, lfcThreshold = 1)
summary(res_lfc1)
plotMA(res_lfc1, ylim = c(-5,5))

# Check a single gene expression in the samples
genecounts <- plotCounts(dds, gene = 'ENSMUST00000005841.15',
                         intgroup = c('group', 'condition'),
                         returnData = TRUE)
geneplot <- ggplot(genecounts, aes(x = condition, y = count, color = group)) + geom_beeswarm(cex = 3) + geom_point(size = 2)
geneplot + stat_compare_means() 
geneplot + ylim(0, 7000)



