
library(tidyverse)
library(Biostrings)
library(ShortRead)
library(Rsubread)
library(RNASeqR)
library(stringr)
library(DESeq2)
library(ggplot2)

directorio <- "C:\\Users\\andoni\\OneDrive\\Escritorio\\TFM\\fastq_tfm\\RNAseq CD4.CD8 jan2023\\MUSehmzR\\soapnuke\\clean"

subcarpetas <- list.files(path = directorio, full.names = TRUE)

#################################### ALIGNMENT ####################################

setwd("C:\\Users\\andoni\\OneDrive\\Escritorio\\TFM")
buildindex(basename = "ncbi_m39", 
           reference= "GCF_000001635.27_GRCm39_genomic.fna", 
           gappedIndex = F,
           indexSplit = T, 
           memory=10000)

for (i in 1:length(subcarpetas)) {
  archivos <- list.files(path = subcarpetas[i+1], pattern = ".fq.gz", full.names = TRUE)
  output_prefix <- sub("_1\\.fq\\.gz$", "", basename(archivos[1]))
  align(index = "ncbi_m39",
        readfile1 = archivos[1],
        readfile2 = archivos[2],
        output_format = "BAM",
        output_file = paste0(subcarpetas[i+1], "/",output_prefix,".bam"), 
        nthreads = 10)
}

################################### COUNTS MATRIX ###################################

setwd("C:\\Users\\andoni\\OneDrive\\Escritorio\\TFM")

archivos_bam <- list.files("./bam_files/", pattern = ".bam$")
archivo_gtf_anot <- list.files(pattern = ".gtf.*$")

fc_SE <- featureCounts(files=paste0("./bam_files/",archivos_bam), annot.ext= archivo_gtf_anot,
                        isGTFAnnotation = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",
                        isPairedEnd = TRUE, requireBothEndsMapped = TRUE, nthreads = 10, 
                        strandSpecific = 0, countMultiMappingReads = TRUE, useMetaFeatures = TRUE)

write.table(fc_SE$counts,file="./raw_read_counts.txt",sep="\t", quote=F, append=F)

################################### DATA FILTERING ###################################

setwd("C:\\Users\\andoni\\OneDrive\\Escritorio\\TFM")

# Import raw counts' table

raw_counts <- read.delim("./raw_read_counts.txt", stringsAsFactors = FALSE)
colnames(raw_counts) <- c("P47", "P48", "P49", "P50", "P51", "P52", "P71", "P72", "P73",
                          "P74", "P75", "P76", "P83", "P84", "P85", "P86", "P87", "P88")

# Convert it into a data frame and check duplicates

class(raw_counts)
raw_counts <- as.data.frame(raw_counts)
rownames(raw_counts) <- toupper(rownames(raw_counts))

sum(duplicated(rownames(raw_counts))) # No duplicates

# Check for NAs

sum(is.na(raw_counts)) # No NAs

# Create condition table for DESeq2 analysis

cond_treatment <- c(rep("untreated",6), rep("treated",6), rep("control", 6))
cond_treatment <- factor(cond_treatment)

cond_cell <- rep(c("CD4","CD8"),9)
cond_cell <- factor(cond_cell)

cond_all <- data.frame(samples = colnames(raw_counts), cell = cond_cell,
                         treatment = cond_treatment)

cond_CD4 <- cond_all[cond_all$cell == "CD4", ]
cond_CD8 <- cond_all[cond_all$cell == "CD8", ]

# Remove IG genes to avoid introducing noise or bias into the analysis

raw_counts <- raw_counts[!str_detect(rownames(raw_counts), "IG[A H K L]{1}[D E G V M]{0,1}\\w+"),]

# Separate TCR genes from the rest of the genes

TCR_genes <- raw_counts[str_detect(rownames(raw_counts), "^TR[A B D G][J D V C]\\w+"),]
raw_counts <- subset(raw_counts, !(row.names(raw_counts) %in% row.names(TCR_genes))) # rest of the genes

# Separate samples into CD4 and CD8 cell types

CD4_counts <- raw_counts[, cond_all$cell == "CD4"]
CD8_counts <- raw_counts[, cond_all$cell == "CD8"]

# Build the DESeq2 object for each group of genes

dds_rest <- DESeqDataSetFromMatrix(countData = as.matrix(raw_counts),
                                   colData = cond_all,
                                   design = ~ cell + treatment)

dds_CD4 <- DESeqDataSetFromMatrix(countData = as.matrix(CD4_counts),
                                  colData = cond_CD4,
                                  design = ~ treatment)

dds_CD8 <- DESeqDataSetFromMatrix(countData = as.matrix(CD8_counts),
                                  colData = cond_CD8,
                                  design = ~ treatment)

# Select only rows with more than 10 counts

dds_rest_filt <- dds_rest[rowSums(counts(dds_rest)) >= 10,] 

dds_CD4_filt <- dds_CD4[rowSums(counts(dds_CD4)) >= 10,]
dds_CD8_filt <- dds_CD8[rowSums(counts(dds_CD8)) >= 10,]

######################## DESeq2 ANALYSIS REST OF THE GENES ######################## 

dds_DGE <- DESeq(dds_rest_filt)
resultsNames(dds_DGE)

dds_DGE_CD4 <- DESeq(dds_CD4_filt)
dds_DGE_CD8 <- DESeq(dds_CD8_filt)

################################### NORMALIZATION ###################################

# Sequencing depth normalization

dds_raw <- estimateSizeFactors(dds_rest_filt)
sizeFactors(dds_raw)

counts_normalized  <- counts(dds_raw, normalized = TRUE) # Sum of the total depth after normalization 
colSums(counts(dds_raw, normalized = TRUE))

# Log2 transformation of read counts

counts_log_normalized <- log2(counts_normalized + 1) # Adding 1 pseudocount

par(mfrow=c(1,2))

boxplot(counts_normalized , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "Untransformed  read  counts", ylab = "Read  counts")

boxplot(counts_log_normalized , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "Log2-transformed  read  counts", ylab = "log2(read  counts)")

library(hexbin)
library(vsn)

#--------------------------------------------------------------------------------------------------

par(mfrow=c(1,2))

SdPlot <- meanSdPlot(counts_normalized, ranks = FALSE, plot = FALSE, lwd = 4)  
SdPlot$gg + ggtitle("Sequencing depth normalized") + ylab("Standard deviation") 

SdPlot_log <- meanSdPlot(counts_log_normalized, ranks = FALSE, plot = FALSE, lwd = 15)
SdPlot_log$gg + ggtitle("Sequencing depth normalized log2(read counts)") + ylab("Standard deviation")

#--------------------------------------------------------------------------------------------------

# Variance stabilizing transformation

vsd <- vst(dds_raw, blind = FALSE)
head(assay(vsd), 3)

# The rlog transformation

rld <- rlog(dds_raw, blind = FALSE)
head(assay(rld), 3)

par(mfrow=c(2,2))

boxplot(counts_normalized , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "Untransformed  read  counts", ylab = "Read  counts")

boxplot(counts_log_normalized , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "Log2-transformed  read  counts", ylab = "log2(read  counts)")

boxplot(assay(vsd) , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "Variance stabilizing transformation of  read  counts", ylab = "Variance of read counts")

boxplot(assay(rld) , notch = TRUE , las = 2, cex.axis = 0.7,
        main = "rlog-transformed  read  counts", ylab = "rlog(read  counts)")

######################## DESeq2 ANALYSIS REST OF THE GENES ######################## 

# Apply a regularized log transformation and obtain a PCA plot

rld <- rlog(dds_DGE, blind = F) 
plotPCA(rld, intgroup = c("cell", "treatment"))

# Extract the results from the comparatives

dds_DGE_treated_vs_control_CD4 <- results(dds_DGE_CD4, contrast = c("treatment", "treated", "control"))
dds_DGE_untreated_vs_control_CD4 <- results(dds_DGE_CD4, contrast = c("treatment", "untreated", "control"))
dds_DGE_untreated_vs_treated_CD4 <- results(dds_DGE_CD4, contrast = c("treatment", "untreated", "treated"))

dds_DGE_treated_vs_control_CD8 <- results(dds_DGE_CD8, contrast = c("treatment", "treated", "control"))
dds_DGE_untreated_vs_control_CD8 <- results(dds_DGE_CD8, contrast = c("treatment", "untreated", "control"))
dds_DGE_untreated_vs_treated_CD8 <- results(dds_DGE_CD8, contrast = c("treatment", "untreated", "treated"))

# Obtain the complete cases

dds_DGE_treated_vs_control_CD4 <- dds_DGE_treated_vs_control_CD4[complete.cases(dds_DGE_treated_vs_control_CD4), ]
dds_DGE_untreated_vs_control_CD4 <- dds_DGE_untreated_vs_control_CD4[complete.cases(dds_DGE_untreated_vs_control_CD4), ]
dds_DGE_untreated_vs_treated_CD4 <- dds_DGE_untreated_vs_treated_CD4[complete.cases(dds_DGE_untreated_vs_treated_CD4), ]

dds_DGE_treated_vs_control_CD8 <- dds_DGE_treated_vs_control_CD8[complete.cases(dds_DGE_treated_vs_control_CD8), ]
dds_DGE_untreated_vs_control_CD8 <- dds_DGE_untreated_vs_control_CD8[complete.cases(dds_DGE_untreated_vs_control_CD8), ]
dds_DGE_untreated_vs_treated_CD8 <- dds_DGE_untreated_vs_treated_CD8[complete.cases(dds_DGE_untreated_vs_treated_CD8), ]

# Plot histogram of p-values frequencies 

par(mfrow=c(2,3))

hist(dds_DGE_treated_vs_control_CD4$pvalue , col = "blue",  xlab = "", border = "white",
     ylab = "Frequency", breaks = 0:40/40, main = "Treated vs control · CD4")
hist(dds_DGE_untreated_vs_control_CD4$pvalue , col = "blue",  xlab = "", border = "white", 
     ylab = "Frequency", breaks = 0:40/40, main = "Untreated vs control · CD4")
hist(dds_DGE_untreated_vs_treated_CD4$pvalue , col = "blue",  xlab = "", border = "white",
     ylab = "Frequency", breaks = 0:40/40, main = "Untreated vs treated · CD4")

hist(dds_DGE_treated_vs_control_CD8$pvalue , col = "blue",  xlab = "", border = "white",
     ylab = "Frequency", breaks = 0:40/40, main = "Treated vs control · CD8")
hist(dds_DGE_untreated_vs_control_CD8$pvalue , col = "blue",  xlab = "", border = "white", 
     ylab = "Frequency", breaks = 0:40/40, main = "Untreated vs control · CD8")
hist(dds_DGE_untreated_vs_treated_CD8$pvalue , col = "blue",  xlab = "", border = "white",
     ylab = "Frequency", breaks = 0:40/40, main = "Untreated vs treated · CD8")

# Plot MA-plots

par(mfrow=c(2,3))

plotMA(dds_DGE_treated_vs_control_CD4, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Treated vs control · CD4")
plotMA(dds_DGE_untreated_vs_control_CD4, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Untreated vs control · CD4")
plotMA(dds_DGE_untreated_vs_treated_CD4, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Untreated vs treated · CD4")

plotMA(dds_DGE_treated_vs_control_CD8, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Treated vs control · CD8")
plotMA(dds_DGE_untreated_vs_control_CD8, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Untreated vs control · CD8")
plotMA(dds_DGE_untreated_vs_treated_CD8, alpha = 0.05, ylim = c(-7,7), colSig = "limegreen", main = "Untreated vs treated · CD8")

# Filter according to log2 fold change

filt_treated_vs_control_CD4 <- subset(dds_DGE_treated_vs_control_CD4, abs(log2FoldChange) >= 1)
filt_treated_vs_control_CD4 <- filt_treated_vs_control_CD4[order(filt_treated_vs_control_CD4$padj), ]
filt_untreated_vs_control_CD4 <- subset(dds_DGE_untreated_vs_control_CD4, abs(log2FoldChange) >= 1)
filt_untreated_vs_control_CD4 <- filt_untreated_vs_control_CD4[order(filt_untreated_vs_control_CD4$padj), ]
filt_untreated_vs_treated_CD4 <- subset(dds_DGE_untreated_vs_treated_CD4, abs(log2FoldChange) >= 1)
filt_untreated_vs_treated_CD4 <- filt_untreated_vs_treated_CD4[order(filt_untreated_vs_treated_CD4$padj), ]

filt_treated_vs_control_CD8 <- subset(dds_DGE_treated_vs_control_CD8, abs(log2FoldChange) >= 1)
filt_treated_vs_control_CD8 <- filt_treated_vs_control_CD8[order(filt_treated_vs_control_CD8$padj), ]
filt_untreated_vs_control_CD8 <- subset(dds_DGE_untreated_vs_control_CD8, abs(log2FoldChange) >= 1)
filt_untreated_vs_control_CD8 <- filt_untreated_vs_control_CD8[order(filt_untreated_vs_control_CD8$padj), ]
filt_untreated_vs_treated_CD8 <- subset(dds_DGE_untreated_vs_treated_CD8, abs(log2FoldChange) >= 1)
filt_untreated_vs_treated_CD8 <- filt_untreated_vs_treated_CD8[order(filt_untreated_vs_treated_CD8$padj), ]

# Plot volcano plots 

library(ggplot2)
library(ggrepel)
library(ggpubr)

volcano_1 <- as.data.frame(filt_treated_vs_control_CD4) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Treated vs control · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

volcano_2 <- as.data.frame(filt_untreated_vs_control_CD4) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Untreated vs control · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

volcano_3 <- as.data.frame(filt_untreated_vs_treated_CD4) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Untreated vs treated · CD4")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

volcano_4 <- as.data.frame(filt_treated_vs_control_CD8) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Treated vs control · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

volcano_5 <- as.data.frame(filt_untreated_vs_control_CD8) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Untreated vs control · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

volcano_6 <- as.data.frame(filt_untreated_vs_treated_CD8) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue"), alpha = 0.2))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Unreated vs treated · CD8")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(alpha = "none")

ggarrange(volcano_1 + geom_text_repel(data = as.data.frame(filt_treated_vs_control_CD4[1:20, ]), 
                                      aes(label = rownames(filt_treated_vs_control_CD4)[1:20]), max.overlaps = 60, size = 2.5), 
          volcano_2 + geom_text_repel(data = as.data.frame(filt_untreated_vs_control_CD4[1:20, ]), 
                                      aes(label = rownames(filt_untreated_vs_control_CD4)[1:20]), max.overlaps = 60, size = 2.5), 
          volcano_3 + geom_text_repel(data = as.data.frame(filt_untreated_vs_treated_CD4[1:20, ]), 
                                      aes(label = rownames(filt_untreated_vs_treated_CD4)[1:20]), max.overlaps = 60, size = 2.5),
          volcano_4 + geom_text_repel(data = as.data.frame(filt_treated_vs_control_CD8[1:20, ]), 
                                      aes(label = rownames(filt_treated_vs_control_CD8)[1:20]), max.overlaps = 60, size = 2.5), 
          volcano_5 + geom_text_repel(data = as.data.frame(filt_untreated_vs_control_CD8[1:20, ]), 
                                      aes(label = rownames(filt_untreated_vs_control_CD8)[1:20]), max.overlaps = 60, size = 2.5), 
          volcano_6 + geom_text_repel(data = as.data.frame(filt_untreated_vs_treated_CD8[1:20, ]), 
                                      aes(label = rownames(filt_untreated_vs_treated_CD8)[1:20]), max.overlaps = 60, size = 2.5),
          nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom", labels = c("A", "B", "C", "D", "E", "G"))

################################### GENE LISTS OBTENTION #################################

# CD4 gene list obtention

filt_treated_vs_control_CD4 <- data.frame(filt_treated_vs_control_CD4[,c("log2FoldChange", "padj")])
filt_treated_vs_control_CD4$gene <- rownames(filt_treated_vs_control_CD4)
rownames(filt_treated_vs_control_CD4) <- NULL

filt_untreated_vs_control_CD4 <- data.frame(filt_untreated_vs_control_CD4[,c("log2FoldChange", "padj")])
filt_untreated_vs_control_CD4$gene <- rownames(filt_untreated_vs_control_CD4)
rownames(filt_untreated_vs_control_CD4) <- NULL

filt_untreated_vs_treated_CD4 <- data.frame(filt_untreated_vs_treated_CD4[,c("log2FoldChange", "padj")])
filt_untreated_vs_treated_CD4$gene <- rownames(filt_untreated_vs_treated_CD4)
rownames(filt_untreated_vs_treated_CD4) <- NULL

CD4_genes <- full_join(filt_treated_vs_control_CD4, filt_untreated_vs_control_CD4, by = "gene")
CD4_genes <- CD4_genes %>% relocate(gene)
CD4_genes <- full_join(CD4_genes, filt_untreated_vs_treated_CD4, by = "gene")
colnames(CD4_genes) <- c("gene", "lfc_treated_vs_control", "padj_treated_vs_control", 
                         "lfc_untreated_vs_control", "padj_untreated_vs_control",
                         "lfc_untreated_vs_treated", "padj_untreated_vs_treated")

# Select genes with a padj < 0.01 and a |lfc| > 1 in at least one of the comparisons

CD4_genes <- CD4_genes[rowSums(CD4_genes[, c("padj_treated_vs_control", "padj_untreated_vs_control", "padj_untreated_vs_treated")] < 0.01, na.rm = T) > 0, ]
CD4_genes <- CD4_genes[rowSums(abs(CD4_genes[, c("lfc_treated_vs_control", "lfc_untreated_vs_control", "lfc_untreated_vs_treated")]) > 1, na.rm = T) > 0, ]

list_untreated_vs_treated_CD4 <- CD4_genes[, c("gene", "lfc_untreated_vs_treated", "padj_untreated_vs_treated")]
list_untreated_vs_treated_CD4 <- list_untreated_vs_treated_CD4[complete.cases(list_untreated_vs_treated_CD4),]
up_list_untreated_vs_treated_CD4 <- list_untreated_vs_treated_CD4[list_untreated_vs_treated_CD4$lfc_untreated_vs_treated > 1,] # Upregulated genes
down_list_untreated_vs_treated_CD4 <- list_untreated_vs_treated_CD4[list_untreated_vs_treated_CD4$lfc_untreated_vs_treated < 1,] # Downregulated genes

# CD8 gene list obtention

filt_treated_vs_control_CD8 <- data.frame(filt_treated_vs_control_CD8[,c("log2FoldChange", "padj")])
filt_treated_vs_control_CD8$gene <- rownames(filt_treated_vs_control_CD8)
rownames(filt_treated_vs_control_CD8) <- NULL

filt_untreated_vs_control_CD8 <- data.frame(filt_untreated_vs_control_CD8[,c("log2FoldChange", "padj")])
filt_untreated_vs_control_CD8$gene <- rownames(filt_untreated_vs_control_CD8)
rownames(filt_untreated_vs_control_CD8) <- NULL

filt_untreated_vs_treated_CD8 <- data.frame(filt_untreated_vs_treated_CD8[,c("log2FoldChange", "padj")])
filt_untreated_vs_treated_CD8$gene <- rownames(filt_untreated_vs_treated_CD8)
rownames(filt_untreated_vs_treated_CD8) <- NULL

CD8_genes <- full_join(filt_treated_vs_control_CD8, filt_untreated_vs_control_CD8, by = "gene")
CD8_genes <- CD8_genes %>% relocate(gene)
CD8_genes <- full_join(CD8_genes, filt_untreated_vs_treated_CD8, by = "gene")
colnames(CD8_genes) <- c("gene", "lfc_treated_vs_control", "padj_treated_vs_control", 
                         "lfc_untreated_vs_control", "padj_untreated_vs_control",
                         "lfc_untreated_vs_treated", "padj_untreated_vs_treated")

# Select genes with a padj < 0.01 and a |lfc| > 1 in at least one of the comparisons

CD8_genes <- CD8_genes[rowSums(CD8_genes[, c("padj_treated_vs_control", "padj_untreated_vs_control", "padj_untreated_vs_treated")] < 0.01, na.rm = T) > 0, ]
CD8_genes <- CD8_genes[rowSums(abs(CD8_genes[, c("lfc_treated_vs_control", "lfc_untreated_vs_control", "lfc_untreated_vs_treated")]) > 1, na.rm = T) > 0, ]

list_untreated_vs_treated_CD8 <- CD8_genes[, c("gene", "lfc_untreated_vs_treated", "padj_untreated_vs_treated")]
list_untreated_vs_treated_CD8 <- list_untreated_vs_treated_CD8[complete.cases(list_untreated_vs_treated_CD8),]
up_list_untreated_vs_treated_CD8 <- list_untreated_vs_treated_CD8[list_untreated_vs_treated_CD8$lfc_untreated_vs_treated > 1,] # Upregulated genes
down_list_untreated_vs_treated_CD8 <- list_untreated_vs_treated_CD8[list_untreated_vs_treated_CD8$lfc_untreated_vs_treated < 1,] # Downregulated genes

# Plot heatmaps with all genes significant in at least one comparison

library(pheatmap)

par(mfrow=c(1,2))

dds_CD4_filt <- dds_CD4_filt[rownames(assay(dds_CD4_filt)) %in% CD4_genes$gene]
rld_CD4 <- rlog(dds_CD4_filt, blind = FALSE)

topVarGenes_CD4 <- order(rowVars(assay(rld_CD4)), decreasing = TRUE)
mat_CD4 <- assay(rld_CD4)[topVarGenes_CD4, ]
mat_CD4 <- mat_CD4 - rowMeans(mat_CD4)

anno_CD4 <- as.data.frame(colData(rld_CD4)[, c("treatment", "cell")])

pheatmap(assay(rld_CD4), clustering_distance_rows = "correlation", clustering_distance_cols = "manhattan", 
         annotation_col = anno_CD4, cellwidth = 20, cellheight = 0.15, scale = "row", 
         fontsize_col = 7, show_rownames = F, cluster_rows = T, cluster_cols = T)

dds_CD8_filt <- dds_CD8_filt[rownames(assay(dds_CD8_filt)) %in% CD8_genes$gene]
rld_CD8 <- rlog(dds_CD8_filt, blind = FALSE)

topVarGenes_CD8 <- order(rowVars(assay(rld_CD8)), decreasing = TRUE)
mat_CD8 <- assay(rld_CD8)[topVarGenes_CD8, ]
mat_CD8 <- mat_CD8 - rowMeans(mat_CD8)

anno_CD8 <- as.data.frame(colData(rld_CD8)[, c("treatment", "cell")])

pheatmap(assay(rld_CD8), clustering_distance_rows = "correlation", clustering_distance_cols = "manhattan", 
         annotation_col = anno_CD8, cellwidth = 20, cellheight = 0.1, scale = "row", 
         fontsize_col = 7, show_rownames = F, cluster_rows = T, cluster_cols = T)

################################### RESULTS ANNOTATION #################################

library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(org.Mm.eg.db)
library(biomaRt)

ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

#----------------------------- Add ensembl gene ID and gene name columns -----------------------------

# CD4 samples

anno_untreated_vs_treated_up_CD4 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                          filters = "external_gene_name",
                                          values = up_list_untreated_vs_treated_CD4$gene,
                                          mart = ensembl)

anno_untreated_vs_treated_down_CD4 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                            filters = "external_gene_name",
                                            values = down_list_untreated_vs_treated_CD4$gene,
                                            mart = ensembl)

# CD8 samples

anno_untreated_vs_treated_up_CD8 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                          filters = "external_gene_name",
                                          values = up_list_untreated_vs_treated_CD8$gene,
                                          mart = ensembl)

anno_untreated_vs_treated_down_CD8 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                            filters = "external_gene_name",
                                            values = down_list_untreated_vs_treated_CD8$gene,
                                            mart = ensembl)

#---------------------------------- Add entrez ID column ---------------------------------- 

# CD4 samples

anno_untreated_vs_treated_up_CD4$entrez_id <- mapIds(org.Mm.eg.db,
                                                    keys = anno_untreated_vs_treated_up_CD4$ensembl_gene_id, 
                                                    column = "ENTREZID",
                                                    keytype = "ENSEMBL",
                                                    multiVals = "first")

anno_untreated_vs_treated_down_CD4$entrez_id <- mapIds(org.Mm.eg.db,
                                                     keys = anno_untreated_vs_treated_down_CD4$ensembl_gene_id, 
                                                     column = "ENTREZID",
                                                     keytype = "ENSEMBL",
                                                     multiVals = "first")

# CD8 samples

anno_untreated_vs_treated_up_CD8$entrez_id <- mapIds(org.Mm.eg.db,
                                                    keys = anno_untreated_vs_treated_up_CD8$ensembl_gene_id, 
                                                    column = "ENTREZID",
                                                    keytype = "ENSEMBL",
                                                    multiVals = "first")

anno_untreated_vs_treated_down_CD8$entrez_id <- mapIds(org.Mm.eg.db,
                                                      keys = anno_untreated_vs_treated_down_CD8$ensembl_gene_id, 
                                                      column = "ENTREZID",
                                                      keytype = "ENSEMBL",
                                                      multiVals = "first")

###################################### ENRICHMENT ANALYSIS ###################################### 

library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(msigdbr)
library(stringr)

# ----------------------------------- msig database and ORA ----------------------------------- 

species <- "Mus musculus"
msig <- list()

msig[["H"]] <- msigdbr(species = species, category = "H")
msig[["C2.KEGG"]] <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG")
msig[["C2.react"]] <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
msig[["GO.BP"]] <- msigdbr(species = species, category = "C5", subcategory = "GO:BP")
msig[["GO.MF"]] <- msigdbr(species = species, category = "C5", subcategory = "GO:MF")

msig.dfs <- list()
for (i in names(msig)){
  msig.dfs[[i]] <- msig[[i]] %>% select(gs_name, gene_symbol, entrez_gene, ensembl_gene) %>% as.data.frame()}

titles <- c("MSIG Hallmarks", "KEGG Pathways", "Reactome Pathways",
            "Biological Processes", "Molecular Functions")

# CD4 cells

# Overrepresented in untreated mice

msig.untreated_vs_treated_up_CD4 <- list()
for (i in names(msig.dfs)) {
  msig.untreated_vs_treated_up_CD4[[i]] <- enricher(gene = anno_untreated_vs_treated_up_CD4$external_gene_name,
                                         TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)}

ora.plots.untreated_vs_treated_up_CD4 <- list()
for (n in  1:length(msig.untreated_vs_treated_up_CD4))  {
  ora.plots.untreated_vs_treated_up_CD4[[n]] <- dotplot(msig.untreated_vs_treated_up_CD4[[n]], 
                                                     showCategory = length(msig.untreated_vs_treated_up_CD4[[n]]$Description),
                                                     title = titles[n])}

ora.plots.untreated_vs_treated_up_CD4[[1]] <- dotplot(msig.untreated_vs_treated_up_CD4[[1]], title = titles[1], font.size = 7)
ora.plots.untreated_vs_treated_up_CD4[[1]]
ora.plots.untreated_vs_treated_up_CD4[[2]] <- dotplot(msig.untreated_vs_treated_up_CD4[[2]], title = titles[2], font.size = 7)
ora.plots.untreated_vs_treated_up_CD4[[2]]
ora.plots.untreated_vs_treated_up_CD4[[3]] <- dotplot(msig.untreated_vs_treated_up_CD4[[3]], title = titles[3], font.size = 7)
ora.plots.untreated_vs_treated_up_CD4[[3]]
ora.plots.untreated_vs_treated_up_CD4[[4]] <- dotplot(msig.untreated_vs_treated_up_CD4[[4]], title = titles[4], font.size = 7)
ora.plots.untreated_vs_treated_up_CD4[[4]]
ora.plots.untreated_vs_treated_up_CD4[[5]] <- dotplot(msig.untreated_vs_treated_up_CD4[[5]], title = titles[5], font.size = 7)
ora.plots.untreated_vs_treated_up_CD4[[5]]

# Overrepresented in treated mice

msig.untreated_vs_treated_down_CD4 <- list()
for (i in names(msig.dfs)) {
  msig.untreated_vs_treated_down_CD4[[i]] <- enricher(gene = anno_untreated_vs_treated_down_CD4$external_gene_name,
                                                    TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)}

ora.plots.untreated_vs_treated_down_CD4 <- list()
for (n in  1:length(msig.untreated_vs_treated_down_CD4))  {
  ora.plots.untreated_vs_treated_down_CD4[[n]] <- dotplot(msig.untreated_vs_treated_down_CD4[[n]], 
                                                        showCategory = length(msig.untreated_vs_treated_down_CD4[[n]]$Description),
                                                        title = titles[n])}

ora.plots.untreated_vs_treated_down_CD4[[1]] <- dotplot(msig.untreated_vs_treated_down_CD4[[1]], title = titles[1], font.size = 7)
ora.plots.untreated_vs_treated_down_CD4[[1]]
ora.plots.untreated_vs_treated_down_CD4[[2]] <- dotplot(msig.untreated_vs_treated_down_CD4[[2]], title = titles[2], font.size = 7)
ora.plots.untreated_vs_treated_down_CD4[[2]]
ora.plots.untreated_vs_treated_down_CD4[[3]] <- dotplot(msig.untreated_vs_treated_down_CD4[[3]], title = titles[3], font.size = 7)
ora.plots.untreated_vs_treated_down_CD4[[3]]
ora.plots.untreated_vs_treated_down_CD4[[4]] <- dotplot(msig.untreated_vs_treated_down_CD4[[4]], title = titles[4], font.size = 7)
ora.plots.untreated_vs_treated_down_CD4[[4]]
ora.plots.untreated_vs_treated_down_CD4[[5]] <- dotplot(msig.untreated_vs_treated_down_CD4[[5]], title = titles[5], font.size = 7)
ora.plots.untreated_vs_treated_down_CD4[[5]]


# CD8 cells

# Overrepresented in untreated mice

msig.untreated_vs_treated_up_CD8 <- list()
for (i in names(msig.dfs)) {
  msig.untreated_vs_treated_up_CD8[[i]] <- enricher(gene = anno_untreated_vs_treated_up_CD8$external_gene_name,
                                                    TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)}

ora.plots.untreated_vs_treated_up_CD8 <- list()
for (n in  1:length(msig.untreated_vs_treated_up_CD8))  {
  ora.plots.untreated_vs_treated_up_CD8[[n]] <- dotplot(msig.untreated_vs_treated_up_CD8[[n]], 
                                                        showCategory = length(msig.untreated_vs_treated_up_CD8[[n]]$Description),
                                                        title = titles[n])}

ora.plots.untreated_vs_treated_up_CD8[[1]] <- dotplot(msig.untreated_vs_treated_up_CD8[[1]], title = titles[1], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_up_CD8[[1]]
ora.plots.untreated_vs_treated_up_CD8[[2]] <- dotplot(msig.untreated_vs_treated_up_CD8[[2]], title = titles[2], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_up_CD8[[2]]
ora.plots.untreated_vs_treated_up_CD8[[3]] <- dotplot(msig.untreated_vs_treated_up_CD8[[3]], title = titles[3], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_up_CD8[[3]]
ora.plots.untreated_vs_treated_up_CD8[[4]] <- dotplot(msig.untreated_vs_treated_up_CD8[[4]], title = titles[4], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_up_CD8[[4]]
ora.plots.untreated_vs_treated_up_CD8[[5]] <- dotplot(msig.untreated_vs_treated_up_CD8[[5]], title = titles[5], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_up_CD8[[5]]

# Overrepresented in treated mice

msig.untreated_vs_treated_down_CD8 <- list()
for (i in names(msig.dfs)) {
  msig.untreated_vs_treated_down_CD8[[i]] <- enricher(gene = anno_untreated_vs_treated_down_CD8$external_gene_name,
                                                      TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)}

ora.plots.untreated_vs_treated_down_CD8 <- list()
for (n in  1:length(msig.untreated_vs_treated_down_CD8))  {
  ora.plots.untreated_vs_treated_down_CD8[[n]] <- dotplot(msig.untreated_vs_treated_down_CD8[[n]], 
                                                          showCategory = length(msig.untreated_vs_treated_down_CD8[[n]]$Description),
                                                          title = titles[n])}

ora.plots.untreated_vs_treated_down_CD8[[1]] <- dotplot(msig.untreated_vs_treated_down_CD8[[1]], title = titles[1], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_down_CD8[[1]]
ora.plots.untreated_vs_treated_down_CD8[[2]] <- dotplot(msig.untreated_vs_treated_down_CD8[[2]], title = titles[2], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_down_CD8[[2]]
ora.plots.untreated_vs_treated_down_CD8[[3]] <- dotplot(msig.untreated_vs_treated_down_CD8[[3]], title = titles[3], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_down_CD8[[3]]
ora.plots.untreated_vs_treated_down_CD8[[4]] <- dotplot(msig.untreated_vs_treated_down_CD8[[4]], title = titles[4], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_down_CD8[[4]]
ora.plots.untreated_vs_treated_down_CD8[[5]] <- dotplot(msig.untreated_vs_treated_down_CD8[[5]], title = titles[5], font.size = 6, width = 0.3)
ora.plots.untreated_vs_treated_down_CD8[[5]]


# Hallmarks graph

ggarrange(ora.plots.untreated_vs_treated_up_CD4[[1]], ora.plots.untreated_vs_treated_down_CD4[[1]], 
          ora.plots.untreated_vs_treated_up_CD8[[1]], ora.plots.untreated_vs_treated_down_CD8[[1]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

# KEGG pathways graph

ggarrange(ora.plots.untreated_vs_treated_up_CD4[[2]], ora.plots.untreated_vs_treated_down_CD4[[2]], 
          ora.plots.untreated_vs_treated_up_CD8[[2]], ora.plots.untreated_vs_treated_down_CD8[[2]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

# Reactome pathways graph

ggarrange(ora.plots.untreated_vs_treated_up_CD4[[3]], ora.plots.untreated_vs_treated_down_CD4[[3]], 
          ora.plots.untreated_vs_treated_up_CD8[[3]], ora.plots.untreated_vs_treated_down_CD8[[3]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))

# Biological processes graph

ggarrange(ora.plots.untreated_vs_treated_up_CD4[[4]], ora.plots.untreated_vs_treated_down_CD4[[4]], 
          ora.plots.untreated_vs_treated_up_CD8[[4]], ora.plots.untreated_vs_treated_down_CD8[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")) 

# Molecular functions graph

ggarrange(ora.plots.untreated_vs_treated_up_CD4[[5]], ora.plots.untreated_vs_treated_down_CD4[[5]], 
          ora.plots.untreated_vs_treated_up_CD8[[5]], ora.plots.untreated_vs_treated_down_CD8[[5]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")) 

# ---------------------------------------- GSEA  ---------------------------------------- 

library(cowplot)

# CD4 samples

lfc_untreated_vs_treated_CD4 <- list_untreated_vs_treated_CD4$lfc_untreated_vs_treated
names(lfc_untreated_vs_treated_CD4) <- unique(anno_untreated_vs_treated_CD4$external_gene_name)
lfc_untreated_vs_treated_CD4 <- sort(lfc_untreated_vs_treated_CD4, decreasing = T)

gsea.untreated_vs_treated_CD4 <- list()
for (i in names(msig.dfs)) {
  gsea.untreated_vs_treated_CD4[[i]] <- GSEA(lfc_untreated_vs_treated_CD4, TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)
}

gsea.plots_untreated_vs_treated_CD4 <- list()
for (n in names(gsea.untreated_vs_treated_CD4)){
  gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[n]] <- gseaplot2(gsea.untreated_vs_treated_CD4[[n]], 
                                                                       1:7, base_size = 7)}

gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[1]]
gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[2]]
gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[3]]
gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[4]]
gsea.plots_untreated_vs_treated_CD4[["gseaplot2"]][[5]]

# CD8 samples

lfc_untreated_vs_treated_CD8 <- list_untreated_vs_treated_CD8$lfc_untreated_vs_treated
names(lfc_untreated_vs_treated_CD8) <- unique(anno_untreated_vs_treated_CD8$external_gene_name)
lfc_untreated_vs_treated_CD8 <- sort(lfc_untreated_vs_treated_CD8, decreasing = T)

gsea.untreated_vs_treated_CD8 <- list()
for (i in names(msig.dfs)) {
  gsea.untreated_vs_treated_CD8[[i]] <- GSEA(lfc_untreated_vs_treated_CD8, TERM2GENE = msig.dfs[[i]], pvalueCutoff = 0.05)
}

gsea.plots_untreated_vs_treated_CD8 <- list()
for (n in names(gsea.untreated_vs_treated_CD8)){
  gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[n]] <- gseaplot2(gsea.untreated_vs_treated_CD8[[n]], 
                                                                       1:7, base_size = 7)}

gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[1]]
gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[2]]
gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[3]]
gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[4]]
gsea.plots_untreated_vs_treated_CD8[["gseaplot2"]][[5]]
