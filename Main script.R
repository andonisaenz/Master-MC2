
library(tidyverse)
library(Rsamtools)
library(Biostrings)
library(ShortRead)
library(Rsubread)
library(RNASeqR)
library(fastqcr)
library(stringr)
library(DESeq2)

directorio <- "C:\\Users\\andon\\OneDrive\\Escritorio\\TFM\\fastq_tfm\\RNAseq CD4.CD8 jan2023\\MUSehmzR\\soapnuke\\clean"

subcarpetas <- list.files(path = directorio, full.names = TRUE)

#################################### ALIGNMENT ####################################

setwd("C:\\Users\\andon\\OneDrive\\Escritorio\\TFM")
buildindex(basename = "ncbi_m39", reference= "GCF_000001635.27_GRCm39_genomic.fna", gappedIndex = F,
           indexSplit = T, memory=10000)

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

setwd("C:\\Users\\andon\\OneDrive\\Escritorio\\TFM")

archivos_bam <- list.files("./bam_files/", pattern = ".bam$")
archivo_gtf_anot <- list.files(pattern = ".gtf.*$")

fc_SE <- featureCounts(files=paste0("./bam_files/",archivos_bam), annot.ext= archivo_gtf_anot,
                        isGTFAnnotation = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",
                        isPairedEnd = TRUE, requireBothEndsMapped = TRUE, nthreads = 10, 
                        strandSpecific = 0, countMultiMappingReads = TRUE, useMetaFeatures = TRUE)

write.table(fc_SE$counts,file="./raw_read_counts.txt",sep="\t", quote=F, append=F)

################################### DATA FILTERING ###################################

setwd("C:\\Users\\andon\\OneDrive\\Escritorio\\TFM")

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
cond_treatment <- factor(cond)

cond_cell <- rep(c("CD4","CD8"),9)
cond_cell <- factor(cond_cell)

cond_table <- data.frame(samples = colnames(raw_counts), cell = cond_cell,
                         treatment = cond_treatment)

# Remove IG genes to avoid introducing noise or bias into the analysis

raw_counts <- raw_counts[!str_detect(rownames(raw_counts), "IG[A H K L]{1}[D E G V M]{0,1}\\w+"),]

# Remove predicted genes from the counts matrix

raw_counts <- raw_counts[!str_detect(rownames(raw_counts), "^GM\\d{5}"),]

# Separate TCR genes from the rest of the genes

TCR_genes <- raw_counts[str_detect(rownames(raw_counts), "^TR[A B D G][J D V C]\\w+"),]
TRAC <- raw_counts[str_detect(rownames(raw_counts), "TRAC$"),]
TRDC <- raw_counts[str_detect(rownames(raw_counts), "TRDC$"),]
TCR_genes <- rbind(TCR_genes, TRAC, TRDC)

raw_counts <- subset(raw_counts, !(row.names(raw_counts) %in% row.names(TCR_genes))) # rest of the genes

# Build the DESeq2 object for each group of genes

dds_TCR <- DESeqDataSetFromMatrix(countData = as.matrix(TCR_genes),
                                  colData = cond_table,
                                  design = ~ cell + treatment)

dds_rest <- DESeqDataSetFromMatrix(countData = as.matrix(raw_counts),
                                   colData = cond_table,
                                   design = ~ cell + treatment)

# Select only rows with more than 10 counts

dds_TCR_filt <- dds_TCR[rowSums(counts(dds_TCR)) >= 10,] 
dds_rest_filt <- dds_rest[rowSums(counts(dds_rest)) >= 10,] 

################################### NORMALIZATION ###################################

# Sequencing depth normalization

dds_raw <- estimateSizeFactors(dds_rest_filt)
sizeFactors(dds_raw)

counts_normalized  <- counts(dds_raw, normalized = TRUE) # Sum of the total depth after normalization 
colSums(counts(dds_raw, normalized=TRUE))

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

SdPlot <- meanSdPlot(counts_normalized, ranks = FALSE, plot = FALSE)  
SdPlot$gg + ggtitle("Sequencing depth normalized") + ylab("Standard deviation")

SdPlot_log <- meanSdPlot(counts_log_normalized, ranks = FALSE, plot = FALSE)
SdPlot_log$gg + ggtitle("Sequencing depth normalized log2(read counts)") + ylab("Standard deviation")

#--------------------------------------------------------------------------------------------------

# Variance stabilizing transformation

vsd <- vst(dds_raw, blind = FALSE)
head(assay(vsd), 3)

# The rlog transformation

rld <- rlog(dds_raw, blind = FALSE)
head(assay(rld), 3)

# --------------- Sample correlation using all the different data transformations ---------------

par(mfrow=c(2,2))

# Raw data normalized by sequencing depth
plot(counts_normalized [,1:2], cex=.1, main = "Normalized by sequencing depth")

# Log^2^ normalization adding 1 pseudocount
plot(counts_log_normalized [,1:2], cex=.1, main = "Normalized log2(read counts)")

# rlog transformed data
rlog_norm_counts <- assay(rld)
plot(rlog_norm_counts[,1:2], cex=.1, main = "rlog transformed", xlim=c(0,18), ylim=c(0,18))

# Variance stabilizing transformation
vsd_norm_counts <- assay(vsd)
plot(vsd_norm_counts[,1:2], cex=.1, main = "Variance stabilizing transformation", xlim=c(0,18), ylim=c(0,18))

######################## DESeq2 ANALYSIS REST OF THE GENES ######################## 

dds_DGE <- DESeq(dds_rest_filt)
resultsNames(dds_DGE)

rld <- rlog(dds_DGE, blind = F) # Regularized log transformation
plotPCA(rld, intgroup = c("cell", "treatment"))

dds_DGE_results <- results(dds_DGE)
head(dds_DGE_results)

filt_DGE_results <- subset(dds_DGE_results, padj < 0.01)
filt_DGE_results <- filt_DGE_results[order(filt_DGE_results$padj), ]

head(filt_DGE_results)

# Histogram of p-values frequencies 

hist(dds_DGE_results$pvalue , col = "blue",  xlab = "", border = "white", 
     ylab = "Frequency", breaks = 0:40/40, main = "Frequencies of p-values")

# Volcano plot

library(ggplot2)

as.data.frame(filt_DGE_results) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue")))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Most significant genes (control VS untreated VS treated)")+
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap

library(pheatmap)

vsd <- vst(dds_rest_filt, blind = FALSE)
head(assay(vsd), 3)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)

anno <- as.data.frame(colData(vsd)[, c("treatment", "cell")])
pheatmap(mat, annotation_col = anno, fontsize_row = 7, cellwidth = 20, cellheight = 7,
         fontsize_col = 7)

############################## DESeq2 ANALYSIS TCR GENES ############################## 

dds_DGE_TCR <- DESeq(dds_TCR_filt)
resultsNames(dds_DGE_TCR)

rld_TCR <- rlog(dds_DGE_TCR, blind = F) 
plotPCA(rld_TCR, intgroup = c("cell", "treatment"))

dds_DGE_TCR_results <- results(dds_DGE_TCR)
head(dds_DGE_TCR_results)

filt_DGE_TCR_results <- subset(dds_DGE_TCR_results, padj < 0.01)
filt_DGE_TCR_results <- filt_DGE_TCR_results[order(filt_DGE_TCR_results$padj), ]

head(filt_DGE_TCR_results)

# Histogram of p-values frequencies 

hist(dds_DGE_TCR_results$pvalue , col = "blue",  xlab = "", border = "white", 
     ylab = "Frequency", breaks = 0:40/40, main = "Frequencies of p-values")

# Volcano plot

library(ggplot2)

as.data.frame(filt_DGE_TCR_results) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = ifelse(log2FoldChange > 0, "Red", "Blue")))+
  scale_color_identity(name = "Gene regulation", labels = c("Downregulated", "Upregulated"),
                       guide = guide_legend(reverse = TRUE)) +
  ggtitle("Most significant genes (control VS untreated VS treated)")+
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap

library(pheatmap)

vsd_TCR <- varianceStabilizingTransformation(dds_TCR_filt)
head(assay(vsd_TCR), 3)

topVarGenes_TCR <- head(order(rowVars(assay(vsd_TCR)), decreasing = TRUE), 30)
mat_TCR <- assay(vsd_TCR)[topVarGenes_TCR, ]
mat_TCR <- mat_TCR - rowMeans(mat_TCR)

anno_TCR <- as.data.frame(colData(vsd_TCR)[, c("treatment", "cell")])
pheatmap(mat_TCR, annotation_col = anno_TCR, fontsize_row = 7, cellwidth = 20, cellheight = 7,
         fontsize_col = 7)

################################### RESULTS ANNOTATION #################################

library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(org.Mm.eg.db)
library(biomaRt)

ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
gene_names <- rownames(filt_DGE_results)

results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "external_gene_name",
                 values = gene_names,
                 mart = ensembl)

filt_DGE_results["ENSEMBL"] <- results[,1]

rownames(filt_DGE_results) <- results[,1]

# Add gene symbol column
filt_DGE_results$SYMBOL <- mapIds(org.Mm.eg.db,
                           keys = gene_names, 
                           column = "SYMBOL",
                           keytype = "SYMBOL",
                           multiVals = "first")

# Add gene name column
filt_DGE_results$GENENAME <- mapIds(org.Mm.eg.db,
                             keys = filt_DGE_results$ENSEMBL, 
                             column = "GENENAME",
                             keytype = "ENSEMBL",
                             multiVals = "first")

# Add functional path column
filt_DGE_results$PATH <- mapIds(org.Mm.eg.db,
                         keys = filt_DGE_results$ENSEMBL, 
                         column = "PATH",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# Add entrez ID column
filt_DGE_results$ENTREZID <- mapIds(org.Mm.eg.db,
                             keys = filt_DGE_results$ENSEMBL, 
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")

################################ CLUSTER PROFILER ################################ 

library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)

OrgDb <- org.Mm.eg.db
genes <- as.character(filt_DGE_results$ENTREZID)

MF_annotation <- enrichGO(gene = genes,
                          OrgDb = OrgDb,
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE)

barplot(MF_annotation, showCategory = 30, font.size = 7)

BP_annotation <- enrichGO(gene = genes,
                          OrgDb = OrgDb,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE)

barplot(BP_annotation, showCategory = 25, font.size = 7)

# Build a gene concept network and an enrichment map for MF GO terms

edox_MF <- setReadable(MF_annotation, 'org.Mm.eg.db', 'ENTREZID')

geneList <- filt_DGE_results$log2FoldChange
names(geneList) <- as.character(unique(filt_DGE_results$ENTREZID))
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

cnetplot(edox_MF, categorySize = "pvalue", foldChange = geneList, 
         node_label_size = 2, cex_label_category = 0.4, cex_label_gene = 0.4)

edox_MF <- pairwise_termsim(edox_MF)
emapplot(edox_MF, cex_label_category = 0.7)

# Build a gene concept network and an enrichment map for BP GO terms

edox_BP <- setReadable(BP_annotation, 'org.Mm.eg.db', 'ENTREZID')

cnetplot(edox_BP, categorySize = "pvalue", foldChange = geneList, 
         node_label_size = 2, cex_label_category = 0.4, cex_label_gene = 0.4)

edox_BP <- pairwise_termsim(edox_BP)
emapplot(edox_BP, cex_label_category = 0.5)

################################### KEGG PATHWAYS ################################### 

library(pathview)

kegg_genes <- enrichKEGG(gene = genes, organism = 'mmu', pvalueCutoff = 0.05)

browseKEGG(kegg_genes, 'mmu04110')
browseKEGG(kegg_genes, 'mmu05132')
browseKEGG(kegg_genes, 'mmu05171')

pathview(gene.data = geneList,
         pathway.id = "mmu04110",
         species = 'mmu',
         limit = list(gene=max(abs(geneList)), cpd=1))
