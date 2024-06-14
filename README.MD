# mRNA Analysis of Innate-like T Cells in Mice Treated with LiNKTcellsofαGalCer/CD1d-NP

#Introduction: This repository contains the code and instructions for performing a full mRNA analysis of innate-like T cells in mice treated with LiNKTcellsofαGalCer/CD1d-NP 
#and untreated controls. The analysis includes various visualization techniques such as MA plots, PCA, and volcano plots, among others. The primary focus is on using Subread 
#for alignment and subsequent analysis of the data.

#Build an mouse index
#For mice I downloaded 19:04.2024: Release M34 (GRCm39)
#https://www.gencodegenes.org/
chmod +x GRCm39.primary_assembly.genome.fa.gz
chmod +x gencode.vM34.primary_assembly.annotation.gtf.gz
subread-buildindex -o /archive/pumla/mouse/indexM/GRCm39 /archive/pumla/mouse/indexM/GRCm39.primary_assembly.genome.fa

# Obtain a fastq file
fastq-dump --stdout -X 2 SRR13893223

# To obtain fast.gz file of a paired end read
module avail
module load sratoolkit/3.0.2
fastq-dump  --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR13893223

# To obtain a fastqc report
module load fastqc
fastqc SRR13893223_pass_1.fastq.gz
fastqc SRR13893223_pass_2.fastq.gz

# To obtain alignment
module load subread
chmod +x subread_align.sh
bash subread_align.sh

#I would already have my subread_align.sh file in filezilla written as:
#!/bin/bash
subread-align -i /archive/pumla/mouse/indexM/GRCm39 -r /archive/pumla/mouse/fastqM/SRR13893223_pass_1.fastq.gz -R /archive/pumla/mouse/fastqM/SRR13893223_pass_2.fastq.gz -t 0 --multiMapping -B 1 -T 16 -o /archive/pumla/mouse/outputsM/SRR13893223.bam --sortReadsByCoordinates 2>&1 | tee -a subread_results_23.txt

# To obtain my counts
#First create a file like below to keep in pn filezilla for example a file called `inkt`;
 
 library(Rsubread)
 library(BiocParallel)
 library(limma)
 library(edgeR)
 library(stringr)


#featurecount####
bams <- c("SRR13893218.bam",
          "SRR13893219.bam",
          "SRR13893220.bam",
          "SRR13893221.bam",
          "SRR13893222.bam",
          "SRR13893223.bam")

counts <- featureCounts(files = bams, annot.ext = "/archive/pumla/mouse/annoM/gencode.vM34.primary_assembly.annotation.gtf.gz", 
                        isGTFAnnotationFile= TRUE,
                        GTF.featureType= c("exon"),
                        GTF.attrType="gene_id",
                        useMetaFeatures=T,
                        countMultiMappingReads=T,
                        isPairedEnd=TRUE,
                        nthreads=24
)

counts_df <- as.data.frame(counts$counts)

write.table(counts$counts, "/archive/pumla/mouse/annoM/countsM/counts.tsv", quote = F, col.names = T, sep = "\t")

#calculate fpkm values
z <- DGEList(counts=counts$counts, genes=counts$annotation[,c("GeneID","Length")])
z<- calcNormFactors(z)
RPKM <-rpkm(z)

write.table(RPKM, "/archive/pumla/mouse/annoM/countsM/FPKM.tsv", quote = F, col.names = T, sep = "\t")

done

#To obtain my counts
cat inkt.r
Rscript inkt.r

# PCA
# Setting up environment 
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Set the working directory
setwd("~/inktset1/PCA")

fpkm_data <- read.csv("FPKM.csv", header = TRUE, stringsAsFactors = FALSE)

# Read design.csv containing SampleID and Treatment
design_data <- read.csv("Design.csv", header = TRUE, stringsAsFactors = FALSE)

# View the contents of fpkm_data
View(fpkm_data)

# View the contents of design_data
View(design_data)

# Display the structure of fpkm_data
str(fpkm_data)

# Display the first few rows of fpkm_data
head(fpkm_data)

# Display the structure of design_data
str(design_data)

# Extract numeric data for PCA (assuming columns 2 to end are numeric)
fpkm_numeric <- fpkm_data[, -1]

# Log2 transformation
fpkm_log <- log2(fpkm_numeric + 1)

# Transpose data for PCA (assuming rows are samples and columns are features)
fpkm_log_t <- t(fpkm_log)

# Perform PCA
pca_result <- prcomp(fpkm_log_t, center = TRUE, scale. = FALSE)

# Extract PCA scores
pc_scores <- as.data.frame(pca_result$x)

# Add Treatment information to PCA scores based on SampleID
pc_scores$SampleID <- rownames(pc_scores)
pc_scores_with_treatment <- merge(pc_scores, design_data, by = "SampleID")

# Create PCA plot
p <- ggplot(pc_scores_with_treatment, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(size = 2) +
  geom_text(size = 3, hjust = -0.2, vjust = 0, color = "black") +
  scale_color_manual(values = c("healthy" = "#1f78b4", "tumour" = "#b2df8a")) +
  labs(x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
       color = "Treatment") +
  theme_minimal()

# Display the plot
print(p)

# to remove print 
ggplot(pc_scores_with_treatment, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(size = 3, hjust = -0.2, vjust = 0, color = "black") +
  scale_color_manual(values = c("healthy" = "#1f78b4", "treated" = "#b2df8a")) +
  labs(x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
       color = "Treatment") +
  theme_minimal()


# To try and view some of the data 
# Calculate the range of PC1 and PC2 scores
pc1_range <- range(pc_scores_with_treatment$PC1, na.rm = TRUE)
pc2_range <- range(pc_scores_with_treatment$PC2, na.rm = TRUE)

# Adjust axis limits with some padding
x_limits <- c(pc1_range[1] - 5, pc1_range[2] + 5)  # Adjust the padding (e.g., -5 and +5)
y_limits <- c(pc2_range[1] - 5, pc2_range[2] + 5)  # Adjust the padding (e.g., -5 and +5)

# Create PCA plot with adjusted axis limits and reduced point size
p <- ggplot(pc_scores_with_treatment, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(size = 2) +  # Reduce point size to 2 (adjust as needed)
  geom_text(size = 3, hjust = -0.2, vjust = 0, color = "black") +
  scale_color_manual(values = c("healthy" = "#1f78b4", "treated" = "#b2df8a")) +
  labs(x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
       color = "Treatment") +
  theme_minimal() +
  scale_x_continuous(limits = x_limits) +  # Adjust x-axis limits
  scale_y_continuous(limits = y_limits)    # Adjust y-axis limits

# Display the updated plot
print(p)

#Adjusting the size and saving

# Create PCA plot
p <- ggplot(pc_scores_with_treatment, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(size = 3, hjust = -0.2, vjust = 0, color = "black") +
  scale_color_manual(values = c("healthy" = "#1f78b4", "tumour" = "#b2df8a")) +
  labs(x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
       color = "Treatment") +
  theme_minimal()

# Save plot with adjusted size
ggsave("pca_plot.png", p, width = 10, height = 8)  # Adjust width and height as needed


.........# in a square
# Create PCA plot with adjusted dimensions
p <- ggplot(pc_scores_with_treatment, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(size = 3, hjust = -0.2, vjust = 0, color = "black") +
  scale_color_manual(values = c("healthy" = "#1f78b4", "tumour" = "#b2df8a")) +
  labs(x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)"),
       color = "Treatment") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(0.1, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm")  # Adjust plot margins for more space
  ) +
  coord_fixed(ylim = c(min(pc_scores_with_treatment$PC2) - 20, max(pc_scores_with_treatment$PC2) + 20))  # Adjust ylim to extend the range of y-axis

# Display the plot
print(p)

# Volcano
# Setting up environment 
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load libraries
library(edgeR)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set the working directory
setwd("~/inktset1/Gene names retrived")

# Load count data
count_data <- read.csv('Counts.csv', header = TRUE, row.names = 1)

# Load sample information
sample_info <- read.csv('Design.csv', header = TRUE, row.names = 1)

# Set up the DGEList object
dge <- DGEList(counts = count_data, group = sample_info$Treatment)

# Filter lowly expressed genes
keep <- rowSums(cpm(dge) > 1) >= 1
dge <- dge[keep,]

# Estimate dispersion
dge <- estimateGLMCommonDisp(dge)
dge <- estimateGLMTrendedDisp(dge)
dge <- estimateGLMTagwiseDisp(dge)

# Perform differential expression analysis
fit <- glmQLFit(dge)
qlf <- glmQLFTest(fit, coef = 2)

# Extract results
edger_result <- topTags(qlf, n = nrow(dge))

# Save the edgeR result
write.table(edger_result$table, "edger_results.csv", sep = ",", col.names = NA)

#go and name ensembl_gene_id on the actual edger_results table csv file

#read csv file
edger_results <- read.csv("edger_results.csv")

# Read the CSV file into a data frame
fpkm_data <- read.csv("FPKM.csv")

#sometimes prior to merging you want to remove the versioned gene id and use unversioned IDs

colnames(edger_results)

# Function to remove versioned part from Ensembl gene IDs
remove_version <- function(gene_ids) {
  pattern <- "^(\\w+)\\.\\d+$"  # Regular expression to match the part before the dot
  unversioned_ids <- sub(pattern, "\\1", gene_ids)  # Extract the part before the dot
  return(unversioned_ids)  # Return the unversioned gene IDs
}

# Apply the function to the ensembl_gene_id column in edger_results
edger_results$ensembl_gene_id <- remove_version(edger_results$ensembl_gene_id)


#checking to make sure both file have Gene id as first column
colnames(edger_results)


#now i want to assign gene names and gene type to my data set

# Load biomaRt library
library(biomaRt)

# Use Ensembl biomart for mouse
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")  # for mouse

# Specify attributes to retrieve
attributes <- c("ensembl_gene_id", "external_gene_name", "description")

# Retrieve gene information for mouse
gene_info_mouse <- getBM(attributes = attributes, mart = ensembl)

# Write the gene_mapping_table to a CSV file
write.csv(gene_info_mouse, file = "gene_info_mouse.csv", row.names = FALSE)

colnames(edger_results)

colnames(gene_info_mouse)

# Assuming edger_results and gene_info_mouse dataframes are available

# Merge based on ensembl_gene_id
merged_data <- merge(edger_results, gene_info_mouse, by = "ensembl_gene_id")

# Write the merged data to a CSV file
write.csv(merged_data, file = "merged_data.csv", row.names = FALSE)

colnames(merged_data)

# Sort merged_data by absolute log fold change (descending order)
merged_data <- merged_data %>% 
  arrange(desc(abs(logFC)))

# Select the top 10 most differentially expressed genes
top_10_genes <- head(merged_data, 10)

# Create Volcano Plot
ggplot(merged_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(merged_data$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(x = "Log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot") +
  theme_minimal() +
  geom_text(data = filter(merged_data, FDR < 0.05), aes(label = external_gene_name), vjust = -0.5, color = "red")

# Filter significant genes
significant_genes <- merged_data %>% filter(FDR < 0.05)

# Select top 10 significant genes
top_10_significant_genes <- head(significant_genes, 10)

# Create Volcano Plot
ggplot(merged_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(merged_data$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text(data = top_10_significant_genes, aes(label = external_gene_name), vjust = -0.5, color = "red") +
  labs(x = "Log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot") +
  theme_minimal()

# Select top 5 significant genes
top_5_significant_genes <- head(significant_genes, 5)

# Create Volcano Plot for Top 5 Significant Genes
ggplot(merged_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(merged_data$FDR < 0.05, "red", "black"), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text(data = top_5_significant_genes, aes(label = external_gene_name), vjust = -0.5, color = "red") +
  labs(x = "Log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot for Top 5 Significant Genes") +
  theme_minimal()

# Heatmap
# Setting up environment 
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects

#load libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

#set the working directory
setwd("~/inktset1/Heatmap")

#load the count data
count_data <- read.csv('Counts.csv',header=TRUE,row.names=1)

#load the fpkm data
fpkm_data <- read.csv('FPKM.csv', header=TRUE, row.names=1)


#check column names
colnames(count_data)
head(count_data, 50)

#load the sample information
sample_info <- read.csv('Design.csv',header=TRUE,row.names=1)

#ensure there are no missing values in your count data
sum(is.na(count_data))
sum(is.na(sample_info))

#set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)

#create a deseq object and import the count data and sample information
dds <- DESeqDataSetFromMatrix(countData = count_data,colData = sample_info, design = ~Treatment)

# Set the reference for the Treatment factor and add "peritumour" level
dds$Treatment <- factor(dds$Treatment, levels = c("healthy", "treated"))

#filter the genes based on count per sample
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

#perform the statistical test(s) to identify differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

#change DESeq Object to R object(dataframe)
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)
head(deseq_result, 50)

# Order the results table by increasing p-value (the lowest p-value first, most statistically significant appearing first)
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered, 50)

#make some queriers
#Is blablah gene differentially expressed
deseq_result["23444555",]

#extract the most differentially expressed genes due to the treatment
#select genes with a significant change in gene expression (adjusted p-value below 0.05)
#And log2fold change <1 and >1

#step 1: filter based on p adjusted values
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)

#step 2: filter based on fold changes. here we will use a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

dim(deseq_result)
dim(filtered)   #those that are differentially expressed

#save the deseq result. we will save both the original data(res) and the filtered
write.csv(deseq_result,'de_results.all.csv')
write.csv(deseq_result,'de_results.filtered.csv')


#visualization 
#Heatmap

# Select top 30 genes with the lowest p-values
top_hits <- row.names(deseq_result_ordered)[1:30]

rld <- rlog(dds, blind=FALSE)
pheatmap(assay(rld)[top_hits,], 
         cluster_rows=TRUE, 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sample_info,
         color=colorRampPalette(brewer.pal(9, "Blues"))(100),
         fontsize_row=8, 
         fontsize_col=8)

# Log transform FPKM data
log_fpkm_data <- log2(fpkm_data + 1)

# Calculate Z-scores
cal_z_score <- function(x) {(x - mean(x)) / sd(x)}
zscore_all <- t(apply(log_fpkm_data, 1, cal_z_score))
zscore_subset <- zscore_all[top_hits,]

# Heatmap of log-transformed FPKM values for the top 30 genes
pheatmap(log_fpkm_data[top_hits,], 
         cluster_rows=TRUE, 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sample_info,
         color=colorRampPalette(brewer.pal(9, "Blues"))(100),
         fontsize_row=8, 
         fontsize_col=8)

# Heatmap of Z-scores for the top 30 genes
pheatmap(zscore_subset, 
         cluster_rows=TRUE, 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sample_info,
         color=colorRampPalette(brewer.pal(9, "RdBu"))(100),
         fontsize_row=8, 
         fontsize_col=8)

#Select Top 30 Genes with the Lowest Adjusted p-Values (padj)
# Order results by padj (adjusted p-value)
deseq_result_ordered <- deseq_result[order(deseq_result$padj),]

# Select top 30 genes with the lowest padj values
top_hits <- row.names(deseq_result_ordered)[1:30]

# Perform regularized log transformation
rld <- rlog(dds, blind=FALSE)

# Heatmap of log-transformed counts for the top 30 genes
pheatmap(assay(rld)[top_hits,], 
         cluster_rows=TRUE, 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sample_info,
         color=colorRampPalette(brewer.pal(9, "Blues"))(100),
         fontsize_row=8, 
         fontsize_col=8)
