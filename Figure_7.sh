#!/bin/bash

# Performs quality control checks on fastq files using FASTQC and aggregates the results using MultiQC.

# Directory containing fastq files
fastqDir="/path/to/fastq/files"

# Output directory for FASTQC results
outputDir="/path/to/fastqc/output"

# Run FASTQC on all fastq files in the specified directory
fastqc -o ${outputDir} ${fastqDir}/*.fastq.gz

# Run MultiQC to aggregate the FASTQC reports
multiqc ${outputDir} -o ${outputDir}



#!/bin/bash

# Trimming adapters and filtering for quality in RNA-seq data using Cutadapt.

fastqFiles='list_fastq.txt'  # List of fastq file names

# Loop through each fastq file listed and perform trimming and quality filtering
for fastq in $(cat $fastqFiles)
do
    echo "Trimming and filtering sample: $fastq"
    cutadapt -j 10 -q 30,30 -m 20 -a ADAPTER_FWD -A ADAPTER_REV -o "${fastq}_trimmed.1.fastq.gz" -p "${fastq}_trimmed.2.fastq.gz" "${fastq}_R1.fastq.gz" "${fastq}_R2.fastq.gz"
    echo "Completed sample: $fastq"
done



#!/bin/bash

# Performs read mapping using HISAT2 and directly pipes the output to Samtools for conversion and sorting.


# List of trimmed fastq files
fastqFiles='list_fastq.txt'

# Loop through each file in the list and perform mapping and post-processing
for fastq in $(cat $fastqFiles)
do
    echo "Mapping and processing sample: $fastq"

    # Mapping reads to the reference genome with HISAT2 --> Use samtools for conversion and sorting
    hisat2 -p 12 -x path_to_hisat2_index -1 "${fastq}_trimmed.1.fastq.gz" -2 "${fastq}_trimmed.2.fastq.gz" | \
    samtools sort -@ 12 -o "${fastq}.sorted.bam"

    echo "DONE Processing sample: $fastq"
done




#!/bin/bash

# Runs FeatureCounts to count features from BAM files.


# Running FeatureCounts to count features from BAM files
featureCounts -p -t exon -g gene_name --extraAttributes gene_id,gene_type -T 10 -a path_to_annotation_file -o FeatureCounts_Output *.sorted.bam > featurecounts.screen-output.log





# Load necessary libraries
library(DESeq2)
library(writexl)
library(ggplot2)
library(org.Hs.eg.db)
library(edgeR)
library(sva)

# Read count matrix and metadata
counts <- as.matrix(read.delim("path/to/count_data.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1))
metadata <- as.data.frame(readxl::read_excel("path/to/metadata.xlsx"))
metadata$condition <- factor(metadata$condition)

# Ensure the samples in metadata align with the columns in the count data
stopifnot(identical(metadata$sample, colnames(counts)))

# Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)

# Filter lowly expressed genes
dds <- dds[rowSums(counts(dds)) >= 10,]

# Run the differential expression analysis
dds <- DESeq(dds)

# Transform count data using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
vsd_transformed_data <- assay(vsd)

# Optional: Save the VST-transformed data
write_xlsx(as.data.frame(vsd_transformed_data), "vsd_transformed_data.xlsx")

# Enhanced Heatmap of the Top 2000 Most Variable Genes
library(pheatmap)

# Selecting the top 2000 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 2000)
heatmap_data <- assay(vsd)[topVarGenes, ]

# Define a color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# Generating the heatmap
pheatmap(heatmap_data,
         color = my_palette,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = metadata)


# PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
plotPCA(vsd, intgroup = "condition")

# Saving PCA results
write_xlsx(as.data.frame(pca_data), "PCA_results.xlsx")




# UMAP Analysis
# It's assumed that you have already installed and loaded the 'umap' package
library(umap)
library(ggplot2)

# Preparing data for UMAP
# Ensure that gene names are not included as part of the data
vsd_data_for_umap <- t(assay(vsd))

# Running UMAP
set.seed(42)  # for reproducibility
umap_results <- umap(vsd_data_for_umap)

# Creating a data frame for plotting
umap_data <- data.frame(UMAP1 = umap_results$layout[,1],
                        UMAP2 = umap_results$layout[,2],
                        condition = metadata$condition)

# Plotting UMAP results
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point() +
  labs(title = "UMAP plot of RNA-seq data", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# Differential Expression Results Visualization
# Assuming 'res' holds the results from DESeq2

# Volcano Plot
# Adjusting the threshold for significance
alpha <- 0.05

# Creating a volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(col = padj < alpha & abs(log2FoldChange) > log2(1.5)), 
             size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Volcano plot of differential expression",
       x = "Log2 fold change",
       y = "-Log10 p-value") +
  theme_minimal() +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed")

# Saving the DESeq2 results
# Convert the DESeq2 results to a dataframe for easier handling
res_df <- as.data.frame(res)

# Optionally, filter and annotate the results for significantly differentially expressed genes
significant_res <- subset(res_df, padj < alpha)

# Save the full results and the filtered significant results to an Excel file
write_xlsx(list(Full_DESeq2_Results = res_df,
                Significant_DESeq2_Results = significant_res), 
           "DESeq2_results.xlsx")
