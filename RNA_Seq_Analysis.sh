#!/bin/bash
# Author: P D Rohini
# Date: 27/10/2024

# Dataset: RNA-seq data (https://drive.google.com/drive/u/1/folders/1wIr_bwhWk3gCb4QjnWZrhxNbbh_GZCPz)

# Check if necessary programs are installed
command -v fastqc >/dev/null 2>&1 || { echo >&2 "FastQC is required but it's not installed. Aborting."; exit 1; }
command -v cutadapt >/dev/null 2>&1 || { echo >&2 "Cutadapt is required but it's not installed. Aborting."; exit 1; }
command -v fastp >/dev/null 2>&1 || { echo >&2 "Fastp is required but it's not installed. Aborting."; exit 1; }
command -v hisat2 >/dev/null 2>&1 || { echo >&2 "HISAT2 is required but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "SAMtools is required but it's not installed. Aborting."; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo >&2 "featureCounts is required but it's not installed. Aborting."; exit 1; }

# FastQC for quality control of raw FASTQ files
echo "Running FastQC on raw FASTQ files..."
fastqc *.fastq

# Install MultiQC if not already installed
pip install multiqc  # Ensure pip is installed and configured
multiqc .  # Generate a summary report

# Cutadapt to trim low-quality reads
echo "Trimming low-quality reads using Cutadapt..."
for file in *.fastq.gz; do
    cutadapt -m 20 -q 20 -o trimmed_"$file" "$file"
done
fastqc trimmed_*.fastq.gz  # Check trimmed files
multiqc .  # Generate report for trimmed files

# Fastp for further processing
echo "Running Fastp for further read processing..."
fastp -i trimmed_Heart_ZT0_1_R1.fastq.gz -I trimmed_Heart_ZT0_1_R2.fastq.gz \
-o cleaned_R1.fastq.gz -O cleaned_R2.fastq.gz \
--trim_front1 5 --trim_tail1 5 --cut_mean_quality 20 \
--detect_adapter_for_pe --length_required 20 --html fastp_report.html
fastqc cleaned_*.fastq.gz  # Check cleaned files
multiqc .  # Report on cleaned files

# HISAT2 index creation
echo "Building HISAT2 index..."
conda install -c bioconda hisat2
hisat2 --version  # Check HISAT2 version
hisat2-build GRCm39.primary_assembly.genome.fa genome_index

# Alignment of reads with genome_index
# Define the input and output directories
input_dir="/media/rohini/a364e659-dcb0-418b-8128-88023676592e/biostate"
output_dir="/media/rohini/san"  # Update this to your desired output directory

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# List of sample pairs for alignment
samples=(
    "cleaned_Heart_ZT0_1_R1.fastq.gz cleaned_Heart_ZT0_1_R2.fastq.gz"
    "cleaned_Heart_ZT0_2_R1.fastq.gz cleaned_Heart_ZT0_2_R2.fastq.gz"
    "cleaned_Heart_ZT12_1_R1.fastq.gz cleaned_Heart_ZT12_1_R2.fastq.gz"
    "cleaned_Heart_ZT12_2_R1.fastq.gz cleaned_Heart_ZT12_2_R2.fastq.gz"
    "cleaned_Liver_ZT0_1_R1.fastq.gz cleaned_Liver_ZT0_1_R2.fastq.gz"
    "cleaned_Liver_ZT0_2_R1.fastq.gz cleaned_Liver_ZT0_2_R2.fastq.gz"
    "cleaned_Liver_ZT12_1_R1.fastq.gz cleaned_Liver_ZT12_1_R2.fastq.gz"
    "cleaned_Liver_ZT12_2_R1.fastq.gz cleaned_Liver_ZT12_2_R2.fastq.gz"
)

# Run HISAT2 and SAMtools for each sample pair
for sample in "${samples[@]}"; do
    # Split the sample pair into R1 and R2
    read -r r1 r2 <<< "$sample"

    # Full paths for R1 and R2
    r1_path="${input_dir}/${r1}"
    r2_path="${input_dir}/${r2}"

    # Extract the base name without the leading "cleaned_" and the file extension
    base_name="${r1##cleaned_}"  # Remove "cleaned_" prefix
    base_name="${base_name%%_R*}"  # Remove "_R*" part

    # Define output file name
    output_file="${output_dir}/${base_name}_sorted.bam"  # Name based on input files

    # Run HISAT2 and SAMtools to align reads and sort
    hisat2 -x "${input_dir}/genome_index" \
        -1 "$r1_path" \
        -2 "$r2_path" \
        --very-sensitive --threads 8 | \
        samtools sort -@ 8 -o "$output_file" --output-fmt bam

    echo "Processed: $r1 and $r2 -> $output_file"
done

# Count features using featureCounts
echo "Counting features using featureCounts..."
featureCounts -T 8 -p -a /media/rohini/san/gencode.vM35.basic.annotation.gtf \
  -o /media/rohini/san/gene_counts_paired.csv /media/rohini/san/*.bam

# Execute R script for Differential Gene Expression Analysis & Functional Enrichment Analysis
Rscript -e "

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)  
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggtree)

#Alignment Statistics
# Sample data
samples <- c("Heart_ZT0_1", "Heart_ZT0_2", "Heart_ZT12_1", "Heart_ZT12_2", 
            "Liver_ZT0_1", "Liver_ZT0_2", "Liver_ZT12_1", "Liver_ZT12_2")

# Assigned reads
assigned_reads <- c(58938673, 48502407, 60904988, 23387384, 37458596, 56691703, 26701998, 19879985)

# Unassigned multi-mapping reads
unassigned_multimapping <- c(39446669, 31901118, 43474042, 14419741, 
                     80114626, 105466664, 59774825, 37543445)

# Calculate total reads
total_reads <- assigned_reads + unassigned_multimapping

# Calculate multimapping percentage
multimapping_percentage <- (unassigned_multimapping / total_reads) * 100

# Create a data frame for the results
results <- data.frame(Sample = samples, Assigned_Reads = assigned_reads, Unassigned_Multimapping = unassigned_multimapping, Total_Reads = total_reads, Multimapping_Percentage = multimapping_percentage)

# Print results
print(results)

# Optional: Plotting the results for better visualization

ggplot(results, aes(x = Sample, y = Multimapping_Percentage)) +
 geom_bar(stat = "identity", fill = "blue") + theme_minimal() + labs(title = "Multimapping Percentage per Sample", x = "Sample", y = "Multimapping Percentage (%)") +
             theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Data import and preprocessing
gene_counts <- read.csv("gene_counts_paired.csv", row.names = 1)
gene_counts <- select(gene_counts, 6:13)

# Get the number of samples
num_samples <- ncol(gene_counts) 

# Create sample metadata based on the number of columns
if (num_samples == 8) {
    sample_metadata <- data.frame(
        row.names = colnames(gene_counts),  
        condition = c("Heart_ZT0_1", "Heart_ZT0_2", "Heart_ZT12_1", "Heart_ZT12_2", 
                      "Liver_ZT0_1", "Liver_ZT0_2", "Liver_ZT12_1", "Liver_ZT12_2")
    )
} else {
    stop("Number of samples does not match the expected count. Please check your gene_counts file.")
}

# Print sample metadata to confirm
print(sample_metadata)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = gene_counts,
    colData = sample_metadata,
    design = ~ tissue + time + tissue:time  
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract normalized counts for reproducibility check
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
normalized_counts <- assay(vsd)

# Scatter plot: Sample correlation
plot(
    normalized_counts[, "Heart_ZT0_1"], normalized_counts[, "Heart_ZT0_2"],
    xlab = "Heart ZT0 Replicate 1", ylab = "Heart ZT0 Replicate 2",
    main = "Scatter Plot: Reproducibility Check"
)

# Compute sample-to-sample distances
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)

# Create heatmap of sample distances
pheatmap(sample_dist_matrix, 
         clustering_distance_rows = sample_dist, 
         clustering_distance_cols = sample_dist, 
         main = "Sample-to-Sample Distance Heatmap")

# Perform PCA using transformed counts
pca_data <- prcomp(t(assay(vsd)))

# Create a data frame for plotting PCA
pca_df <- as.data.frame(pca_data$x)
pca_df$condition <- colData(dds)$condition

# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4) +
    labs(title = "PCA: Overall Pattern of Variation", x = "PC1", y = "PC2") +
    theme_minimal()

# Tissue-specific DEGs
res_tissue <- results(dds, contrast = c("tissue", "Heart", "Liver"))
# Time-specific DEGs
res_time <- results(dds, contrast = c("time", "ZT0", "ZT12"))
# Interaction effect
res_interaction <- results(dds, name = "tissueHeart.timeZT12")

# Create a metadata frame matching sample names
sample_metadata <- data.frame(
    sample = colnames(gene_counts),
    tissue = c("Heart", "Heart", "Heart", "Heart", 
               "Liver", "Liver", "Liver", "Liver"),
    time = c("ZT0", "ZT0", "ZT12", "ZT12", 
             "ZT0", "ZT0", "ZT12", "ZT12")
)

# Ensure tissue and time are treated as factors
sample_metadata$tissue <- factor(sample_metadata$tissue)
sample_metadata$time <- factor(sample_metadata$time)

# Check metadata to ensure it aligns
print(sample_metadata)

# Create DESeq2 dataset again
dds <- DESeqDataSetFromMatrix(
    countData = gene_counts,
    colData = sample_metadata,
    design = ~ tissue + time + tissue:time
)

# Run DESeq2 to fit the model
dds <- DESeq(dds)

# List the names available for contrasts
resultsNames(dds)

# Tissue-specific DEGs
res_tissue <- results(dds, name = "tissue_Liver_vs_Heart")
# Time-specific DEGs
res_time <- results(dds, name = "time_ZT12_vs_ZT0")
# Interaction effect
res_interaction <- results(dds, name = "tissueLiver.timeZT12")

# Check summary of results
summary(res_tissue)
summary(res_time)
summary(res_interaction)

# Paired contrast: Heart vs Liver at ZT0
res_heart_liver_ZT0 <- results(dds, contrast = list(c("tissue_Liver_vs_Heart")))
# Paired contrast: Heart vs Liver at ZT12
res_heart_liver_ZT12 <- results(dds, contrast = list(c("tissueLiver.timeZT12")))

# Check summary of the paired contrasts
summary(res_heart_liver_ZT0)
summary(res_heart_liver_ZT12)

# Volcano plot for tissue-specific DEGs
EnhancedVolcano(res_tissue,
                lab = rownames(res_tissue),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Tissue-Specific DEGs: Heart vs Liver',
                pCutoff = 0.05,
                FCcutoff = 1)

# Volcano plot for time-specific DEGs
EnhancedVolcano(res_time,
                lab = rownames(res_time),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Time-Specific DEGs: ZT0 vs ZT12',
                pCutoff = 0.05,
                FCcutoff = 1)

# Extract significant DEGs from tissue-specific analysis
sig_tissue <- res_tissue[which(res_tissue$padj < 0.05), ]
# Get the gene names and subset the original counts for these genes
sig_genes <- rownames(sig_tissue)
heatmap_counts <- gene_counts[sig_genes, ]

# Normalize counts using variance stabilizing transformation (VST)
vst_counts <- vst(dds, blind = TRUE)

# Heatmap of significant DEGs
pheatmap(assay(vst_counts)[sig_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = sample_metadata,
         main = "Clustering of Significant DEGs")                                                                                                         

# Convert top_counts to a data frame
top_counts_df <- as.data.frame(top_counts)

# Check the structure
str(top_counts_df)

# Move row names to a column named "GeneID"
top_counts_df$GeneID <- rownames(top_counts)
rownames(top_counts_df) <- NULL  # Remove row names

melted_data <- pivot_longer(
    data = top_counts_df, 
    cols = -GeneID,  # All columns except "GeneID"
    names_to = "Sample", 
    values_to = "Expression"
)

head(melted_data)

# Merge melted_data with sample_metadata
melted_data <- merge(melted_data, sample_metadata, by.x = "Sample", by.y = "row.names")

# Plot expression patterns of top DEGs
ggplot(melted_data, aes(x = time, y = Expression, color = tissue, group = tissue)) +
    geom_line() +
    facet_wrap(~ GeneID, scales = "free_y") +
    theme_bw() +
    ggtitle("Expression Patterns of Top 10 DEGs")

# Install necessary libraries
BiocManager::install(c("clusterProfiler", "enrichplot", "ggtree", "org.Mm.eg.db"),  
                     force = TRUE)

# Load additional necessary libraries
library(clusterProfiler)
library(enrichplot)
library(ggtree)
library(org.Mm.eg.db)

# Run DESeq to get results
dds <- DESeq(dds) 
res <- results(dds)

# Step 1: Extract significant genes based on adjusted p-value and log2 fold change
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Step 2: Get the ENSEMBL IDs of significant genes
sig_gene_symbols <- rownames(sig_genes)  
print(head(sig_gene_symbols))  # Verify ENSEMBL IDs

# Step 3: Clean ENSEMBL IDs to remove version numbers
sig_gene_symbols_cleaned <- gsub("\..*$", "", sig_gene_symbols)  # Remove everything after the dot
print(head(sig_gene_symbols_cleaned))  # Check cleaned ENSEMBL IDs

# Step 4: Convert cleaned ENSEMBL IDs to SYMBOLS
converted_genes <- bitr(
    sig_gene_symbols_cleaned, 
    fromType = "ENSEMBL", 
    toType = "SYMBOL", 
    OrgDb = org.Mm.eg.db
)

# Check if conversion was successful
if (nrow(converted_genes) == 0) stop("No valid ENSEMBL IDs were found for conversion")

# Step 5: Generate a named vector for gene symbols with their log2FoldChange values
sig_gene_symbols_vector <- setNames(
    sig_genes[match(converted_genes$ENSEMBL, rownames(sig_genes)), "log2FoldChange"], 
    converted_genes$SYMBOL
)

# Step 6: Perform Gene Ontology (GO) enrichment analysis using clusterProfiler
go_enrichment <- enrichGO(
    gene = names(sig_gene_symbols_vector), 
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL", 
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
)

# Step 7: Check GO enrichment results
print(head(go_enrichment))

# Step 8: Visualize GO enrichment results
dotplot(go_enrichment, showCategory = 20) + ggtitle("GO Enrichment Analysis - Biological Processes")

# Step 1: Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(
    gene = names(sig_gene_symbols_vector),
    organism = 'mmu', 
    pvalueCutoff = 0.05
)

# Step 2: Visualize KEGG enrichment results
barplot(kegg_enrichment, showCategory = 10) + ggtitle("KEGG Pathway Enrichment Analysis")

# Step 3: Plotting significant pathways
pathways <- kegg_enrichment$ID[1:5]  # Get top 5 pathways
for (pathway in pathways) {
    pathview(
        gene.data = sig_gene_symbols_vector, 
        pathway.id = pathway, 
        species = "mmu", 
        out.suffix = pathway
    )
}

# End of script

