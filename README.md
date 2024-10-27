# RNA-Seq Analysis Pipeline

## Overview
This pipeline provides an RNA-Seq analysis workflow, including quality control, trimming, alignment, feature counting, and differential expression analysis.

## Steps
1. Quality Control with FastQC and MultiQC
2. Trimming with Cutadapt
3. Cleaning with Fastp
4. Indexing and Alignment with HISAT2
5. Counting Features with featureCounts
6. Differential Expression Analysis in R

## Requirements
- FastQC
- MultiQC
- Cutadapt
- Fastp
- HISAT2
- SeqKit
- featureCounts
- R (with DESeq2, ggplot2, and pheatmap)

## Usage
Run each step in sequence. Update the input/output directories as needed in each script.
