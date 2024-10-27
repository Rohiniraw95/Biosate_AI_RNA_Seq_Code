#!/bin/bash

# Quality Control with FastQC and MultiQC
echo "Running FastQC..."
fastqc *.fastq
pip install multiqc
multiqc .

# Trimming with Cutadapt
echo "Running Cutadapt..."
for file in *.fastq.gz; do
    cutadapt -m 20 -q 20 -o "trimmed_${file}" "$file"
done
fastqc trimmed_*.fastq.gz
multiqc .

# Additional Trimming with Fastp
echo "Running Fastp on sample trimmed_Heart_ZT0_1..."
fastp -i trimmed_Heart_ZT0_1_R1.fastq.gz -I trimmed_Heart_ZT0_1_R2.fastq.gz \
-o cleaned_R1.fastq.gz -O cleaned_R2.fastq.gz \
--trim_front1 5 --trim_tail1 5 --cut_mean_quality 20 \
--detect_adapter_for_pe --length_required 20 --html fastp_report.html
fastqc cleaned_R*.fastq.gz
multiqc .

# HISAT2 Indexing
echo "Building HISAT2 index..."
conda install -c bioconda hisat2 -y
hisat2 --version
hisat2-build GRCm39.primary_assembly.genome.fa genome_index

# Remove Gaps in Genome Fasta
echo "Removing gaps in genome fasta with SeqKit..."
conda install -c bioconda seqkit -y
seqkit seq --upper-case --remove-gaps GRCm39.primary_assembly.genome.fa.gz -o cleaned_genome.fa
hisat2-build cleaned_genome.fa genome_index

# Alignment of Reads with HISAT2
input_dir="/media/rohini/a364e659-dcb0-418b-8128-88023676592e/biostate"
output_dir="/media/rohini/san"  # Update this to your desired output directory

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# List of Sample Pairs
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
for sample in "${samples[@]}"
do
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

    # Run HISAT2 and SAMtools
    echo "Aligning $r1 and $r2..."
    hisat2 -x "${input_dir}/genome_index" \
        -1 "$r1_path" \
        -2 "$r2_path" \
        --very-sensitive --threads 8 | \
        samtools sort -@ 8 -o "$output_file" --output-fmt bam

    echo "Processed: $r1 and $r2 -> $output_file"
done

# Convert BAM Files to CSV with FeatureCounts
echo "Generating gene count CSV with FeatureCounts..."
featureCounts -T 8 -p -a /media/rohini/san/gencode.vM35.basic.annotation.gtf \
  -o /media/rohini/san/gene_counts_paired.csv /media/rohini/san/*.bam

echo "All tasks completed successfully."

