#!/bin/bash

# This script runs Cutadapt and Trimmomatic on paired-end sequencing data.
# Ensure that you have the correct adapter sequences for your specific library preparation kit and sequencing platform.

# Set your input files, adapters file, and output directories
READ1=input_read1.fastq.gz
READ2=input_read2.fastq.gz
ADAPTERS=adapters.fasta
CUTADAPT_OUTPUT="cutadapt_output" # include path of the output 
TRIMMOMATIC_OUTPUT="trimmomatic_output" # include path of the output

# Trimmomatic settings - adjust as needed
TRIMMOMATIC_SETTINGS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Load the necessary modules in the ASC environment
module load cutadapt trimmomatic

# Create output directories if they do not exist
mkdir -p $CUTADAPT_OUTPUT $TRIMMOMATIC_OUTPUT

# Run Cutadapt to trim adapter sequences from the reads
echo "Running Cutadapt..."
cutadapt -a file:$ADAPTERS -A file:$ADAPTERS -o ${CUTADAPT_OUTPUT}/trimmed_read1.fastq.gz -p ${CUTADAPT_OUTPUT}/trimmed_read2.fastq.gz $READ1 $READ2
echo "Cutadapt completed."

# Run Trimmomatic to further trim and filter reads based on quality and length, using the output from Cutadapt
echo "Running Trimmomatic..."
trimmomatic PE -threads 4 -phred33 \
  ${CUTADAPT_OUTPUT}/trimmed_read1.fastq.gz ${CUTADAPT_OUTPUT}/trimmed_read2.fastq.gz \
  ${TRIMMOMATIC_OUTPUT}/trimmed_paired_read1.fastq.gz ${TRIMMOMATIC_OUTPUT}/trimmed_unpaired_read1.fastq.gz \
  ${TRIMMOMATIC_OUTPUT}/trimmed_paired_read2.fastq.gz ${TRIMMOMATIC_OUTPUT}/trimmed_unpaired_read2.fastq.gz \
  $TRIMMOMATIC_SETTINGS
echo "Trimmomatic completed."

echo "Pipeline finished."
