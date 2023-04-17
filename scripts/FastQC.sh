#!/bin/bash

# This script runs FastQC on raw sequencing data and compresses the results into a tar.gz file.
# It is designed for use on the Alabama Supercomputer Center (ASC) or other high-performance computing environments.

echo "Loading required modules..."
module load fastqc/0.10.1

# Define variables for directories and output
QUALITY_DIR=/scratch/OmarHas/Quality
RAWDATA_DIR=/scratch/OmarHas/RawData
FASTQC_RESULTS_DIR=/home/aubclsb0314/RNA-seq/fastqc
PRE_CLEAN_QUALITY=PreCleanQuality

echo "Creating directories..."
mkdir -p $RAWDATA_DIR
mkdir -p $QUALITY_DIR/$PRE_CLEAN_QUALITY

echo "Running FastQC on raw data files..."
cd $RAWDATA_DIR
fastqc *.fastq --outdir=$QUALITY_DIR/$PRE_CLEAN_QUALITY

echo "Compressing FastQC results into a tar.gz file..."
cd $QUALITY_DIR/$PRE_CLEAN_QUALITY
OUTPUT_NAME=fastqc_results
tar cvzf $OUTPUT_NAME.tar.gz $QUALITY_DIR/$PRE_CLEAN_QUALITY/*

echo "Script completed."
