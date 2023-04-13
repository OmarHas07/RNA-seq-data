#!/bin/bash

### Load the modules you need to use
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load trimmomatic/0.38
module load cutadabt


### Define variable for directories and define array for samples from csv file.
CLEAND=/scratch/OmarHas/AT/Cleaned
WORKDIR=/scratch/OmarHas/AT/RAW_Data
OUTDIR=/scratch/OmarHas/AT/Cleaned
FASTQCRESULTS=/home/aubclsb0314/RNA-seq/fastqc


samples=(`awk -F, '{print $1}' ~/Transcriptomes_RunInfo.csv`)

### Go to the working directoru, copy an adapter file here and loop through the samples

cd $WORKDIR

cp /home/aubclsb0314/RNA-Seq/AdaptersToTrim_All.fa .

for x in ${samples[@]}; do

# Trimming with Trimmomatic
java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar PE -threads 6 -phred33 \
    ${x}_1.fastq ${x}_2.fastq \
    ${x}_1_paired.fastq ${x}_1_unpaired.fastq \
    ${x}_2_paired.fastq ${x}_2_unpaired.fastq \
    ILLUMINACLIP:AdaptersToTrim_All.fa:1:35:20 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:30 MINLEN:25

# Trimming with Cutadapt
cutadapt -a adapter_1 -A adapter_2 -o ${x}_1_paired_trimmed.fastq -p ${x}_2_paired_trimmed.fastq ${x}_1_paired.fastq ${x}_2_paired.fastq -q 20,20 --minimum-length=36

done
# copy the cleand data to another directory. 
cp *_paired_trimmed.fastq $CLEAND

### Run fastqc on the cleaned paired files

fastqc -t 6 *paired_trimmed.fastq

#Copy fastqc output to another directory. 

# Copy FastQC output files to the new directory
cp *_fastqc.zip $FASTQCRESULTS
cp *_fastqc.html $FASTQCRESULTS