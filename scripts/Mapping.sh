#!/bin/bash

# This script performs RNA-Seq data processing and analysis using the Alabama Super Computer (ASC)
# It includes loading necessary modules, creating directories, building a reference index for HISAT2,
# mapping reads to the reference genome, converting SAM to BAM, sorting and indexing BAM files,
# assembling transcripts with StringTie, and preparing read count data for downstream analysis.

echo "Loading required modules..."
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread
module load gffcompare

# Set stack size to unlimited
ulimit -s unlimited

# Define variables for directories and reference genome
WORKDIR=/scratch/OmarHas/AT/RAW_Data
CLEAND=/scratch/OmarHas/AT/Cleaned
REFD=/scratch/OmarHas/mapping/ref
MAPD=/scratch/OmarHas/mapping/Map_Hisat
COUNTSD=/scratch/OmarHas/mapping/counts
RESULTSD=/home/aubclsb0314/RNA-seq/csv_counts
REF=AT

echo "Creating directories..."
mkdir -p $REFD $MAPD $COUNTSD $RESULTSD

echo "Creating reference index for HISAT2..."
cd $REFD
gffread $REF.gff3 -T -o $REF.gtf
extract_splice_sites.py $REF.gtf > $REF.ss
extract_exons.py $REF.gtf > $REF.exon
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fna AT_index

echo "Mapping reads to reference genome..."
cd $CLEAND
ls | grep ".fastq" | cut -d "_" -f 1 | sort | uniq > list
cd $MAPD
mv $CLEAND/list .

# Loop through all unique FASTQ files and map them to the reference genome
echo "Processing FASTQ files..."
while read i; do
  echo "Processing sample $i..."
  hisat2 -p 6 --dta --phred33 -x "$REFD"/AT_index -1 "$CLEAND"/"$i"_1_paired.fastq -2 "$CLEAND"/"$i"_2_paired.fastq -S "$i".sam
  samtools view -@ 6 -bS "$i".sam > "$i".bam
  samtools sort -@ 6 "$i".bam "$i"_sorted
  samtools flagstat "$i"_sorted.bam > "$i"_Stats.txt

  mkdir "$COUNTSD"/"$i"
  stringtie -p 6 -e -B -G "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i" "$MAPD"/"$i"_sorted.bam
done < list

echo "Copying files to results directory..."
cp *.txt $RESULTSD

echo "Preparing read count data for downstream analysis..."
python prepDE.py -i $COUNTSD
cp *.csv $RESULTSD

echo "Script completed."
