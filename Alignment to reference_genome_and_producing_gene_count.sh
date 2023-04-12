#!/bin/bash

#Mapping the raw data to a reference genome
#Loading the packages
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/

#Setting the stack size
ulimit -s unlimited

#Define Variables
WD=/scratch/team5/mapping_test
CLEAND=/scratch/team5/CleanData
REFD=/scratch/team5/mapping_test/ref
REF=GCF_004329235.1_PodMur_1.0_genomic
MAPD=/scratch/team5/mapping_test/Map_Hisat
COUNTSD=/scratch/team5/mapping_test/counts
RESULTSD=/home/aubclsb0314/project/counts

#Creating directories and subdirectories
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

#Creating Index file
cd $REFD
gffread $REF.gff3 -T -o $REF.gtf
extract_splice_sites.py $REF.gtf > $REF.ss
extract_exons.py $REF.gtf > $REF.exon
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fna PodMur_index

#Mapping the raw data to the reference genome
cd $CLEAND
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list
cd $MAPD
mv $CLEAND/list .

#Loop through all the unique fasta files and map them to the reference genome
while read i;
do
hisat2 -p 6 --dta --phred33 -x "$REFD"/PodMur_index -1 "$CLEAND"/"$i"_1_paired.fastq -2 "$CLEAND"/"$i"_2_paired.fastq -S "$i".sam
samtools view -@ 6 -bS "$i".sam > "$i".bam
samtools sort -@ 6 "$i".bam "$i"_sorted
samtools flagstat "$i"_sorted.bam > "$i"_Stats.txt

mkdir "$COUNTSD"/"$i"
stringtie -p 6 -e -B -G "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i" "$MAPD"/"$i"_sorted.bam

done < list

#Copying the files to the results directory
cp *.txt $RESULTSD
python prepDE.py -i $COUNTSD
cp *.csv $RESULTSD
