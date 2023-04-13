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


# I am using Alabama Super Computer (ASC) to run this script so I am just loading the packages;
#if you are using other SC, probably you'll be using the same lines of code, you'll need to check specific insturctions for the SC you're using.
#If you're working on your local computer however, you can use anaconda or pip to install these packages.



#Setting the stack size
ulimit -s unlimited

#Because I am using the ASC I will set the stack size to unlimited which is a indicates the memories we are using that stores the temporary information 
#like variables and function calling.
#However, I don't recommend setting the stack size to unlimited if you're using a local computer or a VM,
#as most operating systems might face stack overflow and shutdown the whole program.




#Define Variables

#Just like what we've done in the cleaning step.For sack of reproductivity, and to make this script easier to use for other projects, we will set variables to everything,
#the reference genome, raw data, directories and subdirectories.

First, Making a working directory
WORKDIR=/scratch/OmarHas/AT/RAW_Data
CLEAND=/scratch/OmarHas/AT/Cleaned
REFD=/scratch/OmarHas/mapping/ref
MAPD=/scratch/OmarHas/mapping/Map_Hisat
COUNTSD=/scratch/OmarHas/mapping/counts
RESULTSD=/home/aubclsb0314/RNA-seq/csv_counts

#making a variable for the name of the reference genome. I downloaded the annotation and fasta files from TAIR, they came with a different name but I changed 
#it to something easier, but it won't matter as long as you're going to use a variable. 

REF=AT

#Creating directories and subdirectories 
# -p option just makes any directory or subdirectory that's missing from the path. It won't hurt even if that directory exists. 

mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

#Creating Index file
#Next big step is creating reference index our HiSat2 mapper
#First, go to the directory where you downloaded the reference genome
#Second, we will need to identify the exons and splicing sites; 
#we will be using two Python scripts that are commonly used in RNA-Seq data analysis pipelines "extract_splice_sites.py" and "extract_exons.py.
#These two scripts are included with thte HiSat2 package.
#In addtion to using gffread which converts the annotation file from gff3 to gft format for HiSat to use

cd $REFD
gffread $REF.gff3 -T -o $REF.gtf
extract_splice_sites.py $REF.gtf > $REF.ss
extract_exons.py $REF.gtf > $REF.exon
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fna AT_index

#Mapping the raw data to the reference genome
#All of the steps above were just preparing the files to be used by HiSat2. Now that we created the index, gft files,
#we will move on to map our data to the reference genome.
#First, move to the directory where our raw data are present

cd $CLEAND

#Second step is creating a list of all of the fastq files that we need to map.
#This syntax catches all of files that has .fastq using grep, and cut the underscore from the name,
#then it only includes one version of the file using uniq and outputs the content into a list.
#Thrid, we will need to move the mapping directory, and move the unique ids from the orignal files to map

ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list
cd $MAPD
mv $CLEAND/list .

#Loop through all the unique fasta files and map them to the reference genome


while read i;
do
hisat2 -p 6 --dta --phred33 -x "$REFD"/AT_index -1 "$CLEAND"/"$i"_1_paired.fastq -2 "$CLEAND"/"$i"_2_paired.fastq -S "$i".sam
samtools view -@ 6 -bS "$i".sam > "$i".bam
samtools sort -@ 6 "$i".bam "$i"_sorted
samtools flagstat "$i"_sorted.bam > "$i"_Stats.txt

mkdir "$COUNTSD"/"$i"
stringtie -p 6 -e -B -G "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i" "$MAPD"/"$i"_sorted.bam

done < list

#Copying the files to the results directory
cp *.txt $RESULTSD

# the PreDE.py is a python script that gather of the read counts from Stringtie into one file which is the csv file that
#we're going to import into R for differential genes expression. 

python prepDE.py -i $COUNTSD
cp *.csv $RESULTSD
