# Purpose of this script 
To analyze RNA-Seq data, the following steps are typically taken:

- Ensuring the **quality** of reads with `Fastqc`
- **Trimming** (post-cleaning) with `trimmomatic`
- **Identify exons and introns
- ** Mapping to a reference genome using `Hisat`
- **aligning the annotation to the RNA transcript 
- Producing **gene counts** with `stringtie `
- Analysis of **differential gene expression** between treatments using `DESeq2` or `edgeR`

# Languages used in this script: 
- Bash 
- Python
- R 

# Number for scripts 

Typicall I will try to produce one script to the cleaning step and another script from identifying exons and introns to a gene count file that can be used in R. 

# Mapping the raw data to a reference genome

For this particular script we will be using Hisat2 to map our raw data to reference genome. My RNA-seq data comes from  on the asexual whiptail lizard (Aspidoselis tesselata) samples; but there is no annotated reference genome for this species. Hisat2 is mapper that requires reference genome; so we are going to use a well-annotated reference genome from a close relative of the whiptail lizard called Podarcis muralis. 

I downloaded the fna (Fasta file), gff3, and gtf files from NCBI; these files are attached in this respiratory, and you can also find them in the link below
https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_004329235.1/

### Loading the packages 

I am using Alabama Super Computer (ASC) to run this script so I am just loading the packages; if you are using other SC, probably you'll be using the same lines of code, you'll need to check specific insturctions for the SC you're using. If you're working on your local computer however, you can use anaconda or pip to install these packages. 

``` 

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

```

### Setting the stack size

Because I am using the ASC I will set the stack size to unlimited which is a indicates the memories we are using that stores the temporary information like variables and function calling. However, I don't recommend setting the stack size to  unlimited if you're using a local computer or a VM, as most operating systems might face stack overflow and shutdown the whole program. 

```
ulimit -s unlimited
``` 

### Define Variables 

For sack of reproductivity, and to make this script easier to use for other projects, we will set variables to everything, the reference genome, raw data, directories and subdirectories. 

First, Making a working directory

```
WD=/scratch/team5/mapping_test  

```

Then, making another variable to the directory where we have our trimmed and cleaned data that was produced from trimmoatic

```
CLEAND=/scratch/team5/CleanData 

```

Path to our reference genome

```
REFD=/scratch/team5/mapping_test/ref 
```

We als going to make a variable for the name of the reference genome files that comes before the .fna, .gtf, and the gff3 extentions, when I downloaded the the data from NCBI it has that name before the extention: GCF_004329235.1_PodMur_1.0_genomic, you can change it to a simpler one; but I will just leave it anyone because we will use the variable instead 

```
REF=GCF_004329235.1_PodMur_1.0_genomic    

```
Then, a directory for the Hisat output 

```
MAPD=/scratch/team5/mapping_test/Map_Hisat  

```

Directory for Stringtie output 

```
COUNTSD=/scratch/team5/mapping_test/counts 

```

And, directory where we will have the data that we are going to look at such as the statisitics of Hisat2 and most importantly csv file that has the read counts of our data

```
RESULTSD=/home/aubclsb0314/project/counts 

```

Next, we will make these paths by creating directories and subdirectories of the variables above, I am going to use mkdir -p, the -p option just create any directory or subdirectory missing from the path. 

```
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD


```


### Creating Index file 


Next big step is creating reference index our HiSat2 mapper

First, go to the directory where you downloaded the reference genome 

``` 
cd $REFD

```

Second, we will need to identify the exons and splicing sites; we will be using two Python scripts that are commonly used in RNA-Seq data analysis pipelines "extract_splice_sites.py" and "extract_exons.py. These two scripts are included with thte HiSat2 package. 

In addtion to using gffread which converts the annotation file from gff3 to gft format for HiSat to use 

``` 
gffread $REF.gff3 -T -o $REF.gtf               
extract_splice_sites.py $REF.gtf > $REF.ss
extract_exons.py $REF.gtf > $REF.exon
```

Finally, creating a HiSat2 index file for the reference genome. 

```
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fna PodMur_index
``` 

PodMur is going to be the name of my index files; you can change this name to the whatever name you want. After these steps you should see some new files in your reference directory: .exon, .ss and the index files. 

### Mapping the raw data to the reference genome

All of the steps above were just preparing the files to be used by HiSat2. Now that we created the index, gft files, we will move on to map our data to the reference genome. 

First, move to the directory where our raw data are present

```
cd $CLEAND
```
**Notice: in the cleaned data directory you should have two files for each sample, 1 for the forward, and 2 for the reverse reads. 
**

Second step is creating a list of all of the fastq files that we need to map. 

``` 
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list  
``` 
This syntax catches all of files that has .fastq using **grep**, and **cut** the underscore from the name, then it only includes one version of the file using **uniq** and outputs the content into a list. 

Thrid, we will need to move the mapping directory, and move the unique ids from the orignal files to map 

```
cd $MAPD
mv $CLEAND/list  .

```

#### Next step is creating a while loop to go through all the unique fasta files and map them to the reference genome 

```
while read i;
do

  hisat2 -p 6 --dta --phred33       \
    -x "$REFD"/PodMur_index       \
    -1 "$CLEAND"/"$i"_1_paired.fastq  -2 "$CLEAND"/"$i"_2_paired.fastq      \
    -S "$i".sam
    
    samtools view -@ 6 -bS "$i".sam > "$i".bam  
     samtools sort -@ 6  "$i".bam    "$i"_sorted
      samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt


    ```


The option -p indicates the number of processors that we are going to use, dta reports aligments for later use of StringTie to count data. Sam view converts the SAM files to BAM file, and Sam sort, sorts the BAM file into sorted BAM file which is going to be used later by StringTie for read counts. Then, we creat some stats for about the aligned reads and out put them into stats files. The syntax above might differen a little bit if you're using a local computer. 

After creating the sorted_bam file we are ready to use StringTie to produce the read count file. 


### Creating the CSV fiel that contains read counts of the genes 

First, to keep things organized we will move to the count data directory

```
mkdir "$COUNTSD"/"$i"

```

Second, we will use stringTie to produce the read count

```
stringtie -p 6 -e -B -G  "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i"   "$MAPD"/"$i"_sorted.bam

```
Finally, we will close the while loop we created

``` 
done>list 
```

### Copying the files to the results directory 

We will then only move the files that we are going to look at such as the stats files and the read counts files to the results directory

```
cp *.txt $RESULTSD

```
That's to move the stats.txt files to the results directory. 

Then, we will use PrepDe.py script which is also include in the StringTie package to extract read count information from the output StringTie, and creates gene and transcript count matrices in CSV format that can be used for downstream analysis such as differential gene expression, or co-expression analysis. 

``` 
python /home/$MyID/class_shared/prepDE.py -i $COUNTSD

```

then copy this file to the result directory

```
cp *.csv $RESULTSD


```
