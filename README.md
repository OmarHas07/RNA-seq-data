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

## loading the packages 

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

Because I am using the ASC I will set the stack size to unlimited which is a indicates the memories we are using that stores the temporary information like variables and function calling. However, I don't recommend setting the stack size to  unlimited if you're using a local computer or a VM, as most operating systems might face stack overflow and shutdown the whole program. 

```
ulimit -s unlimited
``` 

