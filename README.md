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
