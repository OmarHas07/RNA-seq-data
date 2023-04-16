# Load required libraries
library(DESeq2)
library(tidyverse)

# Replace this with the path to your gene read count file
gene_read_count_file <- "gene_read_counts.csv"

# Replace this with the path to your gene lengths file
gene_lengths_file <- "gene_lengths.csv"

# Load gene read count data and gene lengths
read_counts <- read_csv(gene_read_count_file) %>%
  column_to_rownames("gene_id")

gene_lengths <- read_csv(gene_lengths_file) %>%
  column_to_rownames("gene_id")

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = DataFrame(condition = factor(rep("A", ncol(read_counts)))),
                              design = ~ condition)
dds <- DESeq(dds, fitType="local")

# Normalize read counts using DESeq2
normalized_counts <- counts(dds, normalized = TRUE)

# Calculate FPKM values
gene_lengths_kb <- gene_lengths / 1000
total_mapped_reads_millions <- sum(normalized_counts) / 1e6
fpkm <- sweep(normalized_counts, 1, gene_lengths_kb, FUN = "/") / total_mapped_reads_millions

# Save FPKM values as a CSV file
write.csv(fpkm, file = "fpkm_values.csv")
