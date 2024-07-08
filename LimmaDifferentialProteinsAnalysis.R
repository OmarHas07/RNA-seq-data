#### WD 

### 2 Hours analysis 
### numbers of DAPs for each treatment against control (signficant) is much higher than in 48 hours. 

setwd("/Users/omarhasannin/Documents/Research Projects/Chapter 1/Proteome data /DEP")

### lIbraries 
library(edgeR)
library(limma)
library(dplyr)
library(tidyverse)
library(dplyr)
library(tidyverse)
library(igraph)
library(ggraph)
library(readxl)
library(readr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(limma)
### Load Data
abundance_data <- read.csv("Abundances.csv", header = TRUE, row.names = 1)
head(abundance_data)
metadata <- read.csv("MetaDataT2.csv") 
metadata<- metadata %>% 
  dplyr::select(Sample, Treatment)
head(metadata)
# Ensure that the columns in abundance_data are in the same order as the rows in metadata
abundance_data <- abundance_data[, metadata$Sample]

###### Keeping only ProteinIDs that have at least 3 replicates higher than 0.01 in each of the treatments 
# Convert row names to a column 
abundance_data$ProteinID <- rownames(abundance_data)


# Convert metadata$Treatment to a factor to ensure consistent ordering
metadata$Treatment <- factor(metadata$Treatment, levels = unique(metadata$Treatment))

# Initialize a logical vector to store if each protein passes the criterion
protein_passes_criterion <- rep(FALSE, nrow(abundance_data))

# Loop through each treatment
for(treatment in levels(metadata$Treatment)) {
  # Get the columns (samples) that belong to this treatment
  treatment_samples <- metadata$Sample[metadata$Treatment == treatment]
  
  # Count the valid data points for each protein in this treatment
  valid_counts <- rowSums(abundance_data[, treatment_samples, drop = FALSE] > 0.01)
  
  # Update the proteins that pass the criterion
  protein_passes_criterion <- protein_passes_criterion | (valid_counts >= 3)
}

# Filter the abundance data
filtered_abundance_data <- abundance_data[protein_passes_criterion, ]

# Check the resulting filtered data
print(head(filtered_abundance_data)) ### Worked
### it removed 89 protein IDs, now 4159 out of 4246 present in the data. 

### Remove ProteinID column
filtered_abundance_data <- filtered_abundance_data %>%  
  dplyr::select(-ProteinID)

### Saving the filtered data 

write_csv(filtered_abundance_data, "0.01>_validC_>=_3_abundances_2Hours.csv")
######################################################################
##### 


log2_data <- log2(filtered_abundance_data + 1)



# Log-transform the data and apply the voom transformation
#v <- voom(filtered_data, plot=TRUE)

# Convert 'Treatment' to a factor and ensure its levels
metadata$Treatment <- factor(metadata$Treatment)
levels(metadata$Treatment)

# Create the design matrix. Assuming 'Treatment' is the column in your metadata
design <- model.matrix(~ 0 + Treatment, data=metadata)
colnames(design) <- levels(metadata$Treatment)

# Fit the linear model
fit <- lmFit(log2_data, design)



# Now, create the contrast matrix
contrast.matrix <- makeContrasts(
  vsControl_tZ = tZ - Control,
  vsControl_DZ = DZ - Control,
  vsControl_cZ = cZ - Control,
  vsControl_iP = iP - Control,
  vsControl_tZ7G = tZ7G - Control,
  vsControl_tZ9G = tZ9G - Control,
  vsControl_DZ7G = DZ7G - Control,
  vsControl_DZ9G = DZ9G - Control,
  vsControl_iP7G = iP7G - Control,
  vsControl_iP9G = iP9G - Control,
  vsControl_cZ9G = cZ9G - Control,
  levels = design
)




# Apply the contrasts to the linear model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract the results for each contrast
tZ <- topTable(fit2, coef="vsControl_tZ", n=Inf, adjust="BH")
DZ <- topTable(fit2, coef="vsControl_DZ", n=Inf, adjust="BH")
iP <- topTable(fit2, coef="vsControl_iP", n=Inf, adjust="BH")
cZ <- topTable(fit2, coef="vsControl_cZ", n=Inf, adjust="BH")
tZ7G <- topTable(fit2, coef="vsControl_tZ7G", n=Inf, adjust="BH")
tZ9G <- topTable(fit2, coef="vsControl_tZ9G", n=Inf, adjust="BH")
DZ7G <- topTable(fit2, coef="vsControl_DZ7G", n=Inf, adjust="BH")
DZ9G <- topTable(fit2, coef="vsControl_DZ9G", n=Inf, adjust="BH")
iP7G <- topTable(fit2, coef="vsControl_iP7G", n=Inf, adjust="BH")
iP9G <- topTable(fit2, coef="vsControl_iP9G", n=Inf, adjust="BH")
cZ9G <- topTable(fit2, coef="vsControl_cZ9G", n=Inf, adjust="BH")


### Testing 
tZviP <- topTable(fit2, coef="vstZ_iP", n=Inf, adjust="BH") %>% 
  filter(abs(logFC) > 1.5, adj.P.Val < 0.05) ### five genes

tZ7GvtZ <- topTable(fit2, coef="vstZ7G_tZ", n=Inf, adjust="BH") %>% 
  filter(abs(logFC) > 1.5, adj.P.Val < 0.05) ### 2 Genes
########################

# List of all result data frames
result_list <- list(tZ=tZ, DZ=DZ, iP=iP, cZ=cZ, tZ7G=tZ7G, tZ9G=tZ9G, 
                    DZ7G=DZ7G, DZ9G=DZ9G, iP7G=iP7G, iP9G=iP9G, cZ9G=cZ9G)

# Apply the filtering to each data frame
filtered_results <- lapply(result_list, function(df) {
  df %>% 
    filter(abs(logFC) > 1.5, adj.P.Val < 0.05)
})

# Assuming list_of_dfs is your list of dataframes
annotated_dfs <- lapply(filtered_results, function(filtered_results) {
  # Extract gene IDs (row names)
  IDs <- rownames(filtered_results)
  gene_ids_no_version <- gsub("\\..*$", "", IDs)
  
  # Get annotations for gene IDs
  annotations <- select(org.At.tair.db, 
                        keys = gene_ids_no_version, 
                        columns = c("SYMBOL", "GENENAME", "ARACYC"), 
                        keytype = "TAIR")
  
  # Keep only the first gene name per gene ID
  unique_annotation <- annotations %>%
    dplyr::group_by(TAIR) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Merge annotations back with the original dataframe
  # Ensure that the TAIR column and gene_ids_no_version are correctly aligned for merging
  # This step may require adjustments based on how your data and annotations are structured
  # For example, you might need to add TAIR as a column to your dataframe before merging
  filtered_results$TAIR <- gene_ids_no_version # This line might need adjustment
  merged_df <- merge(filtered_results, unique_annotation, by.x = "TAIR", by.y = "TAIR", all.x = TRUE)
  
  return(merged_df)
})


# Now filtered_results is a list of filtered data frames

tZ_filtered <- (annotated_dfs$tZ) 
DZ_filtered <- (annotated_dfs$DZ) 
iP_filtered <- (annotated_dfs$iP) 
cZ_filtered <- (annotated_dfs$cZ) 
tZ7G_filtered <- (annotated_dfs$tZ7G) 
tZ9G_filtered <- (annotated_dfs$tZ9G) 
DZ7G_filtered <- (annotated_dfs$DZ7G) 
DZ9G_filtered <- (annotated_dfs$DZ9G) 
iP7G_filtered <- (annotated_dfs$iP7G) 
iP9G_filtered <- (annotated_dfs$iP9G) 
cZ9G_filtered <- (annotated_dfs$cZ9G) 

tZ_filtered$ProteinIDs <- rownames(tZ_filtered)
DZ_filtered$ProteinIDs <- rownames(DZ_filtered)
iP_filtered$ProteinIDs <- rownames(iP_filtered)
cZ_filtered$ProteinIDs <- rownames(cZ_filtered)
tZ7G_filtered$ProteinIDs <- rownames(tZ7G_filtered)
tZ9G_filtered$ProteinIDs <- rownames(tZ9G_filtered)
DZ7G_filtered$ProteinIDs <- rownames(DZ7G_filtered)
DZ9G_filtered$ProteinIDs <- rownames(DZ9G_filtered)
iP7G_filtered$ProteinIDs <- rownames(iP7G_filtered)
iP9G_filtered$ProteinIDs <- rownames(iP9G_filtered)
cZ9G_filtered$ProteinIDs <- rownames(cZ9G_filtered)





write_csv(tZ_filtered, "tZ_filtered.csv")
write_csv(DZ_filtered, "DZ_filtered.csv")
write_csv(iP_filtered, "iP_filtered.csv")
write_csv(cZ_filtered, "cZ_filtered.csv")
write_csv(tZ7G_filtered, "tZ7G_filtered.csv")
write_csv(tZ9G_filtered, "tZ9G_filtered.csv")
write_csv(DZ7G_filtered, "DZ7G_filtered.csv")
write_csv(DZ9G_filtered, "DZ9G_filtered.csv")
write_csv(iP7G_filtered, "iP7G_filtered.csv")
write_csv(iP9G_filtered, "iP9G_filtered.csv")
write_csv(cZ9G_filtered, "cZ9G_filtered.csv")

#####################################################\\\


## Volcano plot ## testing


library(ggplot2)



# For example, "logFC" for log fold change and "adj.P.Val" for adjusted p-values
tZ_filtered$LogFC <- tZ_filtered$logFC  # Adjust column name 
tZ_filtered$P.Value <- tZ_filtered$adj.P.Val  # Adjust column name 



# Creating the volcano plot
tZ <- ggplot(tZ_filtered, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "tZ 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")
tZ

#### iP
# For example, "logFC" for log fold change and "adj.P.Val" for adjusted p-values
iP_filtered$LogFC <- iP_filtered$logFC  # Adjust column name 
iP_filtered$P.Value <- iP_filtered$adj.P.Val  # Adjust column name 



# Creating the volcano plot
iP <- ggplot(iP_filtered, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "iP 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")

wrap_plots(tZ, iP)

##### DZ

# For example, "logFC" for log fold change and "adj.P.Val" for adjusted p-values
DZ_filtered$LogFC <- DZ_filtered$logFC  # Adjust column name 
DZ_filtered$P.Value <- DZ_filtered$adj.P.Val  # Adjust column name 

DZ <- ggplot(DZ_filtered, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "DZ 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")

#### cZ96

cZ96_filtered$LogFC <- cZ96_filtered$logFC  # Adjust column name 
cZ96_filtered$P.Value <- cZ96_filtered$adj.P.Val  # Adjust column name 

cZ96 <- ggplot(cZ96_filtered, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "cZ96 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")


cZ96_filtered$LogFC <- cZ96_filtered$logFC  # Adjust column name 
cZ96_filtered$P.Value <- cZ96_filtered$adj.P.Val  # Adjust column name 

cZ96 <- ggplot(cZ96_filtered, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "cZ96 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")




wrap_plots(tZ, iP, cZ96, DZ)
ggsave("Volcano Plots BaseForms.png", height = 4.5, width = 8.5, bg = "white")



tZviP


tZviP$LogFC <- tZviP$logFC  # Adjust column name 
tZviP$P.Value <- tZviP$adj.P.Val  # Adjust column name 

tZviP <- ggplot(tZviP, aes(x = LogFC, y = -log10(adj.P.Val))) + 
  geom_point(alpha = 0.4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove background lines
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) +
  labs(title = "tZviP 2 Hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-Value)") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed")



tZviP
