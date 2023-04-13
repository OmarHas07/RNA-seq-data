# library dependencies
library(dplyr)
library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(readr)
library(patchwork)
library(RColorBrewer)
library(viridis)


set.seed(666)

setwd("C:/Users/Omar/OneDrive - Auburn University/AT DATA/RNA-Seq Data/R analysis/R")

# read in TPM expression matrix
Exp_table <- read_csv("dataset.csv", col_types = cols())
head(Exp_table)
dim(Exp_table)
# 33548; 143

# read in metadata file
Metadata <- read_csv("MetadataR.csv")
head(Metadata)
dim(Metadata)
#  142;  7




## go from wide to tidy format
Exp_table_long <- Exp_table %>% 
  pivot_longer(cols = !gene_id, names_to = "library", values_to = "FPKM") %>% 
  mutate(logFPKM = log10(FPKM + 1)) 

head(Exp_table_long)



# Make the PCA input data numeric matrix - trouble with this next line
Exp_table_log_wide <- Exp_table_long %>% 
  select(gene_id, library, logFPKM) %>% 
  pivot_wider(names_from = library, values_from = logFPKM) %>% 
  dplyr::group_by(gene_id)

head(Exp_table_log_wide)

# prcomp() performs PCA for you, given a numeric matrix, which is just the transposed Exp_table_log_wide, but without the gene ID column. as.data.frame(t(summary(my_pca)$importance)) saves the standard deviation and proportion of variance into a data table. In this case, the 1st PC accounts for 43% of the variance in this experiment. The 2nd PC accounts for 10% of the variance.
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))


pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

# To make a PCA plot, we will graph the data stored in my_pca$x, which stores the coordinates of each library in PC space. Let's pull that data out and annotate them (with metadata).
PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Sample = row.names(.)) %>% 
  full_join(Metadata %>% 
              select(Sample, Treatment, Timepoint, Group, Type, TypeTime), by = "Sample")
head(PCA_coord)

# defining factors - stopped here for PCA plots, need to edit ###############
#Then plotting th output from the PCA analysis with ggplot this is PC1 versus PC2 by replicate
PCA_by_Group <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Group), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(n = 8, "Accent")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Group") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

#Print the plot as output to view
PCA_by_Group


PCA_by_treatment <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(n = 3, "Accent")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

#Print the plot as output to view
PCA_by_treatment

# This is PC1 versus PC2 colored by tissue type
PCA_by_timepoint2 <- PCA_coord %>% filter(Timepoint == "T2") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 2 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint2

#saving the plot

wrap_plots(PCA_by_timepoint2)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint48 <- PCA_coord %>% filter(Timepoint == "T48") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 48 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint48

#saving the plot

wrap_plots(PCA_by_timepoint48)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint96 <- PCA_coord %>% filter(Timepoint == "T96") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 96 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint96

#saving the plot

wrap_plots(PCA_by_timepoint96)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint144 <- PCA_coord %>% filter(Timepoint == "T144") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 144 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint144

wrap_plots(PCA_by_timepoint144)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

wrap_plots(PCA_by_timepoint2, PCA_by_timepoint48, PCA_by_timepoint96, PCA_by_timepoint144, nrow = 2, ncol = 2)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

# This is PC1 versus PC2 colored by BaseVConjVControl
PCA_by_BaseVConjVControl <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = TypeTime), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Base vs N-glucosides vs Control") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_BaseVConjVControl

# This is PC1 versus PC2 colored by Type
PCA_by_Type <- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = type), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Type") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_Type

#Saving the plots as PNGs 

wrap_plots(PCA_by_Type)
ggsave("../Results/PCA_by_Type.png", height = 4, width = 8.5, bg = "white")

wrap_plots(PCA_by_BaseVConjVControl)
ggsave("../Results/PCA_by_BaseVConjVControl.png", height = 4, width = 8.5, bg = "white")

wrap_plots(PCA_by_timepoint144)
ggsave("../Results/PCA_by_timepoint144.png", height = 4, width = 8.5, bg = "white")

wrap_plots(PCA_by_treatment)
ggsave("../Results/PCA_by_treatment.png", height = 4, width = 8.5, bg = "white")

wrap_plots(PCA_by_Group)
ggsave("../Results/PCA_by_Group.png", height = 4, width = 8.5, bg = "white")



#Now that we have a better handle on the data we will do the co-expression analysis
#Average across the reps
# We start from the long (tidy) table we made earlier. I also pulled the metadata as well to guide the averaging process. 
# By = c("library"="Run) inside full_join() deals with the fact that the library ID is called library 
# in the long table, but Run in the metadata. 
# group_by() followed by summarise(mean = ...) - this takes each gene, tissue, and dev_stage, and computes the mean. 
# The elegance of a tidyverse based workflow is that you do not have to do loops! 
# You let group_by() do the heavy lifting. This could take a moment. This step is doing a lot of mean calculations.

Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord %>% 
              select(Sample, Treatment, Timepoint, Group), 
            by = c("library"="Sample")) %>% 
  group_by(gene_id,Treatment, Timepoint, Group) %>% 
  summarise(mean.logFPKM= mean(logFPKM)) %>% 
  ungroup()  

head(Exp_table_long_averaged)

# Calculate the z-score once you average across the replicates
# The z-score is the difference from mean over standard deviation and standardizes expression patterns (mean = 1, sd = 1)
# In this step, we are grouping by gene. 
# Tissue-stages with higher expression will have a higher z score and vice versa. 
# Note that this is completely relative to each gene itself. 

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(gene_id) %>% 
  mutate(z.score = (mean.logFPKM - mean(mean.logFPKM))/sd(mean.logFPKM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)



Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_id %in% high_var_genes5000$gene_id
         )
head(Exp_table_long_averaged_z_high_var)


# Gene selection
# The next step is correlating each gene to every other gene. However, we have almost 57k genes in this dataset. 
# The number of correlations scales to the square of number of genes. 
# To make things faster and less cumbersome, we can select only the high variance genes. 
# The underlying rationale is if a gene is expressed at a similar level across all samples, it is unlikely that is specifically involved in the biology in a particular stage or tissue.
# There are multiple ways to selecting for high variance genes, and multiple cutoffs. 
# For example, you can calculate the gene-wise variance of logTPM for all genes, and take the upper third. You can only take genes with a certain expression level (say > 5 tpm across all tissues), then take high variance gene. These are arbitrary. You do you.
# Using log(TPM) reduces the bias towards highly expressed genes
# Then I filtered for top 33% high var genes.

#To get the excel file of top 5000 high variance genes
high_var_genes <- Exp_table_long_averaged_z %>% 
  group_by(gene_id) %>% 
  summarise(var = var(mean.logFPKM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes)
dim(high_var_genes)
# 11172     2

#write_xlsx(high_var_genes,"Library\\CloudStorage\\Box-Box\\CSM Leisner Lab\\Ravneet\\JGI\\Coexpanalysis\\high_var_genes.xlsx")

# The above chunk just listed the high var genes, now we need to filter those out in the long table that contains the z-scores.
# For the sake of this example, let's just take top 5000 genes with highest var as a quick exercise
high_var_genes5000 <- high_var_genes %>% 
  slice_max(order_by = var, n = 5000) 

head(high_var_genes5000)

Exp_table_long_averaged_z_high_var %>% 
  group_by(gene_id) %>% 
  count() %>% 
  nrow()

# A good way to check if you have included enough genes in your analyses is to check if your bait genes are among the top var genes.
#high_var_genes5000 %>% 
  #filter(str_detect(geneIDs, Baits$X2[1]))
# 1 Soltu.DM.05G026370.1 0.101

#high_var_genes5000 %>% 
  #filter(str_detect(geneIDs, Baits$X2[2]))
# 1 Soltu.DM.06G025210.1 0.176

#high_var_genes5000 %>% 
  #filter(str_detect(geneIDs, Baits$X2[3]))
# 1 Soltu.DM.05G024030.1 0.624

#high_var_genes5000 %>% 
  #filter(str_detect(geneIDs, Baits$X2[4]))
# 1 Soltu.DM.05G024040.1 0.561

# none of the other bait genes show up.

## IT isn't clearly if the bait genes are in the top 5000 but You could look at the file "high_var_genes5000" and check by hand

# The %in% operator filters geneIDss that are present in high_var_genes5000$geneIDs, thus retaining only high var genes.

# write_csv(high_var_genes5000, "/Users/amg0167/Desktop/Coexpression/high_var_genes5000.csv")

Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_id %in% high_var_genes5000$gene_id)

head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::count() %>% 
  nrow()
# 5000
Exp_table_long_averaged_z_high_var_count <- ncol(Exp_table_long_averaged_z_high_var) - 1
Exp_table_long_averaged_z_high_var_count # this is solution to above code - 6 rows
#5
all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(gene_id) %>% 
  summarise(var = var(mean.logFPKM)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 

#bait_var <- all_var_and_ranks %>% 
  mutate(geneIDs2 = str_sub(geneIDs, start = 1, end = 19)) %>% 
  filter(geneIDs2 %in% Baits$X2) %>% 
  group_by(geneIDs2) %>% 
  slice_max(n = 1, order_by = var)

#bait_var ## this is empty?

all_var_and_ranks %>% 
  ggplot(aes(x = var, y = rank)) +
  geom_rect( 
    xmax = max(high_var_genes5000$var), 
    xmin = min(high_var_genes5000$var),
    ymax = nrow(all_var_and_ranks),
    ymin = nrow(all_var_and_ranks) - 5000,
    fill = "dodgerblue2", alpha = 0.2
  ) +
  geom_line(size = 1.1) +
  geom_hline(
    data = bait_var, aes(yintercept = rank),
    color = "tomato1", size = 0.8, alpha = 0.5
  ) +
  geom_vline(
    data = bait_var, aes(xintercept = var), 
    color = "tomato1", size = 0.8, alpha = 0.5
  ) + 
  labs(y = "rank",
       x = "var(log10(FPKM))",
       caption = "Blue box = top 5000 high var genes.\nRed lines = bait genes.") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  )


# Gene-wise correlation
# We will use the cor() function in R. 
# But the cor() only take vector or matrix as input, so we need to go from long to wide again.

z_score_wide <- Exp_table_long_averaged_z_high_var %>% 
  select(gene_id, Group, z.score) %>% 
  pivot_wider(names_from = Group, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_id
head(z_score_wide)

# The GroupName column contains info for both treatment,cultivar and tissue, which we can recall using the metadata. 
# After long to wide transformation, the GroupName column now becomes the column name of this wide table. 
# Then we produce the correlation matrix. 
# The underlying math here is R takes each column of a matrix and correlates it to every other columns. 
# To get this to work on our wide table, we remove the geneIDs column, transpose it, and feed it into cor().

cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)
# 5000 5000

# Edge selection
# This is building up to a network analysis, where each gene is node, and each correlation is an edge. I have two ways to do this.
# 1. t distribution approximation
# 2. Empirical determination using rank distribution

# For the t distribution you need to calculate the number of tissue, cultivar and treatment combos
number_of_timepoint_tissue_trt <- ncol(z_score_wide) - 1
number_of_timepoint_tissue_trt
## [1] 48
## OR you can find it this way
PCA_coord %>% 
  dplyr::group_by(Timepoint, Treatment) %>% 
  dplyr::count() %>% 
  nrow()
## [1] 48

# Since the correlation matrix is symmetrical along the diagonal we can eliminate redundant data.
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA

# Now we can compute a t statistic from r and compute a p value using the t distribution. 
# This chunk converts the correlation matrix into a data table. 
# Then it goes from wide to long using pivot_longer(). 
# After that, everything is normal dyplr verbs, such as mutate() and filter(). 
# P values are computed using the t distribution. 
# Depending on the sign of t, the upper of lower tail probability is taken. 
# Finally, the p values are adjusted for multiple comparisons using FDR. 
# This step can take a while. Turning a large wide table to a long table always takes a while. Your computer may not have enough memory to run this step if you put in many genes. In this case we only used 5000 genes, so no problem.

edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_matrix)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_of_timepoint_tissue_trt-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_timepoint_tissue_trt-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_timepoint_tissue_trt-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 
head(edge_table)

# You can look at various adjusted p value cutoffs and the corresponding r value before proceeding. 
# Let's say we just look at positively correlated genes.
edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)
# A tibble: 10 × 6
# from                 to                       r     t p.value    FDR
# <chr>                <chr>                <dbl> <dbl>   <dbl>  <dbl>
#  1 Soltu.DM.03G026800.1 Soltu.DM.11G004640.1 0.486  2.36  0.0149 0.0500

edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)
# A tibble: 10 × 6
#from      to            r     t p.value    FDR
#<chr>     <chr>     <dbl> <dbl>   <dbl>  <dbl>
 # 1 AT1G14880 AT2G25950 0.342  2.47 0.00871 0.0100

# So it might make more sense to look at our bait genes and also the R value as a cutoff
# Empirical determination using bait genes and rank distribution
#I just picked two genes I thought were in the high variance group and were both ZIP factors
#edge_table %>% 
  #filter(str_detect(from, "AT5G56970
#") &
         #  str_detect(to,"AT3G48100")) 

edge_table %>% 
  filter(str_detect(from, "AT5G56970") &
           str_detect(to,"AT1G19050") |
           str_detect(from, "AT1G01050") &
           str_detect(to,"AT5G56970")  ) 


  # A tibble: 1 × 6
#from      to             r     t  p.value      FDR
#<chr>     <chr>      <dbl> <dbl>    <dbl>    <dbl>
 # 1 AT1G01050 AT5G56970 -0.862 -11.5 1.94e-15 6.99e-15


# You can also look at the distribution or R values
# Here I randomly sampled 20k edges and plot a histogram. 
# You can plot the whole edge table, but it will take a lot longer to make the graph. 
# When you sample large enough, it does not change the shape of the distribution. 
### Using a r value of 0. doesn't sem to work with the distribution of your data
# r = red line
edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.8, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


ggsave("../Results/r_histogram.png", height = 3.5, width = 5, bg = "white")

# Find positively correlated genes using the r cutoff of 0.8
edge_table_select <- edge_table %>% 
  filter(r >= 0.85)
dim(edge_table_select)
# 2058435             6

# Before we move forward, we can examine the correlation between two bait genes using a scatter plot.
# Here each dot is a library. You can annotate the libraries using metadata, which is now part of PCA_coord. 

Bait_cor_by_Treatment <- z_score_wide %>% 
  filter(gene_id == "AT5G56970" |
           gene_id == "AT1G19050") %>% 
  select(-gene_id) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Group = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Group") %>% 
  ggplot(aes(x = AT5G56970,
             y = AT1G19050)) +
  geom_point(aes(fill = Treatment), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(12, option = "D")) +
  labs(x = "ZIP4 z score",
       y = "ZIP5 z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
Bait_cor_by_Treatment

Bait_cor_by_treatment <- z_score_wide %>% 
  filter(geneIDs == "Soltu.DM.05G024030.1" |
           geneIDs == "Soltu.DM.05G024040.1") %>% 
  select(-geneIDs) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Group = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Group") %>% 
  ggplot(aes(x = Soltu.DM.05G024030.1,
             y = Soltu.DM.05G024040.1)) +
  geom_point(aes(fill = Treatment), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(2, option = "D")) +
  labs(x = "ZIP4 z score",
       y = "ZIP5 z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
Bait_cor_by_treatment

Bait_cor_by_timepoint <- z_score_wide %>% 
  filter(geneIDs == "Soltu.DM.05G024030.1" |
           geneIDs == "Soltu.DM.05G024040.1") %>% 
  select(-geneIDs) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Group = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Group") %>% 
  ggplot(aes(x = Soltu.DM.05G024030.1,
             y = Soltu.DM.05G024040.1)) +
  geom_point(aes(fill = Timepoint), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(4, option = "D")) +
  labs(x = "ZIP4 z score",
       y = "ZIP5 z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
Bait_cor_by_timepoint

wrap_plots(Bait_cor_by_tissue, Bait_cor_by_treatment,  Bait_cor_by_timepoint, nrow = 2)

# Module detection
# We will be using the Leiden algorithm to detect module, which is a graph based clustering method
# We will be using igraph to do some of the downstream analyses. 
# While you can get Leiden as a standalone package, Leiden is also part of the igraph package. 
# The first thing to do is producing a graph object, also known as a network object.

# To make a graph object, you need a edge table. We already made that, which is edge_table_select, a edge table that we filtered based on some kind of r cutoff. 
# Optionally, we can also provide a node table, which contains information about all the notes present in this network. We can make that.
# We need to two things.
# Non-redundant gene IDs from the edge table.
# Functional annotation, which I downloaded.
# Functional annotation
#funct_anno <- read_delim("DM_1-3_516_R44_potato.v6.1.func_annot_transcripts.txt", 
                         #delim = "\t", escape_double = FALSE, 
                       #  trim_ws = TRUE)

funct_anno <- read_csv("FunctionalAnnotation.csv", col_names = F)

head(funct_anno)

node_table <- data.frame(
  geneIDs = c(edge_table_select$from, edge_table_select$to) %>% unique()
) %>% 
  left_join(funct_anno, by = c("geneIDs"="X1")) %>% 
 rename(functional_annotation = X2)
head(node_table)
dim(node_table)
# 90786    2

write_csv(node_table, "func_anno_node.csv")

# We have 4623 genes in this network, Note that this is less than the 5000 top var genes we put in, because we filtered out some edges.
# Now let's make the network object.
# Note that I selected the directed = F argument, because we made our network using correlation. Correlation is non-directional, because cor(A,B) = cor(B,A).

my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F)
##THIS CODE IS WORKING 
##the error that you was getting before was because you need 
##to reference a unique node data in vertices. 
my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = unique(c(edge_table_select$from,edge_table_select$to)),
  directed = F)

# Graph-based clustering
# cluster_leiden() runs the Leiden algorithm for you. resolution_parameter controls how many clusters you will get. 
# The larger it is, the more clusters. You can play around with the resolution and see what you get. 
modules <- cluster_leiden(my_network, resolution_parameter = 1.5, 
                          objective_function = "modularity")

# use heuristics to find optimal number of modules
# Here I wrote a function to detect module, pull out number of modules that have >= 5 genes, and count number of genes contained in modules that have >= 5 genes. All in one function.

optimize_resolution <- function(network, resolution){
  modules = network %>% 
    cluster_leiden(resolution_parameter = resolution,
                   objective_function = "modularity")
  
  parsed_modules = data.frame(
    geneIDs = names(membership(modules)),
    module = as.vector(membership(modules)) 
  )
  
  num_module_5 = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    nrow() %>% 
    as.numeric()
  
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  c(num_module_5, num_genes_contained)
  
}

# Then I can test a list of resolutions in this function. Let's test a range of resolution from 0.25 to 5, in steps of 0.25.

optimization_results <- purrr::map_dfc(
  .x = seq(from = 0.3, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.3, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)
head(optimization_results)

# We have the results organized into one tidy data table. We can graph it.
Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 3, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 3, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)
ggsave("../Results/Num.Modules.png", height = 4, width = 8.5, bg = "white")

# ggsave("../Results/Optimize_resolution.png", height = 5, width = 3.2, bg ="white")
# ggsave("../Results/Optimize_resolution.png", height = 5, width = 3.2, bg ="white")

# Let's say we move on with module detection using a resolution of 2. Next, we need to link the module membership to the gene IDs.

my_network_modules <- data.frame(
  geneIDs = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(node_table, by = "geneIDs")

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>%
  print (n = 26)
#    module     n
# <dbl> <int>
#   1      5   574
#   2      7   564
#   3      1   498
#   4     10   474
#   5      8   458
#   6      6   415
#   7      2   370
#   8      9   365
#   9     11   186
#   10      3    58
#   11      4    52
#   12     17    23
#   13     12     5

# new:
#  1      4   575
#  2      6   528
#  3      1   487
#  4      8   473
#  5      7   416
#  6      5   408
#  7      3   389
#  8      2   382
#  9     11   186
#  10      9    96
#  11     10    60
#  12     18    24
#  13     14     7
#  14     38     6
#  15     12     5

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))
#
#A tibble: 14 × 2
# Groups:   module [14]
#module     n
#<dbl> <int>
#  1      4 29946

# Moving forward we will only use modules that have 5 or more genes.
module_5 <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules <- my_network_modules %>% 
  filter(module %in% module_5$module)
head(my_network_modules)


# Module quality control
my_network_modules %>% 
  filter( geneIDs == "AT5G62920" |
            geneIDs == "AT3G48100" |
            geneIDs == "AT2G46310" |
            geneIDs == "AT4G29740")
#these are ARR5,6 and CRF5
#all are present in module 1 




# Module treatment correspondence
# the essence of this workflow is simple, so we will use a simple method: peak expression.

#unify columns names first. 

names(Exp_table_long_averaged_z_high_var)
names(Exp_table_long_averaged_z_high_var)[1] <- "geneIDs"

Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_averaged_z_high_var %>% 
  inner_join(my_network_modules, by = "geneIDs")

head(Exp_table_long_averaged_z_high_var_modules)

# Now we can produce summary statistics for each cluster and look at their expression pattern using mean.
modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, Timepoint, Treatment, Group) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)

# Then we look at at which developmental stage and tissue is each module most highly expressed.
module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp

# You can also QC the clusters via a line graph. It will be too much to look at if graph all the modules, so let's just pick 2.
# This code chunk is very long, because a few things:
# I reordered x-axis to reflect the biological time sequence.
# Overlaid the average of clusters.
# Added a color strip at the bottom to annotate stages, which reduces the amount of text on the figure.
#"Control",	"cZ",	"cZ9G",	"DZ",	"DZ7G",	"DZ9G",	"iP","iP7G","iP9G",	"tZ","tZ7G","tZ9G"
######
#Works
#####

#str_detect(Timepoint, "T2") ~ 1,
#str_detect(Timepoint, "T48") ~ 2,
#str_detect(Timepoint, "T96") ~ 3,
#str_detect(Timepoint, "T144") ~ 4,
module_line_plot <- Exp_table_long_averaged_z_high_var_modules %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "DZ") ~ 3,
    str_detect(Treatment, "iP") ~ 4,
    str_detect(Treatment, "cZ") ~ 5,
    str_detect(Treatment, "tZ7G") ~ 6,
    str_detect(Treatment, "tZ9G") ~ 7,
    str_detect(Treatment, "DZ7G") ~ 8,
    str_detect(Treatment, "DZ9G") ~ 9,
    str_detect(Treatment, "iP7G") ~ 10,
    str_detect(Treatment, "iP9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12,
  
  )) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  filter(module == "1" |
           module == "5") %>% 
  ggplot(aes(x = Treatment, y = z.score)) +
  #do it by treatment
  facet_grid(module ~ Timepoint) +
  geom_line(aes(group = geneIDs), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "1" |
               module == "5" ) %>% 
      mutate(order_x = case_when(
        str_detect(Treatment, "Control") ~ 1,
        str_detect(Treatment, "tZ") ~ 2,
        str_detect(Treatment, "DZ") ~ 3,
        str_detect(Treatment, "iP") ~ 4,
        str_detect(Treatment, "cZ") ~ 5,
        str_detect(Treatment, "tZ7G") ~ 6,
        str_detect(Treatment, "tZ9G") ~ 7,
        str_detect(Treatment, "DZ7G") ~ 8,
        str_detect(Treatment, "DZ9G") ~ 9,
        str_detect(Treatment, "iP7G") ~ 10,
        str_detect(Treatment, "iP9G") ~ 11,
        str_detect(Treatment, "cZ9G") ~ 12,
      )) %>% 
      mutate(Treatment = reorder(Treatment, order_x)),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_lines_color_strip <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "DZ") ~ 3,
    str_detect(Treatment, "iP") ~ 4,
    str_detect(Treatment, "cZ") ~ 5,
    str_detect(Treatment, "tZ7G") ~ 6,
    str_detect(Treatment, "tZ9G") ~ 7,
    str_detect(Treatment, "DZ7G") ~ 8,
    str_detect(Treatment, "DZ9G") ~ 9,
    str_detect(Treatment, "iP7G") ~ 10,
    str_detect(Treatment, "iP9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12,
  )) %>% 
  mutate(stage = factor(Treatment, levels = c(
   "Control",	"cZ",	"cZ9G",	"DZ",	"DZ7G",	"DZ9G",	"iP","iP7G","iP9G",	"tZ","tZ7G","tZ9G"
  ))) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(. ~ Timepoint) +
  geom_tile(aes(fill = Treatment)) +
  scale_fill_manual(values = viridis(12, option = "D")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))


 ggsave("../Results/module1and5_line_plots.png", height = 4, width = 8.2, bg = "white")



# Module line plot for specific modules
module_line_plot <- Exp_table_long_averaged_z_high_var_modules %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "cZ") ~ 3,
    str_detect(Treatment, "DZ") ~ 4,
    str_detect(Treatment, "iP") ~ 5,
    str_detect(Treatment, "DZ7G") ~ 6,
    str_detect(Treatment, "DZ9G") ~ 7,
    str_detect(Treatment, "iP7G") ~ 8,
    str_detect(Treatment, "iP9G") ~ 9,
    str_detect(Treatment, "tZ7G") ~ 10,
    str_detect(Treatment, "tZ9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12, )) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  filter(module == "7") %>% 
  ggplot(aes(x = Treatment, y = z.score)) +
  facet_grid(module ~ Timepoint) +
  geom_line(aes(group = geneIDs), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "7") %>% 
      mutate(order_x = case_when(
        str_detect(Treatment, "Control") ~ 1,
        str_detect(Treatment, "tZ") ~ 2,
        str_detect(Treatment, "cZ") ~ 3,
        str_detect(Treatment, "DZ") ~ 4,
        str_detect(Treatment, "iP") ~ 5,
        str_detect(Treatment, "DZ7G") ~ 6,
        str_detect(Treatment, "DZ9G") ~ 7,
        str_detect(Treatment, "iP7G") ~ 8,
        str_detect(Treatment, "iP9G") ~ 9,
        str_detect(Treatment, "tZ7G") ~ 10,
        str_detect(Treatment, "tZ9G") ~ 11,
        str_detect(Treatment, "cZ9G") ~ 12,
      )) %>% 
      mutate(Treatment = reorder(Treatment, order_x)),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))



## heat map - works?

modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, 0.95)

# The 95th percentile of averaged z score is 1.3. We can probably roughly clipped the z-scores at 1.5 or -1.5
modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 1.5 ~ 1.5,
    mean.z < -1.5 ~ -1.5,
    T ~ mean.z
  ))

modules_mean_z <- modules_mean_z %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4,
  )) %>%  
  mutate(Timeopint = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Treatment, start = 1, end = 2),
    T ~ Timepoint
  )) %>% 
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) 

head(modules_mean_z)

#########################################################
module_heatmap <- modules_mean_z %>%
  ggplot(aes(x = Treatment, y = as.factor(module))) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z.clipped), color = "grey80") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), labels = c("< -1.5", "0", "> 1.5")) +
  labs(x = NULL,
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    panel.spacing = unit(0.5, "lines")
  )

heatmap_color_strip1 <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint),
  stringsAsFactors = F
) %>%
  #  filter(Cultivar!= "Root") %>%
  # filter(str_detect(Cultivar, "Flyer|Clark|Loda") == F) %>%
  # filter((Tissue == "Leaf" &
  #           str_detect(Cultivar, "Flyer|Clark|Loda"))==F) %>%
  # filter((str_detect(Tissue, "Pod") &
  #           str_detect(Cultivar, "Flyer|Clark|Loda"))==F) %>%
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>%
  mutate(stage = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Timepoint, start = 1, end = 2),
    T ~ Timepoint
  )) %>%
  mutate(stage = factor(stage, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>%
  
  mutate(Timepoint = reorder(Timepoint, order_x)) %>%
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Treatment)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

heatmap_color_strip2 <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint),
  stringsAsFactors = F
) %>%
  # filter(Tissue != "Leaf") %>%
  # filter(str_detect(Cultivar, "Flyer|Clark|Loda") == F) %>%
  # filter((Tissue == "Pod" &
  #           str_detect(Cultivar, "Flyer|Clark|Loda"))==F) %>%
  # filter((str_detect(Tissue, "Root") &
  #           str_detect(Cultivar, "Flyer|Clark|Loda"))==F) %>%
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>%
  mutate(stage = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Timepoint, start = 1, end = 2),
    T ~ Timepoint
  )) %>%
  mutate(stage = factor(stage, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>%
  mutate(Timepoint = reorder(Timepoint, order_x)) %>%
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Timepoint)) +
  scale_fill_manual(values = viridis(4, option = "D")) +
  labs(fill = "Timepoint") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )
wrap_plots(module_heatmap, heatmap_color_strip1, heatmap_color_strip2,
           nrow = 3, heights = c(1, 0.08, 0.08), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggsave("../Results/module_heatmap.png", height = 4.8, width = 10, bg = "white")
#######

module_peak_exp <- module_peak_exp %>% 
  mutate(order_y = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>%  
  mutate(peak_exp = reorder(Timepoint, order_y)) 

modules_mean_z_reorded <- modules_mean_z %>% 
  full_join(module_peak_exp %>% 
              select(module, peak_exp, order_y), by = c("module")) %>% 
  mutate(module = reorder(module, -order_y))

head(modules_mean_z_reorded)
##start
###
###
###
###
module_heatmap <- modules_mean_z_reorded %>% 
  ggplot(aes(x = Treatment, y = as.factor(module))) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z.clipped), color = "grey80") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), labels = c("< -1.5", "0", "> 1.5")) +
  labs(x = NULL,
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    panel.spacing = unit(0.5, "lines") 
  )

heatmap_color_strip1 <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 

  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>% 
  mutate(Timepoint = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Timepoint, start = 1, end = 2),
    T ~ Timepoint
  )) %>% 
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Treatment)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

heatmap_color_strip2 <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>% 
  mutate(Timepoint = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Timepoint, start = 1, end = 2),
    T ~ Timepoint
  )) %>% 
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Treatment)) +
  scale_fill_manual(values = viridis(12, option = "D")) +
  labs(fill = "Treatment") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )


wrap_plots(module_heatmap, heatmap_color_strip1, heatmap_color_strip2, 
           nrow = 3, heights = c(1, 0.08, 0.08), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggsave("../Results/module_heatmap.svg", height = 4.8, width = 10, bg = "white")
ggsave("../Results/module_heatmap.png", height = 4.8, width = 10, bg = "white")


# Network graph
# pull out neighbors of bait genes

neighbors_of_bait <- c(
  neighbors(my_network, v = "AT3G44880"), # PAO
  neighbors(my_network, v = "AT1G17020"), # SRG1 
  neighbors(my_network, v = "AT4G26150") #  GATA22
  #neighbors(my_network, v = "Soltu.DM.10G025990.1") # seed specific - "oleosin"
) %>% 
  unique()  

length(neighbors_of_bait)
# 169

# We can make a sub-network object. First we subset edges in the network.

subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait) &
           to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

subnetwork_genes <- c(subnetwork_edges$from, subnetwork_edges$to) %>% unique()
length(subnetwork_genes)
dim(subnetwork_edges)
# 592
# 1456    6

# We can constrain the edges such that both the start and end of edges are neighbors of baits. 
# I also filtered for highly correlated neighbors (top 5 edges/node based on r value). 
# We still have 5714 edges and 2195 nodes. Note that the most correlated edges for each bait many have overlaps, 
# so the total number of edges remaining will be less than what you think.

# Then we subset nodes in the network.

subnetwork_nodes <- node_table %>% 
  filter(geneIDs %in% subnetwork_genes) %>% 
  left_join(my_network_modules, by = "geneIDs") %>% 
  left_join(module_peak_exp, by = "module") %>% 
  mutate(module_annotation = case_when(
    str_detect(module, "11|6|12") ~ "T2",
    module == "17" ~ "T48",
    module == "6" ~ "T96",
    T ~ "T144"
  ))

dim(subnetwork_nodes)

# Then make sub-network object from subsetted edges and nodes.

my_subnetwork <- graph_from_data_frame(subnetwork_edges,
                                       vertices = subnetwork_nodes,
                                       directed = F)

my_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "white", shape = 21, size = 2,
                  aes(fill = module_annotation)) + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent")[c(1,3,6)], "grey30"),
                    limits = c("Small.Tuber", "Medium.Tuber", "Big.Tuber")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )




# save gene lists of modules

Module7 <- my_network_modules %>%
  filter(module == "7")

#write_excel_csv(Module7, "Module7_genes.csv", col_names = T)





#

neighbors_of_SP6A <- c(
  neighbors(my_network, v = "Soltu.DM.05G026370.1"), # SP6A
  neighbors(my_network, v = "Soltu.DM.05G024030.1"), # SP5G-A
  neighbors(my_network, v = "Soltu.DM.05G024040.1") # SP5G-B
) %>% 
  unique()  

length(neighbors_of_SP6A)
# 169

my_TFs <- my_network_modules %>% 
  filter(geneIDs %in% names(neighbors_of_SP6A)) %>% 
  filter(str_detect(functional_annotation, "PEBP"))

TF_TPM <- Exp_table_long %>% 
  filter(geneIDs %in% my_TFs$geneIDs) %>% 
  inner_join(PCA_coord, by = c("library"="Sample")) %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "AmbientT") ~ 1,
    str_detect(Treatment, "ElevatedT") ~ 2,
  )) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  ggplot(aes(x = Treatment, y = logTPM)) +
  facet_grid(geneIDs ~ Tissue, scales = "free_y") +
  geom_point(aes(fill = Tissue), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(group = geneIDs), 
               fun = mean, alpha = 0.8, size = 1.1, color = "grey20") +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.background = element_blank()
  )

wrap_plots(TF_TPM, module_lines_color_strip, 
           nrow = 2, heights = c(1, 0.05))



