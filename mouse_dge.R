#!/usr/bin/env Rscript

# DGE Analysis with DESeq2 and edgeR
# Input: gene_counts.txt (featureCounts output)
# Output: significant gene list, visualizations, filtered gene_counts

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("edgeR")

setwd("~/Desktop/")
print(getwd())

library(dplyr)
library(biomaRt)
library(DESeq2)
library(edgeR)

# For visualizations
library(limma) 
library(ggplot2)
library(pheatmap)

# Load data
gene_counts <- read.table("~/Desktop/mouse_data/fcounts_results/gene_counts.txt", header = TRUE, sep = "\t", row.names = 1) 

# Step 1. Clean up column names

# Display column names
print(colnames(gene_counts)) 

# Drop the unnecessary columns
df <- gene_counts[, !(names(gene_counts) %in% c("Chr", "Start", "End", "Strand", "Length"))]  

current_columns <- colnames(df) # current column names

cat("Current columns of the data, gene_counts: \n")
print(current_columns)

# Remove “mouse_data.star_results.renamed.” from the column names
new_columns <- gsub("mouse_data\\.star_results\\.renamed\\.", "", current_columns)

# Remove “.bam” from column names
newer_columns <- gsub("\\.bam", "", new_columns)

# Assign new column names
colnames(df) <- newer_columns

cat("Cleaned up column names and data: \n")
head(df) # glance at dataset

# Step 2. Filter gene_counts data

# Filter for mouse protein-coding genes only
ensembl <- useMart(biomart = "ensembl") # connect to the Ensembl database using biomaRt
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

ensembl_gene_ids <- rownames(df) # Retrieve gene biotypes for your Ensembl IDs
gene_info <- getBM(filters = "ensembl_gene_id_version", 
                   values = ensembl_gene_ids, 
                   attributes = c("ensembl_gene_id_version", "gene_biotype"), 
                   mart = ensembl) 

# Filter for counts of protein-coding genes 
protein_coding_genes <- gene_info %>% filter(gene_biotype == "protein_coding")
protein_coding_counts <- df[rownames(df) %in% protein_coding_genes$ensembl_gene_id, ] 
dim(protein_coding_counts) # [1] 21773 (genes)   10 (samples)

# Quality Control:

# Remove genes with low counts (decreases noise):
# Only genes with at least 5 counts present in at least 25% of samples are included 
# in differential expression analyses.
df <- protein_coding_counts 
min_count <- 5
sample_threshold <- ceiling(0.25 * ncol(df))  # 25% of samples
keep_genes <- rowSums(df >= min_count) > sample_threshold
filtered_counts_df <- df[keep_genes, ]

cat("Filtered gene counts\n")
head(filtered_counts_df)

saveRDS(filtered_counts_df, "mouse_data/filtered_gene_counts.rds")

###################
# DGE with DESeq2 #
###################

print(getwd())
setwd("~/Desktop")

filtered_counts <- readRDS("mouse_data/filtered_gene_counts.rds")

# Note: Unpaired design means mice samples are independent of each other. 

# Step 1: Build metadata file (for an unpaired design) 
num_samples <- 5 # Assuming each condition has 5 samples

sample <- colnames(filtered_counts) 
condition <- factor(c(rep("base", num_samples), rep("drug", num_samples))) # Same order as gene_counts columns 
metadata <- data.frame(sample, condition) # unpaired design
metadata$condition <- relevel(metadata$condition, ref = "base")  # Set the reference condition

# Print the metadata
print(metadata)

# Step 2. Run DESeq2 algorithm 
unpaired_design_filtered <- DESeqDataSetFromMatrix(countData = filtered_counts, 
                                                   colData = metadata, 
                                                   design = ~ condition) # unpaired design
deseq2_results_filtered <- DESeq(unpaired_design_filtered, test = "Wald")
deseq2_results_table_filtered <- results(deseq2_results_filtered)

# Print headers and first 10 rows of DESEq2 results
cat("DESeq2 results\n")
head(deseq2_results_table_filtered)

# Step 3. Create a list of the Significant Genes (via DESeq2)

# Filter the genes which have an adjusted p-value (Benjamini-Hochberg) < 0.05 and log2(fold change) >= 1.00:
sig_genes_deseq2 <- deseq2_results_table_filtered[ !is.na(deseq2_results_table_filtered$padj) & 
                                                     (deseq2_results_table_filtered$padj < 0.05) & 
                                                     (abs(deseq2_results_table_filtered$log2FoldChange) >= 1.00), ] 

# Print number of significant genes (adjusted p-value < 0.05 and log2(Fold Change) >= 1.00)
print("Number of significant genes (adjusted p-value < 0.05 and log2(Fold Change) >= 1.00): ")
nrow(sig_genes_deseq2) 

saveRDS(sig_genes_deseq2, "mouse_data/sig_genes_deseq2.rds") # Save significant genes to a file

deseq2_results_table <- deseq2_results_table_filtered[!is.na(deseq2_results_table_filtered$padj), ]

saveRDS(deseq2_results_table, "mouse_data/res_deseq2.rds") # Save DESeq2 output

##################
# DGE with edgeR #
##################

# Step 1. Build design matrix for edgeR

# Create an DGE list for edgeR from the gene_counts data
dge <- DGEList(counts = filtered_counts) 

# Ensure your metadata has the sample names as row names
rownames(metadata) <- metadata$sample # use metadata from deseq2

# Check if the sample names in dge$samples are in metadata
print(all(rownames(dge$samples) %in% rownames(metadata))) # Expect TRUE

# If all sample names are present, you can add the 'condition' information
# from metadata to dge$samples, ensuring the order matches
dge$samples$condition <- metadata[rownames(dge$samples), "condition"]

# For unpaired design: 
design_edger <- model.matrix(~ condition, data = dge$samples)

# Design matrix column names:
print(colnames(design_edger)) # "coef" = entry on the right (e.g., "conditiondrug")

# Step 2. Filter and normalize data

# Filter data using CPM
min_samples <- 3 # gene must have CPM >= 1 in at least 30% of the samples to be kept
cpm_dge <- cpm(dge) # Compute CPM (Counts Per Million)
keep_cpm <- rowSums(cpm_dge >= 1) >= min_samples
keep_fbe <- filterByExpr(dge, design = design_edger)

# Option 1: Uses “|” for OR (less strict filtering) (Option 2 uses “&” for AND)
# keep_combined_or <- keep_cpm | keep_fbe
# dge <- dge[keep_combined_or, keep.lib.sizes = FALSE]
# dge$samples$lib.size <- colSums(dge$counts) # recalculate library sizes
# dge <- calcNormFactors(dge) # Normalization on dge

# Option 3: Combines both by dropping low information rows:
keep_fbe <- filterByExpr(dge,
                         design   = design_edger,
                         min.count = 10,      # CPM of 1 in a 10 M-read library
                         min.prop  = 0.3)     # 30 % of libraries
dge <- dge[keep_fbe, keep.lib.sizes = FALSE]
dge$samples$lib.size <- colSums(dge$counts) # recalculate library sizes
dge <- calcNormFactors(dge) # Normalization on dge

# Step 3: Run the edgeR algorithm

# Estimate dispersion (related to the variance of the negative binomial distribution) 
dge <- estimateDisp(dge, design = design_edger) # unpaired design

# Use either Model 1 or 2 for modeling the data: 
# Use Model 1 unless Model 1 gives 0 or few results

# Model 1
# fit <- glmQLFit(dge, design_edger)
# qlf <- glmQLFTest(fit, coef = "conditiondrug")
# top_tags_all <- topTags(qlf, n = Inf)$table

# Number of genes with FDR < 0.05 and log2(fold change) >= 1
# cat("Number of genes with FDR < 0.05 and log2(fold change) >= 1: \n")
# sum(top_tags_all$FDR < 0.05 & abs(top_tags_all$logFC) >= 1)
# = 0

# Model 2
fit_lrt <- glmFit(dge, design_edger)
lrt <- glmLRT(fit_lrt, coef = "conditiondrug")
top_tags_lrt <- topTags(lrt, n = Inf)$table

# Number of genes with FDR < 0.05 and log2(fold change) >= 1
cat("Number of genes with FDR < 0.05 and log2(fold change) >= 1: \n")
sum(top_tags_lrt$FDR < 0.05 & abs(top_tags_lrt$logFC) >= 1)
# > 0

# Step 4: Create a list of the Significant Genes (via edgeR) using Model 2
top_tags_2 <- as.data.frame(top_tags_lrt) # top_tags_lrt is Model 2
sig_genes_edger <- subset(top_tags_2, FDR < 0.05 & abs(logFC) >= 1)

# Save the results of edgeR and the significant genes
saveRDS(sig_genes_edger, "mouse_data/sig_genes_edger.rds") # significant genes via edgeR (list of genes)

top_tags_lrt <- top_tags_lrt[!is.na(top_tags_lrt$FDR), ]  
saveRDS(top_tags_lrt, "mouse_data/res_edger.rds") # results of edgeR

# Step 5. Create a robust significant gene list using DESeq2 and edgeR gene lists
library(tibble)   # rownames_to_column()

deseq2_df <- as.data.frame(sig_genes_deseq2) %>% rownames_to_column("gene")
edger_df  <- as.data.frame(sig_genes_edger)  %>% rownames_to_column("gene")
common_sig_genes <- intersect(deseq2_df$gene, edger_df$gene)

print("number of significant genes common to both DESeq2 and edgeR: ")
length(common_sig_genes) 

saveRDS(common_sig_genes, "mouse_data/common_sig_genes.rds")

##################
# Visualizations # Exploratory Data Analysis
##################

# MA plot
# M = log₂ fold change (on the y-axis)
# A = mean expression (on the x-axis, in log scale)

######################
# MA plot for DESeq2 # deseq2_results_filtered
######################

res <- results(deseq2_results_filtered)  
sum(!is.na(res$padj) & res$padj < 0.05)  # 248 significant
table(is.na(res$padj)) # TRUE means low counts

# Plot with custom point size and highlight order
png("mouse_data/DESeq2_MAplot.png", width = 1200, height = 1000, res = 150)
DESeq2::plotMA(res, alpha = 0.05, ylim = c(-5, 5),
               main = "MA plot (DESeq2)", colSig = "red")
legend("topright",
       legend = c("Significant (FDR < 0.05)",
                  "Not significant",
                  "Excluded (low counts)"),
       col = c("red", "grey", "blue"),
       pch = 16, cex = 0.7)
dev.off()

sum(res$log2FoldChange > 5, na.rm = TRUE)   # too high - represented by a triangle
sum(res$log2FoldChange < -5, na.rm = TRUE) # too low - represented by a triangle

#####################
# MA plot for edgeR # qlf/lrt or top_tags_lrt
#####################

library(limma)
library(ggplot2)

# ranked genes
# res_edger <- topTags(qlf, n = Inf)$table # Model 1
top_tags_lrt <- topTags(lrt, n = Inf)$table # Model 2

# calculate logCPM from the original dge object
logCPM_dge <- cpm(dge, log = TRUE)
# res_edger$AvgLogCPM <- rowMeans(logCPM_dge[rownames(res_edger), ]) # Model 1
top_tags_lrt$AvgLogCPM <- rowMeans(logCPM_dge[rownames(top_tags_lrt), ]) # Model 2

# Create the MA plot for edgeR - Model 1
#ma_plot_edger <- ggplot(res_edger, aes(x = AvgLogCPM, y = logFC)) +
#  geom_point(aes(color = ifelse(abs(logFC) >= 1 & FDR < 0.05, "significant", "not significant")), size = 2) +
#  scale_color_manual(values = c("significant" = "red", "not significant" = "grey")) +
#  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = "black") +
#  labs(title = "MA Plot", x = "Average LogCPM", y = "Log2 Fold Change", color = "Significance") +
#  theme_bw() +
#  geom_vline(xintercept = mean(res_edger$AvgLogCPM), linetype = "dotted", color = "blue") 
#
#ggsave("mouse_data/ma_plot_edger.png", plot = ma_plot_edger, width = 8, height = 6, dpi = 300)

# Create the MA plot for edgeR - Model 2
ma_plot_edger_2 <- ggplot(top_tags_lrt, aes(x = AvgLogCPM, y = logFC)) +
  geom_point(aes(color = ifelse(abs(logFC) >= 1 & FDR < 0.05, "significant", "not significant")), size = 0.5) +
  scale_color_manual(values = c("significant" = "red", "not significant" = "grey")) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = "black") +
  labs(title = "MA Plot (edgeR)", x = "Average LogCPM", y = "Log2 Fold Change", color = "Significance") +
  theme_bw() +
  geom_vline(xintercept = mean(top_tags_lrt$AvgLogCPM), linetype = "dotted", color = "blue") 

ggsave("mouse_data/edgeR_MAplot.png", plot = ma_plot_edger_2, width = 8, height = 6, dpi = 300)

#################
### PCA Plots ###
#################

#######################
# PCA Plot for DESeq2 # deseq2_results_filtered
#######################

# 1. Variance stabilizing transformation (VST) or regularized log transformation (rlog)
#    This helps to account for the mean-variance relationship in RNA-Seq data

# deseq2_results_table <-- deseq2_results_filtered with no NAs

vsd <- vst(deseq2_results_filtered, blind = FALSE) # Using VST
# Alternatively, you could use:
# rld <- rlog(dds, blind = FALSE)

# 2. Perform PCA
pcaData_deseq2 <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)

# 3. Create a data frame for plotting
percentVar_deseq2 <- pcaData_deseq2$sdev^2 / sum(pcaData_deseq2$sdev^2)
pcaDf_deseq2 <- data.frame(PC1 = pcaData_deseq2$x[, 1],
                           PC2 = pcaData_deseq2$x[, 2],
                           group = colData(deseq2_results_filtered)$condition) # Assuming your condition column is named 'condition'

# 4. Generate the PCA plot using ggplot2
pcaPlot_deseq2 <- ggplot(pcaDf_deseq2, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar_deseq2[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_deseq2[2] * 100), "% variance")) +
  ggtitle("PCA Plot (DESeq2 - VST)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("mouse_data/DESeq2_PCA_2D_plot.png", plot = pcaPlot_deseq2, width = 8, height = 6, dpi = 300)

##############################
# 3-dim PCA plot for DESeq2 ## pcaData_deseq2
##############################

library(plotly)
library(htmlwidgets)

# 1. Create the 3D scatter plot using plotly
# Calculate the percentage of variance explained by PC3
percentVar_pc3_deseq2 <- pcaData_deseq2$sdev[3]^2 / sum(pcaData_deseq2$sdev^2) * 100

# Create the 3D scatter plot using plotly
pcaPlot_3d_deseq2 <- plot_ly(
  data = pcaDf_deseq2,
  x = ~PC1,
  y = ~PC2,
  z = ~pcaData_deseq2$x[, 3],
  color = ~group,
  size = 4,
  marker = list(symbol = "circle"),
  text = ~paste(
    "Sample:",
    rownames(pcaDf_deseq2),
    "<br>Group:",
    group,
    "<br>PC1:",
    round(PC1, 2),
    "<br>PC2:",
    round(PC2, 2),
    "<br>PC3:",
    round(pcaData_deseq2$x[, 3], 2)
  ),
  hoverinfo = "text",
  type = "scatter3d",  # Explicitly set the plot type
  mode = "markers"     # Explicitly set the mode
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", round(percentVar_deseq2[1] * 100), "% variance")),
      yaxis = list(title = paste0("PC2: ", round(percentVar_deseq2[2] * 100), "% variance")),
      zaxis = list(title = paste0("PC3: ", round(percentVar_pc3_deseq2, 2), "% variance"))
    ),
    title = "3D PCA Plot (DESeq2 - VST)"
  )

saveWidget(pcaPlot_3d_deseq2, "mouse_data/deseq2_PCA_3D.html")

######################
# PCA plot for edgeR # dge
######################
# use dge after normalization and dispersion calculation
v <- voom(dge, design = design_edger, plot = TRUE)

# 2. Perform PCA on the normalized and transformed data
pcaData_edger <- prcomp(t(v$E), center = TRUE, scale. = TRUE)

# 3. Create a data frame for plotting
percentVar_edger <- pcaData_edger$sdev^2 / sum(pcaData_edger$sdev^2)
pcaDf_edger <- data.frame(PC1 = pcaData_edger$x[, 1],
                          PC2 = pcaData_edger$x[, 2],
                          group = dge$samples$condition) # Assuming your condition info is in dgeList$samples$condition

# 4. Generate the PCA plot using ggplot2
pcaPlot_edger <- ggplot(pcaDf_edger, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar_edger[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_edger[2] * 100), "% variance")) +
  ggtitle("PCA Plot (edgeR - voom)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("mouse_data/edgeR_PCA_2D_plot.png", plot = pcaPlot_edger, width = 8, height = 6, dpi = 300)

#########################################
# 3-dimensional plot for PCA with edgeR # pcaData_edger
#########################################

# Assuming the 2-dim plot was completed
# Calculate percentage variance explained by PC3
percentVar_pc3_edger <- pcaData_edger$sdev[3]^2 / sum(pcaData_edger$sdev^2) * 100

# Create data frame for 3D PCA plot
pcaDf_3d_edger <- data.frame(
  PC1 = pcaData_edger$x[, 1],
  PC2 = pcaData_edger$x[, 2],
  PC3 = pcaData_edger$x[, 3],
  group = dge$samples$condition,
  sample = rownames(dge$samples) # Add sample names for hover text
)

# Create 3D plot with plotly
pcaPlot_3d_edger <- plot_ly(
  data = pcaDf_3d_edger,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~group,
  size = 4,
  marker = list(symbol = "circle"),
  text = ~paste("Sample:", sample, "<br>Group:", group,
                "<br>PC1:", round(PC1, 2), "<br>PC2:", round(PC2, 2), "<br>PC3:", round(PC3, 2)),
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", round(percentVar_edger[1] * 100), "% variance")),
      yaxis = list(title = paste0("PC2: ", round(percentVar_edger[2] * 100), "% variance")),
      zaxis = list(title = paste0("PC3: ", round(percentVar_pc3_edger, 2), "% variance"))
    ),
    title = "3D PCA Plot (edgeR - voom)"
  )

# Save the 3D plot
saveWidget(pcaPlot_3d_edger, "mouse_data/edgeR_PCA_3D.html")

###########################################
# Heatmap for the common genes via DESeq2 # unpaired_design_filtered
###########################################

#  Use DESeq2's variance-stabilized transformed data (vsd)
#  Don't use raw counts.
vsd <- vst(unpaired_design_filtered, blind = FALSE) # Assuming 'unpaired_design_filtered' is the DESeqDataSet object
mat_deseq2 <- assay(vsd)[rownames(assay(vsd)) %in% common_sig_genes, ]
mat_scaled_deseq2 <- t(scale(t(mat_deseq2))) # center and scale the data

pheatmap(mat_scaled_deseq2,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE, # Or FALSE if you have too many genes
         show_colnames = TRUE,
         annotation_col = as.data.frame(colData(vsd)), # Sample annotations
         main = "Heatmap of Common Significant Genes (DESeq2 - VST)",
         filename = "mouse_data/common_sig_genes_deseq2_heatmap_ensembl.png") # Saves to a file
# This produces a heatmap with ensembl names for the genes

# The following code produces a heatmap with geneSymbols instead:

# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("AnnotationDbi")

library(org.Mm.eg.db)
library(AnnotationDbi)

current_heatmap_genes_deseq2 <- rownames(mat_scaled_deseq2)

# Remove version numbers from Ensembl IDs
current_heatmap_genes_unversioned_deseq2 <- sub("\\.\\d+$", "", current_heatmap_genes_deseq2)

gene_symbols_mapped_deseq2 <- mapIds(org.Mm.eg.db,
                                     keys = current_heatmap_genes_unversioned_deseq2,
                                     keytype = "ENSEMBL",
                                     column = "SYMBOL",
                                     multiVals = "first") # Keep this for handling multiple mappings

final_gene_names_deseq2 <- ifelse(is.na(gene_symbols_mapped_deseq2),
                                  current_heatmap_genes_deseq2, # Fallback to original versioned Ensembl ID
                                  gene_symbols_mapped_deseq2)

# Ensure the names of this vector match the original versioned Ensembl IDs,
# which allows correct lookup when assigning to rownames.
names(final_gene_names_deseq2) <- current_heatmap_genes_deseq2

# Assign these as row names to the heatmap matrix
rownames(mat_scaled_deseq2) <- final_gene_names_deseq2

pheatmap(mat_scaled_deseq2,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = as.data.frame(colData(vsd)),
         main = "Heatmap of Common Significant Genes - DESeq2",
         fontsize_main = 2,
         margins = c(top = 10, right = 5, bottom = 5, left = 5),
         filename = "mouse_data/common_sig_genes_deseq2_heatmap_symbols.png")

###########################################
# Publication-quality heatmap via DESeq2 ## mat_scaled_deseq2
###########################################
## ---------------------------------------------------------------
## 0.  Preprocessing
## ---------------------------------------------------------------
library(AnnotationDbi)
library(org.Mm.eg.db)      # ↔ your organism

mat_scaled <- mat_scaled_deseq2              # rows = genes, cols = samples

## column annotation
anno_col <- as.data.frame(colData(vsd)) |>
  dplyr::select(sizeFactor, condition, sample)   # keep only needed

## set colors for annotation bars
library(RColorBrewer)
anno_colours <- list(
  sizeFactor = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  condition  = c(base = "#F0E442", drug = "#009E73"),
  sample     = setNames(brewer.pal(n = nrow(anno_col), "Paired"),
                        rownames(anno_col))
)

## ---------------------------------------------------------------
## 1.  HEAT-MAP
## ---------------------------------------------------------------
library(pheatmap)

## pick a balanced red-blue palette
pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

p <- pheatmap(
  mat_scaled,
  color                 = pal,
  cluster_rows          = TRUE,
  cluster_cols          = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method     = "complete",
  show_rownames         = TRUE,
  show_colnames         = TRUE,
  annotation_col        = anno_col,
  annotation_colors     = anno_colours,
  annotation_names_col  = FALSE,     # hide long titles
  fontsize              = 8,         # general text
  fontsize_row          = 6,         # smaller gene labels
  fontsize_col          = 8,
  border_color          = NA,        # no cell borders
  main                  = "Common significant genes (DESeq2, VST)",
  legend                = TRUE
)

## ---------------------------------------------------------------
## 2.  SAVE  (vector PDF + high-res PNG)
## ---------------------------------------------------------------
## vector-quality PDF 
pdf("mouse_data/pub_common_sig_genes_deseq2_heatmap_symbols.pdf", width = 7, height = 4.5)
print(p)
dev.off()

## 300-dpi PNG for slides / email
png("mouse_data/pub_common_sig_genes_deseq2_heatmap_symbols.png",
    width = 2100, height = 1350, res = 300)
print(p)
dev.off()

##########################################
# Publication-quality heatmap via edgeR ## dge
##########################################

#  Use voom-transformed and normalized logCPM values
v <- voom(dge, design_edger, plot = FALSE) # edgeR objects
mat_edger <- v$E[rownames(v$E) %in% common_sig_genes, ]
mat_scaled_edger <- t(scale(t(mat_edger))) #  Center and scale

anno_col_edger_modified <- dge$samples %>% dplyr::select(lib.size, condition)

# This produces a heatmap with ensembl names for the genes
pheatmap(mat_scaled_edger,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = anno_col_edger_modified, # Sample annotations
         main = "Heatmap of Common Significant Genes (edgeR)",
         filename = "mouse_data/common_sig_genes_edger_heatmap_ensembl.png")

# Heatmap for edgeR (using gene symbols (or Ensembl IDs if gene_symbol not found))

current_heatmap_genes_edger <- rownames(mat_scaled_edger)

# Remove version numbers from Ensembl IDs:
current_heatmap_genes_unversioned_edger <- sub("\\.\\d+$", "", current_heatmap_genes_edger)

gene_symbols_mapped_edger <- mapIds(org.Mm.eg.db,
                             keys = current_heatmap_genes_unversioned_edger,
                             keytype = "ENSEMBL",
                             column = "SYMBOL")

final_gene_names_edger <- ifelse(is.na(gene_symbols_mapped_edger),
                                  # Fallback to original versioned Ensembl ID
                                  current_heatmap_genes_edger,                                  			
                                  gene_symbols_mapped_edger)

# Ensure the names of this vector match the original versioned Ensembl IDs, 
# which allows correct lookup when assigning to rownames:
  
names(final_gene_names_edger) <- current_heatmap_genes_edger

rownames(mat_scaled_edger) <- final_gene_names_edger

anno_col_edger_modified <- dge$samples %>% dplyr::select(lib.size, condition)


pheatmap(mat_scaled_edger,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = anno_col_edger_modified, # Sample annotations
         main = "Heatmap of Common Significant Genes (edgeR)",
         filename = "mouse_data/common_sig_genes_edger_heatmap_symbols.png")

###########################
# Volcano Plot for DESeq2 # deseq2_results_table_filtered
###########################

deseq2_volcano_df <- as.data.frame(deseq2_results_table_filtered)
deseq2_volcano_df <- deseq2_volcano_df[complete.cases(deseq2_volcano_df), ]
deseq2_volcano_df$significant <- ifelse(deseq2_volcano_df$padj < 0.05, "Significant", "Not Significant")

logFC_cutoff <- 1 
deseq2_volcano_df$sig_logFC <- ifelse(deseq2_volcano_df$padj < 0.05 & abs(deseq2_volcano_df$log2FoldChange) >= logFC_cutoff, "Sig & LogFC", deseq2_volcano_df$significant)


# Create the Volcano Plot with ggplot2
volcano_plot_deseq2 <- ggplot(deseq2_volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig_logFC), size = 1) + # Color by combined significance
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "grey", "Sig & LogFC" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Significance cutoff line
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") + # logFC cutoff lines
  labs(
    title = "Volcano Plot (DESeq2)",
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted p-value)",
    color = "Significance"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("mouse_data/DESeq2_Volcano_Plot.png", plot = volcano_plot_deseq2, width = 8, height = 6, dpi = 300)

##########################
# Volcano Plot for edgeR # top_tags_lrt
##########################

edgeR_volcano_df <- as.data.frame(top_tags_lrt) 
edgeR_volcano_df <- edgeR_volcano_df[complete.cases(edgeR_volcano_df), ]
edgeR_volcano_df$significant <- ifelse(edgeR_volcano_df$FDR < 0.05, "Significant", "Not Significant")

logFC_cutoff <- 1 # Or your chosen cutoff
edgeR_volcano_df$sig_logFC <- ifelse(edgeR_volcano_df$FDR < 0.05 & abs(edgeR_volcano_df$logFC) >= logFC_cutoff, "Sig & LogFC", edgeR_volcano_df$significant)


# Create the Volcano Plot with ggplot2
volcano_plot_edger <- ggplot(edgeR_volcano_df, aes(x = logFC, y = -log10(FDR))) + # Use logFC and FDR
  geom_point(aes(color = sig_logFC), size = 1) +
  scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "grey", "Sig & LogFC" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot (edgeR)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("mouse_data/edgeR_Volcano_Plot.png", plot = volcano_plot_edger, width = 8, height = 6, dpi = 300)

##########################################
# publication quality of volcano plot: ### deseq2_volcano_df
##########################################

# Load necessary libraries
library(ggplot2)
library(ggrepel) # For non-overlapping text labels
library(org.Mm.eg.db)
library(dplyr)

# Getting gene names
deseq2_volcano_df$ensembl_id <- rownames(deseq2_volcano_df)
deseq2_volcano_df$gene_symbol <- mapIds(
  x = org.Mm.eg.db,
  keys = deseq2_volcano_df$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first" # Or "filter", "list", etc., depending on how you want to handle one-to-many mappings
)
deseq2_volcano_df <- deseq2_volcano_df %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol), ensembl_id, gene_symbol)
  )

# Assuming deseq2_volcano_df is already prepared as per your existing code.
# Ensure 'deseq2_volcano_df' contains:
# - log2FoldChange
# - padj (adjusted p-value)
# - gene_name (or similar identifier for labeling, if desired)
# - sig_logFC (the column you use for coloring, which defines significance categories)

# --- Define cutoffs (from your existing code or chosen values) ---
p_cutoff <- 0.05
logFC_cutoff <- 1 # Example: 2-fold change

# --- Re-categorize 'sig_logFC' for more descriptive colors ---
deseq2_volcano_df <- deseq2_volcano_df %>%
  mutate(
    sig_logFC_category = case_when(
      padj < p_cutoff & log2FoldChange > logFC_cutoff ~ "Upregulated (Significant)",
      padj < p_cutoff & log2FoldChange < -logFC_cutoff ~ "Downregulated (Significant)",
      TRUE ~ "Not Significant"
    )
  )

# --- Count significant genes for annotation ---
num_up <- sum(deseq2_volcano_df$sig_logFC_category == "Upregulated (Significant)", na.rm = TRUE)
num_down <- sum(deseq2_volcano_df$sig_logFC_category == "Downregulated (Significant)", na.rm = TRUE)
num_not_sig <- sum(deseq2_volcano_df$sig_logFC_category == "Not Significant", na.rm = TRUE)

# --- Create the production-quality volcano plot ---
volcano_plot_production <- ggplot(deseq2_volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
  # Add points, colored by the new category
  geom_point(aes(color = sig_logFC_category), size = 1.5, alpha = 0.8) +
  
  # Define manual colors for the categories
  scale_color_manual(
    values = c(
      "Upregulated (Significant)" = "red",
      "Downregulated (Significant)" = "blue",
      "Not Significant" = "grey70" # Use a lighter grey for less emphasis
    ),
    name = "Gene Expression Status" # More descriptive legend title
  ) +
  
  # Significance cutoff line (p-value = 0.05)
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
  # Label the p-value cutoff line
  annotate("text", x = max(deseq2_volcano_df$log2FoldChange, na.rm = TRUE) * 0.95,
           y = -log10(p_cutoff) + 0.1,
           label = paste0("p-value < ", p_cutoff),
           hjust = 1, vjust = 0, size = 3) +
  
  # logFC cutoff lines (e.g., Fold Change > 2 or < 0.5)
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
  # Label the logFC cutoff line for up-regulated
  annotate("text", x = logFC_cutoff + 0.1,
           y = max(-log10(deseq2_volcano_df$padj), na.rm = TRUE) * 0.95,
           label = paste0("log2FC > ", logFC_cutoff),
           hjust = 0, vjust = 1, size = 3) +
  # Label the logFC cutoff line for down-regulated
  annotate("text", x = -logFC_cutoff - 0.1,
           y = max(-log10(deseq2_volcano_df$padj), na.rm = TRUE) * 0.95,
           label = paste0("log2FC < ", -logFC_cutoff),
           hjust = 1, vjust = 1, size = 3) +
  
  # Labels for title, axes, and legend
  labs(
    title = "Volcano Plot: Differential Gene Expression",
    subtitle = paste0(
      "Upregulated: ", num_up,
      " | Downregulated: ", num_down,
      " | Not Significant: ", num_not_sig
    ),
    x = expression(log[2]~"Fold Change"), # Use expression for subscript
    y = expression(-log[10]~"("~Adjusted~italic("p")~"-value)") # Use expression for subscript and italic p
  ) +
  
  # Use a classic theme and customize elements for a clean look
  theme_classic() + # A clean, common theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and bold title
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey50"), # Centered subtitle for counts
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2), # Lighter grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = 0.5), # Ensure axis lines are visible
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Adjust plot margins
  ) +
  
  # Add labels for top N most significant genes
  # You would need a 'gene_name' column in your deseq2_volcano_df
  # Adjust 'top_n' as needed
  geom_text_repel(
    data = deseq2_volcano_df %>%
      filter(sig_logFC_category != "Not Significant") %>%
      arrange(padj) %>%
      head(10), # Label top 10 significant genes
    aes(label = gene_symbol), # Assuming 'gene_symbol' is your gene identifier column
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    min.segment.length = 0.2, # Draw segments if needed
    segment.color = 'grey50'
  )

# Print the plot
print(volcano_plot_production)
ggsave("mouse_data/DESeq2_Volcano_Plot_pub.png", plot = volcano_plot_production, width = 8, height = 6, dpi = 300)

################################################
######## PUBLICATION-QUALITY MA PLOT ########### deseq2_volcano_df
################################################

# Load necessary libraries
library(ggplot2)
library(ggrepel) # For non-overlapping text labels
library(dplyr)   # For %>% pipe operator and mutate/case_when

# --- Define cutoffs (same as for volcano plot for consistency) ---
p_cutoff <- 0.05
logFC_cutoff <- 1 # Example: 2-fold change, adjust as needed

# --- Categorize genes for coloring based on significance and fold change ---
# This ensures that your 'deseq2_volcano_df' has the 'sig_status_ma' column
# It requires 'log2FoldChange', 'padj', and 'baseMean' columns to be present.
deseq2_volcano_df <- deseq2_volcano_df %>%
  mutate(
    sig_status_ma = case_when(
      padj < p_cutoff & log2FoldChange > logFC_cutoff ~ "Upregulated (Significant)",
      padj < p_cutoff & log2FoldChange < -logFC_cutoff ~ "Downregulated (Significant)",
      TRUE ~ "Not Significant"
    )
  )

# --- Count significant genes for plot subtitle annotation ---
num_up_ma <- sum(deseq2_volcano_df$sig_status_ma == "Upregulated (Significant)", na.rm = TRUE)
num_down_ma <- sum(deseq2_volcano_df$sig_status_ma == "Downregulated (Significant)", na.rm = TRUE)
num_not_sig_ma <- sum(deseq2_volcano_df$sig_status_ma == "Not Significant", na.rm = TRUE)

# --- Create the production-quality MA plot ---
ma_plot_production <- ggplot(deseq2_volcano_df, aes(x = log10(baseMean), y = log2FoldChange)) +
  # Add points, colored by their significance status
  geom_point(aes(color = sig_status_ma), size = 1.5, alpha = 0.8) +
  
  # Define manual colors for each category
  scale_color_manual(
    values = c(
      "Upregulated (Significant)" = "red",
      "Downregulated (Significant)" = "blue",
      "Not Significant" = "grey70" # Use a lighter grey for less emphasis
    ),
    name = "Gene Expression Status" # More descriptive legend title
  ) +
  
  # Horizontal line at 0 fold change (no differential expression)
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6) +
  
  # Fold Change cutoff lines (e.g., at log2FC = 1 and -1)
  geom_hline(yintercept = logFC_cutoff, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
  geom_hline(yintercept = -logFC_cutoff, linetype = "dashed", color = "darkgreen", linewidth = 0.6) +
  
  # Annotate the fold change cutoff lines directly on the plot
  annotate("text", x = min(log10(deseq2_volcano_df$baseMean), na.rm = TRUE),
           y = logFC_cutoff + 0.1,
           label = paste0("log2FC = ", logFC_cutoff),
           hjust = 0, vjust = 0, size = 3, color = "darkgreen") +
  annotate("text", x = min(log10(deseq2_volcano_df$baseMean), na.rm = TRUE),
           y = -logFC_cutoff - 0.1,
           label = paste0("log2FC = ", -logFC_cutoff),
           hjust = 0, vjust = 1, size = 3, color = "darkgreen") +
  
  # Labels for title, axes, and legend
  labs(
    title = "       MA Plot: Differential Gene Expression vs. Mean Expression",
    subtitle = paste0(
      "Upregulated: ", num_up_ma,
      " | Downregulated: ", num_down_ma,
      " | Not Significant: ", num_not_sig_ma
    ),
    x = expression(log[10]~"Mean Normalized Counts"), # Formats x-axis label with subscript
    y = expression(log[2]~"Fold Change")             # Formats y-axis label with subscript
  ) +
  
  # Apply a classic theme and customize elements for a clean, publication-ready look
  theme_classic() + # Provides a clean background and axis lines
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Centered, bold title
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey50"), # Centered subtitle for counts
    axis.title = element_text(size = 12, face = "bold"), # Bold axis titles
    axis.text = element_text(size = 10), # Adjust axis text size
    legend.title = element_text(size = 11, face = "bold"), # Bold legend title
    legend.text = element_text(size = 10), # Adjust legend text size
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2), # Lighter major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines for cleaner look
    axis.line = element_line(color = "black", linewidth = 0.5), # Ensure axis lines are visible
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Adjust plot margins
  ) +
  
  # Label the top 10 most significant genes (by adjusted p-value)
  # Uses ggrepel to prevent labels from overlapping
  geom_text_repel(
    data = deseq2_volcano_df %>%
      filter(sig_status_ma != "Not Significant") %>% # Only label significant genes
      arrange(padj) %>% # Sort by adjusted p-value to get the most significant first
      head(10), # Select the top 10 genes
    aes(label = gene_symbol), # Use the 'gene_symbol' column for labels
    size = 3.5, # Adjust label text size
    box.padding = unit(0.4, "lines"), # Padding around text
    point.padding = unit(0.3, "lines"), # Minimum distance from text to point
    min.segment.length = 0.2, # Draw segments to point if needed
    segment.color = 'grey50', # Color of the segment lines
    max.overlaps = Inf # Allows all labels to be shown, even if they would overlap slightly. Use with caution for very dense plots.
  )

# Print the final MA plot
print(ma_plot_production)
ggsave("mouse_data/DESeq2_MA_Plot_pub.png", plot = ma_plot_production, width = 8, height = 6, dpi = 300)
