#!/usr/bin/env Rscript

setwd("~/Desktop/")

library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Load and preprocess data
rmats_data <- read.table("mouse_data/rmats_results/SE.MATS.JCEC.txt", header = TRUE, sep = "\t") # all one line, no “/“ in front of the folder “test_data”

# Define filtering thresholds
fdr_threshold <- 0.05
min_inclusion_level <- 0.1
max_inclusion_level <- 0.9

# Remove missing values
rmats_data_no_na_string <- rmats_data %>% filter(!grepl("NA", IncLevel1) & !grepl("NA", IncLevel2))

# Compute the average inclusion
rmats_data_cleaned <- rmats_data_no_na_string %>% rowwise() %>%
  mutate(AvgInc = mean(as.numeric(unlist(strsplit(paste(IncLevel1, IncLevel2, sep = ","), ","))), na.rm = TRUE)) %>%
  ungroup()

# Filter by FDR
filtered_fdr <- rmats_data_cleaned %>% filter(FDR <= fdr_threshold)

# Filter by average inclusion level: 
filtered_inclusion <- filtered_fdr %>%
  filter(AvgInc >= min_inclusion_level & AvgInc <= max_inclusion_level)

# For troubleshooting
print("Number of events after inclusion level filtering and FDR filtering: ")
nrow(filtered_inclusion) 

df  <- as.data.frame(filtered_inclusion)

# For troubleshooting
head(df)

# Extract PSI values for Group 1 and Group 2 (Control and Treatment, respectively)
# IncLevel1: PSI values for Group 1 (e.g., Control)
# IncLevel2: PSI values for Group 2 (e.g., Treatment)
psi_data <- df %>%
  select(IncLevel1, IncLevel2) %>%
  rowwise() %>%
  mutate(
    Baseline_PSI = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE), # Group1_PSI
    Drug_PSI = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE) # Group2_PSI
  ) %>%
  ungroup()

# Reshape the data for plotting (this means to shift from wide format to long format):
psi_long <- psi_data %>%
  pivot_longer(cols = c(Baseline_PSI, Drug_PSI), names_to = "Condition", values_to = "PSI") %>%
  filter(!is.na(PSI))

# 2. Perform Summary Statistics and Data Distribution
# Generate summary statistics of the inclusion levels (PSI values) across samples 

summary_stats <- psi_long %>%
  group_by(Condition) %>% # need to give an example for group_by
  summarize(
    Mean = mean(PSI, na.rm = TRUE),
    Median = median(PSI, na.rm = TRUE),
    SD = sd(PSI, na.rm = TRUE),
    Min = min(PSI, na.rm = TRUE),
    Max = max(PSI, na.rm = TRUE)
  )

print(summary_stats)

# Condition      Mean   Median  SD    Min       Max
# <chr>          <dbl>  <dbl>   <dbl> <dbl>     <dbl>
# 1 Baseline_PSI 0.661  0.804   0.315 0         1
# 2 Drug_PSI     0.605  0.724   0.291 0.0184    1

######################################################
# 1. Density Plots in Alternative Splicing Analysis ##
# Density Plot for filtered values                  ##
######################################################

#  1) Check the directory:
getwd()

# Make and save density plot:
se_density_plot <- ggplot(psi_long, aes(x = PSI, fill = Condition)) +
    geom_density(alpha = 0.6) +
    theme_minimal() +
    labs(title = "Density Plot of PSI Values", x = "PSI", y = "Density") +
    scale_fill_manual(values = c("#69b3a2", "#404080"))

ggsave("mouse_data/rmats_plots/se_density_plot.png", plot = se_density_plot, width = 6, height = 4, dpi = 300)
  
#################################################
# 2. Boxplots in Alternative Splicing Analysis ##
# Boxplot for PSI Values:                      ##
#################################################
se_boxplot <- ggplot(psi_long, aes(x = Condition, y = PSI, fill = Condition)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Boxplot of PSI Values", x = "Condition", y = "PSI") +
  scale_fill_manual(values = c("#69b3a2","#404080"))

ggsave("mouse_data/rmats_plots/se_boxplot.png", plot = se_boxplot, width = 6, height = 4, dpi = 300) 

##########################################
# 3. Principal Component Analysis (PCA) ##
# PCA plot on PSI Values                ##
##########################################
library(stringr)

# Expand the 4 (number of samples) PSI values per group/condition into individual 
# sample columns:
psi_matrix <- rmats_data_no_na_string %>%
  rowwise() %>%
  mutate(
    Baseline = list(as.numeric(str_split(IncLevel1, ",")[[1]])),
    Drug  = list(as.numeric(str_split(IncLevel2, ",")[[1]]))
  ) %>%
  ungroup()

# Number of samples per group (for example, 5)
n_samples <- 5

psi_df <- psi_matrix %>%
  rowwise() %>%
  mutate(
    # Convert list elements to separate columns
    !!!setNames(
      lapply(1:n_samples, function(i) substitute(Baseline[[i]], list(i = i))),
      paste0("Baseline", 1:n_samples)
    ),
    !!!setNames(
      lapply(1:n_samples, function(i) substitute(Drug[[i]], list(i = i))),
      paste0("Drug", 1:n_samples)
    )
  ) %>%
  ungroup() %>%
  select(starts_with("Baseline"), starts_with("Drug"))

psi_df_clean <- psi_df %>% select(-Baseline, -Drug)

# Transpose to have rows = samples, columns = splice events:
psi_transposed <- as.data.frame(t(psi_df_clean))
colnames(psi_transposed) <- paste0("Event", seq_len(ncol(psi_transposed)))
psi_transposed$Sample <- rownames(psi_transposed) # Add new column (sample type and name)

# Add condition labels to each sample (new column):
psi_transposed$Condition <- rep(c("Baseline", "Drug"), each = n_samples)

# Check the file and remove two columns:
# head(psi_transposed)
pca_input <- psi_transposed %>% select(-Sample, -Condition)

# Run PCA computation (filter out zero variance events):
pca_input_clean <- pca_input[, apply(pca_input, 2, function(col) var(col, na.rm = TRUE) > 0)]
pca_res <- prcomp(pca_input_clean, scale. = TRUE)

# Plot PCA:
pca_df <- as.data.frame(pca_res$x)
pca_df$Condition <- psi_transposed$Condition
pca_df$Sample <- psi_transposed$Sample

# pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  # geom_text(vjust = -1) +
  labs(title = "PCA of PSI Values by Sample")

ggsave("mouse_data/rmats_plots/se_pca_plot.png", plot = pca_plot, width = 6, height = 4, dpi = 300)

#PCA summary:
# print(summary(pca_res)) # number of PCA components = min(number of samples, number of splicing events)

####################
# 4. Volcano Plot ##
####################

library(ggplot2)
library(dplyr)
library(ggrepel)

rmats_data_volcano <- rmats_data %>%
  filter(!grepl("NA", IncLevel1) & !grepl("NA", IncLevel2)) %>%
  mutate(
    PSI_Group1 = sapply(strsplit(IncLevel1, ","), function(x) mean(as.numeric(x), na.rm = TRUE)),
    PSI_Group2 = sapply(strsplit(IncLevel2, ","), function(x) mean(as.numeric(x), na.rm = TRUE)),
    Delta_PSI = PSI_Group2 - PSI_Group1,
    log2_Delta_PSI = log2(PSI_Group2 + 1e-6) - log2(PSI_Group1 + 1e-6),
    # neg_log10_FDR = -log10(FDR)
    neg_log10_FDR = -log10(ifelse(FDR == 0, .Machine$double.xmin, FDR)) # Added ifelse
  )

# Identify top 10 significant events based on FDR
top10_events <- rmats_data_volcano %>% arrange(FDR) %>% head(10)

se_volcano_plot <- ggplot(rmats_data_volcano, aes(x = Delta_PSI, y = neg_log10_FDR)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "#E41A1C")) + # Use "FALSE" and "TRUE" for logical color mapping
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text_repel(data = top10_events, # Use geom_text_repel
                  aes(label = geneSymbol),
                  size = 3, max.overlaps = Inf, # Adjust max.overlaps as needed
                  box.padding = 0.5, point.padding = 0.5,
                  segment.color = 'grey50') +
  theme_minimal() +
  labs(title = "Volcano Plot of PSI Values", x = "Delta PSI", y ="-log10(FDR)")

ggsave("mouse_data/rmats_plots/se_volcano_plot.png", plot = se_volcano_plot, width = 6, height = 4, dpi = 300)

###############
# 5. MA Plot ##
###############
library(ggplot2)
library(dplyr)

# Data Preprocessing
# Calculate PSI values for Group 1 and Group 2
rmats_data_ma <- rmats_data %>% rowwise() %>%
  filter(!grepl("NA", IncLevel1) & !grepl("NA", IncLevel2)) %>%
  mutate( 									PSI_Group1_values = list(as.numeric(unlist(strsplit(IncLevel1, ",")))), PSI_Group2_values = list(as.numeric(unlist(strsplit(IncLevel2, ",")))), PSI_Group1 = mean(PSI_Group1_values, na.rm = TRUE), 			PSI_Group2 = mean(PSI_Group2_values, na.rm = TRUE), 			Mean_PSI = mean(c(PSI_Group1_values, PSI_Group2_values), na.rm = TRUE), log2_FC = log2((PSI_Group2 + 1e-6) / (PSI_Group1 + 1e-6)),	neg_log10_FDR = -log10(FDR),
                   Delta_PSI = PSI_Group2 - PSI_Group1) %>% 
  ungroup()

# Highlight significant events (e.g., FDR < 0.05):
rmats_data_sig <- rmats_data_ma %>% mutate(Significant = ifelse(FDR < 0.05, "Yes", "No"))

# Identify top 10 most significant events:
top10_events_ma <- rmats_data_sig %>% arrange(FDR) %>% head(10)

# MA Plot for all events:
se_ma_plot <- ggplot(rmats_data_sig, aes(x = Mean_PSI, y = log2_FC)) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "#E41A1C")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_text(data = top10_events_ma, aes(label = geneSymbol), size = 3, vjust = -1) +
  theme_minimal() +
  labs(
    title = "MA Plot of PSI Values (Alternative Splicing)",
    x = "Mean PSI",
    y = "Log2 Fold Change in PSI",
    color = "Significant"
  )

ggsave("mouse_data/rmats_plots/se_ma_plot.png", plot = se_ma_plot, width = 6, height = 4, dpi = 300)

# revised plot
# --- MA Plot (Revised) ---
rmats_data_ma_revised <- rmats_data %>%
  filter(!grepl("NA", IncLevel1) & !grepl("NA", IncLevel2)) %>%
  rowwise() %>%
  mutate(
    PSI_Group1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
    PSI_Group2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE),
    # Calculate Average PSI for X-axis (A-value)
    Average_PSI = (PSI_Group1 + PSI_Group2) / 2, # Or use your existing Mean_PSI if it's (G1+G2)/2
    # Calculate Delta PSI for Y-axis (M-value)
    Delta_PSI = PSI_Group2 - PSI_Group1,
    # Handle FDR for significance coloring
    FDR_for_plot = ifelse(FDR == 0, .Machine$double.xmin, FDR)
  ) %>%
  ungroup() %>%
  mutate(Significant = ifelse(FDR_for_plot < 0.05 & abs(Delta_PSI) > 0.1, "Yes", "No")) # Also include Delta PSI threshold for significance

# Identify top 10 most significant events (based on revised data)
top10_events_ma_revised <- rmats_data_ma_revised %>%
  filter(Significant == "Yes") %>% # Filter for significant events first
  arrange(FDR_for_plot) %>%
  head(10)

se_ma_plot_revised <- ggplot(rmats_data_ma_revised, aes(x = Average_PSI, y = Delta_PSI)) + # Changed y to Delta_PSI
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "#E41A1C")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Use ggrepel for text labels if possible
  geom_text_repel(data = top10_events_ma_revised, aes(label = geneSymbol), size = 3, vjust = -1,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50') +
  theme_bw() +
  labs(
    title = "MA Plot of PSI Values (Alternative Splicing)",
    x = "Average PSI",
    y = "Delta PSI (Drug - Baseline)", # Label Y-axis clearly
    color = "Significant"
  )

# print(se_ma_plot_revised)

ggsave("mouse_data/rmats_plots/se_ma_plot_revised.png", plot = se_ma_plot_revised, width = 6, height = 4, dpi = 300)

######################################
# 6. Heatmap of Significant Events: ##
######################################

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize) 

heatmap_significant_events <- rmats_data_ma %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)  # Select top 10 significant events

# Extract PSI matrix for heatmap
heatmap_psi_matrix <- heatmap_significant_events %>%
  select(ID, geneSymbol, exonStart_0base, exonEnd, PSI_Group1, PSI_Group2) %>%
  mutate(identifier = paste(geneSymbol, ID, sep = "_")) %>%
  pivot_longer(
    cols = starts_with("PSI_Group"), 
    names_to = "Group", 
    values_to = "PSI"
  ) %>%
  pivot_wider(
    id_cols = identifier,
    names_from = Group, 
    values_from = PSI) %>%
  rename(PSI_Baseline = PSI_Group1, PSI_Drug = PSI_Group2)

# Set rownames as event IDs
heatmap_psi_matrix <- as.data.frame(heatmap_psi_matrix) 

rownames(heatmap_psi_matrix) <- heatmap_psi_matrix$identifier

heatmap_psi_matrix <- heatmap_psi_matrix %>% select(-identifier) %>% as.matrix() 

# heatmap_psi_matrix <- heatmap_psi_matrix[, !(colnames(heatmap_psi_matrix) %in% c("geneSymbol", "exonStart_0base", "exonEnd", "identifier"))]

# head(heatmap_psi_matrix)

# Z-score normalization
psi_zscore <- t(scale(t(heatmap_psi_matrix)))

col_fun = colorRamp2(c(min(psi_zscore), 0, max(psi_zscore)), c("blue", "white", "red"))

# Define the file path
output_file_path <- "mouse_data/rmats_plots/se_heatmap.png"

# Save the heatmap using a graphics device
png(filename = output_file_path, width = 6, height = 4, units = "in", res = 300)

# Generate Heatmap
draw(Heatmap(
  matrix = psi_zscore,
  name = "Z-score of PSI", # This is the name given to the heatmap itself, and also the default legend title
  col = col_fun, # Use the color mapping function
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Heatmap of Top 10 Significant Splicing Events",
  row_names_gp = gpar(fontsize = 6), # Font size for row names
  column_names_gp = gpar(fontsize = 8), # Font size for column names
  border = TRUE, # Adds a border around the heatmap body
  # Controls for the heatmap legend
  heatmap_legend_param = list(
    title = "Normalized PSI", # This is the direct argument for the legend title
    direction = "vertical", # or "horizontal"
    at = c(floor(min(psi_zscore)), 0, ceiling(max(psi_zscore))), # Explicit breaks for the legend
    labels = c(floor(min(psi_zscore)), "0", ceiling(max(psi_zscore))) # Labels for the breaks
  )
))
dev.off() # Close the graphics device


##########################################
# 7. Pathway Analysis for rMATS Results ##
##########################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(reactome.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)

getwd()

# Load data
rmats_data <- read.table("mouse_data/rmats_results/SE.MATS.JCEC.txt", header = TRUE, sep = "\t")

# Filter for significant events based on FDR < 0.05:
significant_events <- rmats_data %>% filter(FDR < 0.05)

# Extract unique gene symbols:
gene_symbols <- unique(significant_events$geneSymbol)

# --- DIAGNOSTIC STEP 1: Check how many unique gene symbols you have ---
cat("Number of unique gene symbols:", length(gene_symbols), "\n")
if(length(gene_symbols) == 0) {
  stop("No significant gene symbols found. Check your FDR filter or rmats_data.")
}

# Convert gene symbols to Entrez IDs:
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Store the output of bitr in a data frame
gene_entrez_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# --- DIAGNOSTIC STEP 2: Check mapping success ---
cat("Number of gene symbols before mapping:", length(gene_symbols), "\n")
cat("Number of gene symbols mapped to ENTREZID:", nrow(gene_entrez_df), "\n")
cat("Number of unique ENTREZIDs mapped:", length(unique(gene_entrez_df$ENTREZID)), "\n")
cat("First few rows of gene_entrez_df:\n")
print(head(gene_entrez_df)) # Print the data frame itself to see mapped and unmapped SYMBOLs



# Ensure there are valid Entrez IDs
gene_entrez <- na.omit(gene_entrez$ENTREZID)

# --- DIAGNOSTIC STEP 3: Check the final set of ENTREZIDs (now a vector) ---
cat("Number of ENTREZIDs after removing NAs:", length(gene_entrez), "\n")
if(length(gene_entrez) == 0) {
  stop("No ENTREZIDs remaining after mapping and NA removal. Cannot perform enrichment.")
}
cat("First few ENTREZIDs (vector format):\n")
print(head(gene_entrez)) # See what the IDs look like (vector)

#####################################
# Pathway Enrichment Analysis (ORA) #
#####################################

# KEGG Pathway Analysis:
kegg_results <- enrichKEGG(gene = gene_entrez,
                           organism = "mmu",
                           pvalueCutoff = 1,  # Most lenient p-value = 1; default is 0.05
                           qvalueCutoff = 1)  # Most lenient q-value = 1; default is 0.05
if (is.null(kegg_results) || nrow(kegg_results) == 0) {
  message("Still no KEGG results even with pvalueCutoff=1 and qvalueCutoff=1. This points to deeper issues.")
} else {
  message("KEGG results found with lenient cutoffs!")
  print(head(kegg_results))
}

# View top 10 KEGG pathways:
head(kegg_results)

# (Scatterplot) Visualization of KEGG Pathways:
se_kegg <- dotplot(kegg_results, showCategory = 10) +
  ggtitle("Top 10 Enriched KEGG Pathways") +
  theme_bw()

ggsave("mouse_data/rmats_plots/se_kegg.png", plot = se_kegg, width = 6, height = 4, dpi = 300)

##############################
# Reactome Pathway Analysis: #
##############################

# (Scatterplot) visualization of Reactome pathways:

library(ReactomePA)
library(org.Mm.eg.db) # mouse

# Perform Reactome pathway enrichment analysis
reactome_results <- enrichPathway(gene         = gene_entrez,
                                     organism     = "mouse", 
                                     pvalueCutoff = 0.05, # set to 1 if no results to check for errors
                                     pAdjustMethod = "BH", # Benjamini-Hochberg
                                     qvalueCutoff = 0.05, # set to 1 if no results to check for errors
                                     readable = TRUE) 

head(reactome_results)

reactome_results@result$GeneRatio <- parse_ratio(reactome_results@result$GeneRatio)

plot_data <- reactome_results@result %>%
  arrange(p.adjust) %>% # Sort by p.adjust (ascending)
  head(10)

# Create the dotplot manually using ggplot2
se_reactome <- ggplot(plot_data,
                           aes(x = GeneRatio,
                               y = reorder(Description, GeneRatio), # Order by GeneRatio for aesthetic
                               size = Count,
                               color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") + # Adjust colors as desired
  scale_size_continuous(range = c(2, 10), name = "Gene Count") + # Adjust size range
  labs(title = "Reactome Enrichment (ORA)",
       x = "Gene Ratio",
       y = "Pathway Description") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10), # Adjust y-axis text
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# Save the plot
ggsave(filename = "mouse_data/rmats_plots/se_reactome.png",
       plot = se_reactome,
       width = 8, height = 6, dpi = 300)

# GO Analysis:
go_results <- enrichGO(gene = gene_entrez, 
                         OrgDb = org.Mm.eg.db, 
                         keyType = "ENTREZID", 
                         ont = "ALL",  # BP=Biological Process, ”CC"=Cellular Component, "MF"=Molecular Function, or "ALL" 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, # cutoff for adjusted p-values (FDR) 
                         readable = TRUE)

# View the results:
head(go_results)

# Scatterplot: Shows the top 10 enriched GO terms
se_go_plot <- dotplot(go_results, showCategory = 10) +
  ggtitle("Top 10 GO Terms") +
  theme_bw()
ggsave("mouse_data/rmats_plots/se_go_plot.png", plot = se_go_plot, width = 8, height = 6, units = "in", dpi = 300)

# Enrichment map: Visualizes the relationships between enriched GO terms based on shared genes.

# BiocManager::install(“enrichplot")
library(enrichplot)
library(GOSemSim)

mouse_go_data <- godata(annoDb = 'org.Mm.eg.db', ont = "BP") # (only BP, CC, MF not ALL)

go_enrich <- pairwise_termsim(go_results, method = "Wang", semData = mouse_go_data)

# Show a network of the top 15 terms: 
go_enrich_map <- emapplot(go_enrich, showCategory = 15, layout = "kk") +
  theme_bw()
ggsave("mouse_data/rmats_plots/se_go_enrich_map.png", plot = go_enrich_map, width = 8, height = 6, units = "in", dpi = 300)

# Gene-Concept network: Shows which of your input genes are involved in the enriched GO terms.
library(enrichplot)

go_results <- pairwise_termsim(go_results)

gc_network <- emapplot(go_results, showCategory = 10) +
  theme_bw()

ggsave("mouse_data/rmats_plots/se_go_emapplot.png", plot = gc_network, width = 8, height = 6, dpi = 300)

#######################################
# Gene Set Enrichment Analysis (GSEA) #
#######################################
# GSEA ranks genes based on splicing change magnitude (e.g., ΔPSI or log2 fold change).

library(dplyr)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggplot2)

# Get working directory:
getwd()

# Prepare gene list ranked by log2 fold change (ΔPSI):
rmats_data_gsea <- rmats_data %>%
  filter(!grepl("NA", IncLevel1) & !grepl("NA", IncLevel2)) %>%
  mutate(
    PSI_Group1 = sapply(strsplit(IncLevel1, ","), function(x) mean(as.numeric(x), na.rm = TRUE)),
    PSI_Group2 = sapply(strsplit(IncLevel2, ","), function(x) mean(as.numeric(x), na.rm = TRUE)),
    Delta_PSI = PSI_Group2 - PSI_Group1,
    # log2_Delta_PSI = log2(PSI_Group2 + 1e-6) - log2(PSI_Group1 + 1e-6),
    neg_log10_FDR = -log10(FDR),
    log2FC = log2(PSI_Group2 + 1e-6) - log2(PSI_Group1 + 1e-6)
  )

# Deduplicate all genes by max |log2FC|:
rmats_dedup <- rmats_data_gsea %>%
  group_by(geneSymbol) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()

# Map gene symbols to Entrez
gene_map <- bitr(rmats_dedup$geneSymbol, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)

# Merge the two datasets
rmats_merged <- inner_join(rmats_dedup, gene_map, by = c("geneSymbol" = "SYMBOL"))

# Create GSEA gene list
gene_list <- rmats_merged$log2FC
names(gene_list) <- rmats_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA with Reactome Pathways:
# Load pathway gene sets (e.g., Reactome, KEGG)
reactome_pathways <- reactomePathways(names(gene_list))

set.seed(42)  # for reproducibility
gene_list_jittered <- gene_list + runif(length(gene_list), min = -1e-6, max = 1e-6)
gene_list_jittered <- sort(gene_list_jittered, decreasing = TRUE)

gsea_results_reactome <- fgsea(
  pathways = reactome_pathways,
  stats    = gene_list_jittered,
  minSize  = 15,
  maxSize  = 500
)

# View top GSEA pathways
head(gsea_results_reactome)

# Plot top 5 GSEA pathways
library(ggplot2)

top_gsea_reactome <- gsea_results_reactome[order(gsea_results_reactome$padj), ][1:5, ]

gsea_reactome_plot <- ggplot(top_gsea_reactome, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 5 Enriched Pathways by GSEA", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_bw()

ggsave("mouse_data/rmats_plots/se_gsea_reactome_plot.png", plot = gsea_reactome_plot, width = 8, height = 6, units = "in", dpi = 300)

# GSEA with KEGG Pathway Enrichment:

gsea_kegg <- gseKEGG(
    geneList    = gene_list_jittered,
    organism    = "mmu",      # mmu = Mus musculus
    minGSSize   = 15,
    maxGSSize   = 500,
    pvalueCutoff = 0.05 # set to 1 if no results to troubleshoot
  )

# Order and filter the top 5 enriched KEGG pathways by adjusted p-value
gsea_kegg_result <- gsea_kegg@result
top_gsea_kegg <- gsea_kegg_result[order(gsea_kegg_result$p.adjust), ][1:5, ]

# Plot top 5 KEGG GSEA pathways
gsea_kegg_plot <- ggplot(top_gsea_kegg, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top 5 Enriched KEGG Pathways by GSEA",
    x = "KEGG Pathway",
    y = "Normalized Enrichment Score (NES)"
  )

ggsave("mouse_data/rmats_plots/se_gsea_kegg_plot.png", plot = gsea_kegg_plot, width = 8, height = 6, units = "in", dpi = 300)

# GSEA with GO with BP chosen:

gsea_go <- gseGO(
  geneList = gene_list_jittered,
  OrgDb = org.Mm.eg.db,
  ont = "BP",      # BP, CC, MF, or "ALL"
  keyType = "ENTREZID",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

gsea_go_result <- as.data.frame(gsea_go)

# Top 5 GO terms by adjusted p-value
top_gsea_go <- gsea_go_result[order(gsea_go_result$p.adjust), ][1:5, ]

gsea_go_plot <- ggplot(top_gsea_go, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 5 Enriched GO Biological Processes by GSEA",
    x = "GO Term",
    y = "Normalized Enrichment Score (NES)"
  )

ggsave("mouse_data/rmats_plots/se_gsea_go_plot.png", plot = gsea_go_plot, width = 8, height = 6, units = "in", dpi = 300)

############
## Heatmap #
############

# Heatmap for alternative splicing
library(dplyr)
library(tidyr)

library(pheatmap)

# Filter for top 10 significant events based on FDR:
heatmap_significant_events <- rmats_data_ma %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)  # Select top 10 significant events

# Extract PSI matrix for heatmap
heatmap_psi_matrix <- heatmap_significant_events %>%
  select(ID, geneSymbol, exonStart_0base, exonEnd, PSI_Group1, PSI_Group2) %>%
  mutate(identifier = paste(geneSymbol, ID, sep = "_")) %>%
  pivot_longer(
    cols = starts_with("PSI_Group"), 
    names_to = "Group", 
    values_to = "PSI"
  ) %>%
  pivot_wider(
    id_cols = identifier,
    names_from = Group, 
    values_from = PSI) %>%
  rename(PSI_Baseline = PSI_Group1, PSI_Drug = PSI_Group2)

# Set rownames as event IDs
heatmap_psi_matrix <- as.data.frame(heatmap_psi_matrix) 

rownames(heatmap_psi_matrix) <- heatmap_psi_matrix$identifier

heatmap_psi_matrix <- heatmap_psi_matrix %>% select(-identifier)

# heatmap_psi_matrix <- heatmap_psi_matrix[, !(colnames(heatmap_psi_matrix) %in% c("geneSymbol", "exonStart_0base", "exonEnd", "identifier"))]

head(heatmap_psi_matrix)

# Z-score normalization
psi_zscore <- t(scale(t(heatmap_psi_matrix)))

# Define the file path
output_file_path <- "mouse_data/rmats_plots/se_heatmap.png"

# Define legend breaks and labels for the custom title
# Calculate a sensible range for Z-scores to set legend breaks
min_z <- floor(min(psi_zscore))
max_z <- ceiling(max(psi_zscore))
# Example breaks: -2, 0, 2. Adjust based on your actual z-score range.
# Ensure the last label is your desired title.
custom_legend_breaks <- c(min_z, 0, max_z)
# custom_legend_labels <- c(as.character(min_z), "0", paste(max_z, "\nNormalized PSI")) # Use \n for line break
custom_legend_labels <- c("Normalized PSI\n\n", "0", as.character(max_z))

# Save the heatmap using a graphics device
png(filename = output_file_path, width = 6, height = 4, units = "in", res = 300)

# Generate Heatmap
pheatmap(
  psi_zscore,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of Significant Splicing Events",
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  legend_breaks = custom_legend_breaks,
  legend_labels = custom_legend_labels
)
dev.off() # Close the graphics device

# Try complexheatmap

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize) 

heatmap_significant_events <- rmats_data_ma %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)  # Select top 10 significant events

# Extract PSI matrix for heatmap
heatmap_psi_matrix <- heatmap_significant_events %>%
  select(ID, geneSymbol, exonStart_0base, exonEnd, PSI_Group1, PSI_Group2) %>%
  mutate(identifier = paste(geneSymbol, ID, sep = "_")) %>%
  pivot_longer(
    cols = starts_with("PSI_Group"), 
    names_to = "Group", 
    values_to = "PSI"
  ) %>%
  pivot_wider(
    id_cols = identifier,
    names_from = Group, 
    values_from = PSI) %>%
  rename(PSI_Baseline = PSI_Group1, PSI_Drug = PSI_Group2)

# Set rownames as event IDs
heatmap_psi_matrix <- as.data.frame(heatmap_psi_matrix) 

rownames(heatmap_psi_matrix) <- heatmap_psi_matrix$identifier

heatmap_psi_matrix <- heatmap_psi_matrix %>% select(-identifier) %>% as.matrix() 

# heatmap_psi_matrix <- heatmap_psi_matrix[, !(colnames(heatmap_psi_matrix) %in% c("geneSymbol", "exonStart_0base", "exonEnd", "identifier"))]

head(heatmap_psi_matrix)

# Z-score normalization
psi_zscore <- t(scale(t(heatmap_psi_matrix)))

col_fun = colorRamp2(c(min(psi_zscore), 0, max(psi_zscore)), c("blue", "white", "red"))

# Define the file path
output_file_path <- "mouse_data/rmats_plots/se_heatmap.png"

# Save the heatmap using a graphics device
png(filename = output_file_path, width = 6, height = 4, units = "in", res = 300)

# Generate Heatmap
draw(Heatmap(
  matrix = psi_zscore,
  name = "Z-score of PSI", # This is the name given to the heatmap itself, and also the default legend title
  col = col_fun, # Use the color mapping function
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Heatmap of Top 10 Significant Splicing Events",
  row_names_gp = gpar(fontsize = 6), # Font size for row names
  column_names_gp = gpar(fontsize = 8), # Font size for column names
  border = TRUE, # Adds a border around the heatmap body
  # Controls for the heatmap legend
  heatmap_legend_param = list(
    title = "Normalized PSI", # This is the direct argument for the legend title
    direction = "vertical", # or "horizontal"
    at = c(floor(min(psi_zscore)), 0, ceiling(max(psi_zscore))), # Explicit breaks for the legend
    labels = c(floor(min(psi_zscore)), "0", ceiling(max(psi_zscore))) # Labels for the breaks
  )
))
dev.off() # Close the graphics device










