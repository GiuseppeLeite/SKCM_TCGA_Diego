# --------------------------------------------------------------------------------------------------
# Comprehensive RNA-seq Data Analysis Pipeline for TCGA: Normalization, Clustering, DEG Analysis, & Volcano Plot Visualization
# Giuseppe Leite, PhD


# ---------------------------- Load Required Libraries ---------------------------- 
library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(tidyverse)
library(DESeq2)
library(parameters)      
library(pheatmap)
library(RColorBrewer)

# ---------------------------- Step 1: Data Information --------------------------- 
# Export UCSC Xena metadata information to CSV for reference
if (!file.exists("00_tbl_xena_hub_info.csv")) {
  data(XenaData)
  write.csv(XenaData, "00_tbl_xena_hub_info.csv")
  cat("Xena data information exported successfully.\n")
} else {
  cat("Xena data information already exists. Skipping export.\n")
}

# ---------------------- Step 2: Download TCGA Expression Data --------------------- 
# Download RSEM expected counts (log2(count+1)) from the Toil recompute
expression_file <- "./tcga_gene_expected_count.gz"

if (!file.exists(expression_file)) {
  gene_expected_cnt_toil <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
    XenaFilter(filterCohorts = "PANCAN") %>%
    XenaFilter(filterDatasets = "tcga_gene_expected_count")
  
  XenaQuery(gene_expected_cnt_toil) %>%
    XenaDownload(destdir = "./")
  cat("Gene expression data downloaded successfully.\n")
} else {
  cat("Gene expression data already exists. Skipping download.\n")
}

# ---------------------------- Step 3: Clinical Data ------------------------------ 
# Load clinical data and extract sample identifiers
sample_annotation <- read.csv("sample_annotation_cleaned.csv", sep = ";")

# Include "sample" to ensure row identifiers from the expression dataset
filter_expr <- c(sample_annotation$SAMPLE_ID, "sample")

# ---------------------------- Step 4: Expression Data ---------------------------- 
# Efficiently load only selected columns (samples of interest)
expr_subset_by_samp <- fread(expression_file, select = filter_expr)
expr_subset_by_samp <- column_to_rownames(expr_subset_by_samp, var = "sample")

# ---------------------- Step 5: Align Clinical and Expression Data ---------------- 
colnames_selected <- intersect(sample_annotation$SAMPLE_ID, colnames(expr_subset_by_samp))
expr_subset_by_samp <- expr_subset_by_samp[, colnames_selected]
sample_annotation <- sample_annotation[sample_annotation$SAMPLE_ID %in% colnames_selected, ]

# Confirm alignment explicitly
stopifnot(identical(colnames(expr_subset_by_samp), sample_annotation$SAMPLE_ID))

# ---------------------------- Step 6: Gene Mapping ------------------------------- 
# Load and merge gene annotation mappings
if (!file.exists("gencode.v23.annotation.gene.probemap.txt")) {
  stop("Gene mapping file missing. Please download it and place it in the working directory.")
}
gencode_probemap <- fread("gencode.v23.annotation.gene.probemap.txt")
probemap <- fread("zz_gencode.v23.annotation.csv")
probemap <- merge(probemap, gencode_probemap, by = "Gene_Symbol") %>% select(id, Gene_Symbol)

# ----------------------- Step 7: Subset Protein-Coding Genes ---------------------- #
expr_subset_by_samp <- rownames_to_column(expr_subset_by_samp, var = "id")
expr_all <- merge(probemap, expr_subset_by_samp, by = "id")

# Remove duplicates by keeping highest expressed
expr_final <- expr_all %>%
  mutate(mean_expr = rowMeans(select(., -id, -Gene_Symbol), na.rm = TRUE)) %>%
  arrange(Gene_Symbol, desc(mean_expr)) %>%
  distinct(Gene_Symbol, .keep_all = TRUE) %>%
  select(-mean_expr, -id) %>%
  column_to_rownames(var = "Gene_Symbol")

# --------------------- Step 8: Back-transform Expression Data --------------------- 
# Convert from log2(expected_counts + 1) back to expected counts
expr_bt <- round((2^expr_final) - 1, 0)

# remove low count genes 

# Keep genes with count >= 10 in at least 5% of samples
min_count <- 10
min_samples <- ceiling(0.05 * 361)

keep_genes <- rowSums(expr_bt >= min_count) >= min_samples
expr_bt_filtered <- expr_bt[keep_genes, ]

cat("Number of genes before filtering:", nrow(expr_bt), "\n")
cat("Number of genes after filtering:", sum(keep_genes), "\n")

# Confirm final alignment again explicitly
stopifnot(identical(colnames(expr_bt_filtered), sample_annotation$SAMPLE_ID))

write.csv(expr_bt_filtered, "01_expected_cnt_bt.csv")

# ---------------------- Step 9: Normalization with DESeq2 ------------------------- 
# Normalize data (no differential expression testing, design ~1)
col_data <- data.frame(row.names = colnames(expr_bt_filtered))

ds <- DESeqDataSetFromMatrix(countData = expr_bt_filtered, colData = col_data, design = ~ 1)

# Estimates dispersions
ds <- DESeq(ds) # this triggers a warning due to design = ~1, which is expected

vsd <- vst(ds)  # Variance stabilizing transformation

# ------------------- Step 10: Data Preparation for Clustering --------------------- 
# Transpose normalized data for analysis of selected genes
vsd_df <- assay(vsd) %>% as.data.frame() %>% rownames_to_column("gene_symbol")

vsd_transposed <- vsd_df %>% column_to_rownames("gene_symbol") %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample_id")

genes_of_interest <- c("HMGB1", "CALR", "GSDMA", "GSDMB", "GSDMC", "GSDMD", "MLKL", "PRDX3")

gene_data <- select(vsd_transposed, sample_id, all_of(genes_of_interest))

gene_data_scaled <- gene_data %>% mutate(across(-sample_id, ~as.numeric(scale(.x))))

# ------------------------- Step 11: Clustering Analysis --------------------------- 

n_clust <- n_clusters(gene_data_scaled[,-1], standardize = FALSE, package = "all", fast = TRUE)

print(n_clust)
plot(n_clust)

# Run k-means (example with 4 clusters)
k <- 4
km <- kmeans(gene_data_scaled[,-1], centers = k, nstart = 50)
gene_data_scaled$cluster <- factor(km$cluster)

s1 <- plot(n_clust) +
  theme_bw() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black")
  )
print(s1)

# Save the cluster determination plot
dev.copy2pdf(file = "Optimal_Clusters.pdf", width = 5, height = 4)
dev.off()


# (Opcional) salvar
write.csv(gene_data_scaled, "Clustered_Samples_VSD.csv", row.names = FALSE)

# -------------------------- Step 12: PCA and Heatmap ------------------------------ 
# PCA
pca_result <- prcomp(gene_data_scaled[,2:9])


pca_data <- as.data.frame(pca_result$x) %>%
  mutate(cluster = gene_data_scaled$cluster)  # use correct column name

cluster_colors <- c("1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73", "4" = "#D55E00")

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(
    title = "Principal Component Analysis",
    x = "PCA 1", y = "PCA 2",
    color = "Cluster"
  ) +
  theme_bw() +
  scale_color_manual(values = cluster_colors)
print(pca_plot)

# Save PCA Plot as PDF
ggsave("PCA_Clusters.pdf", plot = pca_plot, width = 5, height = 4)

# ---------------------- Step 13: Heatmap Visualization of Clusters ---------------------- 

# Set sample_id as row names to ensure correct alignment
rownames(gene_data_scaled) <- gene_data_scaled$sample_id

# Subset gene expression data, removing sample_id and cluster columns
gene_expr <- gene_data_scaled[, !(colnames(gene_data_scaled) %in% c("sample_id", "cluster"))]

# Create a data frame indicating cluster membership for annotations
annotation_row <- data.frame(cluster = factor(gene_data_scaled$cluster))
rownames(annotation_row) <- gene_data_scaled$sample_id

# Order samples by cluster to clearly visualize grouping patterns
ordered_samples <- order(annotation_row$cluster)
gene_expr <- gene_expr[ordered_samples, ]
annotation_row <- annotation_row[ordered_samples, , drop = FALSE]

# Define a color scale for expression z-scores (blue: low, white: neutral, red: high)
breaks <- seq(-3, 3, length.out = 101)
colors <- colorRampPalette(c("#00008B", "white", "#f70d1a"))(100)

cluster_colors <- c("1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73", "4" = "#D55E00")
annotation_colors <- list(cluster = cluster_colors)


# Generate and display the heatmap (samples clustered manually, genes clustered automatically)
pheatmap_result <- pheatmap(gene_expr,
         annotation_row = annotation_row,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colors,
         breaks = breaks,
         annotation_colors = annotation_colors)

pheatmap_result

# Save PCA Plot as PDF

dev.copy2pdf(file = "heatmap.pdf", width = 5, height = 6)
dev.off()

# ------------------- Step 14: Differential Expression Analysis by Cluster -------------------

# Extract the sample IDs from the count matrix
sample_ids_in_counts <- colnames(expr_bt_filtered)

# Match cluster assignments from gene_data_scaled to the sample IDs
# gene_data_scaled should have columns "sample_id" and "cluster"
cluster_vector <- gene_data_scaled$cluster[match(sample_ids_in_counts, gene_data_scaled$sample_id)]

# Safety check to ensure all samples have a corresponding cluster assignment
if (any(is.na(cluster_vector))) {
  stop("Mismatch between sample IDs in count matrix and cluster assignment!")
}

# Create the colData data frame for DESeq2 using the cluster information
col_data_cluster <- data.frame(cluster = factor(cluster_vector))
rownames(col_data_cluster) <- sample_ids_in_counts

# Create a DESeqDataSet using the count data and the new colData
ds_cluster <- DESeqDataSetFromMatrix(countData = expr_bt_filtered, 
                                     colData = col_data_cluster, 
                                     design = ~ cluster)

# Set cluster "2" as the reference level.
ds_cluster$cluster <- relevel(ds_cluster$cluster, ref = "2")

# Run the DESeq2 analysis
ds_cluster <- DESeq(ds_cluster)


# Compare Cluster 2 vs. Cluster 1.
res_cluster2_vs_cluster1 <- results(ds_cluster, contrast = c("cluster", "1", "2"))
res_cluster2_vs_cluster1$log2FoldChange <- -res_cluster2_vs_cluster1$log2FoldChange
summary(res_cluster2_vs_cluster1)

# Compare Cluster 2 vs. Cluster 3.
res_cluster2_vs_cluster3 <- results(ds_cluster, contrast = c("cluster", "3", "2"))
res_cluster2_vs_cluster3$log2FoldChange <- -res_cluster2_vs_cluster3$log2FoldChange
summary(res_cluster2_vs_cluster3)

# Compare Cluster 2 vs. Cluster 4.
res_cluster2_vs_cluster4 <- results(ds_cluster, contrast = c("cluster", "4", "2"))
res_cluster2_vs_cluster4$log2FoldChange <- -res_cluster2_vs_cluster4$log2FoldChange
summary(res_cluster2_vs_cluster4)

# Save the differential expression results to CSV files.
write.csv(as.data.frame(res_cluster2_vs_cluster1), file = "DE_results_cluster2_vs_cluster1.csv", row.names = TRUE)
write.csv(as.data.frame(res_cluster2_vs_cluster3), file = "DE_results_cluster2_vs_cluster3.csv", row.names = TRUE)
write.csv(as.data.frame(res_cluster2_vs_cluster4), file = "DE_results_cluster2_vs_cluster4.csv", row.names = TRUE)

# ------------------- Step 15: Define a reusable function to generate volcano plots ----------------------------

create_volcano_plot <- function(data, top_genes, title, file_name, fc_cutoff = 0.38, p_cutoff = 0.048) {
  # Create custom key values for coloring the volcano plot based on log2 fold change and significance
  keyvals <- ifelse(
    data$log2FoldChange < -fc_cutoff & data$padj < p_cutoff, "#00008B",  # Down-regulated genes in blue
    ifelse(data$log2FoldChange > fc_cutoff & data$padj < p_cutoff, "#f70d1a", "gray")  # Up-regulated genes in red, non-significant as gray
  )
  keyvals[is.na(keyvals)] <- "gray"
  
  # Define names for the key values for legend annotation
  names(keyvals)[keyvals == "#f70d1a"] <- "Up-regulated genes"
  names(keyvals)[keyvals == "gray"] <- "NS"
  names(keyvals)[keyvals == "#00008B"] <- "Down-regulated genes"
  
  # Generate the volcano plot using EnhancedVolcano
  volcano_plot <- EnhancedVolcano(
    data, 
    lab = rownames(data),           # Label each point using row names (typically gene symbols)
    x = "log2FoldChange", 
    y = "padj",
    title = title,
    selectLab = rownames(top_genes), # Highlight the top genes of interest
    pCutoff = p_cutoff, 
    FCcutoff = fc_cutoff,
    pointSize = 2.0, 
    labSize = 2.9,
    colCustom = keyvals, 
    colAlpha = 4 / 5,
    boxedLabels = TRUE, 
    legendPosition = "none",
    legendLabSize = 9, 
    legendIconSize = 3.0,
    drawConnectors = TRUE, 
    widthConnectors = 0.3,
    colConnectors = "black",
    gridlines.major = TRUE, 
    gridlines.minor = TRUE,
    border = "full", 
    borderWidth = 0.5, 
    borderColour = "black"
  ) +
    ylab("-log10(BH)") +
    theme(
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.ticks = element_line(size = 0.5),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # Save the plot as a PDF file with the specified file name
  ggsave(file_name, plot = volcano_plot, width = 4, height = 5)
  
  # Return the plot for immediate visualization
  return(volcano_plot)
}

# ----------------------- Step 16: Generate Volcano Plots for Each Contrast ---------------------------------------

# Contrast 1: Volcano Plot for Cluster 2 vs Cluster 1
res_cluster2_vs_cluster1_df <- as.data.frame(res_cluster2_vs_cluster1)
cluster2_vs_cluster1_sig <- res_cluster2_vs_cluster1_df %>% filter(padj <= 0.05)
cluster2_vs_cluster1_sig <- cluster2_vs_cluster1_sig[order(cluster2_vs_cluster1_sig$log2FoldChange), ]
top_genes_21 <- rbind(head(cluster2_vs_cluster1_sig, 5), tail(cluster2_vs_cluster1_sig, 5))
volcano_cluster2_vs_cluster1 <- create_volcano_plot(
  data = res_cluster2_vs_cluster1,
  top_genes = top_genes_21,
  title = "Cluster 2 vs Cluster 1",
  file_name = "cluster2_vs_cluster1.pdf"
)
print(volcano_cluster2_vs_cluster1)


# Contrast 2: Volcano Plot for Cluster 2 vs Cluster 3
res_cluster2_vs_cluster3_df <- as.data.frame(res_cluster2_vs_cluster3)
cluster2_vs_cluster3_sig <- res_cluster2_vs_cluster3_df %>% filter(padj <= 0.05)
cluster2_vs_cluster3_sig <- cluster2_vs_cluster3_sig[order(cluster2_vs_cluster3_sig$log2FoldChange), ]
top_genes_23 <- rbind(head(cluster2_vs_cluster3_sig, 5), tail(cluster2_vs_cluster3_sig, 5))
volcano_cluster2_vs_cluster3 <- create_volcano_plot(
  data = res_cluster2_vs_cluster3,
  top_genes = top_genes_23,
  title = "Cluster 2 vs Cluster 3",
  file_name = "cluster2_vs_cluster3.pdf"
)
print(volcano_cluster2_vs_cluster3)


# Contrast 3: Volcano Plot for Cluster 2 vs Cluster 4
res_cluster2_vs_cluster4_df <- as.data.frame(res_cluster2_vs_cluster4)
cluster2_vs_cluster4_sig <- res_cluster2_vs_cluster4_df %>% filter(padj <= 0.05)
cluster2_vs_cluster4_sig <- cluster2_vs_cluster4_sig[order(cluster2_vs_cluster4_sig$log2FoldChange), ]
top_genes_24 <- rbind(head(cluster2_vs_cluster4_sig, 5), tail(cluster2_vs_cluster4_sig, 5))
volcano_cluster2_vs_cluster4 <- create_volcano_plot(
  data = res_cluster2_vs_cluster4,
  top_genes = top_genes_24,
  title = "Cluster 2 vs Cluster 4",
  file_name = "cluster2_vs_cluster4.pdf"
)
print(volcano_cluster2_vs_cluster4)



