# -------------------------------------------------------------------------------------
# Process and Combine GSEA Results for Multiple Cluster Comparisons
#
# This script reads GSEA results files for three contrasts (cluster2_vs_cluster1, 
# cluster2_vs_cluster3, and cluster2_vs_cluster4), cleans the pathway names, filters 
# for significantly enriched pathways (padj < 0.05), and combines the filtered results 
# into a single data frame with an added "contrast" label.
#
# Author: Giuseppe Leite
# Date: 15/04/2025
# -------------------------------------------------------------------------------------

library(tidyverse)
library(stringr)

# Define file paths and their corresponding contrast labels.
file_paths <- c(
  "DE_results_cluster2_vs_cluster1_hallmark_MSigDB.txt",
  "DE_results_cluster2_vs_cluster3_hallmark_MSigDB.txt",
  "DE_results_cluster2_vs_cluster4_hallmark_MSigDB.txt"
)
contrast_labels <- c("cluster2_vs_cluster1", "cluster2_vs_cluster3", "cluster2_vs_cluster4")

# Process and combine GSEA results from all contrasts.
# The map2_dfr() function iterates simultaneously over the file paths and contrast labels,
# reading, cleaning, filtering, and labeling each dataset before combining them into one data frame.
combined_data <- purrr::map2_dfr(file_paths, contrast_labels, function(file, contrast_label) {
  # Read the GSEA results from a tab-delimited file.
  gsea_data <- read.delim(file, sep = "\t")
  
  # Clean pathway names:
  #   1. Remove the "HALLMARK_" prefix.
  #   2. Replace underscores with spaces.
  #   3. Convert text to sentence case (only the first letter capitalized).
  gsea_data$pathway <- gsub("^HALLMARK_", "", gsea_data$pathway)
  gsea_data$pathway <- gsub("_", " ", gsea_data$pathway)
  gsea_data$pathway <- str_to_sentence(gsea_data$pathway)
  
  # Filter for significantly enriched pathways with an adjusted p-value (padj) less than 0.05.
  filtered_data <- gsea_data[gsea_data$padj < 0.05, ]
  
  # Add a new column indicating the contrast label for the current dataset.
  filtered_data$contrast <- contrast_label
  
  # Return the filtered dataset.
  filtered_data
})

# Preview the combined GSEA results data frame.
head(combined_data)
