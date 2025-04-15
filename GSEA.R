# Gene Set Enrichment Analysis (GSEA)
# This method helps to identify gene sets that are significantly enriched
# in a specific biological condition, such as elderly versus adult.

# Load the required libraries for data manipulation, GSEA, and access to the MSigDB database.
library(tidyverse)
library(fgsea)
library(msigdbr)
library(stats)

# --------------------- de_results_cluster2_vs_cluster1 -------------------------
# Read the data file containing the results for the elderly versus adult analysis.
# Here we assume the file contains a list of genes with associated values.
de_results_cluster2_vs_cluster1 <- read.csv("DE_results_cluster2_vs_cluster1.csv", row.names = 1)

# Add a column with gene symbols (from row names) for further processing.
de_results_cluster2_vs_cluster1 <- rownames_to_column(de_results_cluster2_vs_cluster1, var = "gene_symbol")

# Extract the ranking values (stat values) for each gene to be used in GSEA.
ranks <- de_results_cluster2_vs_cluster1 %>%
  select(gene_symbol, stat)

# Order the rankings in decreasing order based on the 'stat' column.
ranks <- ranks[order(ranks$stat, decreasing = TRUE), ]

# Convert the ranking table to a named vector (names = gene symbols, values = stat).
ranks <- deframe(ranks)
head(ranks, 20)  # Display the first 20 genes with their ranks for inspection

# Define the list of Hallmark gene sets (well-defined biological processes) using msigdbr.
# These gene sets are groups of genes associated with specific biological processes or conditions,
# such as inflammation or cell proliferation.
pathways_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = map(genes, compose(as_vector, unname))) %>%
  deframe()

# Perform the enrichment analysis (GSEA) to identify which Hallmark gene sets
# are most associated with the difference between elderly and adult, using 10,000 permutations.
set.seed(123)  # Set a random seed for reproducibility.
fgsea_res <- fgsea(
  pathways = pathways_hallmark,
  stats = ranks,
  nPerm = 10000,
  minSize = 15,
  maxSize = 500
)

# Remove the 'leadingEdge' information from the results (genes directly driving the enrichment)
# to simplify result visualization.
fgsea_res$leadingEdge <- NULL

# Save the Hallmark analysis results to a file for later inspection.
write.table(fgsea_res, file = "DE_results_cluster2_vs_cluster1_hallmark_MSigDB.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# --------------------- de_results_cluster2_vs_cluster3 -------------------------
# Read the data file containing the results for the elderly versus adult analysis.
de_results_cluster2_vs_cluster3 <- read.csv("DE_results_cluster2_vs_cluster3.csv", row.names = 1)

# Add gene symbols as a column.
de_results_cluster2_vs_cluster3 <- rownames_to_column(de_results_cluster2_vs_cluster3, var = "gene_symbol")

# Extract the ranking values (stat values) for each gene.
ranks <- de_results_cluster2_vs_cluster3 %>%
  select(gene_symbol, stat)

# Order the rankings in decreasing order.
ranks <- ranks[order(ranks$stat, decreasing = TRUE), ]

# Convert the ranking table to a named vector.
ranks <- deframe(ranks)
head(ranks, 20)  # Display the top 20 genes with their ranks for inspection

# Define the list of Hallmark gene sets using msigdbr.
pathways_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = map(genes, compose(as_vector, unname))) %>%
  deframe()

# Perform GSEA for the analysis using 10,000 permutations.
set.seed(123)  # Set a random seed for reproducibility.
fgsea_res <- fgsea(
  pathways = pathways_hallmark,
  stats = ranks,
  nPerm = 10000,
  minSize = 15,
  maxSize = 500
)

# Remove the 'leadingEdge' information from the results.
fgsea_res$leadingEdge <- NULL

# Save the Hallmark analysis results to a file.
write.table(fgsea_res, file = "DE_results_cluster2_vs_cluster3_hallmark_MSigDB.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# --------------------- de_results_cluster2_vs_cluster4 -------------------------
# Read the data file containing the results for the elderly versus adult analysis.
de_results_cluster2_vs_cluster4 <- read.csv("DE_results_cluster2_vs_cluster4.csv", row.names = 1)

# Add gene symbols as a column.
de_results_cluster2_vs_cluster4 <- rownames_to_column(de_results_cluster2_vs_cluster4, var = "gene_symbol")

# Extract the ranking values (stat values) for each gene to be used in GSEA.
ranks <- de_results_cluster2_vs_cluster4 %>%
  select(gene_symbol, stat)

# Order the rankings in decreasing order based on 'stat'.
ranks <- ranks[order(ranks$stat, decreasing = TRUE), ]

# Convert the ranking table to a named vector.
ranks <- deframe(ranks)
head(ranks, 20)  # Display the first 20 genes with their ranks for inspection

# Define the list of Hallmark gene sets using msigdbr.
pathways_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = map(genes, compose(as_vector, unname))) %>%
  deframe()

# Perform the GSEA to identify which Hallmark gene sets are most associated with the difference
# between elderly and adult, using 10,000 permutations for better accuracy.
set.seed(123)  # Set a random seed for reproducibility.
fgsea_res <- fgsea(
  pathways = pathways_hallmark,
  stats = ranks,
  nPerm = 10000,
  minSize = 15,
  maxSize = 500
)

# Remove the 'leadingEdge' information from the GSEA results.
fgsea_res$leadingEdge <- NULL

# Save the Hallmark analysis results to a file for later inspection.
write.table(fgsea_res, file = "DE_results_cluster2_vs_cluster4_hallmark_MSigDB.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
