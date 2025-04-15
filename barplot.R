library(ggplot2)
library(scales)
library(stringr)

# --------------------- Graph 1 - Cluster 2 vs Cluster 1 ---------------------
# Load the data from the GSEA results file.
gsea_data_cluster2_vs_cluster1 <- read.delim("DE_results_cluster2_vs_cluster1_hallmark_MSigDB.txt", sep = "\t")

# Clean pathway names:
# - Remove the "HALLMARK_" prefix.
# - Replace underscores with spaces.
# - Capitalize the first letter of each pathway.
gsea_data_cluster2_vs_cluster1$pathway <- gsub("^HALLMARK_", "", gsea_data_cluster2_vs_cluster1$pathway)
gsea_data_cluster2_vs_cluster1$pathway <- gsub("_", " ", gsea_data_cluster2_vs_cluster1$pathway)
gsea_data_cluster2_vs_cluster1$pathway <- str_to_sentence(gsea_data_cluster2_vs_cluster1$pathway)

# Filter the data for pathways with an adjusted p-value (padj) less than 0.05.
filtered_data_cluster2_vs_cluster1 <- gsea_data_cluster2_vs_cluster1[gsea_data_cluster2_vs_cluster1$padj < 0.05, ]

# Create a bar plot for the filtered data.
plot_cluster2_vs_cluster1 <- ggplot(filtered_data_cluster2_vs_cluster1, aes(x = NES, y = reorder(pathway, NES), fill = NES)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.95) + 
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#f70d1a", midpoint = 0) +
  theme_bw() +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Add extra space on the x-axis
  scale_y_discrete(labels = label_wrap(45))

# Display the plot.
print(plot_cluster2_vs_cluster1)

# Save the plot as a PDF.
ggsave("cluster2_vs_cluster1_barplot.pdf", plot = plot_cluster2_vs_cluster1, width = 6, height = 5)



# --------------------- Graph 2 - Cluster 2 vs Cluster 3 ---------------------
# Load the data from the GSEA results file.
gsea_data_cluster2_vs_cluster3 <- read.delim("DE_results_cluster2_vs_cluster3_hallmark_MSigDB.txt", sep = "\t")

# Clean the pathway names.
gsea_data_cluster2_vs_cluster3$pathway <- gsub("^HALLMARK_", "", gsea_data_cluster2_vs_cluster3$pathway)
gsea_data_cluster2_vs_cluster3$pathway <- gsub("_", " ", gsea_data_cluster2_vs_cluster3$pathway)
gsea_data_cluster2_vs_cluster3$pathway <- str_to_sentence(gsea_data_cluster2_vs_cluster3$pathway)

# Filter the data for pathways with padj < 0.05.
filtered_data_cluster2_vs_cluster3 <- gsea_data_cluster2_vs_cluster3[gsea_data_cluster2_vs_cluster3$padj < 0.05, ]

# Create a bar plot for the filtered data.
plot_cluster2_vs_cluster3 <- ggplot(filtered_data_cluster2_vs_cluster3, aes(x = NES, y = reorder(pathway, NES), fill = NES)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.95) +
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#f70d1a", midpoint = 0) +
  theme_bw() +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Add extra space on the x-axis
  scale_y_discrete(labels = label_wrap(45))

# Display the plot.
print(plot_cluster2_vs_cluster3)

# Save the plot as a PDF.
ggsave("cluster2_vs_cluster3_barplot.pdf", plot = plot_cluster2_vs_cluster3, width = 6, height = 5)



# --------------------- Graph 3 - Cluster 2 vs Cluster 4 ---------------------
# Load the data from the GSEA results file.
gsea_data_cluster2_vs_cluster4 <- read.delim("DE_results_cluster2_vs_cluster4_hallmark_MSigDB.txt", sep = "\t")

# Clean the pathway names.
gsea_data_cluster2_vs_cluster4$pathway <- gsub("^HALLMARK_", "", gsea_data_cluster2_vs_cluster4$pathway)
gsea_data_cluster2_vs_cluster4$pathway <- gsub("_", " ", gsea_data_cluster2_vs_cluster4$pathway)
gsea_data_cluster2_vs_cluster4$pathway <- str_to_sentence(gsea_data_cluster2_vs_cluster4$pathway)

# Filter the data for pathways with padj < 0.05.
filtered_data_cluster2_vs_cluster4 <- gsea_data_cluster2_vs_cluster4[gsea_data_cluster2_vs_cluster4$padj < 0.05, ]

# Create a bar plot for the filtered data.
plot_cluster2_vs_cluster4 <- ggplot(filtered_data_cluster2_vs_cluster4, aes(x = NES, y = reorder(pathway, NES), fill = NES)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.95) +
  scale_fill_gradient2(low = "#00008B", mid = "white", high = "#f70d1a", midpoint = 0) +
  theme_bw() +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Add extra space on the x-axis
  scale_y_discrete(labels = label_wrap(45))

# Display the plot.
print(plot_cluster2_vs_cluster4)

# Save the plot as a PDF.
ggsave("cluster2_vs_cluster4_barplot.pdf", plot = plot_cluster2_vs_cluster4, width = 6, height = 5)
