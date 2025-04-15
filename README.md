# TCGA RNA-seq Analysis Pipeline

## Overview

This repository contains an end-to-end RNA-seq data analysis pipeline for The Cancer Genome Atlas (TCGA) RNA-seq data. The pipeline performs multiple steps including:

- **Data Acquisition & Preprocessing:**  
  Downloads TCGA RNA-seq expression data from UCSC Xena, back-transforms log₂ expression counts, maps gene identifiers, and filters protein-coding genes.

- **Normalization & Clustering:**  
  Normalizes data with DESeq2 (using variance stabilizing transformation), clusters samples (with PCA and heatmap visualization), and conducts differential expression analysis.

- **Differential Expression & GSEA:**  
  Performs differential expression analysis (DESeq2) and Gene Set Enrichment Analysis (GSEA) to identify enriched gene sets (e.g., Hallmark pathways) based on the differential expression results.

- **Visualization:**  
  Generates volcano plots, bar plots, and Kaplan–Meier survival curves for data exploration and result presentation.

## Repository Files

- **TCGA_RNAseq_Pipeline.R.R**  
  The main R script that executes the complete RNA-seq analysis pipeline—from data download and preprocessing to clustering, normalization, and differential expression analysis.

- **GSEA.R**  
  An R script to perform Gene Set Enrichment Analysis (GSEA) using the fgsea and msigdbr packages. This script reads differential expression results and identifies significantly enriched gene sets.

- **barplot.R**  
  Contains code to generate bar plots of normalized enrichment scores (NES) for significant pathways identified by GSEA.

- **Kaplan_Meier.R**  
  Script for performing Kaplan–Meier survival analysis.

- **gencode.v23.annotation.gene.probemap.txt**  
  Gene annotation mappings used in the analysis for mapping gene identifiers.

- **zz_gene.protein.coding.csv**  
  A CSV file listing protein-coding genes.

- **sample_annotation_cleaned.csv**  
  Cleaned sample annotation data used for aligning the RNA-seq expression data with clinical metadata.

## Requirements

To run the scripts in this repository, you will need:

- **R Version:** 4.4.2 or later  
- **Packages:** tidyverse, DESeq2, fgsea, msigdbr, EnhancedVolcano, pheatmap, RColorBrewer, parameters, stats, ggplot2, scales, stringr


## How to Use

1. **Run the Main Pipeline:**  
   Execute `TCGA_RNAseq_Pipeline.R.R` to download, preprocess, normalize, cluster, and perform differential expression analysis on the TCGA RNA-seq data.

2. **Perform GSEA Analysis:**  
   Run `GSEA.R` to conduct Gene Set Enrichment Analysis based on the differential expression results. Ensure that the DE results files (e.g., `DE_results_cluster2_vs_cluster1.csv`) are in your working directory.

3. **Generate Visualization Plots:**  
   Run `barplot.R` and `Kaplan_Meier.R` to generate the respective visualization plots for enriched pathways and survival analysis.
