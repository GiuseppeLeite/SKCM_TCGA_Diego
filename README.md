TCGA RNA-seq Analysis Pipeline

Overview

This repository contains an end-to-end RNA-seq data analysis pipeline for The Cancer Genome Atlas (TCGA) RNA-seq data. The pipeline performs multiple steps including:

Data Acquisition & Preprocessing:
Downloads TCGA RNA-seq expression data from UCSC Xena, back-transforms log2 expression counts, maps gene identifiers, and filters protein-coding genes.

Normalization & Clustering:
Normalizes data with DESeq2 (including variance stabilizing transformation), clusters samples (with PCA and heatmap visualization), and conducts differential expression analysis.

Differential Expression & GSEA:
Performs differential expression analysis (DESeq2) and Gene Set Enrichment Analysis (GSEA) to identify enriched gene sets (e.g., Hallmark pathways) based on the differential expression results.

Visualization:
Generates volcano plots, bar plots, and Kaplan–Meier survival curves for data exploration and result presentation.

Repository Files
TCGA_RNAseq_Pipeline.R.R
The main R script that executes the complete RNA-seq analysis pipeline—from data download and preprocessing to clustering, normalization, and differential expression analysis.

GSEA.R
An R script to perform Gene Set Enrichment Analysis (GSEA) using the fgsea and msigdbr packages. This script reads differential expression results and identifies significantly enriched gene sets.

barplot.R
Contains code to generate bar plots of normalized enrichment scores (NES) for significant pathways identified by GSEA.

Kaplan_Meier.R
Script for performing Kaplan–Meier survival analysis.

gencode.v23.annotation.gene.probemap.txt
Gene annotation mappings used in the analysis for mapping gene identifiers.

zz_gencode.v23.annotation.csv
Additional gene annotation information in CSV format.

zz_gene.protein.coding.csv
A CSV file listing protein-coding genes.

sample_annotation_cleaned.csv
Cleaned sample annotation data used for aligning the RNA-seq expression data with clinical metadata.
