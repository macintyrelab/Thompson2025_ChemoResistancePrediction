# This script makes and saves the STAR gene expression data for the TCGA-OV cohort

# Libraries
rm(list = ls(all = TRUE))
library(TCGAbiolinks) # This code will only work using v2.29.6 
library(this.path)

# Paths
anthracycline_dir <- dirname(this.path())
input_data <- file.path(dirname(dirname(anthracycline_dir)), 'Input_Data')
setwd(file.path(input_data, 'TCGA_expression'))

# build query to GDC database
query_rna_seq <- GDCquery(project = 'TCGA-OV', 
                          data.category = 'Transcriptome Profiling', 
                          data.type = 'Gene Expression Quantification', 
                          experimental.strategy = 'RNA-Seq', 
                          workflow.type = 'STAR - Counts')

# download data
GDCdownload(query_rna_seq)

# load data to R
tcga_ov_rnaseq <- GDCprepare(query = query_rna_seq)

# save RNA-Seq data
saveRDS(tcga_ov_rnaseq, file = file.path(input_data, 'TCGA.OV_RNAseq.rds'))
