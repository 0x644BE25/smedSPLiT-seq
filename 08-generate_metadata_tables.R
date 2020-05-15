######################################################
# GENERATE METADATA TABLES
#
# GOAL: Generate global and tissue-subset metadata
# tables including cell barcodes, tissues, global 
# clusters, tissue subclusters, sample, treatment,
# timepoint, RNA_nUMI, RNA_nFeature, and 6K subset 
# occupancy. Also generate tables of global and tissue 
# subset markers. Save as a multi-sheet XLSX file.
######################################################

# ================= IMPORTS ==========================

library(openxlsx)
library(Seurat)
source('./utils.R')

# ================= CELL-LEVEL METADATA ==============

s3 <- readRDS('~/RDS/FULL_150PCs_2020-02-17.Rds')
sixKcells <- readRDS('~/RDS/6K_150PCs_2020-01-09.Rds')
sixKcells <- Cells(sixKcells)

metadata <- s3@meta.data[,c('sample','treatment',
                            'Timepoint','nCount_RNA',
                            'nFeature_RNA','seurat_clusters',
                            'Tissue.Type','subclusters')]

metadata <- cbind(Cells(s3),metadata)

colnames(metadata) <- c('cell_ID','sample','treatment','timepoint',
                        'nUMI','nFeature','seurat_clusters',
                        'tissue_type','subclusters')
metadata$in_6K_subset <- FALSE
metadata[sixKcells,'in_6K_subset'] <- TRUE
metadata$barcode <- unlist(lapply(strsplit(unlist(lapply(strsplit(rownames(metadata), '_'), function(x){ x[2] })), '-'), function(x){ x[1] }))

# ================= MARKER DATA ======================

globalMarkers <- read.csv('~/publication_pipeline/global_cluster_markers.csv',stringsAsFactors=FALSE)
subclusterGlobalMarkers <- read.csv('~/publication_pipeline/subcluster_global_markers.csv',stringsAsFactors=FALSE)
tissueMarkers <- read.csv('~/publication_pipeline/tissue_subset_markers.csv',stringsAsFactors=FALSE)

# ================= WRITE XLSX =======================

openxlsx::write.xlsx(list('cell_metadata'=metadata,
                          'global_cluster_markers'=globalMarkers,
                          'tissue_subset_markers'=tissueMarkers,
                          'subcluster_global_markers'=subclusterGlobalMarkers), 
                     file='~/publication_pipeline/data_tables.xlsx')
