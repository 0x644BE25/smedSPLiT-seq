######################################################
# 150 PCs TISSUE ANALYSIS W/ CLUSTER 1
#
# GOAL: Subset data by tissue type, analyze with BBP's
# determined number of PCs, generate 2D/3D UMAPs, 
# cluster, find markers.
######################################################

# ================= IMPORTS ======================== #

library(ggplot2)  # plotting
library(Cairo)    # better plot rendering on Linux
library(Seurat)   # scRNA-seq analysis
library(future)   # parallelization
library(dplyr)    # data-wrangling
library(here)
here <- here::here  # fix lubridate's masking of here()

source(here('utils.R'))

# ================= PARAMS ===========================

makeSubsets <- FALSE
analyzeSubsets <- TRUE
getMarkers <- TRUE

source(here('user_params.R'))
options(future.globals.maxSize=((1024^3)*maxGB))
annot_file <- "planosphere_annotations.txt"

# ================= SUBSET DATA ==================== #
if (makeSubsets) {
  # iterate over tissue types, create subset, save as Rds
  
  # read in data
  s3 <- readRDS(here(paste0(saveAs,'.Rds')))
  prpr('Data read')
  tissues <- levels(s3$Tissue.Type)
  
  for (t in tissues) {
    cells <- Cells(s3)[which(s3$Tissue.Type==t)]
    print(paste(t, ':', length(cells), 'cells'))
    ss <- CreateSeuratObject(counts=s3[['RNA']]@counts[,cells])
    ss <- AddMetaData(ss, metadata=s3$Sample[cells], col.name='Sample')
    ss <- AddMetaData(ss, metadata=s3$Timepoint[cells], col.name='Timepoint')
    ss <- AddMetaData(ss, metadata=s3$Tissue.Type[cells], col.name='Tissue.Type')
    ss <- AddMetaData(ss, metadata=s3$Treatment[cells], col.name='Treatment')
    ss <- AddMetaData(ss, metadata=s3$seurat_clusters[cells], col.name='global_clusters')
    t <- gsub(" ", "_", t)
    saveRDS(ss, here(paste0(saveAs,'_',t,'.Rds')))
  }
  rm(s3)
  prpr('tissue RDS files saved')
}

# ================= ANALYZE SUBSETS ================ #
if (analyzeSubsets) {
  prpr('Now for Analysis...')
  
  # iterate over subset.Rds files, do analysis, resave
  
  # number of PCs for each tissue type
  # determined by BBP looking at first 55 PCs for each
  nPClist <- c('Cathepsin'=20, 'Epidermis'=25, 'Intestine'=20,
               'Muscle'=20, 'Non-differentiated'=20,'Parenchymal'=35, 'Pharynx'=15,
               'Protonephridia'=20, 'Neural'=40, 'Stem Cells'=15)
  
  tissues <- names(nPClist)[10]
  for (t in tissues) {
    nPCs <- nPClist[t]
    prpr(paste('Analyzing',t,'with',nPCs,'PCs'))
    t <- gsub(" ", "_", t)
    ss <- readRDS(here(paste0(saveAs,'_',t,'.Rds')))
    plan('multiprocess', workers=nCores)
    ss <- SCTransform(ss)
    saveRDS(ss, here(paste0(saveAs,'_',t,'.Rds')))
    ss <- RunPCA(ss, npcs=nPCs, verbose=FALSE)
    ss <- RunUMAP(ss, reduction='pca', dims=1:nPCs, n.components=2, verbose=FALSE)
    ss <- RunUMAP(ss, reduction='pca', dims=1:nPCs, n.components=3, verbose=FALSE, reduction.name='umap3D', reduction.key='umap3D')
    ss <- FindNeighbors(ss, dims=1:nPCs, force.recalc=TRUE, verbose=FALSE)
    ss <- FindClusters(ss, verbose=FALSE)
    p <- DimPlot(ss, reduction='umap', group.by='seurat_clusters') + ggtitle(t)
    saveRDS(ss, here(paste0(saveAs,'_',t,'.Rds')))
    prpr(paste(t,'subset analyzed + saved'))
  }
}

# ================= ANNOTATE MARKERS =============== #
if (getMarkers) {
  prpr('Finding markers...')
  
  allMarkers <- NULL
  tissues <- c('Cathepsin','Epidermis','Intestine','Muscle'
               ,'Non-differentiated','Parenchymal','Pharynx',
               'Protonephridia','Stem Cells')
  for (t in tissues) {
    prpr(paste('Finding markers for',t))
    t <- gsub(" ", "_", t)
    ss <- readRDS(here(paste0(saveAs,'_',t,'.Rds')))
    plan('multiprocess', workers=nCores)
    markers <- FindAllMarkers(ss, assay='SCT')
    markers$SmedID <- substr(row.names(markers),1,12)
    markers$tissue <- t
    markers$cluster <- paste0(t,markers$cluster)
    allMarkers <- rbind(allMarkers, markers)
  }
  
  prpr('Annotating markers')
  # get rid of misleading row names
  row.names(allMarkers) <- 1:nrow(allMarkers)
  saveRDS(allMarkers, here(paste0(saveAs,'_tissue_subset_markers.Rds')))
  
  # import annotations, keep only SmedID + names, remove duplicates
  annotations <- distinct(read.delim(here(annot_file), sep=',', header=TRUE)[,1:2])

  # rename columns
  colnames(annotations) <- c('SmedID','annotations')
  # remove annotations that are just SmedIDs
  annotations <- annotations[which(!startsWith(as.character(annotations$annotations), 'SMED3')),]
  # merge 'em! only keep ones that appear in the list of markers
  annot_markers <- merge(x=allMarkers, y=annotations, by='SmedID', all.x=TRUE, all.y=FALSE)
  # save the file!
  saveRDS(annot_markers, here(paste0(saveAs,'_tissue_subset_markers.Rds')))
  write.csv(annot_markers, file=here(paste0(saveAs,'_tissue_subset_markers.csv')), row.names=FALSE)

  # how much annotation did we get?
  prpr(paste0(100*sum(!is.na(annot_markers$annotations))/length(annot_markers$annotations),"% of markers annotated."))
}
