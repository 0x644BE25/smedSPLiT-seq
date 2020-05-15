######################################################
# FIND MARKERS - PARALLELIZED
#
# GOAL: Find markers for each global cluster, annotate
# with Planosphere's Rosetta Stone transcript mapper
# data, and save as TSV file.
######################################################

# ================= IMPORTS ======================== #

library(Seurat)  # single-cell analysis toolkit
library(future)  # parallelization
library(dplyr)   # data-wrangling
library(plyr)    # data-wrangling
library(gtools)  # mixedsort
library(here)    # directory management
here <- here::here  # fix lubridate's masking of here()

source(here('utils.R'))

# ================= PARAMS ===========================

source(here('user_params.R'))
annotFile <- 'planosphere_annotations.txt'

# ================= DATA =============================

s3 <- readRDS(here(paste0(saveAs,'.Rds')))
s3@active.ident <- s3$seurat_clusters
s3@active.assay <- 'SCT'
prpr('data loaded!')

plan('multiprocess', workers=nCores)
allMarkers <- NULL
clusters <- mixedsort(as.character(unique(s3$seurat_clusters)))
for (cluster in clusters) {
  markers <- FindMarkers(s3, ident.1=cluster, verbose=FALSE)
  markers$cluster <- cluster
  markers$smedID <- rownames(markers)
  allMarkers <- rbind(allMarkers,markers)
  prpr(paste('cluster ',cluster,': ',nrow(markers),' markers found'))
  saveRDS(allMarkers, here(paste0(saveAs,'_global_cluster_markers.Rds')))
}
plan('sequential')

row.names(allMarkers) <- 1:nrow(allMarkers)
saveRDS(allMarkers, here(paste0(saveAs,'_global_cluster_markers.Rds')))
prpr('Done finding markers! Annotating.')

# import annotations and filter
annotations <- distinct(read.csv(here(annotFile), sep=',', header=TRUE)[,1:2])
colnames(annotations) <- c('smedID','annotations')
annotations <- annotations[which(!startsWith(as.character(annotations$annotations), 'SMED3')),]
annot_markers <- merge(x=allMarkers, y=annotations, by='smedID', all.x=TRUE, all.y=FALSE)
# save the file!
write.table(annot_markers, file=here(paste0(saveAs,'_global_cluster_markers.tsv')), quote=FALSE, sep='\t', col.names=NA)
saveRDS(annot_markers, here(paste0(saveAs,'_global_cluster_markers.Rds')))
prpr('Markers saved')

# how much annotation did we get?
prpr(paste0(100*sum(!is.na(annot_markers$annotations))/length(annot_markers$annotations),"% of markers annotated."))
