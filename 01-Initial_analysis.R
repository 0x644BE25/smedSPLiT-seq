######################################################
# INITIAL ANALYSIS
# 
# GOAL: Process scRNA data from the Smed regeneration 
# time-series SPLiT-seq experiment using 150 PCs.
######################################################

# ================= IMPORTS ==========================

library(here)       # directory management
library(Seurat)     # single-cell toolkit
library(future)     # parallelization
library(data.table) # fast file reading
library(Matrix)     # sparse matrices
library(ggplot2)
here <- here::here  # fix lubridate's masking of here()

source(here('utils.R'))

# ================= PARAMS ===========================

doCreateS3 <- FALSE
doSCTransform <- FALSE
doDimReduc <- TRUE
doCluster <- TRUE

source(here('user_params.R'))

# ================= CREATE SEURAT OBJECT =============

if (doCreateS3) {
  exprMatrices <- list.files(here('data/expression_matrices'))
  
  allCounts <- NULL
  for (em in exprMatrices) {
    print(em)
    counts <- fread(here(paste0('data/expression_matrices/',em)))
    smeds <- unlist(c(counts[,1]))
    counts <- Matrix(as.matrix(counts[,2:ncol(counts)]))
    rownames(counts) <- smeds
    allCounts <- cbind(allCounts,counts)
  }
  
  s3 <- CreateSeuratObject(allCounts,project='BBP')
  saveRDS(s3,here(paste0(saveAs,'.Rds')))
  rm(counts,allCounts,smeds)
  
  prpr('Seurat object created')
  
  # ================= ADD METADATA =====================
  
  # Timepoint
  day <- unlist(lapply(strsplit(Cells(s3),'-'), function(x){ return(x[2]) }), use.names=FALSE)
  s3 <- AddMetaData(s3, metadata=day, col.name='Timepoint')
  for (i in c(0:2,4,7)) {
    s3$Timepoint[which(s3$Timepoint==paste0('D',i))] <- paste0('D0',i)
  }
  s3$Timepoint <- factor(s3$Timepoint)
  
  # Treatment (rough)
  txs <- unlist(lapply(strsplit(Cells(s3),'_'), function(x){ return(x[1]) }), use.names=FALSE)
  s3 <- AddMetaData(s3, metadata=txs, col.name='Treatment')
  
  # Sample
  sample <- paste0(as.character(s3$Treatment),'_',as.character(s3$Timepoint))
  s3 <- AddMetaData(s3, metadata=factor(sample), col.name='Sample')
  
  # Make Treatment nice
  txList <- c('WT'='Un-irradiated','1250'='Sub-lethal','6K'='Lethal')
  s3$Treatment[which(s3$Treatment=='WT')] <- 'Un-irrdiated'
  s3$Treatment[which(s3$Treatment=='1250')] <- 'Sub-lethal'
  s3$Treatment[which(s3$Treatment=='6K')] <- 'Lethal'
  s3$Treatment <- factor(s3$Treatment)
  
  saveRDS(s3, here(paste0(saveAs,'.Rds')))
  prpr('Treatment, Timepoint, and Sample metadata added, RDS saved')
}
  
#================= SCTRANSFORM ======================

if (doSCTransform) {  
  s3 <- readRDS(here(paste0(saveAs,'.Rds')))
  
  library(future)     # parallelization
  options(future.globals.maxSize=(1024^3)*maxGB) 
  plan('multiprocess', workers=nCores)
  s3 <- SCTransform(s3)
  
  saveRDS(s3, here(paste0(saveAs,'.Rds')))
  prpr('SCTransform complete, RDS saved')
}
  
  # ================= DIMENSIONAL REDUCTION ============

if (doDimReduc) {
  s3 <- readRDS(here(paste0(saveAs,'.Rds')))
  plan('multiprocess', workers=nCores)
  s3 <- RunPCA(s3, npcs=nPCs)
  ggsave(ElbowPlot(s3, ndims=nPCs, reduction='pca'), file=here(paste0(saveAs,'_elbow_plot.png')))
  s3 <- RunUMAP(s3, dims=1:nPCs)
  s3 <- RunUMAP(s3, dims=1:nPCs, n.components=3, reduction.name='umap3D')
  plan('sequential')
  
  saveRDS(s3, here(paste0(saveAs,'.Rds')))
  prpr('dimensional analysis done, RDS saved')
}
# ================= CLUSTERING =======================

if (doCluster) {
  s3 <- readRDS(here(paste0(saveAs,'.Rds')))
  s3 <- FindNeighbors(s3, force.recalc=TRUE, dims=1:nPCs)
  s3 <- FindClusters(s3)
  
  saveRDS(s3, here(paste0(saveAs,'.Rds')))
  prpr('clustering done, final RDS saved')
}