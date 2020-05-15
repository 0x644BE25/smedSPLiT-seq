######################################################
# DOWNSAMPLE - PROPORTIONAL
#
# GOAL: Downsample dataset while preserving proportion 
# of cells in each cluster in each sample.
######################################################

# ================= IMPORTS ==========================

library(ggplot2)  # plotting
library(Cairo)    # better plot rendering on Linux
library(Seurat)   # scRNA-seq analysis
library(future)   # parallelization
library(dplyr)    # data-wrangling
library(here)
here <- here::here  # fix lubridate's masking of here()

source(here('utils.R'))

# ================= PARAMS ===========================

source(here('user_params.R'))
# max number of cells/sample
cellsPerSample <- 6000

# ================= INIT DATA ========================

# read in data, get max cluster int
s3 <- readRDS(here(paste0(saveAs,'.Rds')))
prpr('Data read')
maxClust <- length(levels(s3$seurat_clusters))-1


# ================= GET CELL COUNTS ==================

cellCounts <- data.frame(matrix(0, nrow=length(unique(s3$Sample)), ncol=maxClust+1))
row.names(cellCounts) <- unique(s3$Sample)
colnames(cellCounts) <- 0:maxClust
cellList <- NULL

i <- 0
# iterate over samples
for (s in unique(s3$Sample)) {
  i <- i + 1
  # get cell & cluster lists for current sample
  currSample <- Cells(s3)[which(s3$Sample==s)]
  clusters <- s3@meta.data$seurat_clusters[which(s3$Sample==s)]
  
  # if we need to down-sample this sample, do
  if (length(currSample) > cellsPerSample){
    for (c in 0:maxClust){
      n <- ceiling(cellsPerSample*sum(clusters==c)/length(clusters))
      curr <- currSample[which(clusters==c)]
      cells <- sample(curr, n, replace=FALSE)
      cellList <- c(cellList, cells)
      cellCounts[i,c+1] <- n
    }
    # else use full sample
  } else {
    cellList <- c(cellList, currSample)
    for (c in 0:maxClust){
      n <- sum(clusters==c)
      cellCounts[i,c+1] <- n
    }
  }
}

# ================= SUBSET & SAVE ====================

tinyData <- subset(s3, cells=unlist(cellList))
tinyData <- DietSeurat(tinyData, counts=FALSE, data=TRUE, assays='SCT')
saveRDS(tinyData, here(paste0(saveAs,'_',cellsPerSample,'_per_sample.Rds')))

