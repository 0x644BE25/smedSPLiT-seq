######################################################
# PREP HEATMAP DATA
#
# GOAL: Create mean expression data at the tissue,
# global cluster, and tissue subcluster level for the
# heatmap ShinyApp.
######################################################

# ================= IMPORTS ==========================

library(Seurat)
library(feather)
library(lubridate)
library(gtools)

# ================= PARAMS ===========================

setwd('~/ShinyApps/heatmaps/')

tissues <- list('Wildtype',
                'Cathepsin', 
                'Epidermis',
                'Intestine',
                'Muscle',
                'Neural',
                'Parenchymal',
                'Protonephridia',
                'Pharynx',
                'Stem_cells',
                'Un-Annotated')

fileList <- paste0('~/RDS/tissue-level/',tissues,'_only.Rds')
names(fileList) <- tissues

days <- c('D00','D01','D02','D04','D07','D10','D14')
samples <- paste0(rep(c('WT','1250','6K'),each=7),'_',days)
tx <- c('Lethal','Sub-lethal','WT')

# ================= METHODS ==========================

# informative progress printouts
time <- Sys.time()
prpr <- function(msg){
  t <- Sys.time()
  print(msg)
  deltaT <- round(time_length(difftime(t,time),unit='minute'),2)
  print(paste(' ',t))
  print(paste0('  (',deltaT,' mins total)'))
}

# generate cluster/subcluster level data
getMeanData <- function(data, unit, name) {
  means <- NULL
  units <- sort(unique(as.character(data@meta.data[,unit])))
  for (i in units){
    print(i)
    curr <- data[['SCT']]@data[,which(data@meta.data[,unit]==i)]
    midpoint <- ceiling(nrow(curr)/2)
    top <- rowMeans(as.matrix(curr[1:midpoint,]))
    bottom <- rowMeans(as.matrix(curr[(midpoint+1):nrow(curr),]))
    means <- rbind(means, c(top, bottom))
  }
  if (length(units)==1){
    means <- data.frame(means)
  }
  
  rownames(means) <- units
  if (unit=='sample' & !name=='WT') {
    means <- means[samples,]
  } else if (unit=='sample' & name=='WT') {
    means <- means[samples[1:7],]
  } else if (unit=='treatment') {
    means <- means[tx,]
  } else if (endsWith(unit, 'clusters')){
    means <- means[mixedsort(rownames(means)),]
  }
  saveRDS(means, paste0('~/RDS/heatmap_data/',name,'_',unit,'.Rds'))
  
  return(TRUE)
}

# ================= GENERATE DATA ====================

data <- readRDS('~/RDS/tissue-level/Wildtype_only.Rds')

for (u in c('seurat_clusters','subclusters','sample', 'Tissue.Type','Timepoint')) {
  if(getMeanData(data, u, 'WT')) {
    prpr(paste('WT',u,'data generated'))
  }
}

for (t in tissues[2:length(tissues)]) {
  data <- readRDS(fileList[[t]])
  data$global_clusters <- droplevels(data$global_clusters)
  prpr(paste0(t,' data loaded...'))
  for (u in c('seurat_clusters','global_clusters','sample','treatment','Timepoint')) {
    if (getMeanData(data, u, t)){
      prpr(paste(t, u,'data generated'))
    }
  }
}
prpr('All done!')

# make gene lists for tissue subsets
for (t in c('ALL',tissues)) {
  print(t)
  if (t=='Wildtype') {
    t <- 'WT'
  }
  data <- readRDS(paste0('~/RDS/heatmap_data/',t,'_Timepoint.Rds'))
  print(head(colnames(data)))
  write_feather(data.frame('smeds'=colnames(data)), paste0('~/ShinyApps/heatmaps/smedList_',t,'.fa'))
}
