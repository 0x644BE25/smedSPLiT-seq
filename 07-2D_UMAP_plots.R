######################################################
# 2D UMAPS FOR PAPER
#
# GOAL: Create 2D UMAP plots with uniform sclaes and
# axes as well as nUMI dot-size legends. Plot cluster
# 1 on the bottom and muscle on the top. Do global
# plot with tissue type and split by tx color by time.
# Do tissue subset plots colored by global cluster and
# split by tx color by time.
######################################################

# ================= IMPORTS ==========================

library(viridis)
library(ggplot2)
library(Cairo)
library(Seurat)
library(magrittr)
library(gdata)
source('~/ShinyApps/colors.R')

# ================= PARAMS ===========================

setwd('~/Analysis/For_BBP_conf/paper')
tissues <- c('Cathepsin','Epidermis','Intestine','Muscle','Neural',
             'Parenchymal','Pharynx','Protonephridia','Stem_Cells', 
             'Un-Annotated')
bigSets <- list('Full'='~/RDS/FULL_150PCs_2019-12-31.Rds',
                '6K'='~/RDS/6K_150PCs_2020-01-02.Rds')

saveAs <- '~/Analysis/For_BBP_conf/paper/'

# custom grid style
noax <- list(
  title="",
  zeroline=FALSE,
  showline=FALSE,
  showticklabels=FALSE,
  showgrid=FALSE,
  showspikes=FALSE
)

wideWidth <- 16
squareWidth <- 12

clust_col <- lab_colors[c(1,3,5,7,9,11,2,4,6,8,10,12)]

# ================= GLOBAL DATASETS ==================
for (na in names(bigSets)){
  print(na)
  s3 <- readRDS(bigSets[[na]])
  
  umap <- Embeddings(s3, reduction='umap')
  curr <- data.frame(x=umap[,1], 
                     y=umap[,2], 
                     nUMI=s3$nUMI, 
                     cluster=s3$seurat_clusters,
                     tissue=s3$Tissue.Type,
                     treatment=s3$treatment,
                     time=s3$Timepoint,
                     stringsAsFactors=FALSE)
  n <- 0
  for (i in c(1,2,10,20,30)) {
    x <- max(umap[,1])
    y <- max(umap[,2])+n/2
    size <- i*500
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 'Muscle','WT','D00')
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 'Muscle','Sub-lethal','D00')
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 'Muscle','Lethal','D00')
    n <- n+1
  }
  curr$nUMI <- as.numeric(curr$nUMI)
  curr$x <- as.numeric(curr$x)
  curr$y <- as.numeric(curr$y)
  
  umiCoeff <- 1/(10^4)
  
  # ================= GLOBAL CLUSTERS ==================
  
  p <- ggplot(data=curr, mapping=aes(x=x,y=y, color=cluster)) + 
    geom_point(size=curr$nUMI*umiCoeff) + theme_void() + scale_color_manual(values=c(clust_col, eighty_nine)) +
    ggtitle('All global clusters')
  ggsave(filename=paste0(saveAs,na,'_global_clusters.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_global_clusters_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # ================= TISSUES ==========================
  
  unan <- curr[which(curr$tissue=='Un-Annotated'),]
  musc <- curr[which(curr$tissue=='Muscle'),]
  rest <- curr[which(!curr$tissue %in% c('Un-Annotated','Muscle')),]
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=curr$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=unan$x, y=unan$y, color=unan$tissue), size=unan$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=rest$x,y=rest$y, color=rest$tissue), size=rest$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=musc$x,y=musc$y, color=musc$tissue),size=musc$nUMI*umiCoeff) +
    scale_color_manual(values=c('transparent',lab_colors[c(1:9,12)])) +
    theme_void() +
    ggtitle('Tissues')
  ggsave(filename=paste0(saveAs,na,'_tissues.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_tissues_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # ================= TIMEPOINT ========================
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$time), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=lab_colors) +
    theme_void() +
    ggtitle('Timepoint')
  ggsave(filename=paste0(saveAs,na,'_timepoint.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_timepoint_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # ================= TREATMENT ========================
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$treatment), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=lab_colors[c(1,4,5)]) +
    theme_void() +
    ggtitle('Treatment')
  ggsave(filename=paste0(saveAs,na,'_treatment.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_treatment_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # =============== SPLIT BY TREATMENT ===============
  
  wt <- curr[which(curr$treatment=='WT'),]
  sl <- curr[which(curr$treatment=='Sub-lethal'),]
  le <- curr[which(curr$treatment=='Lethal'),]
  
  n <- 0
  for (i in c(1,2,10,20,30)) {
    x <- max(umap[,1])
    y <- max(umap[,2])+n/2
    size <- i*500
    sl[nrow(sl)+1,] <- c(x,y, i*500, curr$cluster[1], 'Muscle','Sub-lethal','D00')
    le[nrow(le)+1,] <- c(x,y, i*500, curr$cluster[1], 'Muscle','Lethal','D00')
    n <- n+1
  }
  sl$nUMI <- as.numeric(sl$nUMI)
  sl$x <- as.numeric(sl$x)
  sl$y <- as.numeric(sl$y)
  le$nUMI <- as.numeric(le$nUMI)
  le$x <- as.numeric(le$x)
  le$y <- as.numeric(le$y)
  
  umiCoeff <- 1.5*umiCoeff
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x,y=curr$y, color='hmmm'), size=curr$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=wt$x, y=wt$y, color=wt$time), size=wt$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('WT Timepoint')
  ggsave(filename=paste0(saveAs,na,'_WT_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_WT_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x,y=curr$y, color='hmmm'), size=curr$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=sl$x, y=sl$y, color=sl$time), size=sl$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('Sub-lethal Timepoint')
  ggsave(filename=paste0(saveAs,na,'_Sub-lethal_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_Sub-lethal_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x,y=curr$y, color='hmmm'), size=curr$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=le$x, y=le$y, color=le$time), size=le$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('Lethal Timepoint')
  ggsave(filename=paste0(saveAs,na,'_Lethal_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,na,'_Lethal_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
}

# ================= TISSUE SUBSETS ===================

for (t in tissues){
  print(t)
  s3 <- readRDS(paste0('~/RDS/tissue-level/',t,'_only.Rds'))
  umap <- Embeddings(s3, reduction='umap')
  curr <- data.frame(x=umap[,1], 
                     y=umap[,2], 
                     nUMI=s3$nUMI, 
                     cluster=drop.levels(s3$global_clusters),
                     subcluster=s3$seurat_clusters,
                     treatment=s3$treatment,
                     time=s3$Timepoint,
                     stringsAsFactors=FALSE)
  n <- 0
  for (i in c(1,2,10,20,30)) {
    x <- max(umap[,1])
    y <- max(umap[,2])+n/2
    size <- i*500
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 1, 'WT', 'D00')
    n <- n+1
  }
  curr$nUMI <- as.numeric(curr$nUMI)
  curr$x <- as.numeric(curr$x)
  curr$y <- as.numeric(curr$y)
  
  # =============== GLOBAL CLUSTERS ==================
  
  umiCoeff <- min(4/(nrow(curr)), 1/2000)
  p <- ggplot(data=curr, mapping=aes(x=x,y=y, color=cluster)) + 
    geom_point(size=curr$nUMI*umiCoeff) + theme_void() + scale_color_manual(values=c(clust_col, eighty_nine)) +
    ggtitle(paste(t,'global clusters'))
  p
  ggsave(filename=paste0(saveAs,t,'_global_clusters.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_global_clusters_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # =============== SUBCLUSTERS ======================
  
  p <- ggplot(data=curr, mapping=aes(x=x,y=y, color=subcluster)) + 
    geom_point(size=curr$nUMI*umiCoeff) + theme_void() + scale_color_manual(values=c(clust_col, eighty_nine)) +
    ggtitle(paste(t,'sub-clusters'))
  p
  ggsave(filename=paste0(saveAs,t,'_sub-clusters.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_sub-clusters_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # =============== TIMEPOINT ========================
  
  p <- ggplot(data=NULL) + 
    geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$time), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=lab_colors) +
    theme_void() +
    ggtitle('Timepoint')
  ggsave(filename=paste0(saveAs,t,'_timepoint.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_timepoint_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  # =============== SPLIT BY TREATMENT ===============
  
  n <- 0
  for (i in c(1,2,10,20,30)) {
    x <- max(umap[,1])
    y <- max(umap[,2])+n/2
    size <- i*500
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 1, 'Lethal', 'D00')
    curr[nrow(curr)+1,] <- c(x,y, i*500, curr$cluster[1], 1, 'Sub-lethal', 'D00')
    n <- n+1
  }
  curr$nUMI <- as.numeric(curr$nUMI)
  curr$x <- as.numeric(curr$x)
  curr$y <- as.numeric(curr$y)
  
  wt <- curr[which(curr$treatment=='WT'),]
  sl <- curr[which(curr$treatment=='Sub-lethal'),]
  le <- curr[which(curr$treatment=='Lethal'),]
  
  umiCoeff <- umiCoeff*2
  p <- ggplot(data=NULL) +
    geom_point(mapping=aes(x=wt$x, y=wt$y, color=wt$time), size=wt$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=curr$x, y=curr$y, color='hmm'), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('WT')
  ggsave(filename=paste0(saveAs,t,'_WT_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_WT_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  p <- ggplot(data=NULL) +
    geom_point(mapping=aes(x=sl$x, y=sl$y, color=sl$time), size=sl$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=curr$x, y=curr$y, color='hmm'), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('Sub-lethal')
  ggsave(filename=paste0(saveAs,t,'_Sub-lethal_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_Sub-lethal_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
  
  p <- ggplot(data=NULL) +
    geom_point(mapping=aes(x=le$x, y=le$y, color=le$time), size=le$nUMI*umiCoeff) +
    geom_point(mapping=aes(x=curr$x, y=curr$y, color='hmm'), size=curr$nUMI*umiCoeff) +
    scale_color_manual(values=c(lab_colors[1:7],'transparent')) +
    theme_void() +
    ggtitle('Lethal')
  ggsave(filename=paste0(saveAs,t,'_Lethal_by_time.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
  p <- p + theme(legend.position='none')
  ggsave(filename=paste0(saveAs,t,'_Lethal_by_time_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)
}
