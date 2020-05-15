######################################################
# TISSUE LABEL TRANSFER
#
# GOAL: Use Seurat 3's "label transfer" capability to
# project Fincher's cluster IDs onto our SPLiT-seq
# data and vice-versa.
######################################################

# ================= IMPORTS ==========================

library(gdata)
library(Seurat)
library(plotly)
library(openxlsx)
library(lubridate)
library(viridis)
library(ggplot2)
library(Cairo)
library(magrittr)
source('~/ShinyApps/colors.R')

# ================= PARAMS ===========================

saveAs <- '~/Analysis/For_BBP_conf/paper/label_transfer_'

tissues <- c('Cathepsin','Epidermis','Intestine','Muscle','Neural',
             'Parenchymal','Pharynx','Protonephridia','Stem_Cells', 
             'Non-differentiated')

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
umiCoeff <- 1/(10^4)
ptSize <- 0.1

# ================= METHODS ==========================

time <- Sys.time()

prpr <- function(msg){
  t <- Sys.time()
  cat(paste0(msg,'\n'))
  deltaT <- round(time_length(difftime(t,time),unit='minute'),2)
  cat(paste(' ',t,'\n'))
  cat(paste0('  (',deltaT,' mins total)\n\n'))
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# ================= INIT DATA ========================

fin <- readRDS('~/RDS/Fincher_smed_aligned-2020-02-11.Rds')
bbp <- readRDS('~/RDS/FULL_150PCs_2020-02-17.Rds')
prpr('data loaded')

# ================= LABEL TRANSFER ===================

rds.list <- list(fin, bbp)

reference.list <- rds.list["fin"]
smed.anchors <- FindIntegrationAnchors(object.list=rds.list, dims=1:30)
smed.integrated <- IntegrateData(anchorset=smed.anchors, dims=1:30)
prpr('integrated')

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(smed.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
smed.integrated <- ScaleData(smed.integrated, verbose=FALSE)
prpr('scaled')
smed.integrated <- RunPCA(smed.integrated, npcs=30, verbose=FALSE)
prpr('PCA run')
smed.integrated <- RunUMAP(smed.integrated, reduction="pca", dims=1:30)
prpr('UMAPed')
smed.integrated$experiment <- NA
smed.integrated$experiment[which(startsWith(smed.integrated$orig.ident, 'SRR'))] <- 'Fincher'
smed.integrated$experiment[which(smed.integrated$orig.ident=='SeuratProject')] <- 'BBP'
bbp <- subset(smed.integrated, cells=Cells(bbp))
fin <- subset(smed.integrated, cells=Cells(fin))
saveRDS(smed.integrated,'~/RDS/integrated_2020-02-21.Rds')
prpr('integrated RDS saved')
rm(smed.integrated)

smed.anchors <- FindTransferAnchors(reference=fin, query=bbp, dims=1:30)
prpr('found anchors')
predictions <- TransferData(anchorset=smed.anchors, refdata=fin$fincher.tissue, dims=1:30)
bbp <- AddMetaData(bbp, metadata=predictions)

bbp$prediction.match <- bbp$predicted.id==bbp$Tissue.Type
table(bbp$prediction.match)
saveRDS(bbp, '~/RDS/Full_with_label_transfer.Rds')
prpr('AAAAALLLLLL DONE.')

# ================= ADD RESULTS ======================

bbp <- readRDS('~/RDS/Full_with_label_transfer.Rds')
predictions <- factor(bbp$predicted.id, levels=levels(bbp$Tissue.Type))
bbp <- readRDS('~/RDS/FULL_150PCs_2020-02-17.Rds')
bbp <- AddMetaData(bbp, metadata=predictions, col.name='label.transfer')
saveRDS(bbp, '~/RDS/Full_2020-02-24.Rds')

# ================= PLOTS ============================

umap <- Embeddings(bbp, reduction='umap')
curr <- data.frame(x=umap[,1], 
                   y=umap[,2], 
                   nUMI=bbp$nUMI, 
                   cluster=bbp$seurat_clusters,
                   manual=bbp$Tissue.Type,
                   label=bbp$label.transfer,
                   agree=bbp$agreement,
                   treatment=bbp$treatment,
                   time=bbp$Timepoint,
                   stringsAsFactors=FALSE)

curr$nUMI <- as.numeric(curr$nUMI)
curr$x <- as.numeric(curr$x)
curr$y <- as.numeric(curr$y)

# manual annotation

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=curr$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=unan$x, y=unan$y, color=unan$manual), size=unan$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=rest$x,y=rest$y, color=rest$manual), size=rest$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=musc$x,y=musc$y, color=musc$manual),size=musc$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:5,12,6:9)])) +
  theme_void() +
  ggtitle('Manual annotation')
ggsave(filename=paste0(saveAs,'manual_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'manual_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=ptSize) +
  geom_point(mapping=aes(x=unan$x, y=unan$y, color=unan$manual), size=ptSize) +
  geom_point(mapping=aes(x=rest$x,y=rest$y, color=rest$manual), size=ptSize) +
  geom_point(mapping=aes(x=musc$x,y=musc$y, color=musc$manual),size=ptSize) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:5,12,6:9)])) +
  theme_void() +
  ggtitle('Manual annotation')
ggsave(filename=paste0(saveAs,'manual.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'manual_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# label transfer

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=curr$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=unan$x, y=unan$y, color=unan$label), size=unan$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=rest$x,y=rest$y, color=rest$label), size=rest$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=musc$x,y=musc$y, color=musc$label),size=musc$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Label transfer')
ggsave(filename=paste0(saveAs,'prediction_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'prediction_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=ptSize) +
  geom_point(mapping=aes(x=unan$x, y=unan$y, color=unan$label), size=ptSize) +
  geom_point(mapping=aes(x=rest$x,y=rest$y, color=rest$label), size=ptSize) +
  geom_point(mapping=aes(x=musc$x,y=musc$y, color=musc$label), size=ptSize) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Label transfer')
ggsave(filename=paste0(saveAs,'prediction.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'prediction_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# agreement

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=curr$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$agree), size=curr$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1,8)])) +
  theme_void() +
  ggtitle('Agreement')
ggsave(filename=paste0(saveAs,'agreement_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'agreement_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=ptSize) +
  geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$agree), size=ptSize) +
  scale_color_manual(values=c('transparent',lab_colors[c(1,8)])) +
  theme_void() +
  ggtitle('Agreement')
ggsave(filename=paste0(saveAs,'agreement.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'agreement_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# disagreement manual colors

currD <- curr[which(curr$agree %in% c(FALSE,NA)),]
muscD <- curr[which((curr$agree %in% c(FALSE,NA)) & (curr$manual=='Muscle')),]
unanD <- curr[which((curr$agree %in% c(FALSE,NA)) & (curr$manual=='Non-differentiated')),]
restD <- curr[which((curr$agree %in% c(FALSE,NA)) & (!curr$manual %in% c('Muscle','Non-differentiated'))),]

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=currD$x,y=currD$y, color='#FFFFFF'), size=currD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=unanD$x, y=unanD$y, color=unanD$manual), size=unanD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=restD$x,y=restD$y, color=restD$manual), size=restD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=muscD$x,y=muscD$y, color=muscD$manual),size=muscD$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:5,12,6:9)])) +
  theme_void() +
  ggtitle('Manual annotation (disagreements)')
ggsave(filename=paste0(saveAs,'disagreement_manual_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'disagreement_manual_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=currD$x,y=currD$y, color='#FFFFFF'), size=ptSize) +
  geom_point(mapping=aes(x=unanD$x, y=unanD$y, color=unanD$manual), size=ptSize) +
  geom_point(mapping=aes(x=restD$x,y=restD$y, color=restD$manual), size=ptSize) +
  geom_point(mapping=aes(x=muscD$x,y=muscD$y, color=muscD$manual),size=ptSize) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:5,12,6:9)])) +
  theme_void() +
  ggtitle('Manual annotation (disagreements)')
ggsave(filename=paste0(saveAs,'disagreement_manual.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'disagreement_manual_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# disagreement label transfer colors

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=currD$x,y=currD$y, color='#FFFFFF'), size=currD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=unanD$x, y=unanD$y, color=unanD$label), size=unanD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=restD$x,y=restD$y, color=restD$label), size=restD$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=muscD$x,y=muscD$y, color=muscD$label),size=muscD$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Label transfer (disagreements)')
ggsave(filename=paste0(saveAs,'disagreement_label_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'disagreement_label_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=currD$x,y=currD$y, color='#FFFFFF'), size=ptSize) +
  geom_point(mapping=aes(x=unanD$x, y=unanD$y, color=unanD$label), size=ptSize) +
  geom_point(mapping=aes(x=restD$x,y=restD$y, color=restD$label), size=ptSize) +
  geom_point(mapping=aes(x=muscD$x,y=muscD$y, color=muscD$label),size=ptSize) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Label transfer (disagreements)')
ggsave(filename=paste0(saveAs,'disagreement_label.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'disagreement_label_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# ================= QUANTITATIVE ANALYSIS ============

gc <- mixedsort(levels(bbp$seurat_clusters))
sc <- mixedsort(levels(bbp$subclusters))
ti <- mixedsort(levels(bbp$Tissue.Type))
lt <- mixedsort(levels(drop.levels(bbp$label.transfer)))

# global clusters
gres <- NULL
for (g in gc){
  pred <- bbp$label.transfer[which(bbp$seurat_clusters==g)]
  pred <- table(pred)/length(pred)
  ma <- as.character(unique(bbp$Tissue.Type[which(bbp$seurat_clusters==g)]))
  gres <- rbind(gres,c(paste0(ma,'-',g),c(pred)))
}
rownames(gres) <- gres[,1]
gres <- t(gres[mixedsort(rownames(gres)),lt])

global <- plot_ly(z=as.matrix(gres), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    xaxis=list(
      ticktext=colnames(gres),
      tickvals=0:ncol(gres),
      title='global cluster'
    ),
    yaxis=list(
      ticktext=rownames(gres),
      tickvals=0:nrow(gres),
      title='predicted tissue'),
    autosize=TRUE)
global

# subclusters
sres <- NULL
for (s in sc){
  pred <- bbp$label.transfer[which(bbp$subclusters==s)]
  pred <- table(pred)/length(pred)
  sres <- rbind(sres,c(pred))
}
rownames(sres) <- sc
sres <- t(sres[mixedsort(sc),lt])

sub <- plot_ly(z=as.matrix(sres), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    xaxis=list(
      ticktext=colnames(sres),
      tickvals=0:ncol(sres),
      title='subcluster'
    ),
    yaxis=list(
      ticktext=rownames(sres),
      tickvals=0:nrow(sres),
      title='predicted tissue'),
    autosize=TRUE)
sub

# tissues
tres <- NULL
for (t in ti){
  pred <- bbp$label.transfer[which(bbp$Tissue.Type==t)]
  pred <- table(pred)/length(pred)
  tres <- rbind(tres,c(pred))
}
rownames(tres) <- ti
tres <- t(tres[,lt])

tis <- plot_ly(z=as.matrix(tres), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    xaxis=list(
      ticktext=colnames(tres),
      tickvals=0:ncol(tres),
      title='manual annotation'
    ),
    yaxis=list(
      ticktext=rownames(tres),
      tickvals=0:nrow(tres),
      title='predicted tissue'),
    autosize=TRUE)
tis

# save as multi-sheed xlsx
gres <- cbind(rownames(gres),gres)
colnames(gres) <- c('prediction',colnames(gres)[2:ncol(gres)])
sres <- cbind(rownames(sres),sres)
colnames(sres) <- c('prediction',colnames(sres)[2:ncol(sres)])
tres <- cbind(rownames(tres),tres)
colnames(tres) <- c('prediction',colnames(tres)[2:ncol(tres)])

openxlsx::write.xlsx(x=list(gres,sres,tres), 
                     file='~/Analysis/cluster_correlation/label_transfer_results.xlsx', 
                     sheetName=c('global clusters','subclusters','tissues'),
                     )

# ================= MODE ANALYSIS ====================

modes <- NULL

for(t in ti){
  label.transfer <- mode(as.character(bbp$label.transfer[which(bbp$Tissue.Type==t)]))
  modes <- rbind(modes, c('tissue',t, t, label.transfer))
}

for (g in gc){
  manual.annot <- mode(as.character(bbp$Tissue.Type[which(bbp$seurat_clusters==g)]))
  label.transfer <- mode(as.character(bbp$label.transfer[which(bbp$seurat_clusters==g)]))
  modes <- rbind(modes, c('global',g, manual.annot, label.transfer))
}

for (s in sc){
  manual.annot <- mode(as.character(bbp$Tissue.Type[which(bbp$subclusters==s)]))
  label.transfer <- mode(as.character(bbp$label.transfer[which(bbp$subclusters==s)]))
  modes <- rbind(modes, c('sub',s, manual.annot, label.transfer))
}

colnames(modes) <- c('type','group','manual.annotation','label.transfer')
modes <- data.frame(modes)
modes$match <- as.character(modes$manual.annotation)==as.character(modes$label.transfer)
modes$match[which(modes$manual.annotation=='Non-differentiated')] <- NA

openxlsx::write.xlsx(x=modes, file='~/Analysis/cluster_correlation/label_transfer_summary.xlsx')

# ================= CLUSTER 1 PLOTTING ===============

m <- bbp@meta.data
bbp <- readRDS('~/RDS/tissue-level/Non-differentiated_only.Rds')
pred <- m[Cells(bbp),'label.transfer']
bbp <- AddMetaData(bbp, metadata=pred, col.name='label.transfer')

umap <- Embeddings(bbp, reduction='umap')
curr <- data.frame(x=umap[,1], 
                   y=umap[,2], 
                   nUMI=bbp$nUMI, 
                   label=bbp$label.transfer,
                   treatment=bbp$treatment,
                   time=bbp$Timepoint,
                   stringsAsFactors=FALSE)

curr$nUMI <- as.numeric(curr$nUMI)
curr$x <- as.numeric(curr$x)
curr$y <- as.numeric(curr$y)

umiCoeff <- min(4/(nrow(curr)), 1/2000)

# manual annotation

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=curr$nUMI*umiCoeff) +
  geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$label), size=unan$nUMI*umiCoeff) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Cluster 1 label transfer')
ggsave(filename=paste0(saveAs,'C1_manual_nUMI.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'C1_manual_nUMI_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

p <- ggplot(data=NULL) + 
  geom_point(mapping=aes(x=curr$x,y=curr$y, color='#FFFFFF'), size=.5) +
  geom_point(mapping=aes(x=curr$x, y=curr$y, color=curr$label), size=.5) +
  scale_color_manual(values=c('transparent',lab_colors[c(1:9)])) +
  theme_void() +
  ggtitle('Cluster 1 label transfer')
ggsave(filename=paste0(saveAs,'C1_manual.png'), plot=p, type='cairo-png', width=wideWidth, height=squareWidth)
p <- p + theme(legend.position='none')
ggsave(filename=paste0(saveAs,'C1_manual_noLegend.png'), plot=p, type='cairo-png', width=squareWidth, height=squareWidth)

# ================= BREAKOUT PLOTS ===================

bbp <- readRDS('~/RDS/Full_with_label_transfer.Rds')

sa <- mixedsort(unique(bbp$sample))
gc <- mixedsort(unique(bbp$seurat_clusters))
bySample <- NULL
for (s in sa){
  currS <- NULL
  for (g in gc) {
    currS <- c(currS, mean(bbp$agreement[which(bbp$sample==s & bbp$seurat_clusters==g)]))
  }
  bySample <- rbind(bySample,currS)
}
rownames(bySample) <- sa
colnames(bySample) <- gc

byS <- plot_ly(z=as.matrix(bySample), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    title='proportion of cluster manual annotation matching label transfer',
    xaxis=list(
      ticktext=colnames(bySample),
      tickvals=0:ncol(bySample),
      title='global cluster'
    ),
    yaxis=list(
      ticktext=rownames(bySample),
      tickvals=0:nrow(bySample),
      title='sample'),
    autosize=TRUE)
byS


tx <- mixedsort(unique(bbp$treatment))
byTx <- NULL
for (t in tx){
  currS <- NULL
  for (g in gc) {
    currS <- c(currS, mean(bbp$agreement[which(bbp$treatment==t & bbp$seurat_clusters==g)]))
  }
  byTx <- rbind(byTx,currS)
}
rownames(byTx) <- tx
colnames(byTx) <- gc

byTr <- plot_ly(z=as.matrix(byTx), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    title='proportion of cluster manual annotation matching label transfer',
    xaxis=list(
      ticktext=colnames(byTx),
      tickvals=0:ncol(byTx),
      title='global cluster'
    ),
    yaxis=list(
      ticktext=rownames(byTx),
      tickvals=0:nrow(byTx),
      title='treatment'),
    autosize=TRUE)
byTr

tp <- mixedsort(unique(bbp$Timepoint))
byTp <- NULL
for (t in tp){
  currS <- NULL
  for (g in gc) {
    currS <- c(currS, mean(bbp$agreement[which(bbp$Timepoint==t & bbp$seurat_clusters==g)]))
  }
  byTp <- rbind(byTp,currS)
}
rownames(byTp) <- tp
colnames(byTp) <- gc

byTi <- plot_ly(z=as.matrix(byTp), type='heatmap', colors=hm_colors(256)) %>%
  layout(
    title='proportion of cluster manual annotation matching label transfer',
    xaxis=list(
      ticktext=colnames(byTp),
      tickvals=0:ncol(byTp),
      title='global cluster'
    ),
    yaxis=list(
      ticktext=rownames(byTp),
      tickvals=0:nrow(byTp),
      title='timepoint'),
    autosize=TRUE)
byTi

byGc <- gcTis <- NULL
for (g in gc){
  byGc <- c(byGc, mean(bbp$agreement[which(bbp$seurat_clusters==g)]))
  gcTis <- c(gcTis, unique(as.character(bbp$Tissue.Type[which(bbp$seurat_clusters==g)])))
}
names(byGc) <- gc

all <- rbind(bySample,byTp,byTx,byGc)
rownames(all) <- c(rownames(all)[1:(nrow(all)-1)],'ALL')
colnames(all) <- paste0(gcTis,'-',gc)
all <- all[,mixedsort(colnames(all))]

byAll <- plot_ly(z=as.matrix(all), type='heatmap', colors=lt_colors(256)) %>%
  layout(
    title='proportion of cluster manual annotation matching label transfer',
    xaxis=list(
      ticktext=colnames(all),
      tickvals=0:ncol(all),
      title='global cluster'
    ),
    yaxis=list(
      ticktext=rownames(all),
      tickvals=0:nrow(all),
      title='group'),
    autosize=TRUE)
byAll

