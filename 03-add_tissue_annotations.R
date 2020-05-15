######################################################
# ADD TISSUE ANNOTATIONS TO DATA
#
# GOAL: Add manual global cluster tissue annotations 
# to Seurat3 object.
#####################################################

# ================= IMPORTS ==========================

library(Seurat)
library(openxlsx)
library(here)

here <- here::here  # fix lubridate's masking of here()
source(here('utils.R'))

# ================= PARAMS ===========================

source(here('user_params.R'))

# ================= DATA =============================

s3 <- readRDS(here(paste0(saveAs,'.Rds')))
metadata <- read.xlsx(here('data_tables.xlsx'))
rownames(metadata) <- metadata$cell_ID
tissues <- factor(metadata[Cells(s3),'tissue_type'])
s3 <- AddMetaData(s3,metadata=tissues,col.name='Tissue.Type')

saveRDS(s3,here(paste0(saveAs,'.Rds')))
