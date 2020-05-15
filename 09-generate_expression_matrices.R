library(data.table)

full <- readRDS('~/RDS/FULL_150PCs_2020-02-17.Rds')

for (s in unique(full$sample)){
  cells <- Cells(full)[which(full$sample==s)]
  sub <- full[['RNA']]@counts[,cells]
  print(paste0(s,': ', length(cells), ' cells'))
  sub <- data.table(as.matrix(sub), keep.rownames=TRUE)
  fwrite(x=sub, file=paste0('~/publication_pipeline/data/expression_matrices/',s,'_expression_matrix.csv'), row.names=FALSE, col.names=TRUE)
}

