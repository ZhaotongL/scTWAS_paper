library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(edgeR)
source('../../R/scTWAS_IRLS.R')
gene_info <- fread('../Onek1k/1k1k_gene_GRCh37.txt')
data_dir <- './data/ROSMAP/Gene Expression (snRNAseq - DLPFC, Experiment 2)/processed (March 2024 update)'


sum_sparse <- function(mats) {

  all_rows <- sort(unique(unlist(lapply(mats, rownames))))
  all_cols <- sort(unique(unlist(lapply(mats, colnames))))
  rmap <- setNames(seq_along(all_rows), all_rows)
  cmap <- setNames(seq_along(all_cols), all_cols)

  I <- integer(); J <- integer(); X <- numeric()
  for (m in mats) {
    tM <- as(m, "dgTMatrix")
    I <- c(I, rmap[rownames(m)][tM@i + 1L])
    J <- c(J, cmap[colnames(m)][tM@j + 1L])
    X <- c(X, tM@x)
  }
  sparseMatrix(i = I, j = J, x = X,
               dims = c(length(all_rows), length(all_cols)),
               dimnames = list(all_rows, all_cols))
}

                 
for(ct in c('Microglia', 'Astrocyte','cux2+', 'cux2-','Oligodendrocytes', 'OPCs','Inhibitory Neurons')){
    obj <- readRDS(sprintf('%s/scTransform_by_celltype_afterAgg/%s.rds',data_dir,ct))
    print(ct)
    pb_mat <- GetAssayData(obj, assay = 'RNA', layer = 'counts')
    saveRDS(pb_mat,sprintf('%s/scTransform_by_celltype_afterAgg/%s_RawCount.rds',data_dir,ct))
}

mats <- lapply(paste0(sprintf('%s/scTransform_by_celltype_afterAgg/',data_dir), c('Microglia', 'Astrocyte','cux2+', 'cux2-','Oligodendrocytes', 'OPCs','Inhibitory Neurons'),'_RawCount.rds'), readRDS)

pb_mat <- sum_sparse(mats)
min_count <- which(rowSums(pb_mat)<1000)
pb_mat <- pb_mat[-min_count,]
pb_mat <- pb_mat[order(rowSums(pb_mat),decreasing=T),]

dgList <- DGEList(pb_mat, genes=rownames(pb_mat))
dgList <- calcNormFactors(dgList, method = "TMM")
expr_norm = log_transform(dgList) # same as voom(dgList)$E

# SKIP quantile normalization
# https://support.bioconductor.org/p/77664/#77665
expr_norm_inrt <- matrix(NA, nrow = nrow(expr_norm), ncol = ncol(expr_norm))
for(i in 1:nrow(expr_norm)){
    expr_norm_inrt[i,] <- INT(expr_norm[i,])
}
rownames(expr_norm_inrt) = rownames(expr_norm)
colnames(expr_norm_inrt) = colnames(expr_norm)

PB_df = data.frame(gene_name = rownames(expr_norm_inrt))
PB_df = cbind(PB_df,expr_norm_inrt)
PB_df = merge(PB_df,gene_info %>% dplyr::select(-ensembl_gene_id),by.x='gene_name',by.y='external_gene_name',sort=FALSE)
PB_df = PB_df %>% relocate(gene_name, chromosome_name, start_position, end_position)
fwrite(PB_df, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/Bulk_ANTWAS.tsv',data_dir),sep='\t')

write.table(colnames(PB_df),sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s.header',
					    data_dir, 'Bulk'),quote=F,row.names=F,col.names=F)