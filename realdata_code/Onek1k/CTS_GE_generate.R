library(Seurat)
library(dplyr)
library(data.table)
data_dir <- './data/'
gene_info <- fread('./1k1k_gene_GRCh37.txt')

# read cov
for(ct in c('CD4effCM','Bmem','CD4TGFbStim','CD8all','CD8eff','CD8unknown','DC','MonoNC','NKact','MonoC','CD4all','BimmNaive',
           'NKmat','Plasma')){

    obj <- readRDS(sprintf('%s/ct_rds/%s.rds',data_dir,ct))


    Agg_obj <- AggregateExpression(object=obj, assay = 'RNA', group.by='individual', return.seurat = TRUE)
    colnames(Agg_obj) <- gsub("-", "_", colnames(Agg_obj))
    colnames(Agg_obj) <- gsub("g", "", colnames(Agg_obj))
    meta_data1 <- obj@meta.data %>% select(individual,pool) %>% unique()
    meta_data1 <- meta_data1[match(colnames(Agg_obj),meta_data1$individual),]
    rownames(meta_data1) <- meta_data1$individual
    Agg_obj@meta.data <- meta_data1 # %>% select(-individual)

    data <- SCTransform(object = Agg_obj, return.only.var.genes = FALSE)

    ## scTWAS object
    saveRDS(data, file =  sprintf('%s/scTransform_by_celltype_afterAgg/%s.rds',data_dir,ct))


    ## ANTWAS: TMM+log 
    library(edgeR) # for TMM
    library(sva)
    pb_mat <- GetAssayData(Agg_obj, assay = 'RNA', layer = 'counts')
    min_count <- which(rowSums(pb_mat)<1000)
    pb_mat <- pb_mat[-min_count,]
    pb_mat <- pb_mat[order(rowSums(pb_mat),decreasing=T),]

    dgList <- DGEList(pb_mat, genes=rownames(pb_mat))
    dgList <- calcNormFactors(dgList, method = "TMM")
    expr_norm = log_transform(dgList) # same as voom(dgList)$E
    donor_pool <- (Agg_obj@meta.data)$pool

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
    fwrite(PB_df, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_AN.tsv',data_dir,ct),sep='\t')
    write.table(colnames(PB_df),sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s.header',data_dir,ct),quote=F,row.names=F,col.names=F)

}

### NATWAS: first sctransform then mean aggregate pearson residual
for(ct in c('CD4effCM','Bmem','CD4TGFbStim','CD8all','CD8eff','CD8unknown','DC','MonoNC','NKact','MonoC','CD4all','BimmNaive',
           'NKmat','Plasma')){
    obj <- readRDS(sprintf('%s/ct_rds/%s.rds',data_dir,ct))

    sct_obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 1e6)
    unique_donor <- sort(unique(sct_obj@meta.data$individual))
    pr = matrix(NA, nrow = nrow(sct_obj), ncol = length(unique_donor))
    colnames(pr) = unique_donor; rownames(pr) = rownames(sct_obj)
    for(i in unique_donor){
        ind = which(sct_obj@meta.data$individual == i)
        pr[,i] = rowMeans((sct_obj[['RNA']]$data[,ind,drop=FALSE]))
    }

    pb_filter = fread(sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_AN.tsv',data_dir,ct))
    pr_inrt <- matrix(NA, nrow = nrow(pr), ncol = ncol(pr))
    for(i in 1:nrow(pr)){
        pr_inrt[i,] <- INT(pr[i,])
    }
    rownames(pr_inrt) = rownames(pr)
    colnames(pr_inrt) = colnames(pr)

    PR_df = data.frame(gene_name = rownames(pr_inrt))
    PR_df = cbind(PR_df,pr_inrt)
    PR_df = merge(PR_df,gene_info %>% dplyr::select(-ensembl_gene_id),by.x='gene_name',by.y='external_gene_name',sort=FALSE)
    PR_df = PR_df %>% relocate(gene_name, chromosome_name, start_position, end_position)
    PR_df = PR_df[match(pb_filter$gene_name,PR_df$gene_name),]
    fwrite(PR_df, sprintf('%s/scTransform_by_celltype/Expression_matrices/%s_NA.tsv',data_dir,ct),sep='\t')
}
