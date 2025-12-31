library(Seurat)
library(SeuratDisk) # new, for h5Seurat objects
library(dplyr)
library(data.table)
library(Matrix)
require(SeqArray)

source('../../R/scTWAS_IRLS.R')
gene_info <- fread('../Onek1k/1k1k_gene_GRCh37.txt')

suppressMessages(library("optparse"))
option_list = list(
  make_option("--run_subtype", action="store", default="FALSE", type="character",
	      help="whether to analyze cell subtypes"),
  make_option("--ct", action="store", default="Microglia", type="character",
	      help="which cell type"),
  make_option("--match_with_ct", action="store", default="FALSE", type="character",
	      help="whether to match the gene evauated in subtypes with cell type's"))

opt = parse_args(OptionParser(option_list=option_list))
run_subtype <- (opt$run_subtype == 'TRUE')
cts <- opt$ct
match_with_ct <- (opt$match_with_ct == 'TRUE')
file_suffix <- ifelse(match_with_ct, '_matched', '')

sprintf('%s, run subtype: %s, %s', ct, run_subtype, file_suffix) %>% print

data_dir <- './data/ROSMAP/Gene Expression (snRNAseq - DLPFC, Experiment 2)/processed (March 2024 update)'  
if(run_subtype){
    cell_annotation <- read.csv(sprintf('%s/cell-annotation.n424.csv', data_dir))
    # handle Ex and cux2+, cux2-
    annotation_ct <- ifelse(grepl('cux2', ct), 'Excitatory Neurons', ct)
    # get gene expression data for all subtypes
    sub_cts <- table(cell_annotation$state[cell_annotation$cell.type == annotation_ct]) %>% sort(decreasing=T) %>% names
    print(sprintf('%s subtypes:', ct))
    print(sub_cts)
}
ct_seurat_fn <- c('microglia', 'astrocytes',
        		  'cux2+', 'cux2-',
        		  'oligodendroglia', 'oligodendroglia',
        		  'inhibitory',
        		  'vascular.niche', 'vascular.niche', 'vascular.niche', 'vascular.niche', 'excitatory')
names(ct_seurat_fn) <- c('Microglia', 'Astrocyte',
            			 'cux2+', 'cux2-',
            			 'Oligodendrocytes', 'OPCs',
                         'Inhibitory Neurons',
            			 'Endothelial', 'Fibroblast', 'Pericytes', 'SMC', 'excitatory')
subj_var <- 'individualID'
batch_var <- 'batch_combined'


if(ct == 'excitatory'){
    obj_list <- list()
    obj_list[['cux2-']] <- LoadH5Seurat(sprintf("%s/%s.h5Seurat", data_dir, 'cux2-'), assays='RNA')
    obj_list[['cux2+']] <- LoadH5Seurat(sprintf("%s/%s.h5Seurat", data_dir, 'cux2+'), assays='RNA')
    obj <- merge(obj_list[['cux2-']], obj_list[['cux2+']], add.cell.ids = c("cux2-", "cux2+"))
    print('cux2- and cux2+ combined')
    print(dim(obj))
}else{
    obj <- LoadH5Seurat(sprintf("%s/%s.h5Seurat", data_dir, ct_seurat_fn[ct]), assays='RNA')
}

# handle Ex:
# only part of Ex subtypes are saved in either cux2+/cux2-
# focus on those subtypes
if(grepl('cux2', ct) & run_subtype){
	print('handle Excitatory neurons data')
	# subset cell subtypes for Ex
	covered_state <- unique(cell_annotation$state[cell_annotation$barcode %in% colnames(obj)]) %>% sort()
	sub_cts <- sub_cts[sub_cts %in% covered_state]
	print('overlapped cell states:')
	print(sub_cts)
}

# handle cell types underlying vascular niche:
# the seurat object contain cell types other than those annotated in cell_annotation
# remove those cells when performing cell-type level aggregation
if(ct_seurat_fn[ct] %in% c('vascular.niche', 'oligodendroglia')){
	print('handle vascular niche / oligodendroglia data')
	cell_annotation <- read.csv(sprintf('%s/cell-annotation.n424.csv', data_dir))
	obj$barcode <- rownames(obj@meta.data)
	obj <- subset(obj, subset = barcode %in% cell_annotation$barcode[cell_annotation$cell.type == ct])
}


# -
# remove cells that failed to be mapped to an individual (individualID='NA') or sequenced in a duplicate batch
# following Fujita et al., 2024
# -
print(sprintf('total #cells: %i', ncol(obj)))
# remove cells that were not matched to the 424 subjects
obj <- subset(obj, subset = individualID != 'NA')
cell_annotation <- read.csv(sprintf('%s/cell-annotation.n424.csv', data_dir))
obj$barcode <- rownames(obj@meta.data)
obj <- subset(obj, subset = barcode %in% cell_annotation$barcode)
print(sprintf('#cells matched to 424 subjects: %i', ncol(obj)))
# Combine ...-A and ...-B, e.g. 190403-B4-A and 190403-B4-B as one batch
# "Libraries from four channels were pooled and sequenced on one lane of the Illu- mina HiSeq X "
# As there are a total of eight channels on 10x GEM machine, we assume that -A represents the first four libraries, and -B the second set of four libraries.
obj@meta.data$batch_combined <- sapply(1:nrow(obj@meta.data),
                                       function(i){
                                        x = obj@meta.data$batch[i]
                                        substr(x, 1, nchar(x)-2)}) 
# Identify individualIDs sampled in multiple batches
meta_data1 <- obj@meta.data %>% select('individualID', 'batch_combined') %>% unique()
inds_with_dup_batches <- names(which(table(meta_data1$individualID)>1))
keep_batch_list <- character(length(inds_with_dup_batches)) 
names(keep_batch_list) <- inds_with_dup_batches
for(ind in inds_with_dup_batches){
    keep_batch_list[ind] <- names(which.max(table(obj$batch_combined[obj$individualID == ind])))
}
# Remove the batch with less samples
# following the pre-processing in the original paper
# "Among the remaining 436 specimens, 12 individuals were sequenced twice in distinct batches. After comparing sequencing metrics, one of these duplicates was excluded from further analyses."
keep_cell_inds <- rep(T, ncol(obj))
for(ind in names(keep_batch_list)){
    subset_inds <- (obj$individualID == ind & obj$batch_combined != keep_batch_list[ind])
    keep_cell_inds[subset_inds] <- F
}
obj$keep_cell_inds <- keep_cell_inds
obj <- subset(obj, subset = keep_cell_inds)
print(sprintf('#cells with de-duplicated batch: %i', ncol(obj)))



# -
# match ID between snRNA-seq and WGS data
# -
# load WGS IDs
wgs_samples <- seqVCF_SampID(sprintf('%s/ROSMAP/WGS/DEJ_11898_B01_GRM_WGS_2017-05-15_15.recalibrated_variants.vcf.gz', data_dir))
clinical_covar <- fread(sprintf('%s/ROSMAP_clinical.csv', data_dir))
specimen_covar <- fread(sprintf('%s/ROSMAP_biospecimen_metadata.csv', data_dir))
merged_covar <- merge(clinical_covar, specimen_covar, by = 'individualID', all = TRUE)
merged_covar$newID <- paste0(merged_covar$Study, merged_covar$projid)
# load gene IDs
gene_samples <- unique(obj$individualID)
#colnames(Agg_obj[['RNA']]$counts)
# extract covariate data that covers the gene sample
gene_covar <- merged_covar[merged_covar$individualID %in% gene_samples,]
gene_covar$geno_ID <- rep(NA, nrow(gene_covar))
# some WGS genotype ID are based on Study+projid
gene_covar$geno_ID[gene_covar$newID %in% wgs_samples] <- gene_covar$newID[gene_covar$newID %in% wgs_samples]
# others are based on specimenID
gene_covar$geno_ID[gene_covar$specimenID %in% wgs_samples] <- gene_covar$specimenID[gene_covar$specimenID %in% wgs_samples]
geno_covar_uni <- unique(gene_covar[!is.na(gene_covar$geno_ID), c('individualID', 'geno_ID')])
print(dim(geno_covar_uni)) # some individuals with multiple WGS samples
geno_covar_matched <- geno_covar_uni[match(gene_samples, geno_covar_uni$individualID),] # remove one of the WGS samples for those who have more than one
print(dim(geno_covar_matched))

obj$geno_ID <- geno_covar_matched$geno_ID[match(obj$individualID, geno_covar_matched$individualID)] # match the new geno_IDs to cells
subj_var <- 'geno_ID'

if(ct == 'microglia'){
	# save the matching between WGS ID and snRNAseq ID
	write.table(geno_covar_matched, sprintf('%s/WGS_snRNAseq_sample_matching.txt', data_dir))
	# create a fam file for WGS samples that are matched
	wgs_samples_matched <- wgs_samples[match(geno_covar_matched$geno_ID, wgs_samples)]
	print(length(wgs_samples_matched))
	write.table(data.frame(FID=0, IID=wgs_samples_matched),
		    sprintf('%s/ROSMAP/WGS/plink_files/match_WGS_samples_subset.txt', data_dir), sep = '\t', quote = F, col.names = F, row.names = F)
    }
if(ct %in% c('Endothelial', 'Fibroblast', 'Pericytes', 'SMC')){
    print(dim(obj))
    saveRDS(dim(obj), sprintf('%s/scTransform_by_celltype_afterAgg/dimension_%s.rds', data_dir, ct))
    q()
}



if(!run_subtype){
	sub_cts <- cts
}else{
	sub_ct_sum <- matrix(nrow=length(sub_cts), ncol=3)
	rownames(sub_ct_sum) <- sub_cts
	colnames(sub_ct_sum) <- c('ncells', 'ngenes', 'nsubj')
}

# -
# save pseudo-bulk data by cel ltype
# -
for(sub_ct in sub_cts){
    print('---------------------------')
    print(sprintf('-----------%s----------', sub_ct))
    print('---------------------------')
    # -
    # save scTWAS object
    # -
    if(run_subtype){
        obj_sub <- subset(obj, barcode %in% cell_annotation$barcode[cell_annotation$state == sub_ct])
    }else{
    	obj_sub <- obj
    }
    sprintf('#cells in %s: %i', sub_ct, ncol(obj_sub))
    if(run_subtype) sub_ct_sum[sub_ct, 1] <- ncol(obj_sub)

    Agg_obj <- AggregateExpression(object=obj_sub, assay = 'RNA', group.by=subj_var, return.seurat = TRUE)

    meta_data1 <- obj@meta.data %>% select(all_of(subj_var), batch_combined) %>% unique()
    meta_data1 <- meta_data1[match(colnames(Agg_obj),meta_data1[[subj_var]]),]
    rownames(meta_data1) <- meta_data1[[subj_var]]
    Agg_obj@meta.data <- meta_data1 
    
    ## scTransform normalization
    data <- SCTransform(object = Agg_obj, return.only.var.genes = FALSE)
    
    sprintf('#individual samples with %s cells after aggregation: %i', sub_ct, ncol(data)) %>% print
    
    sub_ct_save <- ifelse(grepl('cux2', ct),
                          sprintf('%s_%s', ct,sub_ct),
                          sub_ct)

    # scTWAS object
    saveRDS(data, file =  sprintf('%s/scTransform_by_celltype_afterAgg/%s.rds',data_dir,sub_ct_save))

    # -
    # AN-TWAS
    # -
    ## TMM+log+batch correct
    library(edgeR) # for TMM
    library(sva)
    pb_mat <- GetAssayData(Agg_obj, assay = 'RNA', layer = 'counts')
    # select genes
    if(match_with_ct){
	    print(sprintf('Match with %s', ct))
	    # use the same genes for subtypes that were selected at the cell type level
        # this is used for microglia subtypes to enable evaluating the same genes as in microglia
        ct_PB_df <- read.table(sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_AN.tsv',data_dir,ct), header = T)
        min_count <- which(!rownames(pb_mat) %in% ct_PB_df$gene_name)
    }else{
	    # otherwise, select genes based on total gene counts
        min_count <- which(rowSums(pb_mat)<1000)
    }
    pb_mat <- pb_mat[-min_count,]
    pb_mat <- pb_mat[order(rowSums(pb_mat),decreasing=T),]

    print(dim(pb_mat))

    dgList <- DGEList(pb_mat, genes=rownames(pb_mat))
    dgList <- calcNormFactors(dgList, method = "TMM")
    expr_norm = log_transform(dgList) # same as voom(dgList)$E
    donor_pool <- (Agg_obj@meta.data)[[batch_var]]

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
    fwrite(PB_df, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_AN%s.tsv',data_dir,sub_ct_save, file_suffix),sep='\t')

        write.table(colnames(PB_df),sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s%s.header',
					    data_dir,sub_ct_save,file_suffix),quote=F,row.names=F,col.names=F)

    if(run_subtype){
	  sub_ct_sum[sub_ct,2] <- nrow(expr_norm)
	  sub_ct_sum[sub_ct,3] <- ncol(expr_norm)
    }

    # -
    # NA-TWAS
    # -
    print('scTWAS')
    print(dim(obj_sub))
    sct_obj <- NormalizeData(obj_sub,normalization.method = "LogNormalize", scale.factor = 1e6)
    unique_donor <- sort(unique(sct_obj@meta.data[[subj_var]]))
    # match with PBINT
    print(all(unique_donor == colnames(expr_norm_inrt)))
    unique_donor <- unique_donor[match(colnames(expr_norm_inrt), unique_donor)]
    print(all(unique_donor == colnames(expr_norm_inrt)))
    pr = matrix(NA, nrow = nrow(sct_obj), ncol = length(unique_donor))
    print(length(unique_donor))
    colnames(pr) = unique_donor; rownames(pr) = rownames(sct_obj)
    n_cells <- numeric(length(unique_donor))
    names(n_cells) <- unique_donor
    for(i in unique_donor){
        ind = which(sct_obj@meta.data[[subj_var]] == i)
        pr[,i] = rowMeans((sct_obj[['RNA']]$data[,ind,drop=FALSE]))
    	n_cells[i] = length(ind)
    }
    if(grepl('cux', sub_ct_save)){
	    saveRDS(pr, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_pr_mean%s.rds',data_dir,sub_ct_save, file_suffix))
	    saveRDS(n_cells, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_pr_n_count%s.rds', data_dir,sub_ct_save,file_suffix))
	    print('save PR mean before INT normalization')
    }
    # # SKIP quantile normalization
    # # https://support.bioconductor.org/p/77664/#77665
    pr_inrt <- matrix(NA, nrow = nrow(pr), ncol = ncol(pr))
    print(anyNA(pr))
    for(i in 1:nrow(pr)){
	    pr_inrt[i,] <- INT(pr[i,])
    }
    rownames(pr_inrt) = rownames(pr)
    colnames(pr_inrt) = colnames(pr)

    PR_df = data.frame(gene_name = rownames(pr_inrt))
    PR_df = cbind(PR_df,pr_inrt)
    PR_df = merge(PR_df,gene_info %>% dplyr::select(-ensembl_gene_id),by.x='gene_name',by.y='external_gene_name',sort=FALSE)
    PR_df = PR_df %>% relocate(gene_name, chromosome_name, start_position, end_position)
    PR_df = PR_df[match(PB_df$gene_name,PR_df$gene_name),]
    fwrite(PR_df, sprintf('%s/scTransform_by_celltype_afterAgg/Expression_matrices/%s_NA%s.tsv',
			  data_dir,sub_ct_save,file_suffix),sep='\t')
}

print('Data saved')

if(run_subtype){
    write.table(sub_ct_sum, sprintf('%s/scTransform_by_celltype_afterAgg/%s_subtypes_summary%s.txt',
				    data_dir,ct,file_suffix))
    print(sprintf('%s/scTransform_by_celltype_afterAgg/%s_subtypes_summary%s.txt',data_dir,ct,file_suffix))
}


