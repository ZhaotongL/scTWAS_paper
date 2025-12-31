library(ggplot2)
library(dplyr)
library(ggplot2)
library(data.table)

run_subtype <- T
n_models <- 4
mode <- '_converge'
resdir <- 'results'
batch_suffix <- '_with_hidden_batch_'
data_dir <- './data'
cts <- c('Microglia', 'Astrocyte', 'Inhibitory Neurons', 'Oligodendrocytes', 'cux2+', 'cux2-', 'OPCs')


summary <- function(cell_type, run_subtype, mode, n_models, batch_suffix, gwas, FDR_Stage1 = F, to_save = T){    
    
    wgtdir = sprintf('%s/output/ROSMAP_Columbia/WEIGHTS', data_dir)
    resdir = sprintf('%s/output/ROSMAP_Columbia/result', data_dir)
    pre = 'scTransform_by_celltype_afterAgg'
    if(!run_subtype){
    	ct_list <- gsub(" ", "_", cell_type) #cell_type
    }else{
    	subtypes <- read.table(sprintf("%s/ROSMAP/Gene Expression (snRNAseq - DLPFC, Experiment 2)/processed (March 2024 update)/scTransform_by_celltype_afterAgg/%s_subtypes_summary.txt", data_dir, cell_type))
    	print(subtypes)
        # focus on cell subtyeps with at least 100 cells
    	ct_list <- rownames(subtypes)[subtypes[,1] > 1000]
    }
    res_df1 = NULL
    res_list <- Stage1_sig_list <- Stage2_sig_list <- list()#lapply(1:length(ct_list), function(i) list())
    #names(res_list) <- names(Stage1_sig_list) <- names(Stage2_sig_list) <- ct_list
    
    for(ct in ct_list){
    	if(grepl('cux2', cell_type) & run_subtype) ct <- sprintf('%s_%s', cell_type, ct)
        #for(gwas in gwas_list){
	    ### new
	    if(n_models == 2){
		    full_method_res_list <- list(PRWGT = list(), PBINT = list())
	    }else{
		    full_method_res_list <- list(PRWGT = list(), PBINT = list(), Eric = list(), PRUNWGT = list())
	    }
	    method_names <- names(full_method_res_list)
    	for(chr in 1:22){
    		for(mn in method_names){
    			full_method_res_list[[mn]][[chr]] <- fread(file = sprintf('%s/%s/%s/%s%s_%s_%i.dat',resdir,pre,ct,gwas, batch_suffix, mn, chr))
        	}
    	}
    	method_res_list <- lapply(method_names, function(mn) do.call(rbind, full_method_res_list[[mn]]))
    	names(method_res_list) <- method_names

    	# which PRWGT results to evaluate
    	method_res_list[['PRWGT']]$MODELCV.R2 <- method_res_list[['PRWGT']][[sprintf('MODELCV.R2%s', mode)]]
    	method_res_list[['PRWGT']]$MODELCV.PV <- method_res_list[['PRWGT']][[sprintf('MODELCV.PV%s', mode)]]
    	method_res_list[['PRWGT']]$TWAS.P <- method_res_list[['PRWGT']][[sprintf('TWAS.P%s', mode)]]

    	# control for multiple hypothesis testing by considering the same set of genes across methods (all predictable genes in Stage 1)
    	gene.all = Reduce(union, lapply(method_res_list, function(mr) (mr %>% filter(!is.na(TWAS.P) & !is.na(MODELCV.R2)))$ID)) 
    	# predictable genes 
    	if(FDR_Stage1){
            # under FDR control
            R2sig_gene_list <- lapply(method_res_list, function(mr){
    					  (mr %>% filter(!is.na(TWAS.P)) %>% filter(p.adjust(MODELCV.PV, method='BH')<0.05))$ID})
        }else{
            # under Bonforroni correction
            R2sig_gene_list <- lapply(method_res_list, function(mr){
                (mr %>% filter(!is.na(TWAS.P)) %>% filter(MODELCV.PV<0.05/length(gene.all)))$ID})
        }
    	names(R2sig_gene_list) <- method_names

    	# TWAS significant genes: predictable
                                        
        stage2_gene.all <- Reduce(union, R2sig_gene_list)
    	TWASsig_gene_list <- lapply(method_names, function(mn){
            (method_res_list[[mn]] %>% filter(ID %in% R2sig_gene_list[[mn]])  %>% filter(!is.na(TWAS.P)) %>% filter(TWAS.P<0.05/length(stage2_gene.all)))$ID})

    	res_df1 = rbind(res_df1, c(ct, gwas, 
    				   sapply(R2sig_gene_list, length),
    				   sapply(TWASsig_gene_list, length)))
        #}
        res_list[[ct]] <- method_res_list
        Stage1_sig_list[[ct]] <- R2sig_gene_list
        Stage2_sig_list[[ct]] <- TWASsig_gene_list
    }

    res_df1 = as.data.frame(res_df1)
    colnames(res_df1) = c('ct','gwas',
                          sapply(method_names, function(mn) sprintf('%s_R2sig_gene', mn)),
            		      sapply(method_names, function(mn) sprintf('%s_TWASsig_gene', mn)))

    if(!FDR_Stage1){
        # summary
        resdir <- 'results'
        result_prefix <- sprintf('%s/%s_subtype_%s_%s', resdir, cell_type, run_subtype, gwas)
        result_suffix <- sprintf('%sn_models_%i', batch_suffix, n_models)
        fn <- sprintf('%s_result_summary_%s.txt', result_prefix, result_suffix)
        write.table(res_df1, fn)
        print(res_df1)

        res_list <- list(res = res_list, Stage1_sig = Stage1_sig_list, Stage2_sig = Stage2_sig_list)
        if(to_save){
            # detailed results
            saveRDS(res_list,
                sprintf('%s_detailed_results_%s.rds', result_prefix, result_suffix))
        }else{
            #return(res_list)
        }
    }else{
        return(list(res = res_list, Stage1_sig = Stage1_sig_list, Stage2_sig = Stage2_sig_list))
    }
}

sum_table <- NULL
sum_res_gwas <- list()
gwas <- 'GCST90027158'
sum_res <- list()
for(ct in cts){
    result_prefix <- sprintf('%s/%s_subtype_%s_%s', resdir, ct, run_subtype, gwas)
    result_suffix <- sprintf('%sn_models_%i', batch_suffix, n_models)
    fn <- sprintf('%s_detailed_results_%s.rds', result_prefix, result_suffix)
    summary(ct, run_subtype, mode, n_models, batch_suffix, gwas)
    sum_res[[ct]] <- readRDS(fn)
    sum_table <- rbind(sum_table, cbind(summary_short(sum_res[[ct]]$res, T), F, ct))
}


