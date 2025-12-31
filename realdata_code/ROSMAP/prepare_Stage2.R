library(Seurat)
library(data.table)
library(dplyr)
library(optparse)

option_list = list(
  make_option("--model_hidden_batch", action="store", default=TRUE, type='character',
              help="whether to model hidden batch effects"),
  make_option("--ct", action="store", default="Microglia", type="character",
              help="which cell type"),
  make_option("--run_subtype", action="store", default="FALSE", type="character",
              help="whether to analyze cell subtypes"),
  make_option("--match_with_ct", action="store", default="FALSE", type="character",
              help="whether to match the gene evauated in subtypes with cell type's"))
opt = parse_args(OptionParser(option_list=option_list))
model_hidden_batch <- (opt$model_hidden_batch == 'TRUE')
run_subtype <- (opt$run_subtype == 'TRUE')
ct <- opt$ct
match_suffix <- ifelse(opt$match_with_ct == 'TRUE', '_matched', '')

print(sprintf('model hidden batch effects: %s', model_hidden_batch))
batch_suffix <- ifelse(model_hidden_batch, '_with_hidden_batch_', '')

data_dir <- './data'

if(!run_subtype){
	ct_list <- ct
}else{
	subtypes <- read.table(sprintf("%s/ROSMAP/Gene Expression (snRNAseq - DLPFC, Experiment 2)/processed (March 2024 update)/scTransform_by_celltype_afterAgg/%s_subtypes_summary.txt", data_dir, ct))
	print(subtypes)
	subtype_names <- rownames(subtypes)
	if(grepl('cux2', ct)){
		subtype_names <- sapply(subtype_names, function(sn) sprintf('%s_%s', ct, sn))
	}	
	ct_list <- subtype_names
}
print(ct_list)

for(ct in ct_list){
  PB_df = fread(sprintf('%s/Expression_matrices/%s_PBINT%s.tsv',data_dir,ct, match_suffix))[,1:4]

  ### Stage 1 prediction model results
  pos_dir <- sprintf('%s/output/ROSMAP_Columbia/WEIGHTS/scTransform_by_celltype_afterAgg', data_dir)

  ct_save <- gsub(" ", "_", ct)
  setwd(sprintf('%s/%s',pos_dir,ct_save))
  for(mn in c('scTWAS', 'ANTWAS', 'NATWAS')){
    if(batch_suffix == '_with_hidden_batch_'){
		  t0 = sprintf('ls *%s%s.wgt.RDat > ../%s%s_%s.pos', batch_suffix, mn, ct_save, batch_suffix, mn)
	  }else if(batch_suffix == ''){
		  # exclude the files that matched with "_with_hidden_batch_"
		  t0 = sprintf("ls -1  *%s%s.wgt.RDat | grep -v '_with_hidden_batch_' > ../%s%s_%s.pos", batch_suffix, mn, ct_save, batch_suffix, mn)
	  }
	  system(t0)
  }

  for(m in c('scTWAS', 'ANTWAS', 'NATWAS')){
    pos = fread(sprintf('%s/%s%s_%s.pos',pos_dir,ct_save,batch_suffix,m),header=FALSE)
    pos$gene <- sub('_.*', '', pos$V1)
    pos = merge(pos,PB_df,by.x='gene',by.y='gene_name')
    pos$P0 =  pos$start_position
    pos$P1 = pos$end_position 
    pos = pos %>% select(WGT=V1,ID=gene,CHR=chromosome_name,P0,P1)
    fwrite(pos, sprintf('%s/%s%s_%s.pos',pos_dir,ct_save,batch_suffix,m), sep='\t')
  }
}
