library(data.table)
library(dplyr)
data_dir <- '~/data/scTransform_by_celltype_afterAgg'
pos_dir <- '~/UKB_TWAS/WEIGHTS/'

ct = 'MonoNC'


PB_df = fread(sprintf('%s/Expression_matrices/%s_PB.tsv',data_dir,ct))[,1:4]


setwd(sprintf('%s/%s',pos_dir,ct))
t0 = sprintf('ls *scTWAS.wgt.RDat > ../%s_scTWAS.pos',ct)
system(t0)
t0 = sprintf('ls *ANTWAS.wgt.RDat > ../%s_ANTWAS.pos',ct)
system(t0)
t0 = sprintf('ls *NATWAS.wgt.RDat > ../%s_NATWAS.pos',ct)
system(t0)

for(m in c('scTWAS','ANTWAS','NATWAS')){
    pos = fread(sprintf('%s/%s_%s.pos',pos_dir,ct,m),header=FALSE)
    pos$gene <- sub('_.*', '', pos$V1)
    pos = merge(pos,PB_df,by.x='gene',by.y='gene_name')
    pos$P0 =  pos$start_position
    pos$P1 = pos$end_position
    pos = pos %>% select(WGT=V1,ID=gene,CHR=chromosome_name,P0,P1)
    fwrite(pos, sprintf('%s/%s_%s.pos',pos_dir,ct,m), sep='\t')
}
