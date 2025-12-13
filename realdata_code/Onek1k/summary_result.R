# library(venn)
library(data.table)
library(dplyr)
# library(ggVennDiagram)
source('/home/panwei/lin00374/TWAS_helper.R')

wgtdir = '~/UKB_TWAS/WEIGHTS/'
resdir = '~/UKB_TWAS/result1/'
gwas_list = c('30000_irnt','30010_irnt','30020_irnt','30030_irnt','30040_irnt','30050_irnt','30060_irnt','30070_irnt','30080_irnt','30090_irnt', '30100_irnt','30110_irnt','30120_irnt','30130_irnt','30140_irnt','30150','30160','30170','30180_irnt','30190_irnt', '30200_irnt','30210_irnt','30220_irnt','30230','30240_irnt','30250_irnt','30260_irnt','30270_irnt','30280_irnt','30290_irnt','30300_irnt','RA2014','SLE2015','Asthma2020')
ct_list = c('CD4all','CD8all','CD8eff','Bmem','BimmNaive','MonoC','MonoNC','NKact','NKmat','DC','Plasma','CD4effCM','CD8unknown') 

res_df1 = NULL
for(ct in ct_list){
    for(gwas in gwas_list){
        skip = FALSE
        ANTWAS_res = tryCatch({
            fread(sprintf('%s/%s/%s/%s_ANTWAS.dat',resdir,ct,gwas))
            }, error = function(err) {
             skip <<- TRUE
            })
        if(skip){next;}
        NATWAS_res = fread(sprintf('%s/%s/%s/%s_NATWAS.dat',resdir,ct,gwas))
        scTWAS_res = tryCatch({
            fread(sprintf('%s/%s/%s/%s_scTWAS.dat',resdir,ct,gwas))
            }, error = function(err) {
             skip <<- TRUE
            })
        if(skip){next;}
        gene.all = union((ANTWAS_res %>% filter(!is.na(TWAS.P)))$ID,
                         (scTWAS_res %>% filter(!is.na(TWAS.P_converge) & CHR!=0))$ID)
        gene.all = union(gene.all, (NATWAS_res %>% filter(!is.na(TWAS.P)))$ID)
        ANTWAS_R2sig_gene = (ANTWAS_res %>% filter(ID %in% gene.all & MODELCV.PV<0.05/length(gene.all)))$ID
        NATWAS_R2sig_gene = (NATWAS_res %>% filter(ID %in% gene.all & MODELCV.PV<0.05/length(gene.all)))$ID
        scTWAS_R2sig_gene = (scTWAS_res %>% filter(ID %in% gene.all & MODELCV.PV_converge<0.05/length(gene.all)))$ID

        # fwrite(cbind(gwas,ANTWAS_res %>% filter(ID %in% ANTWAS_R2sig_gene) %>% select(ID,CHR,P0,P1,MODELCV.R2,TWAS.P)),sprintf('%s/../Stage1/%s_ANTWAS.txt',resdir,ct),sep='\t',append=TRUE)
        # fwrite(cbind(gwas,NATWAS_res %>% filter(ID %in% PR_R2sig_gene) %>% select(ID,CHR,P0,P1,MODELCV.R2,TWAS.P)),sprintf('%s/../Stage1/%s_NATWAS.txt',resdir,ct),sep='\t',append=TRUE)
        # fwrite(cbind(gwas,scTWAS_res %>% filter(ID %in% scTWAS_R2sig_gene) %>% select(ID,CHR,P0,P1,MODELCV.R2_converge,TWAS.P_converge)),sprintf('%s/../Stage1/%s_scTWAS.txt',resdir,ct),sep='\t',append=TRUE)

        R2gene.all = (union(ANTWAS_R2sig_gene,scTWAS_R2sig_gene))
        R2gene.all = union(R2gene.all,NATWAS_R2sig_gene)
        a = ANTWAS_res %>% filter(ID %in% ANTWAS_R2sig_gene)  %>% filter(TWAS.P<0.05/length(R2gene.all))
        ANTWAS_TWASsig_gene = a$ID
        # fwrite(cbind(gwas,a %>% select(ID,CHR,P0,P1,TWAS.P)),sprintf('%s/%s_ANTWAS.txt',resdir,ct),sep='\t',append=TRUE)

        
        a = NATWAS_res %>% filter(ID %in% NATWAS_R2sig_gene)  %>% filter(TWAS.P<0.05/length(R2gene.all))
        NATWAS_TWASsig_gene = a$ID
        # fwrite(cbind(gwas,a %>% select(ID,CHR,P0,P1,TWAS.P)),sprintf('%s/%s_NATWAS.txt',resdir,ct),sep='\t',append=TRUE)

        a = scTWAS_res %>% filter(ID %in% scTWAS_R2sig_gene)  %>% filter(TWAS.P_converge<0.05/length(R2gene.all))
        scTWAS_TWASsig_gene = a$ID
        # fwrite(cbind(gwas,a %>% select(ID,CHR,P0,P1,TWAS.P_converge)),sprintf('%s/%s_scTWAS.txt',resdir,ct),sep='\t',append=TRUE)

        res_df1 = rbind(res_df1, c(ct,gwas,length(scTWAS_R2sig_gene),
                                   length(NATWAS_R2sig_gene),
                                   length(ANTWAS_R2sig_gene),length(scTWAS_TWASsig_gene),
                                   length(NATWAS_TWASsig_gene),
                                   length(ANTWAS_TWASsig_gene)))
    }
}
res_df1 = as.data.frame(res_df1)
colnames(res_df1) = c('ct','gwas','scTWAS_S1','NATWAS_S1','ANTWAS_S1','scTWAS_S2','NATWAS_S2','ANTWAS_S2')
