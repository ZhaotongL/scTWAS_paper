#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH -p agsmall,ag2tb
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --array 1-22
#SBATCH -o ~/pbs/assoc-%j.out


module load plink/2.00-alpha-091019
module load plink
module load R/4.4.0-openblas-rocky8 

datadir="~/data"
posdic="~/UKB_TWAS/WEIGHTS"
lddir="~/data/filter_vcf_r08/bfile/"
OUT="~/UKB_TWAS/result1/"
CHR=${SLURM_ARRAY_TASK_ID}


list_of_GWASs=('30000_irnt' '30010_irnt' '30020_irnt' '30030_irnt')
#list_of_GWASs=('30040_irnt' '30050_irnt' '30060_irnt' '30070_irnt')
#list_of_GWASs=('30080_irnt' '30090_irnt' '30100_irnt')
#list_of_GWASs=('30110_irnt' '30120_irnt' '30130_irnt' '30140_irnt')
#list_of_GWASs=('30180_irnt' '30190_irnt' '30150' '30160')
#list_of_GWASs=('30170' '30200_irnt' '30210_irnt' '30220_irnt')
#list_of_GWASs=('30240_irnt' '30250_irnt' '30260_irnt' '30230')
#list_of_GWASs=('30280_irnt' '30270_irnt' '30290_irnt' '30300_irnt')
#list_of_GWASs=('SLE2015' 'RA2014' 'Asthma2020')
list_of_CTs=('CD4all' 'CD8all' 'CD8eff' 'MonoC' 'MonoNC' 'Bmem' 'BimmNaive' 'NKact' 'NKmat' 'DC' 'Plasma' 'CD4effCM' 'CD8unknown')


for gwas in "${list_of_GWASs[@]}"; do
    sumstat="~/data/GWAS/${gwas}_build37_imputed.txt"
    for PRE in "${list_of_CTs[@]}"; do
    mkdir $OUT/${PRE}/ -p
        
     
    Rscript ~/R/FUSION.assoc_test.R \
       --sumstats $sumstat \
       --weights $posdic/${PRE}_NA.pos \
       --weights_dir $posdic/$PRE/ \
       --ref_ld_chr $lddir/chr \
       --chr $CHR \
       --out $OUT/${PRE}/${gwas}_NA_${CHR}.dat \
       --min_r2pred 0.5 
     
    Rscript ~/R/FUSION.assoc_test.R \
       --sumstats $sumstat \
       --weights $posdic/${PRE}_AN.pos \
       --weights_dir $posdic/$PRE/ \
       --ref_ld_chr $lddir/chr \
       --chr $CHR \
       --out $OUT/${PRE}/${gwas}_AN_${CHR}.dat \
       --min_r2pred 0.5 
      
 
     Rscript ~/R/FUSION.assoc_scTWAS.R \
        --sumstats $sumstat \
        --weights $posdic/${PRE}_scTWAS.pos \
        --weights_dir $posdic/$PRE/ \
        --ref_ld_chr $lddir/chr \
        --chr $CHR \
        --out $OUT/${PRE}/${gwas}_scTWAS_${CHR}.dat \
        --min_r2pred 0.5 
    done
done
        

