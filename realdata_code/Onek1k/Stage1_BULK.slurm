#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH -p agsmall,ag2tb
#SBATCH --array=0-282
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH -o ~/pbs/slurm-%j.out


cd $SLURM_SUBMIT_DIR
module load R/4.4.0-openblas-rocky8 

module load plink/2.00-alpha-091019
PLINK="~/plink"

datadir="~/data"
NR=${SLURM_ARRAY_TASK_ID}
gwas="~/data/GWAS/E4_DM1_build37_imputed.txt"


mkdir -p ~/UKB_TWAS/expression_matrices/tmp/BULK

QCOVAR="~/data/OneK1K_covariates_PCs_PFs_14_cell_types/tmp/CD4all_peer_factors.tsv.qcovar"
HEADER="~/data/scTransform_by_celltype_afterAgg/Expression_matrices/Bulk.header"

for ((row_number = $NR*50 + 1; row_number <= ($NR+1)*50; row_number++))
do
    # Read the i-th row and save it in PARAM_i variable
    PARAM=$(sed -n "${row_number}p" ${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/Bulk_INT.tsv)
    GNAME=`echo $PARAM | awk '{ print $1 }'`
    CHR=`echo $PARAM | awk '{ print $2 }'` 
    P0=`echo $PARAM | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
    P1=`echo $PARAM | awk '{ print $4 + 500e3 }'`
    GENO="~/data/filter_vcf_r08/chr${CHR}.dose.filtered.R2_0.8.vcf.gz"

    OUT="~/UKB_TWAS/expression_matrices/tmp/BULK/$GNAME"
    
    # Extract the locus around this gene
    rm -f $OUT.bed
    plink2 --allow-no-sex -silent --vcf $GENO --chr $CHR --from-bp $P0 --to-bp $P1 --hwe 1e-6 --maf 0.01 --make-bed --out ${OUT}_0 ## ADD QC FILTER!!
    Rscript ~/intersect_SNP.R \
    --bfile ${OUT}_0 --out ${OUT}.snp --gwas ${gwas}
    plink2 --allow-no-sex -silent --bfile ${OUT}_0 --extract ${OUT}.snp --make-bed --out ${OUT}
    rm -f ${OUT}_0.*
    if [ ! -f $OUT.bed ]; then
    continue
    fi
    
    # Convert the expression matrix to a PLINK format phenotype
    echo $GNAME $CHR $P0 $P1
    echo $PARAM | tr ' ' '\n' | paste $HEADER $HEADER - | awk '{print "0", $2, $3}' | tail -n+5 > ${OUT}_INT.pheno

    # Run FUSION
    FINAL_OUT="~/UKB_TWAS/WEIGHTS/BULK_revision1/$GNAME" 

    Rscript ~/R/FusionStage1.R \
    --bfile $OUT --tmp ${OUT}.tmp --out ${FINAL_OUT}_BULK --verbose 0 --PATH_plink $PLINK  --qcovar $QCOVAR --scale_pheno --pheno ${OUT}_INT.pheno
    
done
