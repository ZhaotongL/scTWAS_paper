#!/bin/bash
#SBATCH --job-name=Stage1
#SBATCH --partition=su_lab,encore
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -o outs/Stage1-%A_%a.out
#SBATCH --array 0-330

ct='Bulk'
model_hidden_batch=FALSE
match_suffix=FALSE

if [[ "$model_hidden_batch" == "TRUE" ]]; then
    QCOVAR_suffix='_with_hidden_batch'
else
    QCOVAR_suffix=''
fi

if [[ "$match_suffix" == "TRUE" ]]; then
    match_suffix='_matched'
else
    match_suffix=''
fi

echo "This is job number $SLURM_ARRAY_TASK_ID"

conda activate scTWAS_submit

echo "covariate file suffix: $QCOVAR_suffix"

date

GCTA="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"
PLINK="software/plink2"
FUSION="software/FUSION/fusion_twas-master/"

PRE=${ct}
AGG="scTransform_by_celltype_afterAgg"
datadir="./data"
SCTOBJ="${datadir}/scTransform_by_celltype_afterAgg/${PRE}.rds"
NR=${SLURM_ARRAY_TASK_ID}
gwas="${datadir}/GWAS/AD/Bellenguez/GCST90027158_build37_imputed.txt"

outdir="${datadir}/output/ROSMAP_Columbia"
mkdir -p ${outdir}/expression_matrices/tmp/$PRE_save
mkdir -p ${outdir}/WEIGHTS/$AGG/$PRE_save/

QCOVAR="/sulab/users/csu30/projects/single_cell_TWAS/output/ROSMAP_Columbia/covariates${QCOVAR_suffix}.txt"
HEADER="${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/${PRE}${match_suffix}.header"  # 05/22: fixed a bug

for ((row_number = $NR*50 + 1; row_number <= ($NR+1)*50; row_number++))
do
    # Read the i-th row and save it in PARAM_i variable
    PARAM=$(sed -n "${row_number}p" "${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/${PRE}_PBINT${match_suffix}.tsv")
    GNAME=`echo $PARAM | awk '{ print $1 }'`
    CHR=`echo $PARAM | awk '{ print $2 }'` 
    P0=`echo $PARAM | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
    P1=`echo $PARAM | awk '{ print $4 + 500e3 }'` 
    GENO="${datadir}/ROSMAP/WGS/plink_files/chr${CHR}"
    
    OUT="${outdir}/expression_matrices/tmp/$PRE_save/$GNAME"
    
    # skip if already computed
    FINAL_OUT="${outdir}/WEIGHTS/${AGG}/${PRE_save}/${GNAME}${QCOVAR_suffix}"

    gOUT="${outdir}/expression_matrices/tmp/Microglia/$GNAME"
    
    if [ ! -f "${gOUT}.bed" ]; then
        echo "Generating genotype file for $GNAME"
        echo $gOUT
        rm -f "$gOUT.bed"
        ${PLINK} --allow-no-sex -silent --bfile $GENO --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out ${gOUT}_0
        Rscript ../pipeline/intersect_SNP.R \
        --bfile ${gOUT}_0 --out ${gOUT}.snp --gwas ${gwas}
        ${PLINK} --allow-no-sex -silent --bfile ${gOUT}_0 --extract ${gOUT}.snp --make-bed --out ${gOUT}
        rm -f "${gOUT}_0."*
        #if [ ! -f $gOUT.bed ]; then
        #    continue
        #fi
    fi
    

    # Convert the expression matrix to a PLINK format phenotype
    echo $GNAME $CHR $P0 $P1
    echo $PARAM | tr ' ' '\n' | paste "$HEADER" "$HEADER" - | awk '{print "0", $2, $3}' | tail -n+5 > ${OUT}_PBINT.pheno
    
    
    # ANTWAS
     Rscript ../../R/FusionStage1.R \
    --bfile $gOUT --tmp ${OUT}.tmp --out ${FINAL_OUT}_PBINT --verbose 0 --PATH_plink $PLINK --models ind_PBINT  --qcovar $QCOVAR --scale_pheno --pheno "${OUT}_ANTWAS.pheno"

     fi 
done

date
