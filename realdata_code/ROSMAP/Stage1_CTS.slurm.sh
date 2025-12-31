#!/bin/bash
#SBATCH --job-name=Stage1
#SBATCH --partition=su_lab,encore
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -o outs/Stage1-%A_%a.out

ct=$1
model_hidden_batch=$2

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

conda activate scTWAS_new

echo "covariate file suffix: $QCOVAR_suffix"

date

GCTA="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"
PLINK="software/plink2"
FUSION="software/FUSION/fusion_twas-master/"

# concatenate two words if necessary 
PRE=${ct}
PRE_save=$(echo "$ct" | sed 's/ /_/g')
#PRE_save=$PRE
echo $PRE_save

AGG="scTransform_by_celltype_afterAgg"
datadir="./data"
SCTOBJ="${datadir}/scTransform_by_celltype_afterAgg/${PRE}.rds"
NR=${SLURM_ARRAY_TASK_ID}
gwas="${datadir}/GWAS/AD/Bellenguez/GCST90027158_build37_imputed.txt"

outdir="${datadir}/output/ROSMAP_Columbia"
mkdir -p ${outdir}/expression_matrices/tmp/$PRE_save
mkdir -p ${outdir}/WEIGHTS/$AGG/$PRE_save/

QCOVAR="${datadir}/output/ROSMAP_Columbia/covariates${QCOVAR_suffix}.txt"
HEADER="${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/${PRE}${match_suffix}.header"  # 05/22: fixed a bug

for ((row_number = $NR*50 + 1; row_number <= ($NR+1)*50; row_number++))
do
    # Read the i-th row and save it in PARAM_i variable
    PARAM=$(sed -n "${row_number}p" "${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/${PRE}_AN${match_suffix}.tsv")
    GNAME=`echo $PARAM | awk '{ print $1 }'`
    CHR=`echo $PARAM | awk '{ print $2 }'` 
    P0=`echo $PARAM | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
    P1=`echo $PARAM | awk '{ print $4 + 500e3 }'` 
    GENO="${datadir}/ROSMAP/WGS/plink_files/chr${CHR}"
    
    OUT="${outdir}/expression_matrices/tmp/$PRE_save/$GNAME"
    
    # skip if already computed
    FINAL_OUT="${outdir}/WEIGHTS/${AGG}/${PRE_save}/${GNAME}${QCOVAR_suffix}"

    if [ -f "${FINAL_OUT}_NATWAS.wgt.RDat" ]; then
            echo "File already exists. Proceed to the next gene."
            echo "${FINAL_OUT}_NATWAS.wgt.RDat"
    else

    
    # Extract the locus around this gene
    rm -f "$OUT.bed"
    ${PLINK} --allow-no-sex -silent --bfile $GENO --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out ${OUT}_0
    Rscript ../pipeline/intersect_SNP.R \
    --bfile ${OUT}_0 --out ${OUT}.snp --gwas ${gwas}
    ${PLINK} --allow-no-sex -silent --bfile ${OUT}_0 --extract ${OUT}.snp --make-bed --out ${OUT}
    rm -f "${OUT}_0."*
    if [ ! -f $OUT.bed ]; then
    continue
    fi

    # Convert the expression matrix to a PLINK format phenotype
    echo $GNAME $CHR $P0 $P1
    echo $PARAM | tr ' ' '\n' | paste "$HEADER" "$HEADER" - | awk '{print "0", $2, $3}' | tail -n+5 > ${OUT}_ANTWAS.pheno

    ### save NATWAS as PLINK format ###
    PARAM=$(sed -n "${row_number}p" "${datadir}/scTransform_by_celltype_afterAgg/Expression_matrices/${PRE}_NA${match_suffix}.tsv")
    GNAME=`echo $PARAM | awk '{ print $1 }'`
    CHR=`echo $PARAM | awk '{ print $2 }'`
    P0=`echo $PARAM | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
    P1=`echo $PARAM | awk '{ print $4 + 500e3 }'`
    GENO="${datadir}/ROSMAP/WGS/plink_files/chr${CHR}"

    OUT="${outdir}/expression_matrices/tmp/$PRE_save/$GNAME"

    echo $GNAME $CHR $P0 $P1
    echo $PARAM | tr ' ' '\n' | paste "$HEADER" "$HEADER" - | awk '{print "0", $2, $3}' | tail -n+5 > ${OUT}_NATWAS.pheno

    gOUT="${outdir}/expression_matrices/tmp/Microglia/$GNAME"
    
    # -
    # Run different cell-type-specific TWAS methods
    # -

    # Run scTWAS
    Rscript ../../R/scTWAS_IRLS.R --bfile $gOUT --tmp ${OUT}.tmp --out ${FINAL_OUT}_scTWAS --verbose 0  --PATH_plink $PLINK --qcovar $QCOVAR --sctobj "$SCTOBJ" --gene $GNAME --niter 10

    # Run ANTWAS with FUSION
     Rscript ../../R/FusionStage1.R \
    --bfile $gOUT --tmp ${OUT}.tmp --out ${FINAL_OUT}_ANTWAS --verbose 0 --PATH_plink $PLINK --models ind_PBINT  --qcovar $QCOVAR --scale_pheno --pheno "${OUT}_ANTWAS.pheno"

     # Run NATWAS with FUSION
     Rscript ../../R/FusionStage1.R \
    --bfile $gOUT --tmp ${OUT}.tmp --out ${FINAL_OUT}_NATWAS --verbose 0 --PATH_plink $PLINK --models ind_PBINT  --qcovar $QCOVAR --scale_pheno --pheno ${OUT}_NATWAS.pheno
     fi 
done

date
