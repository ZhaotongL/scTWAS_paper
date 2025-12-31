#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH -p su_lab,encore
#SBATCH --mem=25G
#SBATCH --array 1-22
#SBATCH -o outs/Stage2_%A_%a.out

conda activate scTWAS_new

date 

batch_suffix=''
ct=Bulk
run_subtype=FALSE

FUSION="${datadir}/software/FUSION/fusion_twas-master/"

WGT="scTransform_by_celltype_afterAgg"
posdic="${datadir}/output/ROSMAP_Columbia/WEIGHTS"
lddir="${datadir}/ROSMAP/WGS/plink_files"
OUT="${datadir}/output/ROSMAP_Columbia/result"
CHR=${SLURM_ARRAY_TASK_ID}

subtypedir="/sulab/data/ROSMAP/Gene Expression (snRNAseq - DLPFC, Experiment 2)/processed (March 2024 update)/scTransform_by_celltype_afterAgg"

list_of_GWASs=('GCST90027158')
if [ $run_subtype == 'TRUE' ]
then 
	list_of_CTs=$(awk 'NR>1 && $2 > 1000 {print $1}' "${subtypedir}/${ct}_subtypes_summary.txt" | sed 's/"//g'  | tr -d '"')
	list_of_CTs=($list_of_CTs)
else
	list_of_CTs=$ct
	list_of_CTs=$(echo "$ct" | sed 's/ /_/g')
fi

gwas='GCST90027158'
trait='AD/Bellenguez'
sumstat="${datadir}/GWAS/${trait}/${gwas}_build37_imputed.txt"

echo $gwas
echo $trait
date

for PRE in ${list_of_CTs[@]}; do

    mkdir $OUT/$WGT/${PRE}/ -p
  
    echo $PRE    
     
 	# ANTWAS
    Rscript ../../R/FUSION.assoc_test.R \
        --sumstats $sumstat \
        --weights $posdic/$WGT/${PRE}${batch_suffix}_ANTWAS.pos \
        --weights_dir $posdic/$WGT/$PRE/ \
        --ref_ld_chr $lddir/chr \
        --chr $CHR \
        --out $OUT/$WGT/${PRE}/${gwas}${batch_suffix}_ANTWAS_${CHR}.dat \
        --min_r2pred 0.5 
        
	date    
done
