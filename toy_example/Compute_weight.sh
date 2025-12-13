PLINK="YOUR/PATH/TO/plink"


OUT="gene"
SCTOBJ="./sct.rds"

FINAL_OUT1="./RPL13"
GNAME="RPL13"

## Run scTWAS to train STAGE 1 GReX model
Rscript ./scTWAS_IRLS.R \
     --bfile $OUT --tmp ${OUT}.tmp --out ${FINAL_OUT1}_PRWGT --verbose 0  --PATH_plink $PLINK --sctobj ${SCTOBJ} --gene $GNAME --niter 3
 

