# ==== TODO
# * Make sure BLUP/BSLMM weights are being scaled properly based on MAF
suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))
suppressMessages(library('Seurat'))
suppressMessages(library('caret'))


option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--sctobj", action="store", default=NA, type='character',
              help="Path to molecular phenotype Seurat SCT object [optional]"),
  make_option("--gene", action="store", default=NA, type='character',
              help="Gene to be extracted if sctobj is provided"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--qcovar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--batch", action="store", default=NA, type='character',
              help="Path to batch covariates (PLINK format) [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--niter", action="store", default=10, type='double',
              help="How many iterations to update weights, 1 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--scale_pheno", action="store_true", default=FALSE,
              help="Scale phenotype or not [default: %default]") ## NOT USE 

)

opt = parse_args(OptionParser(option_list=option_list))
hsq = 0.1
hsq.pv = NA
if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

## sctransform helper
extract_sct_data <- function(sct_obj, gene){
    #sct_PR = GetAssayData(sct_obj, assay = 'SCT', layer = 'scale.data')
    SCTres = SCTResults(sct_obj,slot='feature.attributes')
    SCTres = SCTres[match(gene,rownames(SCTres)),]
    cell_attr <- SCTResults(sct_obj,slot='cell.attributes')
    mu_mle <- exp(SCTres$`(Intercept)`)
    seq_depths <- as.numeric(cell_attr$umi)
    mu_mat <- outer(mu_mle, seq_depths)
    sigma_sq_mat <- mu_mat + mu_mat^2 / SCTres$theta
    counts_sub <- sct_obj[['RNA']]$counts[match(rownames(SCTres['(Intercept)']),rownames(sct_obj[['RNA']]$counts)),]
    min_var <- SCTResults(sct_obj,slot='arguments')$set_min_var
    sigma_sq_mat_th <- sigma_sq_mat
    sigma_sq_mat_th[sigma_sq_mat_th < min_var] <- min_var
    sct_res_th <- (counts_sub - mu_mat) / sqrt(sigma_sq_mat_th)
    names(sigma_sq_mat_th) = names(sct_res_th) = names(counts_sub)
    names(seq_depths) = names(mu_mat) = names(sct_res_th)

    return(list(raw_counts = counts_sub, weights = 1 / sqrt(sigma_sq_mat_th), mu_mle = mu_mle,
               seq_depths = seq_depths, pearson_res = sct_res_th, mu_mat = mu_mat,
               theta = SCTres$theta, min_var =  min_var))
}

# --- PREDICTION MODELS


# Elastic Net/LASSO
weights.enet = function( genos , pheno , covar = NULL, intercept= T, alpha=0.5 ) {
	eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
	# remove monomorphics
	sds = apply( genos  , 2 , sd )
	keep = sds != 0 & !is.na(sds)
    X = cbind(covar,genos[,keep])
    if(is.null(covar)){
        penalty.v = rep(1,ncol(genos[,keep]))
    }else{
        penalty.v = c(rep(0,ncol(covar)),rep(1,ncol(genos[,keep])))
    }
	enet = cv.glmnet( x=as.matrix(X) , y=pheno , alpha=alpha , nfold=5 , intercept=intercept , standardize=F, penalty.factor = penalty.v )
    if(!is.null(covar)){
        	eff.wgt[ keep ] = coef( enet , s = "lambda.min")[-(1:(ncol(covar)+1))]
            intercept = coef( enet , s = "lambda.min")[2] ## only for PRWGT model
    }else{
        	eff.wgt[ keep ] = coef( enet , s = "lambda.min")[-(1:(0+1))]
            intercept = NULL
    }
	return( list(eff.wgt=eff.wgt,intercept=intercept,enet=enet) )
}

# --- CLEANUP
cleanup = function() {
	if ( ! opt$noclean ) {
		arg = paste("rm -f " , opt$tmp , "*", sep='')
		system(arg)
	}
}

# Perform i/o checks here:
files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$batch) ) files = c(files,opt$batch)

for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}


fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

if(!is.na(opt$sctobj)){
    bulk_obj <- readRDS(opt$sctobj)
    sct_res <- extract_sct_data(bulk_obj,opt$gene)
    pheno = data.frame(FID=0,IID=names(sct_res$pearson_res))
    pheno = cbind(pheno,pr=as.numeric(sct_res$pearson_res),rc=as.numeric(sct_res$raw_counts),
             w = as.numeric(sct_res$weights), s = as.numeric(sct_res$seq_depths), mu = as.numeric(sct_res$mu_mat),
             rc_w = as.numeric(sct_res$raw_counts) * as.numeric(sct_res$weights),
             rc_s = (as.numeric(sct_res$raw_counts)/as.numeric(sct_res$seq_depths)),
             theta = sct_res$theta,     # updated to return theta
             min_var = sct_res$min_var)
    m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
    m.keep = !is.na(m)
    fam = fam[m.keep,]
    m = m[m.keep]
    pheno = pheno[m,]
    pheno_sct = pheno
}
## pheno and fam mathced

# Load in the covariates if needed
pheno_covar = pheno[,1:3]
if ( !is.na(opt$qcovar) ) {
	qcovar = ( read.table(opt$qcovar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(qcovar)-2,"quantitative covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(qcovar[,1],qcovar[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	qcovar = qcovar[m,]
    pheno_covar = cbind(pheno_covar,qcovar[,-(1:2)])
    ## MINOR UPDATE
    pheno_covar[,-(1:3)] <- scale(pheno_covar[,-(1:3)])
    }
## pheno, pheno_covar and fam mathced

batch = NULL
if ( !is.na(opt$batch) ) {
	batch.df = ( read.table(opt$batch,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(batch.df)-2,"batch covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(batch.df[,1],batch.df[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
    pheno_covar = pheno_covar[m.keep,]
	m = m[m.keep]
	batch.df = batch.df[m,]
    ## 
    batch = batch.df
    }
## pheno, pheno_covar, fam and batch mathced

pheno.file = paste(opt$tmp,".pheno",sep='')
write.table(pheno[,c(1,2,8)],quote=F,row.names=F,col.names=F,file=pheno.file)

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)


# read in genotypes
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2
sds = apply(genos$bed,2,sd)
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = genos$fam[,c(1,2,6)]


m = match( paste(pheno[,1],pheno[,2]) , paste(pheno_covar[,1],pheno_covar[,2]) )
if(ncol(pheno_covar)>3){
    pheno_covar = pheno_covar[m,]
    batch = batch[m,]
    ncov = ncol(pheno_covar) - 3
}else{
    pheno_covar = NULL
    ncov = 0
}

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}

if(!is.na(opt$sctobj)){
    m = match( paste(pheno[,1],pheno[,2]) , paste(pheno_sct[,1],pheno_sct[,2]) )
    pheno_sct = pheno_sct[m,]
}


N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- CROSSVALIDATION ANALYSES
set.seed(1) 
cv.performance = matrix(NA,nrow=2,ncol=opt$niter)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = 1:opt$niter

cv.all = pheno
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=opt$crossval,labels=FALSE)

cv.calls = cv.calls1 = matrix(NA,nrow=N,ncol=1)


wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=opt$niter)
colnames(wgt.matrix) = 1:opt$niter
rownames(wgt.matrix) = genos$bim[,2]
intercept = numeric(opt$niter)


for(j in 1:opt$niter){
    new_mean = numeric(N)
    set.seed(123)
    for ( i in 1:opt$crossval ) {
#         if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
        indx = which(folds==i,arr.ind=TRUE)
        cv.train = cv.all[-indx,]
        geno.train = genos$bed[ cv.sample[ -indx ],]
        nasnps.train = apply( is.na(scale(geno.train)) , 2 , sum )
        if ( sum(nasnps.train) != 0 ) {
            geno.train[,nasnps.train != 0] = 0
        }
        covar.train = pheno_covar[cv.sample[ -indx ], -(1:3) ]
        batch.train = batch[cv.sample[ -indx ], ]
        if(!is.na(opt$sctobj)){sct.train = pheno_sct[cv.sample[ -indx ], ]}
        if(is.null(batch.train)){
            covar.train.m = cbind(1*sct.train$s*sct.train$w, covar.train*sct.train$s*sct.train$w)
        }else{
            covar.train.m = cbind(1*sct.train$s*sct.train$w, covar.train*sct.train$s*sct.train$w, batch.train * sct.train$w)
        }
        pred.res = weights.enet( genos = geno.train * sct.train$s * sct.train$w, pheno = as.matrix(sct.train$rc_w) , 
                                covar = covar.train.m, alpha=0.5, intercept=F)
        pred.wgt = pred.res$eff.wgt
        intercept.m = pred.res$intercept
        pred.wgt[ is.na(pred.wgt) ] = 0
        cv.calls[ indx , 1 ] = (genos$bed[ cv.sample[ indx ] , ] * pheno_sct[cv.sample[ indx ], ]$s )  %*% pred.wgt 
        cv.calls1[ indx , 1 ] = (genos$bed[ cv.sample[ indx ] , ] * pheno_sct[cv.sample[ indx ], ]$s * pheno_sct[cv.sample[ indx ], ]$w)  %*% pred.wgt 
        cv.all[indx, 3] = pheno_sct[cv.sample[ indx ], ]$rc_w - intercept.m  * pheno_sct[cv.sample[ indx ], ]$s * pheno_sct[cv.sample[ indx ], ]$w
        if(ncov>0){
        new_mean[indx] = intercept.m  + genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt + as.matrix(pheno_covar[cv.sample[ indx ],-(1:3)]) %*% coef(pred.res$enet,s = "lambda.min")[3:(3+ncov-1)]
            }else{
         new_mean[indx] = intercept.m  + genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt    
        }
    }		
    # compute rsq + P-value for each iteration
	if ( !is.na(sd(cv.calls[,1])) && sd(cv.calls[,1]) != 0 ) {
       reg = summary(lm( pheno_sct[cv.sample, ]$rc ~ pheno_sct[cv.sample, ]$s  + cv.calls[,1]  - 1, weights =  pheno_sct[cv.sample, ]$w^2))
        reg = summary(lm( pheno_sct[cv.sample, ]$rc - reg$coefficients[1,1] * pheno_sct[cv.sample, ]$s  ~ cv.calls[,1]  - 1, weights =  pheno_sct[cv.sample, ]$w^2))
		cv.performance[ 1, j ] = reg$adj.r.sq
		cv.performance[ 2, j ] = reg$coef[1,4]
	} else {
		cv.performance[ 1:2, j ] = NA
	}   
    
    set.seed(1)
    if(is.null(batch)){
        covar.m = cbind(1*pheno_sct$s*pheno_sct$w, pheno_covar[,-(1:3)]*pheno_sct$s*pheno_sct$w)
    }else{
        covar.m = cbind(1*pheno_sct$s*pheno_sct$w, pheno_covar[,-(1:3)]*pheno_sct$s*pheno_sct$w, batch * pheno_sct$w)
    }
    pred.wgt = weights.enet( genos = genos$bed * pheno_sct$s * pheno_sct$w, pheno = as.matrix(pheno_sct$rc_w) , 
                                covar = covar.m, alpha=0.5, intercept=F)
    wgt.matrix[,j] = pred.wgt$eff.wgt
    intercept[j] = pred.wgt$intercept

    new_mu <- pheno_sct[cv.sample, ]$s * new_mean
    new_var <- new_mu + new_mu^2 / pheno_sct$theta[1] 
    new_var[new_var < pheno_sct$min_var[1]] <- pheno_sct$min_var[1] 
    new_w <- 1 / sqrt(new_var)
    pheno_sct$w[cv.sample] <- new_w
    pheno_sct$rc_w <- pheno_sct$rc * pheno_sct$w

}
    

if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')



# save weights, rsq, p-value for each model, and hsq to output
snps = genos$bim
save( wgt.matrix , snps , cv.performance , hsq, hsq.pv, N.tot , intercept,file = paste( opt$out , ".wgt.RDat" , sep='' ) )
# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()

