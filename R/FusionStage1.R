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
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--scale_pheno", action="store_true", default=TRUE,
              help="Scale phenotype or not [default: %default]"),			  

)

opt = parse_args(OptionParser(option_list=option_list))
hsq = 0.1
hsq.pv = NA
models = unique( c(unlist(strsplit(opt$models,','))) )
M = length(models)


if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
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

weights.glmnet = function( genos , pheno , seq_depth, covar = NULL, intercept= T, alpha=0.5 ) {
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
    glmnet.control(mxitnr=500)
	enet = cv.glmnet( x=as.matrix(X) , y=pheno , alpha=alpha , nfold=5, intercept=intercept , standardize=F, penalty.factor = penalty.v, family = "poisson", offset=log(seq_depth))
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

# ---

fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

# Make/fetch the phenotype file
if ( !is.na(opt$pheno) ) {
	pheno.file = opt$pheno
	pheno = read.table(pheno.file,as.is=T)
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	m = m[m.keep]
	pheno = pheno[m,]
} else {
	pheno.file = paste(opt$tmp,".pheno",sep='')
	pheno = fam[,c(1,2,6)]
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}


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
    pheno_covar[,-(1:3)] <- scale(pheno_covar[,-(1:3)])
}

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


if(ncol(pheno_covar)>3){
    covar.file = paste(opt$tmp,".covar",sep='')
#    write.table(pheno_covar_tmp[,-3],quote=F,row.names=F,col.names=F,file=covar.file)
}

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
if(opt$scale_pheno){
    pheno[,3] = scale(pheno[,3])
}
m = match( paste(pheno[,1],pheno[,2]) , paste(pheno_covar[,1],pheno_covar[,2]) )
if(ncol(pheno_covar)>3){
    pheno_covar = pheno_covar[m,]
    batch = batch[m,]
}else{
    pheno_covar = batch = NULL
}

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}


N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- CROSSVALIDATION ANALYSES
set.seed(1) 
cv.performance = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = models

if ( opt$crossval <= 1 ) {
if ( opt$verbose >= 1 ) cat("### Skipping cross-validation\n")
} else {
if ( opt$verbose >= 1 ) cat("### Performing",opt$crossval,"fold cross-validation\n")
cv.all = pheno
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=opt$crossval,labels=FALSE)

cv.calls = matrix(NA,nrow=N,ncol=M)

for ( i in 1:opt$crossval ) {
	if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
	indx = which(folds==i,arr.ind=TRUE)
	cv.train = cv.all[-indx,]
    geno.train = genos$bed[ cv.sample[ -indx ],]
    nasnps.train = apply( is.na(scale(geno.train)) , 2 , sum )
    if ( sum(nasnps.train) != 0 ) {
        geno.train[,nasnps.train != 0] = 0
    }
    covar.train = pheno_covar[cv.sample[ -indx ], -(1:3) ]
    batch.train = batch[cv.sample[ -indx ], ]
	for ( mod in 1:M ) {
        if(is.null(batch.train)){
            covar.train.m = covar.train
        }else{
            covar.train.m = cbind(covar.train,batch.train)
        }
        pred.wgt = weights.enet( genos = geno.train , pheno = as.matrix(cv.train[,3]) , 
                                covar = covar.train.m, alpha=0.5, intercept=T)$eff.wgt
        pred.wgt[ is.na(pred.wgt) ] = 0
        cv.calls[ indx , mod ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
	}
}

# compute rsq + P-value for each model
for ( mod in 1:M ) {
	if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) {
		reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] ))
		cv.performance[ 1, mod ] = reg$adj.r.sq
		cv.performance[ 2, mod ] = reg$coef[2,4]
	} else {
		cv.performance[ 1, mod ] = NA
		cv.performance[ 2, mod ] = NA
	}
}
if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')
}

# --- FULL ANALYSES
if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

# call models to get weights
wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[,2]
intercept = NULL
for ( mod in 1:M ) {
    if(is.null(batch)){
        covar.m = pheno_covar[,-(1:3)]
    }else{
        covar.m = cbind(pheno_covar[,-(1:3)],batch)
    }
    pred.wgt = weights.enet( genos = genos$bed , pheno = as.matrix(pheno[,3]) , 
                                covar = covar.m, alpha=0.5, intercept=T)$eff.wgt
    wgt.matrix[,mod] = pred.wgt
    
}

# save weights, rsq, p-value for each model, and hsq to output
snps = genos$bim
save( wgt.matrix , snps , cv.performance , hsq, hsq.pv, N.tot , intercept,file = paste( opt$out , ".wgt.RDat" , sep='' ) )
# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()

