suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics (must have SNP and Z column headers) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--weights", action="store", default=NA, type='character',
              help="File listing molecular weight RDat files (must have columns WGT,ID,CHR,P0,P1) [required]"),
  make_option("--weights_dir", action="store", default=NA, type='character',
              help="Path to directory where weight files (WGT column) are stored [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--force_model", action="store", default=NA, type='character',
              help="Force specific predictive model to be used, no flag (default) means select most significant cross-val. Options: blup,lasso,top1,enet"),
  make_option("--caviar", action="store_true", default=FALSE,
              help="Generate eCAVIAR-format (Z,LD) files for fine-mapping [default off]"),
  make_option("--jlim", action="store_true", default=FALSE,
              help="NOT IMPLEMENTED: Compute JLIM statistic [Chun et al Nat Genet 2017].\nRequires jlimR library installed. [default: %default]"),			  
  make_option("--max_impute", action="store", default=0.5 , type='double',
              help="Maximum fraction of SNPs allowed to be missing per gene (will be imputed using LD). [default: %default]"),			  
  make_option("--min_r2pred", action="store", default=0.7 , type='double',
              help="Minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores. [default: %default]"),			  
  make_option("--perm", action="store", default=0, type='integer',
              help="Maximum number of permutations to perform for each feature test [default: 0/off]"),
  make_option("--perm_minp", action="store", default=0.05, type='double',
              help="Minimum p-value for which to initiate permutation test, if --perm flag present [default: %default]"),    		  
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
  make_option("--coloc_P", action="store", default=NA, type='double',
              help="P-value below which to compute COLOC statistic [Giambartolomei et al PLoS Genet 2013]\nRequires coloc library installed and --GWASN flag. [default NA/off]"),
  make_option("--GWASN", action="store", default=NA, type='integer',
              help="Total GWAS/sumstats sample size for inference of standard GWAS effect size."),
  make_option("--PANELN", action="store", default=NA, type='character',
              help="File listing sample size for each panel for inference of standard QTL effect size, cross-referenced against 'PANEL' column in weights file")      
)

opt = parse_args(OptionParser(option_list=option_list))

allele.qc = function(a1,a2,ref1,ref2) {
        a1 = toupper(a1)
        a2 = toupper(a2)
        ref1 = toupper(ref1)
        ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

# Load in summary stats
sumstat = read.table(opt$sumstats,head=T,as.is=T)

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = read.table(opt$weights,head=T,as.is=T)
wgtlist = wgtlist[ as.character(wgtlist$CHR) == as.character(opt$chr) , ]
chr = unique(wgtlist$CHR)

N = nrow(wgtlist)
out.tbl = data.frame( "PANEL" = rep(NA,N) , "FILE" = character(N) , "ID" = character(N) , "CHR" = numeric(N) , "P0" = character(N) , "P1" = character(N) ,"HSQ" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) , "EQTL.ID" = character(N) , "EQTL.R2" = numeric(N) , "EQTL.Z" = numeric(N) , "EQTL.GWAS.Z" = numeric(N) , "NSNP" = numeric(N) , 
                     "NWGT" = numeric(N) , "MODELCV.R2" = character(N) , "MODELCV.PV" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , 
                     "NWGT_max" = numeric(N) , "MODELCV.R2_max" = character(N) , "MODELCV.PV_max" = character(N) , "TWAS.Z_max" = numeric(N) , "TWAS.P_max" = numeric(N) , 
                     "NWGT_converge" = numeric(N) , "MODELCV.R2_converge" = character(N) , "MODELCV.PV_converge" = character(N) , "TWAS.Z_converge" = numeric(N) , "TWAS.P_converge" = numeric(N) , converge = numeric(N),
                     stringsAsFactors=FALSE )

# Load in reference data
genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")
sumstat = sumstat[sumstat$SNP %in% genos$bim$V2,]
genos$bim$V7 = ifelse(genos$bim$V5<genos$bim$V6,paste0(genos$bim$V2,':',genos$bim$V5,':',genos$bim$V6),
                      paste0(genos$bim$V2,':',genos$bim$V6,':',genos$bim$V5))
sumstat$V7 = ifelse(sumstat$A1<sumstat$A2,paste0(sumstat$SNP,':',sumstat$A1,':',sumstat$A2),
                      paste0(sumstat$SNP,':',sumstat$A2,':',sumstat$A1))
# Match summary data to input, record NA where summary data is missing
m = match( genos$bim$V7 , sumstat$V7 )
sum.missing = is.na(m)
sumstat = sumstat[m,]
sumstat$SNP = genos$bim[,2]
sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
	genos$bim = genos$bim[qc$keep,]
	genos$bed = genos$bed[,qc$keep]
	sumstat = sumstat[qc$keep,]
}

# TODO: WARNING if too many NAs in summary stats

FAIL.ctr = 0

## For each wgt file:
for ( w in 1:nrow(wgtlist) ) {
	#cat( unlist(wgtlist[w,]) , '\n' )
	# Load weights
	wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
	load(wgt.file)
    colnames(cv.performance) = colnames(wgt.matrix)
	# Remove NAs (these should not be here)
	wgt.matrix[is.na(wgt.matrix)] = 0
    
    snps$V7 = ifelse(snps$V5<snps$V6,paste0(snps$V2,':',snps$V5,':',snps$V6),
                      paste0(snps$V2,':',snps$V6,':',snps$V5))

    if(is.na(abs(cv.performance[1,10] - cv.performance[1,9])<=1e-3)){
        converge = -9
    }else if(abs(cv.performance[1,10] - cv.performance[1,9])<=1e-3){
        converge = 1
    }else{
        converge = 0
    }
	
	# Match up the SNPs and weights
	m = match( snps$V7 , genos$bim$V7 )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	wgt.matrix = wgt.matrix[m.keep,,drop=F]
	cur.genos = scale(genos$bed[,m[m.keep]])
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

	cur.FAIL = FALSE

	# Match up the SNPs and the summary stats
	m = match(cur.bim$V7 , sumstat$V7)
	cur.Z = sumstat$Z[m]

	# which rows have rsq
	row.rsq = which(rownames(cv.performance)=='rsq')
	# which rows have p-values
	row.pval = which(rownames(cv.performance)=='pval')
	
    if ( all(is.na(cv.performance[1,])) ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
		FAIL.ctr = FAIL.ctr + 1
		next
	}
    
    model_list = c('max','first','converge')
	# Identify the best model


    IMPUTE = FALSE
    # Compute LD matrix
        if ( length(cur.Z) == 0 ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
            cur.FAIL = TRUE
            out.tbl$NSNP[w] = NA
        } else if ( !cur.FAIL ) {
            cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	
            out.tbl$NSNP[w] = nrow(cur.LD)
            cur.miss = is.na(cur.Z)
            # Impute missing Z-scores
            if ( sum(cur.miss) != 0 ) {
                if ( sum(!cur.miss) == 0 ) {
                    cat( "WARNING : " , unlist(wgtlist[w,]) , "had no overlapping GWAS Z-scores, skipping this gene\n")
                    cur.FAIL = TRUE
                } else if ( mean(cur.miss) > opt$max_impute ) {
                    cat( "WARNING : " , unlist(wgtlist[w,]) , "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
                    cur.FAIL = TRUE
                } else {
                    IMPUTE = TRUE
                    cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                    cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                    cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                    cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

                    all.r2pred = rep(1,length(cur.Z))
                    all.r2pred[ cur.miss ] = cur.r2pred
                    if ( sum(is.na(all.r2pred)) != 0 ) {
                        cat( "WARNING : " , unlist(wgtlist[w,]) , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
                        cur.FAIL = TRUE
                    } 
                }
            }
        }
        
    for(mod in model_list){
        if(mod=='first'){
            mod.best = 1
        }else if(mod=='converge'){
            mod.best = 10
        }else if(mod=='max'){
            mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
        }
        if(is.na(mod.best)){
            cur.FAIL = TRUE
        }
        if(!cur.FAIL){
            if(is.na(cv.performance[1,mod.best])){
                cur.FAIL = TRUE
            }
        if ( sum(wgt.matrix[, mod.best] != 0) == 0 ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , names(cv.performance)[ mod.best ] , "had", length(cur.Z) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
            cur.FAIL = TRUE
            }
        if(!cur.FAIL & IMPUTE){
            if ( mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) < opt$min_r2pred ) {
                cat( "WARNING : " , unlist(wgtlist[w,]) , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) , "at expression weight SNPs, skipping this gene.\n")
                cur.FAIL = TRUE
            }
        }

        }




		
		if ( !cur.FAIL ) {
			# Compute TWAS Z-score
			cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
			cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
					
			if ( cur.twasr2pred > 0 ) {
				cur.twas = cur.twasz / sqrt(cur.twasr2pred)
			} else {
				cur.FAIL=T
				cat( "WARNING : " , unlist(wgtlist[w,]) , " had zero predictive accuracy, try a different model.\n")
			}
		}
	

	# populate the output
	if ( sum(names(wgtlist) == "PANEL") == 1 ) out.tbl$PANEL[w] = wgtlist$PANEL[w]
	out.tbl$FILE[w] = wgt.file
	out.tbl$CHR[w] = wgtlist$CHR[w]
	out.tbl$P0[w] = wgtlist$P0[w]
	out.tbl$P1[w] = wgtlist$P1[w]
	out.tbl$ID[w] = wgtlist$ID[w]
	if ( exists("hsq") ) {
		out.tbl$HSQ[w] = hsq[1]
	}
        
	eqtlmod = colnames(wgt.matrix) == "top1"
	topeqtl = which.max( wgt.matrix[,eqtlmod]^2 )
	
	if ( cur.FAIL || sum(eqtlmod) == 0 || length(topeqtl) == 0 || is.na(topeqtl) ) {
		out.tbl$EQTL.ID[w] = NA
		out.tbl$EQTL.R2[w] = NA
		out.tbl$EQTL.Z[w] =  NA
		out.tbl$EQTL.GWAS.Z[w] = NA
	} else {
		out.tbl$EQTL.ID[w] = rownames(wgt.matrix)[topeqtl]
		out.tbl$EQTL.R2[w] = cv.performance[1,eqtlmod]
		out.tbl$EQTL.Z[w] = wgt.matrix[ topeqtl , eqtlmod ]
		out.tbl$EQTL.GWAS.Z[w] = cur.Z[ topeqtl ]
	}
	
	topgwas = which.max( cur.Z^2 )
	if ( !cur.FAIL && length(topgwas) != 0 && !is.na(topgwas) ) {
		out.tbl$BEST.GWAS.ID[w] = snps[ topgwas , 2 ]
		out.tbl$BEST.GWAS.Z[w] = cur.Z[ topgwas ]
	} else {
		out.tbl$BEST.GWAS.ID[w] = NA
		out.tbl$BEST.GWAS.Z[w] = NA
	}

    if(mod=='first'){
        out.tbl$MODELCV.R2[w] = paste(format(cv.performance[row.rsq,mod.best],digits=2,trim=T),collapse=',')
        out.tbl$MODELCV.PV[w] = paste(format(cv.performance[row.pval,mod.best],digits=2,trim=T),collapse=',')
        if ( !cur.FAIL ) {
            out.tbl$NWGT[w] = sum( wgt.matrix[,mod.best] != 0 )
            out.tbl$TWAS.Z[w] = cur.twas
            out.tbl$TWAS.P[w] = 2*(pnorm( abs(out.tbl$TWAS.Z[w]) , lower.tail=F))
        } else {
            out.tbl$TWAS.Z[w] = NA
            out.tbl$TWAS.P[w] = NA
        }
    }
    if(mod=='max'){
        out.tbl$MODELCV.R2_max[w] = paste(format(unique(cv.performance[row.rsq,mod.best]),digits=2,trim=T),collapse=',')
        out.tbl$MODELCV.PV_max[w] = paste(format(unique(cv.performance[row.pval,mod.best]),digits=2,trim=T),collapse=',')
        if ( !cur.FAIL ) {
            out.tbl$NWGT_max[w] = sum( wgt.matrix[,mod.best] != 0 )
            out.tbl$TWAS.Z_max[w] = cur.twas
            out.tbl$TWAS.P_max[w] = 2*(pnorm( abs(out.tbl$TWAS.Z_max[w]) , lower.tail=F))
        } else {
            out.tbl$TWAS.Z_max[w] = NA
            out.tbl$TWAS.P_max[w] = NA
        }
    }
    if(mod=='converge'){
        out.tbl$MODELCV.R2_converge[w] = paste(format(unique(cv.performance[row.rsq,mod.best]),digits=2,trim=T),collapse=',')
        out.tbl$MODELCV.PV_converge[w] = paste(format(unique(cv.performance[row.pval,mod.best]),digits=2,trim=T),collapse=',')
        out.tbl$converge[w]= converge
        if ( !cur.FAIL ) {
            out.tbl$NWGT_converge[w] = sum( wgt.matrix[,mod.best] != 0 )
            out.tbl$TWAS.Z_converge[w] = cur.twas
            out.tbl$TWAS.P_converge[w] = 2*(pnorm( abs(out.tbl$TWAS.Z_converge[w]) , lower.tail=F))
        } else {
            out.tbl$TWAS.Z_converge[w] = NA
            out.tbl$TWAS.P_converge[w] = NA
        }
    }


	if ( cur.FAIL ) FAIL.ctr = FAIL.ctr + 1
    # reset cur.FAIL per mod
    cur.FAIL = FALSE 
    }
}

cat("Analysis completed.\n")
cat("NOTE:",FAIL.ctr,"/",nrow(wgtlist),"genes were skipped\n")
if ( FAIL.ctr / nrow(wgtlist) > 0.1 ) {
cat("If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)\n")
cat("Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [https://github.com/bogdanlab/fizi]\n")
}

# WRITE MHC TO SEPARATE FILE
mhc = as.numeric(out.tbl$CHR) == 6 & as.numeric(out.tbl$P0) > 26e6 & as.numeric(out.tbl$P1) < 34e6

out.tbl$P0 = apply( as.matrix(out.tbl$P0) , 1 , toString )
out.tbl$P1 = apply( as.matrix(out.tbl$P1) , 1 , toString )

if ( sum( mhc ) > 0 ) {
	cat("Results in the MHC are written to",paste(opt$out,".MHC",sep=''),", evaluate with caution due to complex LD structure\n")
	write.table( format( out.tbl[mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=paste(opt$out,".MHC",sep='') )
}
write.table( format( out.tbl[!mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=opt$out )
