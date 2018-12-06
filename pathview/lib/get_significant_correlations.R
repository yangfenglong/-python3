#!/NJPROJ1/MICRO/share/software/R-3.3.3/bin/R
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option("--mat1", help="The relative abundance, default is ko_diff_relative.xls"),
  make_option("--mat2", help="The relative abundance, default is phenotype.xls"),
  make_option("--t1", action="store_true", default=FALSE, help="trans the mat1 table [default %default]"),
  make_option("--t2", action="store_true", default=FALSE, help="trans the mat2 table [default %default]"),
  make_option(c("-m","--method"), default="spearman", help="Alternatives: 'pearson' and 'kendall' [default %default]"), 
  make_option("--r", default=0.6, help="The threshold value of correlation coefficent [default %default]"), 
  make_option("--p", default=0.05, help="The threshold value of adjusted pvalues [default %default]"),
  make_option("--outdir", default="./", help="The output file dir [default %default]"),
  make_option("--prefix", default="ko_phenotype", help="prefix of the output corr_file name"),
  make_option("--adjust", default="fdr", help="Adjustment for multiple tests, Alternatives: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'none' [default %default]")
)

opt <- parse_args(
  OptionParser(
	option_list=option_list,
    usage="/NJPROJ1/MICRO/share/software/R-3.3.3/bin/Rscript %prog [options] \n",
	description = "Description\n
	caculate the correlations for 2 matrix use corr.test {psych} and reports probability values.
	select and write the significant correlation table ( r >= r_threshold and p <= p_threshold)", 
	epilogue = "Example: 
	/NJPROJ1/MICRO/share/software/R-3.3.3/bin/Rscript %prog \\
	--mat1 ko_diff_relative.xls --tans --mat2 phenotype.xls \\
	-m spearman --r 0.6 --p 0.05 --prefix ko_pheno

	Output files are: ko_pheno_qsigr.csv, ko_pheno_qsigp.csv
	\nAuthor(s): Yangfenglong\n"
))

prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
if(is.null(opt$mat1) || is.null(opt$mat2)){
  system(paste("Rscript ",prog," -h"))
  quit("no")
}

####################################################################################################
library(psych)
library(tidyr,lib.loc="/NJPROJ2/MICRO/PROJ/yangfenglong/software/Rlib")
library(dplyr,lib.loc="/NJPROJ2/MICRO/PROJ/yangfenglong/software/Rlib")
outdir <- paste(opt$outdir,"/",sep="")

mat1 <- read.table(opt$mat1,header = T,row.names = 1,sep="\t")
mat2 <- read.table(opt$mat2,header = T,row.names = 1,sep="\t")
if(opt$t1) mat1 <- t(mat1)
if(opt$t2) mat2 <- t(mat2)
cor <- corr.test(mat1,mat2,use="pairwise",method=opt$method,adjust=opt$adjust,alpha=.05)

#rvalue
corr <- as.data.frame(cor$r)
corr$genus1 <- rownames(corr)
total_corr <- corr %>% gather(genus2, r, -genus1)

#pvalue
corp <- as.data.frame(cor$p)
corp$genus1 <- rownames(corp)
total_corp <- corp %>% gather(genus2, p, -genus1)

#filter r and p according to threshold
total_data <- data.frame(genus2=total_corr$genus2,genus1=total_corr$genus1,r=total_corr$r,p=total_corp$p)
select_corr <- filter(total_data, abs(r)>=opt$r & p<opt$p)
qsigr <- select_corr[,c(1:3)] %>% spread(genus2,r)
qsigp <- select_corr[,c(1,2,4)] %>% spread(genus2,p)
colnames(qsigr)[1] <- "" 
colnames(qsigp)[1] <- ""

#write.table
write.csv(cor$r, file=paste(outdir, opt$prefix,"_r.csv",sep=""),quote = F,row.names=T)
write.csv(cor$p, file=paste(outdir, opt$prefix,"_p.csv",sep=""),quote = F,row.names=T)
write.csv(qsigr, file=paste(outdir, opt$prefix, "_qsigr.csv", sep=""), quote=F, row.names=F)
write.csv(qsigp, file=paste(outdir, opt$prefix, "_qsigp.csv", sep=""), quote=F, row.names=F)
