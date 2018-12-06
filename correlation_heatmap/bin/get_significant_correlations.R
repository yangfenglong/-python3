#!/NJPROJ1/MICRO/share/software/R-3.3.3/bin/R
############################################################
# 
# author: yangfenglong
# date: 2018年 06月 12日 星期二 17:10:54 CST
# descr: caculator r and p for 2 matrix correlations 
# 
############################################################
corr.get_rpsig <- function(
    mat1,mat2,
	method="spearman",padj="fdr",
	rcut=0.6,pcut=0.05,
	outdir="./",prefix="mt1_vs_mt2"
){
    library(psych)
    library(tidyr,lib.loc="/NJPROJ2/MICRO/PROJ/yangfenglong/software/Rlib")
    library(dplyr,lib.loc="/NJPROJ2/MICRO/PROJ/yangfenglong/software/Rlib")
    outdir <- paste(outdir,"/",sep="")
    
    cor <- corr.test(mat1,mat2,use="pairwise",method=method,adjust=padj,alpha=.05)
    
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
    select_corr <- filter(total_data, abs(r)>=rcut & p<pcut)
    qsigr <- select_corr[,c(1:3)] %>% spread(genus2,r)
	qsig_r <- corr[as.character(qsigr[,1]),colnames(qsigr)[-1]] 
    qsigp <- corp[as.character(qsigr[,1]),colnames(qsigr)[-1]]
    colnames(qsigr)[1]<-"" 
    write.csv(cor$r,file=paste(outdir,prefix,"_r.csv",sep=""),quote = F,row.names=T)
    write.csv(cor$p,file=paste(outdir,prefix,"_p.csv",sep=""),quote = F,row.names=T)
    write.csv(qsigr,file=paste(outdir,prefix,"_qsigr.csv",sep=""),quote = F,row.names=F)
    cor$qsigr <- qsig_r
	cor$qsigp <- qsigp
	
	return(cor)
}
