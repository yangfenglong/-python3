suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option("--mat1", default=NULL, help="The relative abundance ,default is ko_diff_relative.xls"),
  make_option("--mat2", default=NULL, help="The relative abundance ,default is phenotype.xls"),
  make_option("--mat3", default=NULL, help="The relative abundance ,default is metabotype.xls,mat1-vs-mat3"),
  make_option("--t1", action="store_true", default=FALSE, help="trans the mat1 table [default %default]"),
  make_option("--t2", action="store_true", default=FALSE, help="trans the mat2 table [default %default]"),
  make_option("--t3", action="store_true", default=FALSE, help="trans the mat3 table [default %default]"),
  make_option("--outdir", default="./", help="The output dirctory,  [default %default]"),
  make_option("--fontsize_row", default=12, help="The fontsize of row names,  [default %default]"),
  make_option("--fontsize_col", default=12, help="The fontsize of col names,  [default %default]"),
  make_option("--fontsize_grid", default=10, help="The fontsize of grid annotations,  [default %default]"),
  make_option("--row_names_side", default="left", help="should the row names be put on the left or right of the heatmap,  [default %default]"),
  make_option("--col_names_side", default="top", help="should the column names be put on the top or bottom of the heatmap, [default %default]"),
  make_option("--method",  default="spearman", help="the method of calculate correlation,The alternatives  are 'pearson' and 'kendall' [default %default]"),
  make_option("--r", default=0.6, help="The threshold value of correlation coefficent [default %default]"),
  make_option("--p", default=0.05, help="The threshold value of adjusted pvalues [default %default]"),
  make_option("--prefix", default="ko_phenotype", help="prefix of the output vs_file name"),
  make_option("--prefix2", default="ko_metabotype", help="prefix of the output vs2_file name"),
  make_option("--adjust", default="fdr", help="The method of adjustment for multiple tests,The alternatives  are 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'none' [default %default]")
)
opt<-parse_args(OptionParser(usage="%prog [options]\n", option_list=option_list))
prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
if(is.null(opt$mat1)){
  system(paste("/NJPROJ2/MICRO/PROJ/yangfenglong/software/miniconda3/lib/R/bin/Rscript ",prog," -h"))
  quit("no")
}
#################################################################################################

outdir <- paste(opt$outdir,"/",sep="")
if(!file.exists("opt$outdir")){dir.create(path=opt$outdir,recursive = TRUE)} 
ko <- read.table(opt$mat1,header = T,row.names = 1,sep="\t")
if(opt$t1) ko <- t(ko)
ko<- ko[,colSums(ko)!=0]
phe <- read.table(opt$mat2,header = T,row.names = 1, sep="\t")
if(opt$t2) phe <- t(phe)
phe <- phe[,colSums(phe)!=0]
source("/NJPROJ2/MICRO/PROJ/yangfenglong/yanfa/dheatmap/bin/get_significant_correlations.R")
cor<- corr.get_rpsig(ko,phe,
    method=opt$method,padj=opt$adjust,
    rcut=opt$r,pcut=opt$p,
    outdir=opt$outdir,prefix=opt$prefix)
#cor is a list with several tables
#cor$r       cor$t       cor$se      cor$adjust  cor$ci      cor$Call    cor$qsigp   
#cor$n       cor$p       cor$sef     cor$sym     cor$ci.adj  cor$qsigr
mat <-  cor$qsigp
library(ComplexHeatmap)
library(circlize)
color <- colorRamp2(seq(min(cor$qsigr), max(cor$qsigr), length = 3), c("darkblue", "white", "darkred"), space = "RGB")
anno.mat <- matrix(ifelse(mat <= 0.01, "**", ifelse(mat <= 0.05,"*"," ")),nrow = nrow(mat))
ht1 = Heatmap(cor$qsigr,  col = color, name=opt$prefix, 
	row_title = NA, column_title = NA, 
	row_names_side=opt$row_names_side, column_names_side= opt$col_names_side, 
	row_names_gp = gpar(fontsize = opt$fontsize_row),
	column_names_gp = gpar(fontsize = opt$fontsize_col),
	show_row_dend = FALSE, show_column_dend = FALSE,
	na_col="white",
	cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(anno.mat[i,j], x, y, gp = gpar(fontsize = opt$fontsize_grid))
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
    }
)
svg(file=paste(outdir,opt$prefix,"heatmap.svg",sep=""))
ht1
dev.off()

if(!is.null(opt$mat3)){
    metabo <- read.table(opt$mat3,header = T,row.names = 1,sep="\t")
	if(opt$t3) metabo <- t(metabo)
	metabo <- metabo[,colSums(metabo)!=0]
    cor2<- corr.get_rpsig(ko,metabo,
        method=opt$method,padj=opt$adjust,
        rcut=opt$r,pcut=opt$p,
        outdir=opt$outdir,prefix=opt$prefix2)
    mat2 <- cor2$p[row.names(cor$qsigr),colSums(cor2$r)!=0] 
    anno.mat2 <- matrix(ifelse(mat2 <= 0.01, "**", ifelse(mat2 <= 0.05 ,"*"," ")),nrow = nrow(mat2))
	color2 <- colorRamp2(seq(min(cor$qsigr), max(cor$qsigr), length = 3), c("green", "white", "red"), space = "RGB")
    qsigr2 <- cor2$r[row.names(cor$qsigr),colSums(cor2$r)!=0]
	ht2 = Heatmap(qsigr2,col = color2, name = opt$prefix2,
        row_title = NA, column_title = NA,
        row_names_side=opt$row_names_side, column_names_side= opt$col_names_side,
		show_row_names = FALSE,
		row_names_gp = gpar(fontsize = opt$fontsize_row),
		column_names_gp = gpar(fontsize = opt$fontsize_col),
    	show_row_dend = FALSE, show_column_dend = FALSE,
        na_col="white",
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(anno.mat2[i,j], x, y, gp = gpar(fontsize = opt$fontsize_grid))
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        }
    )

    # combine heatmaps
    mat2
    ht_list = ht1 + ht2
    svg(file=paste(outdir,opt$prefix,"_combine_heatmap.svg",sep=""))
    #ht_list
	draw(ht_list, row_title = NULL, column_title = NULL, column_title_side = "top")
    dev.off()
}
