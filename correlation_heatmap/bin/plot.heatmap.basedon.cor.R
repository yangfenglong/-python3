suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option("--qsigr1", default=NULL, help="The corr.test qsignificant filtered r matrix,caculated from /NJPROJ2/MICRO/PROJ/yangfenglong/yanfa/pathview/lib/get_significant_correlations.R, default is *_qsigr.csv"),
  make_option("--corr1", default=NULL, help="The corr1 matrix ,default is *_r.csv"),
  make_option("--qsigp1", default=NULL, help="The corr.test qsignificant filtered p matrix ,default is *_qsigp.csv"),
  make_option("--t1", action="store_true", default=FALSE, help="trans the qsigr1 and qsigp1 table [default %default]"),
  make_option("--corr2", default=NULL, help="The corr2 matrix, default is  default is *_r.csv"),
  make_option("--corp2", default=NULL, help="The corp2 matrix,default is *_p.csv"),
  make_option("--t2", action="store_true", default=FALSE, help="trans the corr2 and corp2 table [default %default]"),
  make_option("--anno1", default=NULL, help="The annotation list,eg: sample\tgroup\tgroup2.. like"),
  make_option("--anno2", default=NULL, help="The annotation list,eg: sample\tgroup\tgroup2.. like"),
  make_option("--top", default=35, help="show top 35 species or kos [default %default]"),
  make_option("--top2", default=35, help="show top 35 metabolites [default %default]"),
  make_option("--cluster_rows", action="store_true", default=FALSE, help="whether make cluster on rows [default %default]"),
  make_option("--show_row_dend", action="store_true", default=FALSE, help="whether show row clusters [default %default]"),
  make_option("--cluster_columns", action="store_true", default=FALSE, help="whether make cluster on columns. [default %default]"),
  make_option("--show_column_dend", action="store_true", default=FALSE, help="whether show column clusters [default %default]"),
  make_option("--outdir", default="./", help="The output dirctory,  [default %default]"),
  make_option("--fontsize_row", default=12, help="The fontsize of row names,  [default %default]"),
  make_option("--fontsize_col", default=12, help="The fontsize of col names,  [default %default]"),
  make_option("--fontsize_grid", default=10, help="The fontsize of grid annotations,  [default %default]"),
  make_option("--row_names_side", default="left", help="should the row names be put on the left or right of the heatmap,  [default %default]"),
  make_option("--col_names_side", default="top", help="should the column names be put on the top or bottom of the heatmap, [default %default]"),
  make_option("--prefix", default="ko_phenotype", help="prefix of the output vs_file name"),
  make_option("--prefix2", default="ko_metabotype", help="prefix of the output vs2_file name")
)
opt<-parse_args(OptionParser(usage="%prog [options]\n", option_list=option_list))
prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
if(is.null(opt$qsigr1) || is.null(opt$qsigp1)){
  system(paste("/NJPROJ2/MICRO/PROJ/yangfenglong/software/miniconda3/lib/R/bin/Rscript ",prog," -h"))
  quit("no")
}
#################################################################################################

outdir <- paste(opt$outdir,"/",sep="")
if(!file.exists("opt$outdir")){dir.create(path=opt$outdir,recursive = TRUE)} 
qsigr <- read.csv(opt$qsigr1,header = T,row.names = 1)
mat <- read.csv(opt$qsigp1,header = T,row.names = 1)
corr1 <- read.csv(opt$corr1,header = T,row.names = 1)
if(opt$t1){
    qsigr <- t(qsigr)
    mat <- t(mat)
}
library(ComplexHeatmap)
library(circlize)
library(stringr)
#set annotation
if(!is.null( opt$anno1)){ha_row <- read.table(opt$anno1,sep="\t",row.names=1,head=T);ha_row<-rowAnnotation(ha_row)}else{ha_row<-NULL}
if(!is.null(opt$anno2)){
    ha_col <- read.table(opt$anno2,sep="\t",row.names=1,head=T)
	ha_col<- HeatmapAnnotation(ha_col)
  if(opt$col_names_side=="top"){
    bottom_annotation=ha_col
    top_annotation=new("HeatmapAnnotation")
  }else if(opt$col_names_side=="bottom"){
    top_annotation=ha_col
    bottom_annotation=new("HeatmapAnnotation")
  }
}else{
    top_annotation=new("HeatmapAnnotation")
	bottom_annotation=new("HeatmapAnnotation")
}  

if(!is.null(opt$top)){
    qsigr <- corr1[row.names(qsigr),colnames(qsigr)]
    qsigr <- qsigr[order(rowSums(abs(qsigr),na.rm=T),decreasing = T)[1:opt$top],]
    mat <- mat[rownames(qsigr),]
}else{
    qsigr <- corr1[row.names(qsigr),colnames(qsigr)] 
    mat <- mat[row.names(qsigr),colnames(qsigr)]
}


color <- colorRamp2(seq(min(qsigr,na.rm=T), max(qsigr,na.rm=T), length = 3), c("darkblue", "white", "darkred"), space = "RGB")
anno.mat <- matrix(ifelse(is.na(mat),"", ifelse(mat <= 0.01, "**", ifelse(mat <= 0.05,"*", ""))),nrow = nrow(mat))
namelen_row<-max(str_length(row.names(qsigr)))
namelen_col<-max(str_length(colnames(qsigr)))
ht1 = Heatmap(qsigr,  col = color, name=opt$prefix, 
    row_title = NA, column_title = NA, 
    row_names_side=opt$row_names_side, column_names_side= opt$col_names_side, 
    row_names_gp = gpar(fontsize = opt$fontsize_row),
    column_names_gp = gpar(fontsize = opt$fontsize_col),
    top_annotation = top_annotation, bottom_annotation=bottom_annotation,
	cluster_rows = opt$cluster_rows, cluster_columns=opt$cluster_columns,
    show_row_dend = opt$show_row_dend, show_column_dend = opt$show_column_dend,
    na_col="white",
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(anno.mat[i,j], x, y, gp = gpar(fontsize = opt$fontsize_grid))
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
    }
)
if(!is.null(ha_row)) {ht1<- ht1+ha_row}
svg(file=paste(outdir,opt$prefix,"_heatmap.svg",sep=""),width = ncol(mat)/5+2+namelen_row/30, height = nrow(mat)/5+2+namelen_col/30)
ht1
dev.off()

if(!is.null(opt$corr2) & !is.null(opt$corp2)){
    corr2 <- read.csv(opt$corr2,header = T,row.names = 1)
	mat2 <- read.csv(opt$corp2,header = T,row.names = 1)
	if(opt$t2){
	  corr2 <- t(corr2)
	  mat2 <- t(mat2)
	}
	if(!is.null(opt$top2)){
		corr2 <- corr2[row.names(qsigr),order(colSums(abs(corr2),na.rm=T),decreasing = T)[1:opt$top2]]
	    mat2 <- mat2[row.names(qsigr),colnames(corr2)]  
	}else{
	    corr2 <- corr2[row.names(qsigr),]
		mat2 <- mat2[row.names(qsigr),]
	}
    anno.mat2 <- matrix(ifelse(is.na(mat2),"", ifelse(mat2 <= 0.01, "**", ifelse(mat2 <= 0.05 ,"*", " "))),nrow = nrow(mat2))
    color2 <- colorRamp2(seq(min(corr2,na.rm=T), max(corr2,na.rm=T), length = 3), c("green", "white", "red"), space = "RGB")
    namelen2_row<-max(str_length(row.names(corr2)))
	namelen2_col<-max(str_length(colnames(corr2)))
    ht2 = Heatmap(corr2,col = color2, name = opt$prefix2,
        row_title = NA, column_title = NA,
        row_names_side=opt$row_names_side, column_names_side= opt$col_names_side,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = opt$fontsize_row),
        column_names_gp = gpar(fontsize = opt$fontsize_col),
	    top_annotation = top_annotation, bottom_annotation=bottom_annotation,  
		cluster_columns=opt$cluster_columns,
        show_row_dend = FALSE, show_column_dend = opt$show_column_dend,
        na_col="white",
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(anno.mat2[i,j], x, y, gp = gpar(fontsize = opt$fontsize_grid))
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        }
    )

    # combine heatmaps
    ifelse(is.null(ha_row),c(ht_list <- ht1 + ht2), c(ht_list<-ht1 + ht2 + ha_row)) # ifelse is vectorized and does not provide special cases for non-vectorized conditions
    svg(file=paste(outdir,opt$prefix,"_combine_heatmap.svg",sep=""),width = ncol(mat)/5+ncol(mat2)/5+2+max(namelen2_row,namelen_row)/30, height = nrow(mat)/5+2+max(namelen2_col,namelen_col)/30)
    #ht_list
    draw(ht_list, row_title = NULL, column_title = NULL, column_title_side = "top")
    dev.off()
}
