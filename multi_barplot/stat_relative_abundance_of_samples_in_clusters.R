args <- commandArgs(T)
if(length(args) < 2){
    cat("usage: <input: clusters.csv> <output: relative.xls>\n")
    cat ("Example: /NJPROJ2/MICRO/PROJ/yangfenglong/software/miniconda3/lib/R/bin/Rscript /NJPROJ2/MICRO/PROJ/yangfenglong/yanfa/barplot/stat_relative_abundance_of_samples_in_clusters.R MultiSamplesClusterDisBar/clusters.csv relative.xls\n")
    quit("no")
}

library(tidyverse)      
a <-read.csv(args[1],head=T)
a$Barcode <- str_split(a[,1],'-',simplify = TRUE)[,2]

# stat counts of each sample in each cluster
a %>% 
    group_by(Cluster) %>% 
	count(Barcode) %>% 
	spread(key=Barcode,value=n) -> b

df <- column_to_rownames(b,colnames(b)[1])

#caculate the relative percent of each sample in each cluster
rwsum <- rowSums(df,na.rm=T)
df_r <- df/rwsum

#print out the relative table with heads
tb <- as.tibble(df_r) %>%
    add_column(.,Samples=rownames(df_r), .before = 1)

write_tsv(tb,args[2])
