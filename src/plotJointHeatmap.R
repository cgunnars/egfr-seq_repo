source('./src/utilities.R')
library(argparse)
prepEnv()
library(RColorBrewer)
parser <- ArgumentParser()

parser$add_argument('-i', type='character', nargs=1)
parser$add_argument('-c', type='character', nargs='+')
parser$add_argument('-s', type='character', nargs='+')
parser$add_argument('-d', type='character', nargs='+')

parser$add_argument('-o', type='character', nargs=1)

args  <- parser$parse_args()
input <- args$i
c     <- args$c
d     <- args$d #other results to plot as annotation variable
s     <- args$s
o     <- args$o

dds <- readRDS(input)
group <- as.character(design(dds))[-1]

basename    <- basename(file_path_sans_ext(input))

de          <- lapply(c, function(x) rownames(read.csv(glue('./data/DE_results/{basename}_{x}_DE.csv'), row.names=1, header=T)))
de_combined <- Reduce(union, de)


de_other     <- lapply(d, function(x) read.csv(x, row.names=1, header=T))


de_other    <- lapply(de_other, function(x) x[ ,'log2FoldChange']) %>% as.data.frame(row.names=Reduce(union, lapply(de_other, function(x) rownames(x))))
colnames(de_other) <- d
ann_row <- de_other[de_combined, , drop=F]

vsd         <- assays(dds)[['vsd']]


## TODO FIX ANNOTATION COLORS
pdf(glue('./fig/heatmap/{o}.pdf'), height=length(de_combined)/6)
hm <- pheatmap(t(scale(t(assay(vsd[de_combined, vsd[[group]] %in% s])))), annotation_row=ann_row)
               
draw(hm)
dev.off()
