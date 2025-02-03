source('./src/utilities.R')
prepEnv()
library(argparse)

parser <- ArgumentParser()
parser$add_argument('-i', type='character', nargs=1)
parser$add_argument('-c', type='character', nargs='+')
parser$add_argument('-s', type='character', nargs='+')
parser$add_argument('-o', type='character', nargs=1)

args  <- parser$parse_args()
input <- args$i
c     <- args$c
s     <- args$s
o     <- args$o

dds <- readRDS(input)
group <- as.character(design(dds))[-1]

basename    <- basename(file_path_sans_ext(input))
de          <- lapply(c, function(x) rownames(read.csv(glue('./data/DE_results/{basename}_{x}_DE.csv'), row.names=1, header=T)))
de_combined <- Reduce(union, de)
ann_row     <- lapply(de, function(x) as.numeric(de_combined %in% x)) %>% as.data.frame()
colnames(ann_row) <- c
rownames(ann_row) <- de_combined

vsd <- assays(dds)[['vsd']]

print(head(vsd[de_combined, vsd[[group]] %in% s]))
print(head(ann_row))
pdf(glue('./fig/heatmap/{basename}_{o}.pdf'), height = length(de_combined) / 6)
hm <- pheatmap(t(scale(t(assay(vsd[de_combined, vsd[[group]] %in% s])))), annotation_row=ann_row)
draw(hm)
dev.off()
