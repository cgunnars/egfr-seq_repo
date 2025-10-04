
source('./src/utilities.R')
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #path to dds object
parser$add_argument('-g', nargs=1, type='character') #path to gene list
parser$add_argument('-o', nargs=1, type='character') #where to save
parser$add_argument('-c', nargs='+', type='character') # conditions

dds      <- readRDS(args$i)
genefile <- args$g
conditions <- args$c

genelist <- readtxt(genefile)

hm <- plotBasicHeatmap(genelist, dds, conditions)

pdf(file=glue('{out}'), getHeight(degs), getWidth(conditions))
draw(hm)
dev.off()
