source('./src/utilities.R')
library(argparse)
library(dplyr)
library(ggplot2)
library(glue)
library(ggrepel)

parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character')
parser$add_argument('-d', nargs=1, type='character')
parser$add_argument('-o', nargs=1, type='character')
args <- parser$parse_args()

i <- read.csv(args$i, row.names=1, header=T) 
d <- read.csv(args$d, row.names=1, header=T)
o <- args$o


degs   <- i[abs(i$log2FoldChange) > 1 & i$padj < 0.05 & !is.na(i$padj), ]
degs_2 <- d[abs(d$log2FoldChange) > 1 & d$padj < 0.05 & !is.na(d$padj), ]

print(i['ctpC', ])
print(d['ctpC', ])
common <- Reduce(intersect, list(rownames(i), rownames(d)))

df <- cbind(
            i[common,'log2FoldChange'], d[common, 'log2FoldChange'], 
            common %in% rownames(degs), common %in% rownames(degs_2)) %>% as.data.frame()
rownames(df) <- common
colnames(df) <- c('fc1', 'fc2', 'is_de', 'is_de2')
df['Label']  <- common
df[!(common %in% rownames(degs)),'Label'] <- ""
ggplot(df, aes(x=fc2, y=fc1, label=Label)) + 
                                geom_point(aes(alpha=is_de, col=is_de2), size=1) +
                                geom_hline(yintercept=0.7, linetype='dotted') +
                                geom_hline(yintercept=-0.7, linetype='dotted') +
                                geom_vline(xintercept=0.7, linetype='dotted') +
                                geom_vline(xintercept=-0.7, linetype='dotted') +
                                geom_text_repel(min.segment.length = 0, max.overlaps=100) + 
                                theme_classic()

ggsave(glue('./fig/biplot/{o}.pdf'))
