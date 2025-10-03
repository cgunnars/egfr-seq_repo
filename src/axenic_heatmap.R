
## for a given combined DEG table and two plot labels (lab1, lab2),
#  generate plots of FCs against each other (sfig6e)
plotDEGs_comparison <- function(combined, lab1, lab2, outname, fcthresh_low=0.7) {
    
    biplot <- ggplot(data=combined, aes(x=FC_2, y=FC_1, col=DE_1)) + 
               geom_point(alpha=0.5, aes(shape=notDE_2)) + 
               geom_hline(yintercept=fcthresh_low, linetype='dashed') + 
               geom_vline(xintercept=fcthresh_low, linetype='dashed') + 
               geom_hline(yintercept=-fcthresh_low, linetype='dashed') + 
               geom_vline(xintercept=-fcthresh_low, linetype='dashed') +
               labs(y=glue('log2FoldChange_{lab1}'), x=glue('log2FoldChange_{lab2}'),
                    col=glue('{lab1} DEG'), shape=glue('unlikely {lab2} DEG')) + theme_classic() 
    
    DEG_bycat <- list(rownames(combined[combined$DE_1 & combined$notDE_2, ]),
                      rownames(combined[combined$DE_1 & !combined$DE_2, ]),
                      rownames(combined[combined$DE_1 & combined$DE_2, ])) 
    names(DEG_bycat) <- c(glue('DE, unlikely {lab2} effect'), glue('DE, not an {lab2} DE'), glue('DE, {lab2} effect'))
    venn <- plot(euler(DEG_bycat, labels=names(DEG_bycat)), quantities=T)     
    p <- plot_grid(biplot, venn)
    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/axenic_heatmap/venn-biplot_{outname}.{ext}'), p, width=8, height=4, create.dir=T))
}
source('src/utilities.R')
prepEnv()
library(scales)


parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #intracellular experiment name
parser$add_argument('-a', nargs=1, type='character') #axenic experiment name
parser$add_argument('-c', nargs=1, type='character') #intracellular exp condition 
parser$add_argument('-r', nargs=1, type='character') #intracellular vehicle condition
parser$add_argument('-g', nargs=1, type='character') #intracellular group name

data_dir = './data/DE_results'
fig_dir  = './fig/axenic_heatmap'

args    <- parser$parse_args()
exp_i   <- args$i
exp_a   <- args$a
cond_i  <- args$c
ctrl_i  <- args$v
group_i <- args$g

dds_i <- readRDS(glue('./{data_dir}/{exp_i}.Rds'))

combined_ia  <- read.csv(glue('{data_dir}/combined/combined_intraaxenic_{cond_i}_{exp_a}.csv'), row.names=1)
hm_ia        <- plotAxenicHeatmap(dds_i, combined_ia, group_i, conditions=c(cond_i, ctrl_i, 'phago_4h'))
pdf(file=glue('{fig_dir}/axenic_heatmap_{cond_i}_{exp_a}.pdf'), width=8, height=10)
draw(hm_ia)
dev.off()

plotDEGs_comparison(combined_ia, 'intracellular', 'axenic', glue('{cond_i}_{exp_a}'))
