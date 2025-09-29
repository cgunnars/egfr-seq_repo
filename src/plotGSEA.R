source('src/utilities.R')
prepEnv()
library(scales)


plotGSEA <- function(fc, gene_lists, outname) {
    gsea_results <- GSEA(geneList = fc, TERM2GENE = gene_lists, verbose = F, 
                         eps=1e-10, minGSSize = 1, pvalueCutoff = 1)
                          
    write.csv(gsea_results@result, glue('./data/enrich/gsea_{outname}.csv'))
    genesets <- unique(gene_lists$term) 
    figs     <- lapply(genesets, function(geneset_plot)
                       gseaplot(gsea_results, geneSetID = geneset_plot, 
                            	by = "runningScore", title = geneset_plot) +
                       geom_text(label=ifelse(gsea_results@result[gsea_results@result$ID == geneset_plot, 'p.adjust'] < 0.01,
                                           paste0(sprintf('NES %.2f \n p=', 
                                                          gsea_results@result[gsea_results@result$ID == geneset_plot, 'NES']), 
                                              scientific(gsea_results@result[gsea_results@result$ID == geneset_plot, 'p.adjust'], 3)),
                                           sprintf('NES %.2f \n p=%.2f',
                                                    gsea_results@result[gsea_results@result$ID == geneset_plot, 'NES'],
                                                    gsea_results@result[gsea_results@result$ID == geneset_plot, 'p.adjust'])
                                           ), x=500, y=0.9, check_overlap=T)+ 
                    scale_y_continuous(limits = c(-1, 1))
                )
    figs      <- plot_grid(plotlist=figs, ncol = length(figs))
    lapply(c('pdf', 'svg'), function(ext) ggsave(glue('./fig/gsea/gsea_{outname}.{ext}'), figs, width=4*length(genesets), height=4))
    
}

parser <- ArgumentParser() #
parser$add_argument('-i', nargs=1, type='character') # combined FC csv
parser$add_argument('-g', nargs='+', type='character') # names of gene sets to plot enrichment
parser$add_argument('-o', nargs=1, type='character')

args <- parser$parse_args()
combined <- read.csv(args$i, row.names=1)
gene_lists <- lapply(args$g, function(i) data.frame(term=i, gene_list=readtxt(glue('./data/gene_lists/{i}.txt'))))
gene_lists <- do.call('rbind', gene_lists)


list_1 <- combined$log2FoldChange
names(list_1) <- rownames(combined)
list_1 <- sort(list_1, decreasing=T)
plotGSEA(list_1, gene_lists, glue('{args$o}'))


