getDE <- function(dds, comparison, filestem, p_thresh=0.05, fc_thresh=1) {
    # pull from design and comparison list to generate the contrast
    print(comparison)
    group      <- unlist(as.character(design(dds))[-1])
    comp <- unlist(str_split(comparison, pattern='_vs_'))
    case <- comp[1]
    ref  <- comp[2]

    #write results tables
    results    <- results(dds, alpha=0.05, contrast=c(group,case,ref))
    results_de <- results %>% subset(log2FoldChange > fc_thresh | log2FoldChange < -fc_thresh) %>%
                              subset(padj < p_thresh)
    results_de   <- results_de[order(results_de$log2FoldChange), ] %>% data.frame
    results_up   <- row.names(results_de[results_de$log2FoldChange > 0, ])
    results_down <- row.names(results_de[results_de$log2FoldChange < 0, ])
    write.csv(results,       glue('./data/DE_results/{filestem}_{comparison}_full.csv'))
    write.csv(results_de,    glue('./data/DE_results/{filestem}_{comparison}_DE.csv'))
    writeLines(results_up,   glue('./data/DE_results/{filestem}_{comparison}_up.txt'))
    writeLines(results_down, glue('./data/DE_results/{filestem}_{comparison}_down.txt'))

    p <- EnhancedVolcano(results, lab=rownames(results), 
                         x='log2FoldChange', y='pvalue', 
                         FCcutoff=1, pCutoff=0.05, pCutoffCol='padj',
                         labSize=5, pointSize=1, 
                         col=c('grey30', 'grey30', 'grey30', 'red2'))
    ggsave(plot=p, filename=glue('./fig/DE_results/{filestem}_volcano_{comparison}.pdf'))
}

source('src/utilities.R')
prepEnv()
library(argparse)

parser <- ArgumentParser()
parser$add_argument('-i', '--input_file', type='character', nargs=1, help='DESeq2 dataset object')
parser$add_argument('-c', '--comparisons_file', type='character', nargs=1, help='txt file of comparisons')

args        <- parser$parse_args()
file        <- args$i
comparisons <- readtxt(args$c)

dds <- readRDS(file)

group_comparison <- as.character(design(dds))[-1]
if ('Donor' %in% colnames(colData(dds))){
    collapse_by        <- paste0(group_comparison, '_Donor')
    dds[[collapse_by]] <- paste0(dds[[group_comparison]], '_', dds[['Donor']]) 
    dds_coll           <- collapseReplicates(dds, dds[[collapse_by]], dds[['Replicate']])
    dds                <- dds_coll
}

vsd <- vst(dds)
dds <- DESeq(dds)
assays(dds)[['vsd']] <- vsd


filestem = basename(file_path_sans_ext(file))
saveRDS(dds, file=glue('./data/DE_results/{filestem}.Rds'))

lapply(comparisons, function(c) getDE(dds, c, filestem))




