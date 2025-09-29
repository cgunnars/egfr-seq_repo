getDE <- function(dds, comparison, filestem, p_thresh=0.05, fc_thresh=1) {
    # pull from design and comparison list to generate the contrast
    group      <- unlist(as.character(design(dds))[-1])
    comp <- unlist(str_split(comparison, pattern='_vs_'))
    case <- comp[1]
    ref  <- comp[2]
    
    dds[[group]] <- relevel(dds[[group]], ref=ref)
    dds        <- DESeq(dds, quiet=T)
    ## TODO APEGLM
    res.ape            <- lfcShrink(dds=dds, coef=glue('{group}_{case}_vs_{ref}'), 
                                    type='apeglm')#, lfcThreshold=fc_thresh)
    res.apethresh      <- lfcShrink(dds=dds, coef=glue('{group}_{case}_vs_{ref}'),
                                    type='apeglm', lfcThreshold=fc_thresh)
    resLA              <- results(dds, lfcThreshold=fc_thresh, altHypothesis='lessAbs', 
                                  name=glue('{group}_{case}_vs_{ref}'))
    res.ape['padj_LA'] <- resLA['padj'] 
    res.ape['svalue'] <- res.apethresh['svalue']
    #write results tables
    results    <- results(dds, alpha=0.05, contrast=c(group,case,ref))
    results_de <- res.ape %>% subset(log2FoldChange > fc_thresh | log2FoldChange < -fc_thresh) %>%
                                  subset(padj < p_thresh)

    results_de   <- results_de[order(results_de$log2FoldChange), ] %>% data.frame
    results_up   <- row.names(results_de[results_de$log2FoldChange > 0, ])
    results_down <- row.names(results_de[results_de$log2FoldChange < 0, ])
    write.csv(res.ape,       glue('./data/DE_results/{filestem}_{comparison}_full.csv'))
    write.csv(results_de,    glue('./data/DE_results/{filestem}_{comparison}_DE.csv'))
    writeLines(results_up,   glue('./data/DE_results/{filestem}_{comparison}_up.txt'))
    writeLines(results_down, glue('./data/DE_results/{filestem}_{comparison}_down.txt'))

    p <- EnhancedVolcano(results, lab=rownames(results), 
                         x='log2FoldChange', y='pvalue', 
                         FCcutoff=1, pCutoff=0.05, pCutoffCol='padj',
                         labSize=5, pointSize=1, 
                         col=c('grey30', 'grey30', 'grey30', 'red2')) + theme_classic()
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




