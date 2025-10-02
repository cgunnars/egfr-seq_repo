getDE <- function(dds, comparison, filestem, p_thresh=0.05, fc_thresh=1, fc_thresh_low=0.7) {
    # pull from design and comparison list to generate the contrast
    group <- unlist(as.character(design(dds))[-1])
    comp <- unlist(str_split(comparison, pattern='_vs_'))
    case <- comp[1]
    ref  <- comp[2]
    
    dds[[group]] <- relevel(dds[[group]], ref=ref)
    dds        <- DESeq(dds, quiet=T)
    ## TODO APEGLM
    res.ape            <- lfcShrink(dds=dds, coef=glue('{group}_{case}_vs_{ref}'), 
                                    type='apeglm')#, lfcThreshold=fc_thresh)
    resLA              <- results(dds, lfcThreshold=fc_thresh_low, altHypothesis='lessAbs', 
                                  name=glue('{group}_{case}_vs_{ref}'))
    res.ape['padj_LA'] <- resLA['padj'] 
    res.ape$DE        <- as.logical(abs(res.ape[, 'log2FoldChange']) > fc_thresh & res.ape[, 'padj'] < p_thresh)
    res.ape$DE[is.na(res.ape$DE)] <- F #set NA p values to neither DE nor not DE
    
    res.ape$notDE <- resLA[, 'padj'] < p_thresh
    res.ape$notDE[is.na(res.ape$notDE)] <- F

    res.ape$category <- 'DE not excluded'
    res.ape[res.ape$DE, 'category'] <- 'DE'
    res.ape[res.ape$notDE, 'category'] <- 'DE unlikely'
    res.ape[!res.ape$DE & !res.ape$notDE & abs(res.ape$log2FoldChange) > fc_thresh_low & (!is.na(res.ape$padj) & res.ape$padj < p_thresh), 'category'] <- 'DE likely'
    
    #write results tables
    results_de <- res.ape %>% subset(log2FoldChange > fc_thresh | log2FoldChange < -fc_thresh) %>%
                                  subset(padj < p_thresh)

    results_de   <- results_de[order(results_de$log2FoldChange), ] %>% data.frame
    results_up   <- row.names(results_de[results_de$log2FoldChange > 0, ])
    results_down <- row.names(results_de[results_de$log2FoldChange < 0, ])
    write.csv(res.ape,       glue('./data/DE_results/{filestem}_{comparison}_full.csv'))
    write.csv(results_de,    glue('./data/DE_results/{filestem}_{comparison}_DE.csv'))
    writeLines(results_up,   glue('./data/DE_results/{filestem}_{comparison}_up.txt'))
    writeLines(results_down, glue('./data/DE_results/{filestem}_{comparison}_down.txt'))
    return(res.ape)
}

plotVolcano <- function(results, comparison, filestem, p_thresh=0.05, fc_thresh=1){
    p <- EnhancedVolcano(results, lab=rownames(results), 
                         x='log2FoldChange', y='pvalue', 
                         FCcutoff=fc_thresh, pCutoff=p_thresh, pCutoffCol='padj',
                         labSize=5, pointSize=1, 
                         col=c('grey30', 'grey30', 'grey30', 'red2')) + theme_classic()
    ggsave(plot=p, filename=glue('./fig/DE_results/{filestem}_volcano_{comparison}.pdf'), create.dir=T)
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


vsd <- vst(dds)
dds <- DESeq(dds)
assays(dds)[['vsd']] <- vsd


filestem = basename(file_path_sans_ext(file))
saveRDS(dds, file=glue('./data/DE_results/{filestem}.Rds'))

results <- lapply(comparisons, function(c) getDE(dds, c, filestem))
lapply(seq(comparisons), function(i) plotVolcano(results[[i]], comparisons[[i]], filestem))




