## for a given dds, 
## provide a list of metadata categories to split sample names into, e.g. sara_d1_A => Drug Day Donor  
makeGeneDf <- function(dds, sample_meta, 
		       drug_levels=c('pel','gef','sara','lap','var','DMSO','phago','tb'), 
		       day_levels=c('d4', 'd3', 'd2','d1','4h','0h'),
		       dose_levels=c('25', '5')) {
    vsd <- assays(dds)[['vsd']]
    genedf <- assay(vsd) %>% data.frame() %>% rownames_to_column(var='row') %>% 
              pivot_longer(., cols=-row) %>% data.frame()
    meta <-   str_split_fixed(genedf$name, pattern='_', length(sample_meta)) %>% data.frame()
    colnames(meta) <- sample_meta

    for (n in sample_meta) {
    	genedf[[n]] <- meta[[n]]
    	if (n == 'Drug') {
	    drug_levels <- drug_levels[drug_levels %in% genedf[[n]]]
	    genedf[[n]] <- factor(genedf[[n]], levels=drug_levels)
	} else if (n == 'Day') {
	    day_levels <- day_levels[day_levels %in% genedf[[n]]]
	    genedf[[n]] <- factor(genedf[[n]], levels=day_levels)
	} else if (n == 'Dose') {
	    drug_levels <- drug_levels[drug_levels %in% genedf[[n]]]
	    genedf[[n]] <- factor(genedf[[n]], levels=dose_levels)
	}
    } 
    return(genedf)
}

library(argparse)
library(tidyverse)
library(glue)
library(DESeq2)


parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character')
parser$add_argument('-s', nargs='+', type='character')

wd   <- './data/DE_results/'
args <- parser$parse_args()
exp <- args$i
sample_meta <- args$s

dds <- readRDS(glue('{wd}/{exp}.Rds'))

df  <- makeGeneDf(dds, sample_meta)

print(head(df))

write.csv(df, glue('./data/clean_dds/{exp}_df.csv'))
