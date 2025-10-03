source('./src/utilities.R')
library(tidyverse)
library(plyr)
library(glue)
library(stats)
library(scales)
library(argparse)


## runs hypergeometric enrichment on gene list using iModulon DB data 
#  in: gene list (Rv or gene name OK) 
#  optional: specific different background size, default behavior is to consider all iModulon-mappable genes
#  returns: list containing table of iModulon counts to plot, results table for the GSEA, iModulon mappings for gene_list 
enrich_imod <- function(gene_list, n_bg=NA) {
    # gene_mod table is tidy w/ iModulon index to gene mapping
    gene_mod  <- read.csv('./data/source_data/imodulon_data/gene_presence_list.csv')
    imod_anno <- read.csv('./data/source_data/imodulon_data/iM_table.csv')
    rv2gene   <- read.csv('./data/source_data/imodulon_data/gene_info.csv')
    # rename iModulon index to iModulon name using iM table
    gene_mod$iModulon_name <- mapvalues(gene_mod$iModulon, from=imod_anno$k, to=imod_anno$name)
    # rename Rv numbers to gene names
    gene_mod$locus         <- mapvalues(gene_mod$Gene, 
                                        from=rv2gene$X, to=rv2gene$gene_name)



    # iModulon table just for DEG list, note that not all genes are represented in a modulon
    mod_shared <- gene_mod[gene_mod$locus %in% gene_list, ]
    # for statistical testing, calculate representation in each modulon, and set singly represented modulons to 'single gene' for plotting
    table <- mod_shared %>% dplyr::count(iModulon_name)
    table$condition = 'condition'
    table[table$n == 1, 'iModulon_name'] <- 'single_gene'

    
    # filter out poorly represented modulons before testing
    table_filt <- table %>% filter(n>1)

    # for each iModulon, calculate hypergeometric enrichment
    imods_test <- unique(table_filt$iModulon_name)
    n_table    <- length(rownames(mod_shared)) #number of genes in the test set (DEGs that have an iModulon mapping)
    if(is.na(n_bg)){
        n_bg       <- length(rownames(gene_mod)) # depends on hypothesis -- here all iModulon-mappable genes
    }
    imods_p    <- lapply(imods_test, function(j) { 
                         n_imodtable <- table_filt[table_filt$iModulon_name == j, 'n'] #in test set, number of iModulon members represented
                         n_imodtotal <- imod_anno[imod_anno$name == j, 'n_genes'] # total number of iModulon members
                         p_imod <- phyper(n_imodtable - 1, n_imodtotal,
                                          n_bg - n_imodtotal, n_table, 
                                          lower.tail=F)
                         return(p_imod)
                         })

    genelist_bymod        <- mod_shared %>% group_by(iModulon_name) %>% group_map(~ .x$locus) %>% as.list()
    names(genelist_bymod) <- mod_shared %>% group_by(iModulon_name) %>% group_map(~ .y) %>% unlist()
    gsea_df    <- data.frame(name=imods_test, p=unlist(imods_p), n=table_filt$n)#
    gsea_df$padj <- p.adjust(gsea_df$p, method='BH')
    gsea_df$genes <- genelist_bymod[imods_test]

    return(list(table, gsea_df, gene_mod))
}

# takes in frequency table for iModulons, enrichment test results, and name for output file
# plots a stacked bar graph, annotated with pvalues for significantly enriched iModulons
plotEnrich <- function(tables, gsea_dfs, outname) {

    tables <- lapply(seq(tables), 
                 function(i) {
                     x <- tables[[i]]
                     x['condition'] <- comp_names[i]

                     p <- gsea_dfs[[i]]$p
                     
                     length(p) <- nrow(x)
                     names         <- c(gsea_dfs[[i]]$name, rep('single gene', nrow(x) - nrow(gsea_dfs[[i]]) ))
                     p <- data.frame('iModulon_name' = names, p)
                     x <- merge(x,p, by.x='iModulon_name', all.x=T) 
                     return(x)
                 })
    if (length(tables) > 1) {
    	table  <- do.call('rbind', tables)
    } else {
    	table <- tables[[1]]
    }
    ggplot(data = table, aes(x=condition, y=n, fill=iModulon_name)) +
        geom_bar(position='fill', stat='identity') +
        scale_fill_brewer(palette='Paired') +
        geom_text(aes(label=ifelse(p < 0.05, scientific(p, digits = 3), "")),
                  position=position_fill(vjust = .5)) +
        theme_classic()
    ggsave(glue('./fig/relative_heatmap/{outname}_allcomps_iModulon.pdf'),
           width=4 * length(tables), height=3, dpi=300, units='in', create.dir=T)

}

## for a given iModulon:gene mapping and a list of genes,
## return a dataframe with gene = index, mapping each gene to iModulon category (2 columns to account for multi-map)
## handles unmapped, multi-mapping genes, flattens iModulon categories that don't map to more than one gene
makeModAnno <- function(gene_mod, gene_list) {
	all_mods     <- gene_mod[gene_mod$locus %in% gene_list, c('locus', 'iModulon_name')]

	## identify genes that don't map and genes that map to an iModulon that is represented once
	single_genes <- all_mods %>% dplyr::count(iModulon_name) %>% filter(n==1) 
	mods_keep    <- all_mods[!all_mods$iModulon_name %in% single_genes$iModulon, ]

	locus_single <- all_mods[all_mods$iModulon_name %in% single_genes$iModulon_name, 'locus']
	mods_single  <- data.frame('locus' = locus_single, 
				   'iModulon_name'=rep('single_gene', length(locus_single)))
	locus_unmapped <- gene_list[!gene_list %in% gene_mod$locus]	
	mods_unmapped  <- data.frame('locus'=locus_unmapped, 
				     'iModulon_name'=rep('unmapped', length(locus_unmapped)))
	mod_anno <- do.call('rbind', list(mods_keep, mods_single, mods_unmapped)) 
		
	# handle duplicate rows --> into two columns for plotting
	nondup     <- mod_anno[!duplicated(mod_anno$locus),]
	duplicates <- mod_anno[duplicated(mod_anno$locus), ]
	mod_anno <- join(nondup, duplicates, by='locus')
	mod_anno[is.na(mod_anno)] <- 'unmapped'
	
	# rename columns and set index for plotting
	rownames(mod_anno) <- mod_anno$locus
	mod_anno <- mod_anno[, c(2,3)]
	colnames(mod_anno) <- c('iModulon_name', 'iModulon_name_2')
	return(mod_anno[gene_list,])
}

prepEnv()

parser <- ArgumentParser()
parser$add_argument('-m', nargs=1, type='character')
parser$add_argument('-c', nargs=1, type='character') # comparison
parser$add_argument('-n', nargs=1, type='integer')
parser$add_argument('-e', nargs='+', type='character')
args <- parser$parse_args()

data_dir   = './data/DE_results'
comparison = args$c
mode       = args$m
n_bg 	   = args$n
if (is.null(n_bg)) { n_bg <- NA}
exp        = args$e

if (mode == 'intraaxenic') { ## in this mode, calculate iModulon enrichment, with and without shared 
	exp_intra    <- exp[[1]]
	exp_axenic   <- exp[[2]]
	comp <- read.csv(glue('{data_dir}/combined/combined_{mode}_{comparison}_{exp_axenic}.csv'), row.names=1)

	intra         <- rownames(comp[comp$DE_1, ])
	unique        <- readtxt(glue('{data_dir}/{exp_intra}_{comparison}_unique.txt'))
	axenic        <- rownames(comp[comp$DE_2, ])
	axenic_likely <- rownames(comp[comp$DE_2 | comp$category %in% 'DE 1, DE 2 likely', ])
	ia            <- rownames(comp[comp$DE_1 & !comp$DE_2, ])
	ia_strict     <- rownames(comp[comp$category == 'DE 1, DE 2 unlikely', ])
	ia_medium     <- rownames(comp[comp$category %in% c('DE 1, DE 2 not excluded', 'DE 1, DE 2 unlikely'), ])
	unique_exclude <- unique[!unique %in% axenic_likely]

	comp_names <- c('intra', 'intra_unique', 'intra_exclude-axenic', 'axenic', 'axenic_likely', 
			'intra_strict-exclude', 'intra_med-exclude', 'intra_unique_exclude')
	gene_lists <- list(intra, unique, ia, axenic, axenic_likely, 
			   ia_strict, ia_medium, unique_exclude)
	n_bg = c(NA, 320, NA, NA, 
		 NA, NA, NA, 320)
} else if (mode == 'drugs') {
	exp_intra <- exp
	likely_shared <- readtxt(glue('{data_dir}/{exp_intra}_likely_shared.txt'))
	pel           <- readtxt(glue('{data_dir}/{exp_intra}_pel_d1_unique.txt'))
	gef           <- readtxt(glue('{data_dir}/{exp_intra}_gef_d1_unique.txt'))
	sara          <- readtxt(glue('{data_dir}/{exp_intra}_sara_d1_unique.txt'))
	
	comp_names    <- c('likely_shared', 'pel_d1_unique', 'gef_d1_unique', 'sara_d1_unique')
	gene_lists <- list(likely_shared, pel, gef, sara)
	n_bg       <- rep(n_bg, length(gene_lists))
} else if (mode == 'single') {
	exp   <- exp
	degs <- read.csv(glue('{data_dir}/{exp}_{comparison}_vs_DMSO_d1_DE.csv'), row.names=1) %>% row.names()
	gene_lists <- list(degs)
	comp_names <- c(comparison)
	n_bg <- c(n_bg)
}
enrich_out  <- lapply(seq(gene_lists), function(i) enrich_imod(gene_lists[[i]], n_bg=n_bg[i]))

tables      <- lapply(enrich_out, function(i) i[[1]])
gsea_dfs    <- lapply(enrich_out, function(i) i[[2]])
gene_mods <- lapply(enrich_out, function(i) i[[3]])

plotEnrich(tables, gsea_dfs, glue('{comparison}_{mode}')) 


lapply(seq(length(gsea_dfs)), 
       function(i) write.csv(as.matrix(gsea_dfs[[i]]), 
			     glue('./data/enrich/{comparison}_{mode}_{comp_names[[i]]}.csv')
			     )
      )


## ALSO PLOT OVERALL EXPRESSION OF IMODS
if (mode == 'drugs') {
	gene_mod <- gene_mods[[1]]
	group = 'Drug_Day'
	conditions = c('pel_d1', 'gef_d1', 'sara_d1')
	all_degs <- lapply(conditions,
                   	   function(x) read.csv(glue('{data_dir}/{exp_intra}_{x}_vs_DMSO_d1_DE.csv'), row.names=1) %>% rownames())
	all_degs <- unique(unlist(all_degs))
	all_fc   <- lapply(conditions,
			   function(x) {df <- read.csv(glue('{data_dir}/{exp_intra}_{x}_vs_DMSO_d1_full.csv', row.names=1))
			   	       rownames(df) <- df$X
			   	       return(df)
			   	       }
			  )
	dds <- readRDS(glue('{data_dir}/{exp_intra}.Rds'))

	### plot overall expression of differentially expressed genes that map to each iModulon
	imod_plots = c('Rv1776c+WhiB4', 'MarR', 'SigC', 'HPT-2b Induced', 'Rv0681', 'Mce1R')
	for (plot in imod_plots) {
    		genes_plot <- gene_mod[gene_mod$iModulon_name %in% plot, 'locus']
    		genes_plot <- genes_plot[genes_plot %in% all_degs]
		
		fc_plot    <- lapply(all_fc, 
				     function(x) 
					     return( x[genes_plot, 'log2FoldChange'])
			            
				    )
		fc_plot    <- do.call('cbind', fc_plot)
		colnames(fc_plot)   <- lapply(conditions, function(x) glue('FC_{x}'))
		rownames(fc_plot)   <- genes_plot
		fc_plot    <- fc_plot %>% as.data.frame()
		
		anno_cols <- list(FC_sara_d1 = getFCColors(),
				  FC_gef_d1  = getFCColors(),
				  FC_pel_d1  = getFCColors())

    		hm <- plotBasicHeatmap(genes_plot, dds, group, c(conditions, 'DMSO_d1'))
		hm <- rowAnnotation(df = fc_plot, col=anno_cols) + hm

    		pdf(glue('./fig/relative_heatmap/hm_pathogen_iModulon_{plot}.pdf'))
    		draw(hm)
    		dev.off()
	}

	## 
	mod_colors <- getiModColors()

	mod_annos  <- lapply(gene_lists, function(x) makeModAnno(gene_mod, x))

	lapply(seq(gene_lists), function(i) {
	       		hm <- plotBasicHeatmap(gene_lists[[i]], dds, group, c(conditions, 'DMSO_d1'))
			hm <- rowAnnotation(df=mod_annos[[i]], 
				   	    col = list(iModulon_name = mod_colors, iModulon_name_2 = mod_colors)) + hm
			pdf(file=glue('./fig/unique-shared_heatmap/{exp_intra}_{comp_names[i]}.pdf'), 
			    height=getHeight(gene_lists[[i]]))
			draw(hm)
			dev.off()	
				})
}

