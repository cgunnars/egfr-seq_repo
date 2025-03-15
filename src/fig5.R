getDEGs_time <- function(cond, exp, refs) {
  up_genes <- lapply(seq(3), 
                    function(i) {readLines(
                    file(glue('{wd}/{exp}_{cond}_d{i}_vs_{refs[i]}_up.txt'))
                   )
                   }
                   )
  names(up_genes) <- lapply(seq(3), function(i) glue('{cond}_d{i}_vs_{refs[i]}'))
  down_genes <- lapply(seq(3), 
                        function(i) {readLines(
                        file(glue('{wd}/{exp}_{cond}_d{i}_vs_{refs[i]}_down.txt'))
                    )
                    }
                    )
  names(down_genes) = names(up_genes)

  num_up     <- unlist(lapply(up_genes, function(i) length(i)))
  num_down   <- unlist(lapply(down_genes, function(i) length(i)))

  df <- data.frame(day = seq(3), n_de = num_up, condition = cond, direction = 'up', ref=refs)
  df <- rbind(df, 
              data.frame(day = seq(3), n_de = num_down, condition = cond, direction = 'down', ref=refs))
  return(list(up_genes=up_genes, down_genes=down_genes, df=df))
}

plotDEGs_time <- function(exp) {
    data_pel <- getDEGs_time(refs = c('phago_4h', 'pel_d1', 'pel_d2'), exp = exp, cond = 'pel')
    data_DMSO <- getDEGs_time(refs = c('phago_4h', 'DMSO_d1', 'DMSO_d2'), exp=exp, cond = 'DMSO')
    data_DMSOvspel <- getDEGs_time(refs = c('DMSO_d1', 'DMSO_d2', 'DMSO_d3'), exp=exp, cond='pel')
    all_data <- rbind(data_pel$df, data_DMSO$df)#, data_DMSO)

    fig_5a <- ggplot(data = all_data, aes(x=day, y=n_de, col=condition)) + 
              geom_line(aes(linetype=direction)) + 
              labs(x='Day', y='Number of DEGs') +
              scale_x_continuous(labels = function(x) round(as.numeric(x))) +
              scale_color_hue(labels = c('vehicle (DMSO)', 'pelitinib')) +
              theme_classic()
    fig_5b <- ggplot(data = data_DMSOvspel$df, aes(x=day, y=n_de)) + 
              geom_line(aes(linetype=direction)) + 
              labs(x='Day', y='Number of DEGs') +
              scale_x_continuous(labels = function(x) round(as.numeric(x))) +
              theme_classic()
    
    fig_5ab <- plot_grid(fig_5a, fig_5b, labels=c('A', 'B'))
    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/fig5/fig5ab.{ext}'), plot=fig_5ab, width=7, height=2.5))
}

getDEGs_method <- function(dds, coef, fc_thresh=c(0.3, 0.5, 0.7, 0.8, 0.9, 1)) {
    # test methods for calling fold changes, filter only significant results
    resape <- lapply(fc_thresh, 
                     function (i) 
                        lfcShrink(dds = dds, coef = coef, 
                                  type='apeglm', lfcThreshold = i) %>% 
                                  subset(svalue < 0.005))
    resLA  <- lapply(fc_thresh, 
                     function(i)
                        results(dds, lfcThreshold = i, 
                                altHypothesis = 'lessAbs', 
                                name = coef) %>% 
                        subset(padj < 0.05))
    resGA  <- lapply(fc_thresh, function(i)
                        results(dds, lfcThreshold = i, 
                                altHypothesis = 'greaterAbs', 
                                name = coef) %>% 
                        subset(padj < 0.05))
    res    <- lapply(fc_thresh, function(i)
                        results(dds, name = coef) %>% 
                        subset(padj < 0.05) %>% 
                        subset(log2FoldChange > i | log2FoldChange < -i))
    
    resallshr <- lfcShrink(dds = dds, coef = coef, type='apeglm')
    resshr <- lapply(fc_thresh, function (i) 
                        resallshr %>% 
                        subset(padj < 0.05) %>% 
                        subset(log2FoldChange > i | log2FoldChange < -i))

    # number of DEGS per method
    nape   <- unlist(lapply(resape, function(i) dim(i)[[1]]))
    nLA    <- unlist(lapply(resLA,  function(i) dim(i)[[1]]))
    nGA    <- unlist(lapply(resGA,  function(i) dim(i)[[1]]))
    nreg   <- unlist(lapply(res,    function(i) dim(i)[[1]]))
    nshr   <- unlist(lapply(resshr, function(i) dim(i)[[1]]))
    nde_df <- do.call('rbind', 
                      list(
                        data.frame(fc_thresh=fc_thresh, n=nape, 
                                   method='apeGLM shrunk, post-hoc above threshold'),
                        data.frame(fc_thresh=fc_thresh, n=nreg, 
                                   method='unshrunk, post-hoc above threshold'),
                        data.frame(fc_thresh=fc_thresh, n=nshr, 
                                   method='apeGLM shrunk, test above threshold'),
                        data.frame(fc_thresh=fc_thresh, n=nGA,
                                   method='Wald test above threshold')))
    
    # number of unspecified for best method
    category_levels <- c('unspecified', 
                         'apeGLM shrunk, post-hoc above threshold', 
                         'Wald test below threshold')
    ngenes <- dim(dds)[1]
    unspec <- ngenes - (nshr + nLA)
    nunspec_df <- do.call('rbind', list(
                  data.frame(fc_thresh=fc_thresh, n=nshr, category='apeGLM shrunk, post-hoc above threshold'),
                  data.frame(fc_thresh=fc_thresh, n=nLA, category='Wald test below threshold'),
                  data.frame(fc_thresh=fc_thresh, n=unspec, category='unspecified')))
    nunspec_df$category <- factor(nunspec_df$category, levels=category_levels) 
                               

    # categorize genes by unspecified, above, or below
    allgenes <- rownames(dds)
    unspec_genes    <- lapply(seq(length(resshr)), 
                           function(i) {
                               allgenes[!allgenes %in% union(rownames(resshr[[i]]), rownames(resLA[[i]]))]
                           })

    unspecgenes_df  <- lapply(seq(length(unspec_genes)), function(i) { 
                                do.call('rbind', list(
                                data.frame(fc_thresh=fc_thresh[i],
                                           category='unspecified',
                                           baseMean=resallshr[unspec_genes[[i]], 'baseMean'],
                                           log2FoldChange=abs(resallshr[unspec_genes[[i]], 'log2FoldChange']),
                                           padj=resallshr[unspec_genes[[i]], 'padj']),
                                data.frame(fc_thresh=fc_thresh[i],
                                           category='apeGLM shrunk, post-hoc above threshold',
                                           baseMean=resshr[[i]][, 'baseMean'],
                                           log2FoldChange=abs(resshr[[i]][, 'log2FoldChange']),
                                           padj=resshr[[i]][, 'padj']),
                                data.frame(fc_thresh=fc_thresh[i],
                                           category='Wald test below threshold',
                                           baseMean=resLA[[i]][, 'baseMean'],
                                           log2FoldChange=abs(resLA[[i]][, 'log2FoldChange']),
                                           padj=resLA[[i]][, 'padj'])
                                 ))
                           })
     unspecgenes_df <- do.call('rbind', unspecgenes_df)
     unspecgenes_df$fc_thresh <- factor(unspecgenes_df$fc_thresh)
     unspecgenes_df$category <- factor(unspecgenes_df$category, 
                                       levels=category_levels)

     results <- list(nde_df, nunspec_df, unspecgenes_df)
     names(results) <- c('nde_df', 'nunspec_df', 'unspecgenes_df')
     
     return(results)
}

plotDEGs_method <- function(dds, coef, fc_thresh) {
    results <- getDEGs_method(dds, coef, fc_thresh)
    nde_df  <- results$nde_df
    sfig_6a <- ggplot(data=nde_df, aes(x=fc_thresh, y=n, col=method)) + geom_line() + 
               labs(x='Log2 fold change threshold', y='Number of DEGs') + theme_classic()
    unspecgenes_df <- results$unspecgenes_df
    lines = data.frame(n=seq(length(fc_thresh)), fc_thresh=fc_thresh, fc=fc_thresh)

    sfig_6b <- ggplot(unspecgenes_df, aes(x=baseMean, y=log2FoldChange, col=category)) + 
               geom_point(alpha=0.3, size=0.5) + scale_x_log10() + 
               geom_abline(data=lines, aes(intercept=fc, slope=0))

    # Use vars() to supply variables from the dataset:
    sfig_6b <- sfig_6b + facet_grid(cols = vars(fc_thresh)) + theme_classic()
     
    nunspec_df <- results$nunspec_df
    sfig_6c <- ggplot(data = nunspec_df, aes(x=fc_thresh, y=n, fill=category)) + geom_area(alpha=0.6) + theme_classic()

    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/sfig6/sfig6a.{ext}'), sfig_6a))
    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/sfig6/sfig6b.{ext}'), sfig_6b))
    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/sfig6/sfig6c.{ext}'), sfig_6c))
}

# for a comparison (coef), test whether is DE (> abs(log2fcthresh_high & < pthresh) or 
#                          test whether is NOT DE (< fcthresh_low by Wald's test)
getDEGs_category <- function(dds, coef, pthresh = 0.05, fcthresh_high=1, fcthresh_low=0.7) {
    res_shr <- lfcShrink(dds, coef, type='apeglm')
    res_LA  <- results(dds, lfcThreshold = fcthresh_low, altHypothesis = 'lessAbs', 
                       name = coef)
    res_shr$DE <- as.logical(abs(res_shr[, 'log2FoldChange']) > fcthresh_high & res_shr[, 'padj'] < pthresh)
    res_shr$DE[is.na(res_shr$DE)] <- F
    res_shr$notDE <- res_LA[, 'padj'] < 0.05
    res_shr$notDE[is.na(res_shr$notDE)] <- F

    return(res_shr)
}
getDEGs_comparison <- function(dds1, dds2, coef1, coef2, fcthresh_low=0.7) {
    res_1 <- getDEGs_category(dds1, coef1, fcthresh_low=fcthresh_low)
    res_2 <- getDEGs_category(dds2, coef2, fcthresh_low=fcthresh_low)

    combined_rownames <- intersect(rownames(res_1), rownames(res_2))
    combined_colnames <- c('log2FoldChange', 'DE', 'notDE', 'baseMean')
    combined <- cbind(res_1[combined_rownames,
                            combined_colnames], 
                      res_2[combined_rownames,
                            combined_colnames])
    colnames(combined) <- c('FC_1', 'DE_1', 'notDE_1', 'baseMean_1', 
                            'FC_2', 'DE_2', 'notDE_2', 'baseMean_2')

    combined$category <- 'mixed evidence'
    combined[combined$DE_1 & !combined$notDE_2, 'category'] <- 'DE 1, DE 2 not excluded'
    combined[combined$DE_1 & combined$FC_2 > fcthresh_low & !combined$DE_2, 'category'] <- 'DE 1, DE 2 likely'
    combined[combined$DE_1 & combined$notDE_2, 'category']  <- 'DE 1, DE 2 unlikely'                          
    combined[combined$DE_2 & combined$notDE_1, 'category']  <- 'DE 2, DE 1 unlikely'
    combined[combined$DE_2 & !combined$notDE_1, 'category'] <- 'DE 2, DE 1 not excluded'
    combined[combined$DE_2 & combined$FC_1 > fcthresh_low & !combined$DE_1, 'category'] <- 'DE 2, DE 1 likely'
    combined[combined$DE_1 & combined$DE_2, 'category']     <- 'DE 1, DE 2'
    combined[combined$notDE_1 & combined$notDE_2, 'category'] <- 'low change'

    return(combined)
}

plotDEGs_comparison <- function(combined, lab1, lab2, fcthresh_low=0.7) {
    sfig_6e <- ggplot(data=combined, aes(x=FC_2, y=FC_1, col=DE_1)) + 
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
    sfig_6f <- plot(euler(DEG_bycat, labels=names(DEG_bycat)), quantities=T)     
    p <- plot_grid(sfig_6e, sfig_6f, labels=c('E', 'F'))
    lapply(c('pdf', 'png'), function(ext) ggsave(glue('./fig/sfig6/sfig6ef.{ext}'), p, width=8, height=4))
}

plotAxenicHeatmap <- function(dds, combined, group, conditions, mode='axenic') {
    set.seed(3)
    vsd  <- assays(dds)[['vsd']]

    col_fun = colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "magenta"))
    if (mode == 'axenic') {
        degs <- rownames(combined[combined$DE_1, ])
    
        annotation_row <- data.frame(category=combined[degs, 'category'], 
                                     fc_axenic=(combined[degs, 'FC_2']), row.names=degs)

        col = list(fc_axenic = col_fun,
               category = hue_pal()(4))
        names(col$category) <- c('DE 1, DE 2', 'DE 1, DE 2 likely', 'DE 1, DE 2 unlikely', 'DE 1, DE 2 not excluded')
    } else {
        degs <- rownames(combined[combined$DE_1 | combined$DE_2, ])
        annotation_row <- data.frame(category=combined[degs,'category'],
                                     fc_1=combined[degs, 'FC_1'], 
                                     fc_2=combined[degs, 'FC_2'],
                                     row.names=degs)
        col = list(fc_1 = col_fun, fc_2 = col_fun, 
                   category = c('#F8766D', '#FF61C3', '#DB72FB', '#00B9E3', '#619CFF', '#00C19F', '#00BA38'))

        names(col$category) <- c('DE 1, DE 2', 'DE 1, DE 2 likely', 'DE 2, DE 1 likely', 'DE 1, DE 2 not excluded', 
                                 'DE 1, DE 2 unlikely', 'DE 2, DE 1 not excluded', 'DE 2, DE 1 unlikely') 
    }

    mat_subset <- vsd[degs, vsd[[group]] %in% conditions]
    scaled_mat <- t(scale(t(assay(mat_subset))))
    column_split = mat_subset[[group]]
    hm <- Heatmap(scaled_mat, name='row z(vsd)', column_split=column_split, row_km = 3, row_km_repeats=50,
                  left_annotation = rowAnnotation(df = annotation_row, col=col)
                  )
    
    
    return(hm)
}


library(glue)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(apeglm)
library(dplyr)
library(eulerr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(scales)

wd = './data/DE_results'
exp_intra = '20240502_pel-timecourse-6donor_pathogen'
dds_intra   <- readRDS(glue('{wd}/{exp_intra}.Rds'))
coefs_intra <- lapply(c('pel', 'gef', 'sara'), function(i) glue('Drug_Day_{i}_d1_vs_DMSO_d1'))

exp_axenic   = '20240816_24h-axenic_pathogen'
dds_axenic   <- readRDS(glue('{wd}/{exp_axenic}.Rds'))
coefs_axenic <- lapply(c('pel', 'gef'), function(i) glue('Drug_{i}_vs_DMSO'))


#plotDEGs_time(exp_intra)
#plotDEGs_method(dds_intra, 'Drug_Day_pel_d1_vs_DMSO_d1', c(0.3, 0.5, 0.7, 0.8, 0.9, 1))


combined_pel_intra_axenic <- getDEGs_comparison(dds_intra, dds_axenic, coefs_intra[[1]], coefs_axenic[[1]])
combined_gef_intra_axenic <- getDEGs_comparison(dds_intra, dds_axenic, coefs_intra[[2]], coefs_axenic[[2]])

plotDEGs_comparison(combined_pel_intra_axenic, 'intracellular', 'axenic')

hm_pel <- plotAxenicHeatmap(dds_intra, combined_pel_intra_axenic, 'Drug_Day', conditions=c('pel_d1', 'DMSO_d1', 'phago_4h'))
png(file='./fig/fig5/fig5d.png', width=8, height=10, units='in', res=480)#, width=8*units(1, 'in'), height=10*units(1,'in'))
draw(hm_pel)
dev.off()

hm_gef <- plotAxenicHeatmap(dds_intra, combined_gef_intra_axenic, 'Drug_Day', conditions=c('gef_d1', 'DMSO_d1', 'phago_4h'))
png(file='./fig/fig5/fig5e.png', width=8, height=10, units='in', res=480)
draw(hm_gef)
dev.off()

combined_pelgef <- getDEGs_comparison(dds_intra, dds_intra, coefs_intra[[1]], coefs_intra[[2]])
hm_pelgef <- plotAxenicHeatmap(dds_intra, combined_pelgef, 'Drug_Day', conditions=c('pel_d1', 'DMSO_d1', 'gef_d1'), mode='joint')
png(file='./fig/fig5/fig5f.png', width=8, height=15, units='in', res=480)
draw(hm_pelgef)
dev.off()
