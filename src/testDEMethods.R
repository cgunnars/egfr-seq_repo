## helper function for plotDEGs_method
#  for a given DESeq2 object dds and comparison coef (in format GROUP_COMPARISON_vs_CONTROL),
#  calculate the number of DEGs against a range of log2fc thresholds (fc_thresh)
#  methods -- resape apeglm,   test against the threshold
#             resLA  Wald,     test  below the threshold
#             resGA  Wald,     test  above the threshold
#             res    unshrunk, post-hoc threshold    
#             resshr shrunk,   post-hoc threshold
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
                                   method='apeGLM shrunk, test above threshold'),
                        data.frame(fc_thresh=fc_thresh, n=nreg, 
                                   method='unshrunk, post-hoc above threshold'),
                        data.frame(fc_thresh=fc_thresh, n=nshr, 
                                   method='apeGLM shrunk, post-hoc above threshold'),
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

# for a given DESeq2 object dds, comparison coef (GROUP_COND_vs_CONTROL), and character vector of absfc_thresholds
# generate plots of the number of positively identifed DEGs per method (sfig6a)
#                                 statistical categorization vs expression level (sfig6b)
#                                 total number of unspecified genes vs fc threshold (sfig6c)
plotDEGs_method <- function(dds, coef, fc_thresh, outname) {
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

    ggsave(glue('./fig/deg-method/sfig6a_{outname}.pdf'), sfig_6a, width=6, height=3)
    ggsave(glue('./fig/deg-method/sfig6b_{outname}.pdf'), sfig_6b, width=15, height=3)
    ggsave(glue('./fig/deg-method/sfig6c_{outname}.pdf'), sfig_6c, width=6, height=3)
}

source('./src/utilities.R')
prepEnv()


parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character')
parser$add_argument('-c', nargs=1, type='character')
parser$add_argument('-o', nargs=1, type='character')

args    <- parser$parse_args()
exp     <- args$i
coef    <- args$c
outname <- args$o

dds <- readRDS(exp)
plotDEGs_method(dds, coef, c(0.3, 0.5, 0.7, 0.8, 0.9, 1), outname)
