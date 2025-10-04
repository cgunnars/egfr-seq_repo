calcPca <- function(vsd, ntop=500) {
    vars <- rowVars(assay(vsd))
    topgenes <- order(vars, decreasing=T)[seq(1,ntop)]
    
    pca <- prcomp(t(assay(vsd)[topgenes,]))
    return(pca)
}

calcLoadings <- function(pca) {
    loadings <- pca$rotation
    pc1 <- loadings[, 'PC1'] %>% sort(., decreasing=F)
    pc2 <- loadings[, 'PC2'] %>% sort(., decreasing=F)
}

pcaPlot <- function(vsd, color, shape, name, ntop=500, npcs=3) {
    pca <- calcPca(vsd, ntop=ntop)
    percentVar <- round(pca$sdev^2 / sum(pca$sdev^2) * 100)

    all_pcnames <- paste0('PC', c(1:3))
    d <- data.frame(V1=pca$x[,all_pcnames[1]],
                    V2=pca$x[,all_pcnames[2]],
                    V3=pca$x[,all_pcnames[3]],
                    color=vsd[[color]], shape=vsd[[shape]], 
                    name=colnames(vsd), colData(vsd))
    colnames(d)[1:3] <- all_pcnames
    attr(d, 'percentVar') <- round(percentVar[1:3])
    for (pcsToUse in list(c(1,2), c(1,3), c(2,3))) {
        pcnames <- all_pcnames[pcsToUse]
        index_1 <- pcsToUse[1]
        index_2 <- pcsToUse[2] 
        p <- ggplot(data=d, aes_string(x=pcnames[1], y=pcnames[2], color=color, shape=shape)) +
             geom_point(size=5) + 
             xlab(glue('PC{pcsToUse[1]}: {percentVar[index_1]}% variance')) +
             ylab(glue('PC{pcsToUse[2]}: {percentVar[index_2]}% variance')) +
             scale_color_manual(values=cols) + theme_classic()
        ggsave(glue('./fig/abx-dose/pca_{name}_{pcsToUse[1]}_{pcsToUse[2]}.pdf'), p, width=5, height=5, create.dir=T)
    }

}


calcDE    <- function(dds, group, cond, ref) {
    dds[[group]] <- relevel(dds[[group]], ref=ref)
    dds          <- DESeq(dds, quiet=T)
    res          <- lfcShrink(dds=dds, coef=glue('Drug_Dose_{cond}_vs_{ref}'), type='apeglm')

    return(res)
}

calcAbxDE <- function(dds, drug, dose) {
    cond = glue('{drug}_{dose}')
    group = 'Drug_Dose'
    refs <- paste0(c('pel', 'gef'), glue('_{dose}'))
    
    res  <- lapply(refs, function(ref) calcDE(dds, group, cond, ref))

    de   <- lapply(res, function(r) r %>% subset(., abs(log2FoldChange) > 1 & padj < 0.05) %>% rownames())

    
    
    names(de) <- c(glue('{drug}_{dose}_vs_pel_{dose}'), glue('{drug}_{dose}_vs_gef_{dose}'))

    return(de)
}

calcDoseDE <- function(dds, drug) {
    cond = glue('{drug}_25')
    group = 'Drug_Dose'
    ref  = glue('{drug}_5')

    res <- calcDE(dds, group, cond, ref)

    de  <- res %>% subset(., abs(log2FoldChange) > 1 & padj < 0.05) %>% rownames()

    names(de) <- c(glue('{cond}_vs_{ref}'))
    return(de)
}


library(glue)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ComplexHeatmap)
library(eulerr)
library('UpSetR')
source('./src/utilities.R')
wd = './data/DE_results'
exp = '20240808_4h-abx_pathogen'
dds <- readRDS(glue('{wd}/{exp}.Rds'))


vsd <- assays(dds)[['vsd']]
vsd$Dose <- factor(vsd$Dose, levels=c('5', '25'))
cols <- getDrugColors()

# Calculate lap_25 vs gef_25 and pel_25
de_lap_25_tables <- lapply(c('lap_25_vs_pel_25', 'lap_25_vs_gef_25'),
		           function(i) read.csv(glue('{wd}/{exp}_{i}_DE.csv', row.names=1))) 
de_lap_25        <- lapply(de_lap_25_tables, function(i) i$X)

de_var_25_tables <- lapply(c('var_25_vs_pel_25', 'var_25_vs_gef_25'),
		           function(i) read.csv(glue('{wd}/{exp}_{i}_DE.csv', row.names=1))) 
de_var_25        <- lapply(de_var_25_tables, function(i) i$X)


de_dose_tables <- lapply(c('pel','gef','var','lap'),
			function(i) read.csv(glue('{wd}/{exp}_{i}_25_vs_{i}_5_DE.csv'), row.names=1))
de_dose        <- lapply(de_dose_tables, function(i) rownames(i))
names(de_dose) <- c('pel', 'gef', 'var', 'lap')


## Lap/var upset plots
var_genes <- list(de_var_25[[1]], de_var_25[[2]], 
                  unname(unlist((de_dose$var)))
                )

names(var_genes) <- c('var_25_vs_pel_25', 'var_25_vs_gef_25', 
                      'var_25_vs_var_5')

pl <- upset(fromList(var_genes), order.by='freq')
pdf('~/egfr-seq_repo/fig/abx-dose/var_upset.pdf', height=5, width=5)
pl
dev.off()

lap_genes <- list(de_lap_25[[1]], de_lap_25[[2]], 
                  unname(unlist((de_dose$lap)))
                  )
names(lap_genes) <- c('lap_25_vs_pel_25', 'lap_25_vs_gef_25', 
                      'lap_25_vs_lap_5')
pl <- upset(fromList(lap_genes), order.by='freq')
pdf('~/egfr-seq_repo/fig/abx-dose/lap_upset.pdf', height=5, width=5)
pl
dev.off()

pl <- upset(fromList(c(var_genes[1:2], lap_genes[1:2])), order.by='freq')
pdf('~/egfr-seq_repo/fig/abx-dose/lap-var_upset.pdf', height=7, width=5)
pl
dev.off()

pl <- upset(fromList(de_dose), order.by='freq')
pdf('~/egfr-seq_repo/fig/abx-dose/dose_upset.pdf', height=7, width=5)
pl 
dev.off()

## HEATMAPS LAP + VAR
lysine_genes <- c('lat', 'pcd', 'Rv3292', 'usfY')
napm_genes   <- c('Rv0045c', 'ino1', 'Rv0047c')
mmps4_genes  <- c('mmpS4', 'mmpL4', 'Rv0452')
pyr_genes    <- c('pyrR', 'pyrB', 'pyrC', 'Rv1382', 'carA', 'carB')

var_genes_plot <- c(lysine_genes, napm_genes, mmps4_genes, pyr_genes)

row_anno <- rowAnnotation(de_pel = var_genes_plot %in% var_genes$var_25_vs_pel_25, 
                          de_gef = var_genes_plot %in% var_genes$var_25_vs_gef_25,
                          de_dose = var_genes_plot %in% var_genes$var_25_vs_var_5,
			  col = list(de_pel  = c('FALSE'='white', 'TRUE'='black'),
				     de_gef  = c('FALSE'='white', 'TRUE'='black'),
				     de_dose = c('FALSE'='white', 'TRUE'='black')))
    
set.seed(3)
mat_subset <- vsd[var_genes_plot, ]
scaled_mat <- t(scale(t(assay(mat_subset))))

column_split = mat_subset$Drug_Dose
row_split    = c(rep('Lysine\ndegradation', length(lysine_genes)), 
                 rep('NapM', length(napm_genes)),
                 rep('mmpS4', length(mmps4_genes)),
                 rep('Pyrimidine\nsalvage', length(pyr_genes)))
hm <- Heatmap(scaled_mat, name='row z(vsd)', column_split=column_split, row_split=row_split, left_annotation=row_anno)
pdf('~/egfr-seq_repo/fig/abx-dose/var-heatmap.pdf', height=5)
draw(hm)
dev.off()


rpl_genes <- rownames(dds)[grep('^rpl', rownames(dds))]
rps_genes <- rownames(dds)[grep('^rps', rownames(dds))]
ribogenes = c(rpl_genes, rps_genes)

atpgenes <- rownames(dds)[grep('^(ndh|atp|nuo|sdh|qcr|cta|cyd)', rownames(dds))]
lexa <- c('lexA', 'recA', 'ruvA', 'ruvB', 'pafA', 'pafB', 'pafC', 'dnaE2', 'ung', 'lhr', 'xthA')
oxidative <- c('katG', 'furA', 'cysD', 'cysN', 'cysO', 'mec', 'mshA', 'mshD', 'mshC')
inib <- c('iniB', 'iniA', 'iniC','mutA', 'mutB')

set.seed(3)
hm<-Heatmap(t(scale(t(assay(vsd[c(atpgenes, ribogenes, lexa, oxidative, inib),])))),
        row_split=c(rep('ETC genes',length(atpgenes)), 
                    rep('Ribosomal genes' ,length(ribogenes)),
                    rep('DNA damage',length(lexa)),
                    rep('Oxidative stress',length(oxidative)),
                    rep('iniBAC',length(inib))),
        column_split=vsd$Drug_Dose, name= 'row z(vsd)')
pdf('~/egfr-seq_repo/fig/abx-dose/lap-heatmap.pdf', height=18)
draw(hm)
dev.off()

hm_peldose <- plotBasicHeatmap(de_dose$pel, dds, 
			       c('pel_25', 'pel_5', 'gef_25', 'gef_5', 'lap_25', 'lap_5', 'var_25', 'var_5'))
pdf(file='./fig/abx-dose/pel-dose-heatmap.pdf', width=8, height=25)#, units='in', res=480)
draw(hm_peldose)
dev.off()

hm_gefdose <- plotBasicHeatmap(de_dose$gef, dds, 
			       c('pel_25', 'pel_5', 'gef_25', 'gef_5', 'lap_25', 'lap_5', 'var_25', 'var_5'))
pdf(file='./fig/abx-dose/gef-dose-heatmap.pdf', width=8, height=25)#, units='in', res=480)
draw(hm_gefdose)
dev.off()


## PCA plots with and without lapatinib 
p_all <- pcaPlot(vsd, color='Drug', shape='Dose', name='all')
p_nolap <- pcaPlot(vsd[, !vsd$Drug == 'lap'], color='Drug', shape='Dose', name='nolap')

pca_nolap <- calcPca(vsd[, !vsd$Drug == 'lap'])
loadings <- pca_nolap$rotation
pc1      <- loadings[,'PC1'] %>% sort(., decreasing=F)
pc2      <- loadings[,'PC2'] %>% sort(., decreasing=F)
## Plot top 20 PC loadings (in each direction)
top_loadings   <- c(pc1[1:10], pc1[491:500])
pc_loadings_df <- top_loadings %>% as.data.frame()
colnames(pc_loadings_df) <- c('PC1')
pc_loadings_df$gene <- factor(names(top_loadings), levels=names(top_loadings))
pc_loadings_df$gef_DEG <- rownames(pc_loadings_df) %in% unname(unlist((de_dose$gef)))
pc_loadings_df$var_DEG <- rownames(pc_loadings_df) %in% unname(unlist((de_dose$var)))


p_loadings <- ggplot(pc_loadings_df, aes(y=PC1,x=gene,fill=gef_DEG)) + geom_bar(stat='identity') + theme_classic() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave('./fig/abx-dose/non-lap_loadings.pdf', p_loadings, width=4.5, height=2, create.dir=T)





