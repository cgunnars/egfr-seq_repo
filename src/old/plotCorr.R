# if thresholding a self-correlation, the diagonal will be == 1 and the rho cutoff won't work well
threshold_cor <- function(rho_cor, p_cor, rho_cutoff = 0.5, p_cutoff = 10e-7) {
    return(rho_cor[(rowMin(p_cor)  < p_cutoff & rowMax(abs(rho_cor))  > rho_cutoff),
                   (colMins(p_cor) < p_cutoff & colMaxs(abs(rho_cor)) > rho_cutoff)]) 
}

source('src/utilities.R')
prepEnv()
library(argparse)

parser <- ArgumentParser()
parser$add_argument('-i', '--input_files', type='character', nargs=2, help='DESeq2 dataset objects')
parser$add_argument('-c', type='character', nargs='+', help='txt file for deg comparisons')
parser$add_argument('-m', '--match_on', type='character', nargs='+', help='metadata to match on')
parser$add_argument('-d', type='character', nargs='+', help='txt file for deg comparisons')

args        <- parser$parse_args()
files       <- args$i
match_group <- args$m
c1          <- args$c
c2          <- args$d

#print(args$m)

dds1 <- readRDS(files[[1]])
dds2 <- readRDS(files[[2]])

# e.g. match on Drug Day Donor Replicate
#      design is Drug_Day or Drug
#      make sure Drug_Day_Donor is in both metadatas
if (length(match_group) > 0) {
    match_values1 <- lapply(match_group, function(x) colData(dds1)[[x]])
    match_values2 <- lapply(match_group, function(x) colData(dds2)[[x]]) 
    match_group   <- paste(match_group, collapse='_')
    if (!match_group %in% colnames(colData(dds1))) {
        colData(dds1)[match_group] <- Reduce(function(x, y) paste0(x, '_', y),
                                             match_values1) 
    }
    if (!match_group %in% colnames(colData(dds2))) {
        colData(dds2)[match_group] <- Reduce(function(x, y) paste0(x, '_', y),
                                             match_values2)
    }
    common <- intersect(dds1[[match_group]], dds2[[match_group]])
    dds1 <- dds1[, dds1[[match_group]] %in% common]
    dds1 <- dds1[, order(dds1[[match_group]])]
    dds2 <- dds2[, dds2[[match_group]] %in% common]
    dds2 <- dds2[, order(dds2[[match_group]])]

}


vsd1 <- assays(dds1)[['vsd']]
vsd2 <- assays(dds2)[['vsd']]

## Select DEGs of interest
basename1 <- basename(file_path_sans_ext(files[[1]]))
de_1 <- lapply(c1, function(x) read.csv(glue('./data/DE_results/{basename1}_{x}_DE.csv'), row.names=1, header=T))
de_names1 <- lapply(de_1, function(x) rownames(x))
genes_1 <- Reduce(union, de_names1)

basename2 <- basename(file_path_sans_ext(files[[2]]))
de_2 <- lapply(c2, function(x) read.csv(glue('./data/DE_results/{basename2}_{x}_DE.csv'), row.names=1, header=T))
de_names2 <- lapply(de_2, function(x) rownames(x))
genes_2 <- Reduce(union, de_names2)

res <- calc_cor(vsd1[genes_1, ], vsd2[genes_2, ])

rho <- res$rho
p   <- res$p

write.csv(x=rho, file=glue('./data/corr_results/{basename1}_{basename2}_rho.csv'))
write.csv(x=p,   file=glue('./data/corr_results/{basename1}_{basename2}_p.csv'))


rho_thresh <- threshold_cor(rho, p)
write.csv(x=rho_thresh, file=glue('./data/corr_results/{basename1}_{basename2}_rho-thresh.csv'))


# annotate heatmap according to which genes are degs
ann1 <- lapply(de_names1, function(x) {as.numeric(colnames(rho_thresh) %in% x)}) %>% as.data.frame()
colnames(ann1) <- c1
rownames(ann1) <- colnames(rho_thresh)

ann2 <- lapply(de_names2, function(x) {as.numeric(rownames(rho_thresh) %in% x)}) %>% as.data.frame()
colnames(ann2) <- c2
rownames(ann2) <- rownames(rho_thresh)


color = colorRampPalette(c("navy", "white", "red"))(5)
hm <- pheatmap(rho_thresh, annotation_col = ann1, annotation_row = ann2, 
               breaks = c(-1, -.7, -0.5, 0.5, 0.7, 1), color = color, 
               legend_breaks=c(-.7, 0, 0.7), legend_labels = c('<-0.7', '0', '>0.7'))

pdf(glue('./fig/corr/{basename1}_{basename2}_rho-thresh.pdf'), 
    width = ncol(rho_thresh) / 6, height = nrow(rho_thresh) / 6) 

draw(hm)
dev.off()


