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
print(args$c)
print(args$d)

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
print(head(vsd1))

## Select DEGs of interest
basename <- basename(file_path_sans_ext(files[[1]]))
print(basename)
de_1 <- lapply(c1, function(x) rownames(read.csv(glue('./data/DE_results/{basename}_{x}_DE.csv'), row.names=1, header=T))) 
genes_1 <- Reduce(union, de_1)

de_2 <- lapply(c1, function(x) rownames(read.csv(glue('./data/DE_results/{basename}_{x}_DE.csv'), row.names=1, header=T)))
genes_2 <- Reduce(union, de_1)
## Read in information about DEGs
#down_all1 <- lapply(files1, function(x) readtxt(x)) %>% Reduce(union, .)
#genes_1   <- union(up_all1, down_all1)

#up_all2   <- lapply(files2, function(x) readtxt(x)) %>% Reduce(union, .)
#down_all2 <- lapply(files2, function(x) readtxt(x)) %>% Reduce(union, .)
#genes_2   <- union(up_all2, down_all2)

print(head(genes_1))
#print(head(vsd1[genes_1, ]))
cor <- calc_cor(vsd1[genes_1, ], vsd2[genes_2, ])
print(max(cor[[1]]))





