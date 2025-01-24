getDE <- function(dds, comparison) {
    #write output full DESeq2 datatable
    #write output select DEG datatable
    #write output list of up genes
    #write output list of down genes
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
dds <- DESeq2(dds)

for(c in comparisons){
    #get DE and write to file
    #make volcano plots and write to file
}




