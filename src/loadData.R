library(argparse)
library(tools)
source('src/utilities.R')

prepEnv()

parser <- ArgumentParser(description='File i/o for loading counts data')
parser$add_argument('-i', '--input_file', type='character', nargs=1, help='input counts data')
parser$add_argument('-l', '--locus_file', nargs=1, default='./data/gene_info/H37Rv_gene-lengths.csv')
parser$add_argument('--pathogen', action='store_true', default=F, help='is pathogen data?')
parser$add_argument('-s', '--sample_info', type='character', nargs='+', help='sample info by category, e.g. Drug Dose Replicate') 
parser$add_argument('-d', '--design', type='character', nargs='+', help='design variables to combine')
parser$add_argument('-m', '--meta', type='character', required=F, nargs='+', help='additional annotations')

args        <- parser$parse_args()
file        <- args$i
is_pathogen <- args$pathogen
sample_list <- args$s
design      <- args$d
locus_file  <- args$l
meta        <- args$m

dds      <- makeDDS(file, sample_list, design, is_pathogen, locus_file)

filebase <- basename(file_path_sans_ext(file))
outstem = './data/raw_dds/'
outext  = '.Rds'

# meta data is in format Name Value, split every other, then combine into named vector
if(length(meta) > 0) {
    meta_cols <- meta[seq(1,length(meta),2)]
    meta_vals <- meta[seq(2,length(meta),2)]
    names(meta_vals) <- meta_cols
    colData(dds) <- cbind(colData(dds), t(meta_vals))
}

saveRDS(dds, glue('{outstem}{filebase}{outext}'))

