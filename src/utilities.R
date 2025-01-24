prepEnv <- function() {
    options(max.print=999999)
    library(argparse)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(glue)
    library(pheatmap)
    library(cowplot)
    library(ggplot2)
    library(tools)
    library(EnhancedVolcano)
    library(DESeq2)
    library("vsn")
    library(RNAseqQC)
    library(ComplexHeatmap)
}

loadHost <- function(host_file) {
    read_data <- readRDS(host_file)
    return(read_data)
}

loadLocus <- function(locus_file) {
    locus_df <- read.csv(locus_file, row.names='Locus')
    return(locus_df)
}

loadPathogen <- function(pathogen_file, locus_file) {
    pathogen_data <- read.csv(pathogen_file, sep='\t', row.names='gene_id')
    #change from rv number to gene names
    locus_df      <- loadLocus(locus_file)
    gene_names    <- locus_df[row.names(pathogen_data), 'Name']
    # where na mapping, revert to original row names in data file 
    gene_names[is.na(gene_names)] <- row.names(pathogen_data)[is.na(gene_names)]
    row.names(pathogen_data)      <- gene_names

    read_data <- pathogen_data %>% select(ends_with('NumReads'))
    colnames(read_data) <- sub('_NumReads', '', colnames(read_data))
    tpm_data  <- pathogen_data %>% select(ends_with('TPM'))
    colnames(tpm_data) <- sub('_TPM', '', colnames(tpm_data))

    return(list(read_data, tpm_data))
}

#inputs: sample_list, categories for sample data
#        design, formula -- would be Drug_Day, Drug_Dose, or Drug
makeDDS <- function(file, sample_list, design, is_pathogen=T, locus_file=NA) {
    if(is_pathogen) {
        pathogen_info <- loadPathogen(file, locus_file)
        read_data <- pathogen_info[[1]]
        tpm_data  <- pathogen_info[[2]]
    } else {
        read_data <- loadHost(file)
        tpm_data  <- NA
    }
    
    # drop genes and samples with no reads 
    read_data <- read_data[rowSums(read_data, na.rm=T) > 0, colSums(read_data, na.rm=T) > 0 ]
    
    samples             <- data.frame(sample=colnames(read_data))
    samples_df          <- samples %>% separate(sample, sample_list)
    row.names(samples_df) <- colnames(read_data)

    # for design, make super-category separated by "_", e.g. Drug_Dose pel_5
    design_values <- lapply(design, function(x) samples_df[[x]])
    design_name   <- paste(design, collapse='_')
    samples_df[design_name] <- Reduce(function(x, y) paste0(x, '_', y), 
                                      design_values)
    

    dds <- DESeqDataSetFromMatrix(countData = read_data,
                                  colData   = samples_df,
                                  design    = as.formula(glue('~ {design_name}')))
    if(!sum(is.na(tpm_data))) {
        tpm_data <- tpm_data[rowSums(tpm_data, na.rm=T) > 0, colSums(tpm_data, na.rm=T) > 0 ]
        assays(dds)[['tpm']] <- tpm_data
    }
    return(dds)
}

# input: two scaled matrices
calc_cor <- function(vsd1, vsd2) {
    rho_cor <- matrix(nrow=nrow(vsd1), ncol(vsd2))
    p_cor   <- matrix(nrow=nrow(vsd1), ncol(vsd2))

    rho_cor <- apply(assay(vsd1), 1, function(col_mat1){
                 apply(assay(vsd2), 1, function(col2, col1){
                    cor.test(col2, col1)$estimate
                 }, col1=col_mat1)
               })

    p_cor   <- apply(assay(vsd1), 1, function(col_mat1){
                 apply(assay(vsd2), 1, function(col2, col1){
                    cor.test(col2, col1)$p.value
                 }, col1=col_mat1)
               })
    return(list(rho_cor, p_cor))
}

readtxt <- function(file) {
    con <- file(file, open="r")
    txt <- readLines(con)
    close(con)
    return(txt)
}
