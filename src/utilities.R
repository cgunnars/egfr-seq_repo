prepEnv <- function() {
    options(max.print=999999)
    library(argparse,quietly=T)
    library(dplyr,quietly=T)
    library(tidyr,quietly=T)
    library(tibble,quietly=T)
    library(stringr,quietly=T)
    library(glue,quietly=T)
    library(eulerr,quietly=T)
    library(pheatmap,quietly=T)
    library(cowplot,quietly=T)
    library(ggplot2,quietly=T)
    library(tools,quietly=T)
    library(EnhancedVolcano,quietly=T)
    library(DESeq2,quietly=T)
    library(apeglm,quietly=T)
    library("vsn",quietly=T)
    library(RNAseqQC,quietly=T)
    library(ComplexHeatmap,quietly=T)
    library(circlize,quietly=T)
    library(RColorBrewer,quietly=T)
    library(scales,quietly=T)
    library(clusterProfiler,quietly=T)
    library(enrichplot,quietly=T)
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


remove_pseudo <- function(genes) {
    genes <- grep(pattern = "^{^RP[0-9]|^CTD-|^AC0[0-9]}", genes, 
                  value=T, invert=T)
    return(genes)
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
    res     <- list(rho_cor, p_cor)
    names(res) <- c('rho', 'p')
    return(res)
}

readtxt <- function(file) {
    con <- file(file, open="r")
    txt <- readLines(con)
    close(con)
    return(txt)
}



# for timecourse with n_timepoints,
# return list of up+down genes + a plottable dataframe of n genes vs time
# additional input: test condition (cond) + references to compare
#                   experiment name (assume standard file structure)
##  helper function for plotDEGs_time

#   references must be provided to specify "adjacent" conditions along the timecourse
#   e.g. pel_d3 can be matched to DMSO_d3 or pel_d2
#   pre: DE analysis has already been run, and txt files w/ DEGs in format wd/exp_cond_d1_vs_ref_d1_up.txt exist
#        timepoints in _dN format with 1-day intervals between
getDEGs_time <- function(wd='./data/DE_results', cond, exp, refs, n_timepoints=3) {
  up_genes <- lapply(seq(length(refs)), 
                    function(i) {readLines(
                    file(glue('{wd}/{exp}_{cond}_d{i}_vs_{refs[i]}_up.txt'))
                   )
                   }
                   )
  names(up_genes) <- lapply(seq(length(refs)), function(i) glue('{cond}_d{i}_vs_{refs[i]}'))
  down_genes <- lapply(seq(length(refs)), 
                        function(i) {readLines(
                        file(glue('{wd}/{exp}_{cond}_d{i}_vs_{refs[i]}_down.txt'))
                    )
                    }
                    )
  names(down_genes) = names(up_genes)

  if (n_timepoints > 1) {
      num_up     <- unlist(lapply(up_genes, function(i) length(i)))
      num_down   <- unlist(lapply(down_genes, function(i) length(i)))

      df <- data.frame(day = seq(n_timepoints), n_de = num_up, condition = cond, direction = 'up', ref=refs)
      df <- rbind(df, 
              data.frame(day = seq(n_timepoints), n_de = num_down, condition = cond, direction = 'down', ref=refs))
      
      out = list(up_genes=up_genes, down_genes=down_genes, df=df)
  } else {
      out = list(up_genes=up_genes, down_genes=down_genes)
  }
  return(out)
}


# for a comparison (coef), test whether is DE (> abs(log2fcthresh_high & < pthresh) or 
#                          test whether is NOT DE (< fcthresh_low by Wald's test)
getDEGs_category <- function(dds, coef, group=NA, pthresh = 0.05, fcthresh_high=1, fcthresh_low=0.7) {
    # set the reference condition if needed (default reference = DMSO)
    if (!coef %in% resultsNames(dds)) {
        ref_split <- coef %>% str_replace(., glue('{group}_'), '') %>% str_split(., '_vs_')
        ref       <- ref_split[[1]][2]
        dds[[group]] <- relevel(factor(dds[[group]]), ref=ref)
        dds <- DESeq(dds)
    }
    # test DE and not DE
    res_shr <- lfcShrink(dds, coef, type='apeglm')
    res_LA  <- results(dds, lfcThreshold = fcthresh_low, altHypothesis = 'lessAbs', 
                       name = coef)
    res_shr$DE <- as.logical(abs(res_shr[, 'log2FoldChange']) > fcthresh_high & res_shr[, 'padj'] < pthresh)
    res_shr$DE[is.na(res_shr$DE)] <- F #set NA p values to neither DE nor not DE
    res_shr$notDE <- res_LA[, 'padj'] < pthresh
    res_shr$notDE[is.na(res_shr$notDE)] <- F

    return(res_shr)
}

# for all DEGs for a two given conditions (coef1 and coef2),
# statistically categorize all genes

## modes -- joint, within the same dds object, compare coef1 and coef2 (e.g. drug A vs drug B)
#           single, within two separate dds objects, compare coef1 and coef1 (e.g. axenic vs intra)
getDEGs_comparison <- function(dds1, dds2, coef1, coef2, group1=NA, group2=NA, fcthresh_low=0.7, mode='single') {
    res_1 <- getDEGs_category(dds1, coef1, group=group1, fcthresh_low=fcthresh_low)
    res_2 <- getDEGs_category(dds2, coef2, group=group2, fcthresh_low=fcthresh_low)
    if (mode == 'joint') {
        ## assumes that dds1 and dds2 are the same thing....
        coef_split1 <- coef1 %>% str_replace(., glue('{group1}_'), '') %>% str_split(., '_vs_')
        coef_split1 <- coef_split1[[1]][1] 
        coef_split2 <- coef2 %>% str_replace(., glue('{group1}_'), '') %>% str_split(., '_vs_')
        coef_split2 <- coef_split2[[1]][1]

        coef_12     <- glue('{group}_{coef_split1}_vs_{coef_split2}')  
        res_12   <- getDEGs_category(dds1, coef_12, group=group1, fcthresh_low=fcthresh_low)

        res_list <- list(res_1, res_2, res_12)
        resnames <- list('1', '2', '12')
    } else {
        res_list <- list(res_1, res_2)
        resnames <- list('1', '2')
    }
    combined_rownames <- Reduce(intersect, lapply(res_list, function(i) rownames(i)))
    combined_colnames <- c('log2FoldChange', 'DE', 'notDE', 'baseMean', 'padj')
    combined <- do.call(cbind,
                        lapply(res_list, function(i) i[combined_rownames, combined_colnames])
                )
    colnames(combined) <- unlist(lapply(resnames, 
                                 function(i) 
                                     c(glue('FC_{i}'), glue('DE_{i}'), glue('notDE_{i}'), glue('baseMean_{i}'), glue('padj_{i}'))
                                 ))
    combined$category <- 'mixed evidence'
    combined[combined$DE_1 & !combined$notDE_2, 'category'] <- 'DE 1, DE 2 not excluded'
    combined[combined$DE_1 & abs(combined$FC_2) > fcthresh_low & !combined$DE_2, 'category'] <- 'DE 1, DE 2 likely'
    combined[combined$DE_1 & combined$notDE_2, 'category']  <- 'DE 1, DE 2 unlikely'                          
    combined[combined$DE_2 & combined$notDE_1, 'category']  <- 'DE 2, DE 1 unlikely'
    combined[combined$DE_2 & !combined$notDE_1, 'category'] <- 'DE 2, DE 1 not excluded'
    combined[combined$DE_2 & abs(combined$FC_1) > fcthresh_low & !combined$DE_1, 'category'] <- 'DE 2, DE 1 likely'
    combined[combined$DE_1 & combined$DE_2, 'category']     <- 'DE 1, DE 2'
    combined[combined$notDE_1 & combined$notDE_2, 'category'] <- 'low change'

    # for joint mode, use differential expression information to triage DEGs 
    if (mode == 'joint') { 
        combined[combined$DE_1 & !combined$notDE_2 & !combined$DE_2 & combined$DE_12, 'category'] <- 'DE 1, DE 2 unlikely'
        combined[combined$DE_2 & !combined$notDE_1 & !combined$DE_1 & combined$DE_12, 'category'] <- 'DE 2, DE 1 unlikely'
    } 
    return(combined)
}

compareDEGs <- function(stem1, stem2, coef1, coef2, mode='axenic', fcthresh_low=0.7) {
    res_1 <- read.csv(glue('./data/DE_results/{stem1}_{coef1}_full.csv'), row.names=1)
    res_2 <- read.csv(glue('./data/DE_results/{stem2}_{coef2}_full.csv'), row.names=1)
    if (mode == 'joint' & stem1 == stem2) {
        ## assumes that dds1 and dds2 are the same thing....
        coef_split1 <- coef1 %>% str_split(., '_vs_')
        coef_split1 <- coef_split1[[1]][1] 
        coef_split2 <- coef2 %>% str_split(., '_vs_')
        coef_split2 <- coef_split2[[1]][1]

        ## we require at least one of these to exist
        coef_12 <- glue('{coef_split1}_vs_{coef_split2}')
        file_12 <- glue('./data/DE_results/{stem1}_{coef_12}_full.csv')
        coef_21 <- glue('{coef_split2}_vs_{coef_split1}')
        file_21 <- glue('./data/DE_results/{stem1}_{coef_21}_full.csv')
        if (file.exists(file_12)) {  
            res_12  <- read.csv(file_12, row.names=1)
        } else if (file.exists(file_21)) {
            res_21  <- read.csv(file_21, row.names=1)
            res_21$log2FoldChange <- -res_21$log2FoldChange
            res_12 <- res_21
        } 
    
        res_list <- list(res_1, res_2, res_12)
        resnames <- list('1', '2', '12')
    } else {
        res_list <- list(res_1, res_2)
        resnames <- list('1', '2')
    }
    combined_rownames <- Reduce(intersect, lapply(res_list, function(i) rownames(i)))
    combined_colnames <- c('log2FoldChange', 'DE', 'notDE', 'baseMean', 'padj')
    combined <- do.call(cbind,
                        lapply(res_list, function(i) i[combined_rownames, combined_colnames])
                )
    colnames(combined) <- unlist(lapply(resnames, 
                                 function(i) 
                                     c(glue('FC_{i}'), glue('DE_{i}'), glue('notDE_{i}'), glue('baseMean_{i}'), glue('padj_{i}'))
                                 ))
    combined$category <- 'mixed evidence'
    combined[combined$DE_1 & !combined$notDE_2, 'category'] <- 'DE 1, DE 2 not excluded'
    combined[combined$DE_1 & abs(combined$FC_2) > fcthresh_low & !combined$DE_2, 'category'] <- 'DE 1, DE 2 likely'
    combined[combined$DE_1 & combined$notDE_2, 'category']  <- 'DE 1, DE 2 unlikely'                          
    combined[combined$DE_2 & combined$notDE_1, 'category']  <- 'DE 2, DE 1 unlikely'
    combined[combined$DE_2 & !combined$notDE_1, 'category'] <- 'DE 2, DE 1 not excluded'
    combined[combined$DE_2 & abs(combined$FC_1) > fcthresh_low & !combined$DE_1, 'category'] <- 'DE 2, DE 1 likely'
    combined[combined$DE_1 & combined$DE_2, 'category']     <- 'DE 1, DE 2'
    combined[combined$notDE_1 & combined$notDE_2, 'category'] <- 'low change'

    # for joint mode, use differential expression information to triage DEGs 
    if (mode == 'joint') { 
        combined[combined$DE_1 & !combined$notDE_2 & !combined$DE_2 & combined$DE_12, 'category'] <- 'DE 1, DE 2 unlikely'
        combined[combined$DE_2 & !combined$notDE_1 & !combined$DE_1 & combined$DE_12, 'category'] <- 'DE 2, DE 1 unlikely'
    } 
    return(combined)
}

### PLOTTING
plotAxenicHeatmap <- function(dds, combined, group, conditions, mode='axenic') {
    set.seed(3)
    vsd  <- assays(dds)[['vsd']]
    
    col_fun = getFCColors()
   
    if (mode == 'axenic') {
        degs <- rownames(combined[combined$DE_1, ])
    
        annotation_row <- data.frame(category=combined[degs, 'category'], 
                                     fc_axenic=(combined[degs, 'FC_2']), row.names=degs)
        col = list(fc_axenic = col_fun,
               category = hue_pal()(4))

        names(col$category) <- c('DE 1, DE 2', 'DE 1, DE 2 likely', 
                                 'DE 1, DE 2 unlikely', 'DE 1, DE 2 not excluded')
        gap = c(1,5,5)
    } else if (mode == 'joint') {
        degs <- rownames(combined[combined$DE_1 | combined$DE_2, ])
        annotation_row <- data.frame(
                                     fc_1=combined[degs, 'FC_1'], 
                                     fc_2=combined[degs, 'FC_2'],
                                     category=combined[degs,'category'],
                                     row.names=degs)
        col = list(fc_1 = col_fun, fc_2 = col_fun, 
                   category = c('#A000C0', '#DF86C4', 
                                '#9B7BEA', '#FFD9D0', 
                                '#CCDFFF', '#F8766D', '#619CFF'))

        names(col$category) <- c('DE 1, DE 2', 'DE 1, DE 2 likely', 
                                 'DE 2, DE 1 likely', 'DE 1, DE 2 not excluded', 
                                 'DE 2, DE 1 not excluded', 'DE 1, DE 2 unlikely', 
                                 'DE 2, DE 1 unlikely') 
        col$category <- as.factor(col$category) #, levels=names(col$category))
        gap = c(1, 5, 1, 5)
    } else if (mode == 'ref') {
        degs <- rownames(combined[combined$DE_1, ])
        annotation_row <- data.frame(
                                     fc_1 = combined[degs, 'FC_1'],
                                     fc_2 = combined[degs, 'FC_2'],
                                     fc_3 = combined[degs, 'FC_3'],
                                     category  =  combined[degs, 'category'],
                                     category2 =  combined[degs, 'category2'])
        col = list(fc_1 = col_fun, fc_2 = col_fun, fc_3 = col_fun, 
                   category = getCategoryColors(), 
                   category2 = getCategoryColors())
        gap = c(4,1,1,1,1,4)
    } else if (mode == 'single') {
        degs <- rownames(combined[combined$DE_1, ])
        
    }
    mat_subset <- vsd[degs, vsd[[group]] %in% conditions]
    scaled_mat <- t(scale(t(assay(mat_subset))))
    mat_subset[[group]] <- factor(mat_subset[[group]], levels=conditions)
    column_split = factor(mat_subset[[group]], levels=conditions)
    if (mode == 'single'){
        hm <- Heatmap(scaled_mat, name='row z(vsd)', column_split=column_split, 
                      row_km = 3, row_km_repeats=50)
    } else {
        hm <- Heatmap(scaled_mat, name='row z(vsd)', column_dend_reorder = FALSE, 
                      column_split=column_split, row_km = 3, row_km_repeats=50,
                      left_annotation = rowAnnotation(df = annotation_row, col=col, 
                                                      gap=unit(gap, "mm"))
                      )
    }
    
    return(hm)
}

getHeight <- function(degs) {
    return(length(degs) / 5 + 4)
}

plotBasicHeatmap <- function(degs, dds, group, conditions, cluster_columns=TRUE, row_km=1) {
    set.seed(3)
    vsd <- assays(dds)[['vsd']]

    mat_subset <- vsd[degs, vsd[[group]] %in% conditions]
    scaled_mat <- t(scale(t(assay(mat_subset))))
    mat_subset[[group]] <- factor(mat_subset[[group]], levels=conditions)
    column_split = factor(mat_subset[[group]], levels=conditions)
    hm <- Heatmap(scaled_mat, name='row z(vsd)', column_split=column_split, cluster_columns=cluster_columns, row_km=row_km)
    return(hm)
}

getDrugColors <- function() {
    cols <- c('pel' = '#FC5A8D', 'pelitinib' = '#FC5A8D',  'pel_d1'  = '#FC5A8D', 
              'DMSO'= '#808080',  'DMSO_d1' = '#808080',
              'gef' = '#87de87', 'gefitinib' = '#87de87', 'gef_d1'  = '#87de87', 
              'lap'='#4169e1' ,   'var'='#6F2DA8', 'rifampicin'='#FFFFFF',
              'sara' = '#C00000', 'saracatinib' = '#C00000', 'sara_d1' = '#C00000')
}

getFCColors <- function() {
    col_fun = colorRamp2(c(-2, 0, 2), c('#018571', 'white', '#A16928'))
    return(col_fun)
}

getCategoryColors <- function(mode='ref') {
    cols        <- c('#A000C0', '#DF86C4', 
                     '#FFD9D0', '#F8766D')
    names(cols) <- c('DE 1, DE 2', 'DE 1, DE 2 likely', 
                     'DE 1, DE 2 not excluded', 'DE 1, DE 2 unlikely')
    return(cols)
}

getiModColors <- function() {
    cols <- c(brewer.pal(n=12, "Paired"), '#FFFFFF')
    names(cols) <- c('Fumarate Reductase', 'HPT-2b Induced', 'Lsr2', 'MarR', 'Mce1R', 'Rv0681', 'Rv1776c+WhiB4','Rv2488c', 'SigC', 'single_gene', 'VirS', 'Zur', 'unmapped')
    return(cols)
}
