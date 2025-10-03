source('src/utilities.R')
library(assertthat)
prepEnv()


parser <- ArgumentParser()
parser$add_argument('-i', nargs='+', type='character') #experiment name(s), max 2, "reference" should be listed first
parser$add_argument('-r', nargs='+', type='character') #ctrl condition(s), max 2
parser$add_argument('-c', nargs='+', type='character') #exp conditions, "reference" should be listed first
args <- parser$parse_args()

exp   = args$i
ctrl  = args$r
conds = args$c

data_dir = './data/DE_results'

if (length(exp) == 1) { #if only one experiment is provided, everything else should match
    assert_that(length(ctrl) == 1)
    ref <- conds[[1]]
    conds_other <- conds[2:length(conds)]
    coef_ref   <- glue('{ref}_vs_{ctrl}')
    coef_other <- lapply(conds_other, function(cond) glue('{cond}_vs_{ctrl}')) 

    # for each DEG in the reference condition, 
    # calculate statistical categories that gene (DE, not DE etc) for tested conditions
    combined_list <- lapply(seq(coef_other), function(i) {
                                df <- compareDEGs(exp, exp, coef_ref, coef_other[[i]], mode='joint')
                                coef_n  <- i+1
                                cols_ro <- c(unlist(lapply(c('1', glue('{coef_n}'), glue('1{coef_n}')), 
                                                    function(j) { 
                                                        c(glue('FC_{j}'), 
                                                          glue('DE_{j}'), 
                                                          glue('notDE_{j}'), 
                                                          glue('baseMean_{j}'), 
                                                          glue('padj_{j}'))}
                                                    )     ), glue('category{i}'))
                                colnames(df) <- cols_ro
                                return(df)
                           })
    combined_ref <- do.call(cbind, combined_list)
    combined_ref <- combined_ref[, !duplicated(colnames(combined_ref))]
    print(head(combined_ref))

    outname <- paste(conds, collapse = '_')

    print(outname) 
    write.csv(combined_ref, glue('{data_dir}/combined/{exp}_{outname}.csv'))
} else {
    assert_that(length(exp) == 2)
    assert_that(length(ctrl) == 2)
    assert_that(length(conds) == 2)

    exp_1 <- exp[[1]]
    exp_2 <- exp[[2]]
    coefs <- lapply(seq(ctrl), function(i) glue('{conds[i]}_vs_{ctrl[i]}'))

    combined_ia <- compareDEGs(exp_1, exp_2, coefs[[1]], coefs[[2]], mode='axenic')
    write.csv(combined_ia, glue('{data_dir}/combined/combined_intraaxenic_{conds[1]}_{exp[2]}.csv'))
}
