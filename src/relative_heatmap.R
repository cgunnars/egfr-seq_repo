source('src/utilities.R')
prepEnv()


parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #experiment name
parser$add_argument('-r', nargs=1, type='character') #reference experimental condition
parser$add_argument('-c', nargs=3, type='character') #exp condition 1
parser$add_argument('-g', nargs=1, type='character') #group
args <- parser$parse_args()

group = args$g
ctrl  = args$r
conds = args$c

ref   = conds[[1]] 
cond1 = conds[[2]]
cond2 = conds[[3]]
conditions = c(ctrl, ref, cond1, cond2)
coef_ref = glue('{group}_{ref}_vs_{ctrl}')
coef_1   = glue('{group}_{cond1}_vs_{ctrl}')
coef_2   = glue('{group}_{cond2}_vs_{ctrl}')


experiment = args$i
fig_dir = './fig/relative_heatmap'
data_dir = './data/DE_results' 

dds_intra <- readRDS(glue('./{data_dir}/{experiment}.Rds'))


# for each DEG in the reference condition, 
# calculate statistical categories that gene (DE, not DE etc) for tested conditions
combined_r1 <- getDEGs_comparison(dds_intra, dds_intra, coef_ref, coef_1, group1=group, group2=group, mode='joint')
cols_r1 <- c(unlist(lapply(c('1', '2', '12'), 
                    function(i) { 
                        c(glue('FC_{i}'), 
                          glue('DE_{i}'), 
                          glue('notDE_{i}'), 
                          glue('baseMean_{i}'), 
                          glue('padj_{i}'))}
                    )     ), 'category')
colnames(combined_r1) <- cols_r1

combined_r2 <- getDEGs_comparison(dds_intra, dds_intra, coef_ref, coef_2, group1=group, group2=group, mode='joint')
cols_r2 <- c(unlist(lapply(c('1', '3', '13'), 
                    function(i) { 
                        c(glue('FC_{i}'), 
                          glue('DE_{i}'), 
                          glue('notDE_{i}'), 
                          glue('baseMean_{i}'), 
                          glue('padj_{i}'))}
                    )     ), 'category2')
colnames(combined_r2) <- cols_r2

# bind comparisons, omit first five columns to avoid duplicating reference data
combined_ref12 <- cbind(combined_r1, combined_r2[6:ncol(combined_r2)])
write.csv(combined_ref12, glue('{data_dir}/combined/{experiment}_{ref}_{cond1}_{cond2}.csv'))

hm_ref12 <- plotAxenicHeatmap(dds_intra, combined_ref12, group, conditions=conditions, mode='ref')

fig_height = 5 + length(combined_ref12[combined_ref12$DE_1, ]) / 2.5 
fig_width  = length(conditions) * 2 

pdf(file=glue('{fig_dir}/relative_heatmap_{experiment}_{ref}_{cond1}_{cond2}.pdf'), 
    width=fig_width, height=fig_height)#, units='in', res=480)
draw(hm_ref12)
dev.off()
