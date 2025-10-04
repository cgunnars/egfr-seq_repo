# for three conditions relative to the same shared ctrl in the same experiment,
# plot a heatmap of DEGs in condition 1 + DE status in other two conditions 
# i.e. mappings of "shared DEG with condition 2?" 

source('src/utilities.R')
prepEnv()


parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #experiment name
parser$add_argument('-r', nargs=1, type='character') #ctrl condition
parser$add_argument('-c', nargs=3, type='character') #exp conditions, "reference" should be listed first
args <- parser$parse_args()

ctrl  = args$r
conds = args$c

ref   = conds[[1]] 
cond1 = conds[[2]]
cond2 = conds[[3]]
conditions = c(ctrl, ref, cond1, cond2)

experiment = args$i
fig_dir = './fig/relative_heatmap'
data_dir = './data/DE_results' 

dds_intra <- readRDS(glue('./{data_dir}/{experiment}.Rds'))

combined_ref12 <- read.csv(glue('{data_dir}/combined/{experiment}_{ref}_{cond1}_{cond2}.csv'), row.names=1)

hm_ref12 <- plotAxenicHeatmap(dds_intra, combined_ref12, conditions=conditions, mode='ref')

fig_height = 5 + length(combined_ref12[combined_ref12$DE_1, ]) / 2.5 
fig_width  = length(conditions) * 2 

pdf(file=glue('{fig_dir}/relative_heatmap_{experiment}_{ref}_{cond1}_{cond2}.pdf'), 
    width=fig_width, height=fig_height)#, units='in', res=480)
draw(hm_ref12)
dev.off()
