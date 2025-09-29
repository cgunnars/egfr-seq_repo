## script for fig4 ghi
#  generates plots for 
## req:
#  combined_1 - 3, DEG data tables relative to other conditions 

source('./src/utilities.R')
prepEnv()

# in: dataframe degs which specifies whether gene is DE / not DE / likely DE  in other conditions
#     string refname which designates the primary condition
# out: dataframe degs containing DE / not DE / likely DE counts for each other condition 
# pre: degs must have columns called category and category2, which specify whether DEG is DE in two other conditions
getCatN <- function(degs, refname) {
    category1 <- degs %>% dplyr::count(category,dir, sort=T)
    category1 <- category1 %>% pivot_longer(cols='category')
    category2 <- degs %>% dplyr::count(category2,dir , sort=T)
    category2 <- category2 %>% pivot_longer(cols='category2')
    
    degs <- rbind(category1, category2)
    degs['condition'] <- refname
    return(degs)
}

# add column to dataframe degs to specify direction of regulation 
addDir <- function(degs) {
    degs['dir'] <- (degs$FC_1 > 0)
    return(degs)
}



parser <- ArgumentParser() ## TODO: extensibility for n conditions 
parser$add_argument('-e', type='character', nargs=1) # experiment name
parser$add_argument('-a', type='character', nargs=1) # condition 1
parser$add_argument('-b', type='character', nargs=1) # condition 2
parser$add_argument('-c', type='character', nargs=1) # condition 3
parser$add_argument('-g', type='character', nargs=1) # group (aka. Drug, Day, etc)
parser$add_argument('-v', type='character', nargs=1) # control group, e.g. vehicle

args <- parser$parse_args()

data_dir   = './data/DE_results'
experiment = args$e #'20240502_pel-timecourse-6donor_pathogen'
group      = args$g

ctrl       = args$v #'DMSO_d1'
ref1       = args$a #'pel_d1'
ref2       = args$b #'gef_d1'
ref3       = args$c #'sara_d1'
refs_all   = c(ref1, ref2, ref3)
conditions = c(ref1, ref2, ref3, ctrl)


# read in categorization tables for each condition's DEGs
combined_1 <- read.csv(glue('{data_dir}/combined/{experiment}_{ref1}_{ref2}_{ref3}.csv'))
combined_2 <- read.csv(glue('{data_dir}/combined/{experiment}_{ref2}_{ref1}_{ref3}.csv'))
combined_3 <- read.csv(glue('{data_dir}/combined/{experiment}_{ref3}_{ref1}_{ref2}.csv'))
combined_all <- list(combined_1, combined_2, combined_3)


degs_all     <- lapply(combined_all,
                       function(x) return(x[x$DE_1,]))
degs_all     <- lapply(degs_all, function(x) addDir(x))
cat_all      <- lapply(seq(3),
                       function(i) getCatN(degs_all[[i]],
                                           refs_all[[i]]))
cat_df       <- do.call('rbind', cat_all)
deg_df       <- do.call('rbind', degs_all)

# set factor levels for plotting 
cat_df$condition <- factor(cat_df$condition, levels=c(ref1, ref2, ref3))
cat_df$dir       <- factor(cat_df$dir, levels=c('TRUE','FALSE'))



## fig 2g
n_df <- cat_df %>% group_by(condition, dir) %>% summarize(n=sum(n)/2) # divide by two to account for double entries
ggplot(data = cat_df, aes(x=name, y=n, fill=value)) + 
    geom_bar(position='fill', stat='identity') + 
    scale_fill_manual(values=getCategoryColors(), 
                      aesthetics='fill') +
    facet_grid(vars(dir), vars(condition)) +
    theme_classic()
ggsave(glue('./fig/combined_bar/{experiment}_fill.pdf'),
        width=6, height=3, dpi=300, units='in')

ggplot(data = n_df, aes(x=condition, y=n)) +
       geom_bar(position='identity', stat='identity') +
       facet_grid(vars(dir)) + theme_classic()
ggsave(glue('./fig/combined_bar/{experiment}_n.pdf'), 
       width=6, height=3, dpi=300, units='in')



degs_1  <- degs_all[[1]]
degs_2  <- degs_all[[2]]
degs_3 <- degs_all[[3]]

## for each gene, categorize as condition-specific or shared
likely = c('DE 1, DE 2 likely', 'DE 1, DE 2')
unlikely = 'DE 1, DE 2 unlikely'
low_info = 'DE 1, DE 2 not excluded'
unlikely_all = c(unlikely, low_info) #allow a gene to be likely exclusive if it's unlikely in one and low info in the other
degs_1likely <- unique(c(degs_1$X, 
                         degs_2[degs_2$category %in% likely, 'X'], 
                         degs_3[degs_3$category %in% likely, 'X'])
                      )
degs_1unique <- degs_1[((degs_1$category  == unlikely & degs_1$category2 %in% unlikely_all) |
                        (degs_1$category2 == unlikely & degs_1$category  %in% unlikely_all)), 'X']
degs_2likely <- unique(c(degs_2$X,
                         degs_1[degs_1$category %in% likely, 'X'],
                         degs_3[degs_3$category2 %in% likely, 'X']))
degs_2unique <- degs_2[((degs_2$category  == unlikely & degs_2$category2 %in% unlikely_all) | 
                        (degs_2$category2 == unlikely & degs_2$category  %in% unlikely_all)), 'X']
degs_3likely <- unique(c(degs_3$X,
                         degs_1[degs_1$category2 %in% likely, 'X'],
                         degs_2[degs_2$category2 %in% likely, 'X']))
degs_3unique <- degs_3[((degs_3$category  == unlikely & degs_3$category2 %in% unlikely_all) |
                        (degs_3$category2 == unlikely & degs_3$category  %in% unlikely_all)), 'X']

all_deglikely <- list(degs_1likely,
                      degs_2likely,
                      degs_3likely)
all_deg <-       list(degs_1likely,
                      degs_2likely,
                      degs_3likely,
                      degs_1unique,
                      degs_2unique,
                      degs_3unique)
shared = Reduce(intersect, all_deglikely)
writeLines(shared, glue('./data/DE_results/{experiment}_likely_shared.txt'))
writeLines(degs_1unique, glue('./data/DE_results/{experiment}_{ref1}_unique.txt'))
writeLines(degs_2unique, glue('./data/DE_results/{experiment}_{ref2}_unique.txt'))
writeLines(degs_3unique, glue('./data/DE_results/{experiment}_{ref3}_unique.txt'))
names(all_deg) <- c(refs_all, 'unique1', 'unique2', 'unique3')

drug_colors <- getDrugColors()
unique_colors <- c('unique1'='#FFFFFF', 'unique2'='#FFFFFF', 'unique3'='#FFFFFF')

venn_colors <- c(drug_colors, unique_colors)
venn_colors <- venn_colors[names(all_deglikely)]
## TODO: fix venn colors... or don't
shared_venn <- plot(euler(all_deg, labels=names(all_deglikely), shape='ellipse'), fills=venn_colors, quantities=T)
ggsave(glue('./fig/combined_bar/{experiment}_vennlikely.pdf'), shared_venn, width=8, height=8)


## draw unique and shared heatmaps
dds <- readRDS(glue('./data/DE_results/{experiment}.Rds'))
hm_unique <- plotBasicHeatmap(degs_1unique, dds, group, conditions)
pdf(file=glue('./fig/unique-shared_heatmap/{experiment}_{ref1}.pdf'), 
    height=getHeight(degs_1unique))#, units='in', res=480)
draw(hm_unique)
dev.off()


hm_unique <- plotBasicHeatmap(degs_2unique, dds, group, conditions)
pdf(file=glue('./fig/unique-shared_heatmap/{experiment}_{ref2}.pdf'), 
    height=getHeight(degs_2unique))#, units='in', res=480)
draw(hm_unique)
dev.off()


hm_unique <- plotBasicHeatmap(degs_3unique, dds, group, conditions)
pdf(file=glue('./fig/unique-shared_heatmap/{experiment}_{ref3}.pdf'),
    height=getHeight(degs_3unique))#, units='in', res=480)
draw(hm_unique)
dev.off()


hm_shared <- plotBasicHeatmap(Reduce(intersect,all_deglikely), dds, group, conditions)
pdf(file=glue('./fig/unique-shared_heatmap/{experiment}_shared.pdf'),
    height=getHeight(Reduce(intersect,all_deglikely)))#, units='in', res=480)
draw(hm_shared)
dev.off()


## reassemble truth table and write to file
rownames <- combined_all[[1]]$X 
data <- lapply(seq(3), function(i) combined_all[[i]][, c('FC_1', 'DE_1', 'notDE_1')])
data <- do.call('cbind', data)
rownames(data) <- rownames

colnames = lapply(refs_all, function(i) c(glue('FC_{i}'), glue('DE_{i}'), glue('notDE_{i}'))) %>% unlist(.)


colnames(data) <- colnames
write.csv(data, glue('./data/DE_results/combined/combined_{experiment}_full-category.csv'))
