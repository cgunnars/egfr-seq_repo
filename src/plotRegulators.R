# generalize to pel
makeGeneDf <- function(dds) {
    vsd <- assays(dds)[['vsd']]
    genedf <- assay(vsd) %>% data.frame() %>% rownames_to_column(var='row') %>% 
              pivot_longer(., cols=-row) %>% data.frame()
    meta <-   str_split_fixed(genedf$name, pattern='_', 3) %>% data.frame()
    colnames(meta) <- c('Drug', 'Day', 'Donor')

    genedf$Drug <- meta$Drug
    genedf$Drug <- factor(genedf$Drug, levels=c('pel','gef','sara','DMSO', 'phago', 'tb'))
    genedf$Day  <- meta$Day
    genedf$Donor <- meta$Donor
    genedf$Drug_Donor <- paste(meta$Drug, meta$Donor, sep='_')
    genedf$Donor_Day  <- as.factor(paste(meta$Donor, meta$Day, sep='_'))
    return(genedf)
}


plotUpregulation <- function(genedf, de) { 
    cols <- getDrugColors()
    p <- ggplot(data=genedf, 
                aes(x=Drug, y=value, col=Drug)) + 
         geom_boxplot(width=0.1, position=position_dodge(0.5)) + 
         geom_point(position=position_dodge(0.5)) +
         ylab('vsd') +
         scale_color_manual(values=cols) + 
         geom_line(aes(group=interaction(Donor, Day)), 
                   col='black', position=position_dodge(0)) + 
	 geom_text(data = de, inherit.aes=F,
		   aes(x=1.5, y=value,
		       label=ifelse(padj < 0.05 & abs(log2FoldChange) > 0.7, 
		       		    sprintf('FC %.2f\np %s', 
			       	            signif(log2FoldChange,2), scientific(padj,digits=2)),
			       	    ""))) +
	 theme_classic()
    g <- p + facet_wrap(vars(row), ncol=9)
    return(g)
}

source('./src/utilities.R')
prepEnv()

## TODO: extensibility by setting group argument here
parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #experiment name
parser$add_argument('-c', nargs=1, type='character') #cond name
parser$add_argument('-r', nargs=1, type='character')
#parser$add_argument('-g', nargs=1, type='character')

args  <- parser$parse_args()
exp   <- args$i
cond  <- args$c
ref   <- args$r
#group <- args$g


wd = './data/DE_results'
dds   <- readRDS(glue('{wd}/{exp}.Rds'))

sigma = c('sigA', 'sigB', 'sigC', 'sigD', 'sigE', 'sigF', 'sigG', 'sigH', 'sigK')

tcs   = c('mprA', 'mprB', 'kdpD', 'kdpE', 'tcrX', 'tcrY', 'mtrA', 'mtrB',
          'phoP', 'phoR', 'senX3', 'regX3',
          'prrA', 'prrB', 'pdtaR', 'pdtaS', 'trcR', 'trcS', 'narS', 'narL',
          'devR', 'devS', 'dosT')

day   = 'd1'
de    <- read.csv(glue('{wd}/{exp}_{cond}_{day}_vs_{ref}_{day}_full.csv'))
colnames(de) <- c('row', colnames(de)[2:length(colnames(de))])

genedf   <- makeGeneDf(dds)
sigma_df <- genedf[genedf$row %in% sigma & 
                   genedf$Drug %in% c(cond, ref) &
                   genedf$Day == day, ]
tcs_df   <- genedf[genedf$row %in% tcs & 
                   genedf$Drug %in% c(cond, ref) &
                   genedf$Day == day, ]

de_sigma <- de[de$row %in% sigma, ]
de_sigma$value <- min(sigma_df$value)+.5

de_tcs <- de[de$row %in% tcs, ]
de_tcs$value <- min(tcs_df$value)+.5
	   
pa <- plotUpregulation(sigma_df, de_sigma)
pb <- plotUpregulation(tcs_df, de_tcs)

a  <- ggsave(glue('./fig/regulators/sigma_{cond}.pdf'), pa,
                  width=9, height=2.5, create.dir=T)
                                                  
b  <- ggsave(glue('./fig/regulators/tcs_{cond}.pdf'), pb, 
      	     width=9, height=7.5, create.dir=T)
