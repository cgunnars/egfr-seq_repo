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


sigma = c('sigA', 'sigB', 'sigC', 'sigD', 'sigE', 'sigF', 'sigG', 'sigH', 'sigK')

tcs   = c('mprA', 'mprB', 'kdpD', 'kdpE', 'tcrX', 'tcrY', 'mtrA', 'mtrB',
          'phoP', 'phoR', 'senX3', 'regX3',
          'prrA', 'prrB', 'pdtaR', 'pdtaS', 'trcR', 'trcS', 'narS', 'narL',
          'devR', 'devS', 'dosT')

day   = 'd1'
de    <- read.csv(glue('data/DE_results/{exp}_{cond}_{day}_vs_{ref}_{day}_full.csv'))
colnames(de) <- c('row', colnames(de)[2:length(colnames(de))])

genedf   <- read.csv(glue('data/clean_dds/{exp}_df.csv'))
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

a  <- ggsave(glue('./fig/regulators/sigma_{exp}_{cond}.pdf'), pa,
                  width=9, height=2.5, create.dir=T)
                                                  
b  <- ggsave(glue('./fig/regulators/tcs_{exp}_{cond}.pdf'), pb, 
      	     width=9, height=7.5, create.dir=T)
