source('./src/utilities.R')

prepEnv()

parser <- ArgumentParser()
parser$add_argument('-i', nargs=1, type='character') #experiment name
parser$add_argument('-c', nargs=1, type='character') #cond name

data_dir = './data/DE_results'

args <- parser$parse_args()
exp  <- args$i
cond <- args$c

ref   = 'DMSO_d1'
dds <- readRDS(glue('{data_dir}/{exp}.Rds')) 
degs <- read.csv(glue('{data_dir}/{exp}_{cond}_vs_{ref}_DE.csv'), row.names=1) %>% row.names()

all_fc   <- lapply(c(cond, ref),
		   function(x) {df <- read.csv(glue('{data_dir}/{exp}_{x}_vs_phago_4h_full.csv', row.names=1))
		   	       rownames(df) <- df$X
		   	       return(df)
		   	       }
		  )
fc_cond  <- all_fc[[1]]
fc_dmso  <- all_fc[[2]]
rows     <- degs
fc_cond <- cbind(fc_cond[rows, 'log2FoldChange'], fc_dmso[rows, 'log2FoldChange'], 
		 ifelse(abs(fc_cond[rows, 'log2FoldChange']) > 1 & abs(fc_dmso[rows,'log2FoldChange']) > 1, 1, 0)) 
rownames(fc_cond) <- rows
colnames(fc_cond) <- c(glue('FC_{cond}'), glue('FC_{ref}'), 'shared_DE')
fc_cond <- fc_cond %>% as.data.frame()

anno_cols <- list(getFCColors(),
		  getFCColors(),
		  c('0' = 'white', '1'='black'))
names(anno_cols) <- colnames(fc_cond)
		  
anno_df <- fc_cond 
hm <- plotBasicHeatmap(degs, dds, c('phago_4h', cond, ref), cluster_columns=F, row_km=4)
hm <- rowAnnotation(df = anno_df, col=anno_cols) + hm
	      

pdf(glue('./fig/time-dependent/phago_{exp}_{cond}_hm.pdf'), width=7, height=getHeight(degs))
draw(hm)
dev.off()




