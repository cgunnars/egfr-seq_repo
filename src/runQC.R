plotQC <- function(d, basename, ext='pdf') {
  fig1 <- plot_total_counts(d)
  fig2 <- plot_library_complexity(d)
  fig3 <- plot_gene_detection(d)
 
  vsd <- vst(d)
  design <- as.character(design(dds))[-1]
  anno_vars = str_split(design, pattern='_')[[1]]
  set.seed(3)

  pdf((glue('{basename}_heatmap.pdf')), width = 8, height = 8)
  fig4 <- plot_sample_clustering(vsd, anno_vars = anno_vars, distance = "euclidean")
  draw(fig4)
  dev.off()
  
  plots <- list(fig1, fig2, fig3)

  p <- cowplot::plot_grid(plotlist = plots, ncol = length(plots))
  cowplot::save_plot(file=glue('{basename}.{ext}'), plot=p, base_asp = length(plots))
}


source('./src/utilities.R')
prepEnv()

parser <- ArgumentParser()
parser$add_argument('-i', '--input_file', type='character', nargs=1, help='DESeq2 dataset object')

args   <- parser$parse_args()
file   <- args$i

dds <- readRDS(file)

fig_stem <- './fig/QC/'
fig_name <- basename(file_path_sans_ext(file))

plotQC(dds, glue('{fig_stem}{fig_name}_pre'))



sample_threshold = ncol(dds) / 2
count_threshold = 3
# take counts data -> cpm
cpm <- assay(dds) / (colSums(assay(dds)) / 1e6)
# require gene counts to be > 5 
cpm_data_threshold = rowSums(cpm > count_threshold, na.rm=T) >= sample_threshold

dds_filt <- dds[cpm_data_threshold, ]


plotQC(dds_filt, glue('{fig_stem}{fig_name}_post'))



