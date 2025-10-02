library(ggplot2)
library(ggridges)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(readxl)

### Re-work data frames to join Adaptation, proliferation and MKI67 data
Adaptive_signature_df <- readRDS("./data/source_data/sci-plex/adaptive_resistance_upgenes_zscored.RDS")

Adaptive_signature_df_BT112 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT112",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT112) <- sapply(colnames(Adaptive_signature_df_BT112), function(x){stringr::str_split(x,pattern = "BT112_")[[1]][2]})
row.names(Adaptive_signature_df_BT112) <- paste0(row.names(Adaptive_signature_df_BT112),"_BT112")

Adaptive_signature_df_BT228 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT228",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT228) <- sapply(colnames(Adaptive_signature_df_BT228), function(x){stringr::str_split(x,pattern = "BT228_")[[1]][2]})
row.names(Adaptive_signature_df_BT228) <- paste0(row.names(Adaptive_signature_df_BT228),"_BT228")

Adaptive_signature_df_BT333 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT333",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT333) <- sapply(colnames(Adaptive_signature_df_BT333), function(x){stringr::str_split(x,pattern = "BT333_")[[1]][2]})
row.names(Adaptive_signature_df_BT333) <- paste0(row.names(Adaptive_signature_df_BT333),"_BT333")

Adaptive_signature_df_joint <- rbind(Adaptive_signature_df_BT112,Adaptive_signature_df_BT228,Adaptive_signature_df_BT333)

hmcols <- colorRampPalette(c("blue","white","red"))(35)

paletteLength <- 35

myBreaks <- c(seq(min(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)]), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)])/paletteLength, max(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)]), length.out=floor(paletteLength/2)))

drug_annotations <- read.csv("./data/source_data/sci-plex/Annotations_final_RG.csv")
drug_annotations <- replace(drug_annotations, drug_annotations=='', "Unknown")
drug_annotations$Reversibility <- sapply(drug_annotations$reversible, function(x){
  if(x == "Yes")return("Reversible")
  if(x == "No")return("Covalent")
  return("Unknown")
})


ctrl_df <- read_excel("./data/lux_data/fig1/clean-data.xlsx")
ctrl_df <- ctrl_df %>% filter(day==5) %>% group_by(condition, donor) %>% summarize(mean_mfi = mean(`norm mfi`))
ctrl_df <- ctrl_df %>% group_by(condition) %>% summarize(median_mfi = median(mean_mfi))

ctrl_df$log2_ctrl <- round(log2(ctrl_df$median_mfi) * 2) / 2

drug_annotations$log2_ctrl <- 0

drug_annotations[tolower(drug_annotations$drug) %in% ctrl_df$condition, 'log2_ctrl']  
print(head(drug_annotations))


ann_col <- list(Reversibility = c("Covalent" = "#FF0022","Reversible" = "#011627","Unknown" = "white"),
                Class = c("acetamide" = RColorBrewer::brewer.pal(9,"Set1")[1], 
                          "aminobenzimidazole" = RColorBrewer::brewer.pal(9,"Set1")[2],
                          "aminoethylamide" = RColorBrewer::brewer.pal(9,"Set1")[3], 
                          "aminonucleoside" = RColorBrewer::brewer.pal(9,"Set1")[4],      
                          "aminopyrazine"  = RColorBrewer::brewer.pal(9,"Set1")[5], 
                          "benzoimidazole" = RColorBrewer::brewer.pal(9,"Set1")[6], 
                          "dicarboxylic acid" = RColorBrewer::brewer.pal(9,"Set1")[7], 
                          "EGFR-activator"  = RColorBrewer::brewer.pal(9,"Set1")[8], 
                          "indole" = RColorBrewer::brewer.pal(9,"Set1")[9],
                          "monoclonal antibody" = RColorBrewer::brewer.pal(8,"Set2")[1], 
                          "naturally-derived" = RColorBrewer::brewer.pal(8,"Set2")[2],
                          "PROTAC" = RColorBrewer::brewer.pal(8,"Set2")[3],             
                          "pyridine" = RColorBrewer::brewer.pal(8,"Set2")[4],           
                          "pyrimidine" = RColorBrewer::brewer.pal(8,"Set2")[5],          
                          "quinazolinamine" = RColorBrewer::brewer.pal(8,"Set2")[6], 
                          "quinazoline" = RColorBrewer::brewer.pal(8,"Set2")[7],
                          "quinoline" = RColorBrewer::brewer.pal(8,"Set2")[8], 
                          "sulfanomide" = RColorBrewer::brewer.pal(3,"Paired")[1],
                          "tyrphostin" = RColorBrewer::brewer.pal(3,"Paired")[2]))
set.seed(3)

pheatmap::pheatmap(Adaptive_signature_df_joint,
                   clustering_method = "ward.D2",
                   cluster_rows = FALSE,
                   gaps_row = c(5,10),
                   treeheight_col = 10,
                   color = hmcols,
                   breaks = myBreaks,
                   border_color = NA,
                   fontsize = 12,
                   cutree_cols = 5,
                   annotation_col = data.frame(row.names = drug_annotations$drug,
                                               Reversibility = drug_annotations$Reversibility,
                                               Class = drug_annotations$class_broad),
                   annotation_colors = ann_col,
                   file = "fig/chem_info/Adaptive_signature_heatmap.pdf",
                   width = 14,
                   height = 6)
