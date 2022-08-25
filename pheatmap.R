df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Liyanage Isoplexis Codeplex121520_raw_data (1).csv")
library(pheatmap)

df_num <- as.matrix(df[,2:24])
rownames(df_num) <- df$`Donor`
df_num_scale = scale(df_num) # scale the data to a distribution with mean as 0 and standard deviation as 1. (x-mean)/sd

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
  }

annotation <- data.frame(sample = rep(c("HD", "ICU","RD"), c(17,20,27)))
row.names(annotation) <- rownames(df_num_scale)
sample       <- c("forestgreen","red2","cornflowerblue")
names(sample ) <- c("HD", "ICU","RD")
anno_colors <- list(sample = sample)

#col.pal <- RColorBrewer::brewer.pal(6, "YlGnBu")

heatplot<-pheatmap(t(df_num_scale), 
                   #color=col.pal,
                   
                   annotation=annotation, 
                   annotation_colors = anno_colors,
                   cluster_cols  = FALSE,
                   border_color= NA,
                   annotation_legend = FALSE,
                   fontsize_col = 7,
                   fontsize_row = 8
         
         #show_colnames=FALSE
         )

#ggsave(heatmap,  filename = "heatmap.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")
cowplot::ggsave2(paste0("heatmap" , ".png"), heatplot, width = 9, height = 6, units = 'in' ,type = "cairo", dpi = 300)

heatplot<-pheatmap(t(log2(df_num_scale)), 
                   
                   annotation=annotation, 
                   annotation_colors = anno_colors,
                   cluster_cols  = FALSE,
                   border_color= NA,
                   annotation_legend = FALSE
                   #show_colnames=FALSE
                   )

