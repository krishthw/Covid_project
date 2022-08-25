library(pheatmap)
#Heat maps for HD_ICU_RD
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df2<-aggregate(df[,7:length(colnames(df))], list(df$`Group`), mean) # mean of groups
df_mean <- as.data.frame(df2[,2:length(colnames(df2))]) # get markes
rownames(df_mean) <- df2$`Group.1` # for HD-RD-ICU data only

# rescaleed means

# 1.normalizing data (min-max scale) --> x(scaled)= x- min(x)) /(max(x)-min(x))
library(scales)
#_____________________________________________________________
rescale.many <- function(dat, column.nos) { 
  nms <- names(dat) 
  for(col in column.nos) { 
    name <- paste(nms[col],".rescaled", sep = "") 
    dat[name] <- rescale(dat[,col]) 
  } 
  cat(paste("Rescaled ", length(column.nos), " variable(s)n")) 
  dat 
} 
#_____________________________________________________________

df_mean_norm<-rescale.many(df_mean,c(1:length(colnames(df_mean)))) # minmaxsclae for whole dataset
df_mean_norm<-df_mean_norm[,(length(colnames(df_mean))+1):length(colnames(df_mean_norm))] # get normalized columns
colnames(df_mean_norm) <- gsub(".rescaled","",colnames(df_mean_norm)) # remove unneccessary stuff from column names

nheatplot<-pheatmap(t(df_mean_norm),
                    cellwidth = 10,
                    cellheight = 10,
                    cluster_cols  = FALSE,
                    cluster_rows = TRUE,
                    legend_breaks=(0:1),
                    legend_labels =c("row min","row max"),
                    color=colorRampPalette(c("#0000CD","white","#FF8C00"))(200))
cowplot::ggsave2(paste0("CovidFinalfeb/Plasma/heatmaps/heatmap_minmax_plasma" , ".png"), nheatplot, width = 6, height = 9, units = 'in' ,type = "cairo", dpi = 300)



#Heat maps for HD_ICU_RD
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Cytokines Anna's data Updated 12.25.2020_Krish.csv")
df<-as.data.frame(df) # save data table as dataframe

df<-df[df[, "Group"] == 'ICU',]
df2<-aggregate(df[,7:length(colnames(df))], list(df$`Diabetes`), mean) # mean of groups
df_mean <- as.data.frame(df2[,2:length(colnames(df2))]) # get markes
rownames(df_mean) <- list('Non Diabetes','Diabetes') # for Diabetes only

df_mean_norm<-rescale.many(df_mean,c(1:length(colnames(df_mean)))) # minmaxsclae for whole dataset
df_mean_norm<-df_mean_norm[,(length(colnames(df_mean))+1):length(colnames(df_mean_norm))] # get normalized columns
colnames(df_mean_norm) <- gsub(".rescaled","",colnames(df_mean_norm)) # remove unneccessary stuff from column names

nheatplot<-pheatmap(t(df_mean_norm),
                    cellwidth = 10,
                    cellheight = 10,
                    cluster_cols  = FALSE,
                    cluster_rows = TRUE,
                    legend_breaks=(0:1),
                    legend_labels =c("row min","row max"),
                    color=colorRampPalette(c("#0000CD","white","#FF8C00"))(200))
cowplot::ggsave2(paste0("heatmap_minmax_anna" , ".png"), nheatplot, width = 6, height = 9, units = 'in' ,type = "cairo", dpi = 300)
















# in case of normalize before get the mean

df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Cytokines Anna's data Updated 12.25.2020_Krish.csv")
df<-as.data.frame(df) # save data table as dataframe


# 1.normalizing data (min-max scale) --> x(scaled)= x- min(x)) /(max(x)-min(x))
library(scales)
#_____________________________________________________________
rescale.many <- function(dat, column.nos) { 
  nms <- names(dat) 
  for(col in column.nos) { 
    name <- paste(nms[col],".rescaled", sep = "") 
    dat[name] <- rescale(dat[,col]) 
  } 
  cat(paste("Rescaled ", length(column.nos), " variable(s)n")) 
  dat 
} 
#_____________________________________________________________

df_norm<-rescale.many(df,c(7:length(colnames(df)))) # minmaxsclae for whole dataset
df_norm<-df_norm[,c(2,(length(colnames(df))+1):length(colnames(df_norm)))] # get normalized columns
colnames(df_norm) <- gsub(".rescaled","",colnames(df_norm)) # remove unneccessary stuff from column names

df_norm_mean<-aggregate(df_norm[,2:length(colnames(df_norm))], list(df_norm$`Group`), mean) # anna

df2 <- as.data.frame(df_norm_mean[,2:length(colnames(df_norm_mean))])

rownames(df2) <- df_norm_mean$`Group.1` # for HD-RD-ICU data only
#rownames(df_mean) <- list('Non Diabetes','Diabetes') # for Diabetes only

nheatplot<-pheatmap(t(df2),
                    cellwidth = 10,
                    cellheight = 10,
                    cluster_cols  = FALSE,
                    cluster_rows = T,
                    legend_breaks=(0:1),
                    legend_labels =c("row min","row max"),
                    color=colorRampPalette(c("#0000CD","white","#FF8C00"))(200))
