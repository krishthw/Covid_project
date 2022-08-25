df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Day 1 IN2Color Coded- December 3rd  2020heatmap .csv")
# ----------------------------------------------------------------------------
# 1.Data Preperation------------------------- for three sets -----------------
# ----------------------------------------------------------------------------
# remove all zero columns
zero<-df[, colSums(df[,-1]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df<-df %>% select(Sample, names(nonzero))

#Look for duplicated columns
## function to test unique columns =========================================
test.unique <- function(df) {  
  length1 <- length(colnames(df))
  length2 <- length(unique(colnames(df)))        
  if (length1 - length2 > 0 ) {
    
    print(paste("There are", length1 - length2, " duplicates", sep=" "))
  }     
}
#===========================================================================
test.unique(df) 

duplicated(colnames(df)) # if duplicated names present output "TRUE"

# it is good to know mean and sd for each group
dfmeans<-df%>%
  group_by(Sample) %>%
  summarise_all(funs(sd(., na.rm=TRUE)))

df_sub<- df %>% dplyr:: select(ends_with('CXCR2P'))

df_num <- as.matrix(df_sub)
rownames(df_num) <- df$`Sample`
df_num_scale = scale(df_num) # scale the data to a distribution with mean as 0 and standard deviation as 1. (x-mean)/sd


annotation <- data.frame(Group = rep(c("HD", "ICU","RD"), c(16,20,31)))
row.names(annotation) <- rownames(df_num_scale)
Group       <- c("forestgreen","red2","cornflowerblue")
names(Group ) <- c("HD", "ICU","RD")
anno_colors <- list(Group = Group)

#col.pal <- RColorBrewer::brewer.pal(6, "YlGnBu")

heatplot<-pheatmap(t(df_num_scale), 
                   #color=col.pal,
                   
                   annotation=annotation, 
                   annotation_colors = anno_colors,
                   cluster_cols  = FALSE,
                   border_color= NA,
                   annotation_legend = FALSE,
                   fontsize_col = 7,
                   fontsize_row = 8,
                   
                   
                   #show_colnames=FALSE
)

#ggsave(heatmap,  filename = "heatmap.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")
cowplot::ggsave2(paste0("heatmap_CXCR2P_raw" , ".png"), heatplot, width = 9, height = 6, units = 'in' ,type = "cairo", dpi = 300)

CXCR2P<-data.frame(rowSums(df_num))
df1 <- cbind(sample = rownames(CXCR2P), CXCR2P)
df2 <- cbind(sample = rownames(annotation),annotation)
df3<-inner_join(df2,df1)
colnames(df3)[3] <- "total_CXCR2P"

write.csv(df3,'total_CXCR2P.csv')
df3<-df3[,-1]
dfsig<-df3


