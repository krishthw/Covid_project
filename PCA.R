library(cowplot)
library(grid)

library(dplyr)
library(ggplot2)
library(DescTools)
library(rstatix)
library(tibble)
library(broom)
library(ggpubr)

# data file
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/IN2-Cytokins-Plasma-Combined- December 23rd  2020(Cleaned).csv")

# __________________________________________________________________________________________
# 1.Data Preperation 
# __________________________________________________________________________________________
#Look for duplicated columns
duplicated(colnames(df)) # if duplicated names present output "TRUE"

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

# remove unnecessary characters
colnames(df) <- gsub(" ","",colnames(df))
colnames(df) <- gsub("-","",colnames(df))
colnames(df) <- gsub("(","",colnames(df),fixed=T)
colnames(df) <- gsub(")","",colnames(df),fixed=T)

# remove all zero columns
zero<-df[, colSums(df[,-c(1:5)]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df<-df %>% select(Group, names(nonzero))

# get the columns for PCA
# innate panel
df %>% select(Eosinophils,Neutrophil,Monocytes,CM,ITM,NCM,CD56BRTCD16BRT,CD56BRTCD16N,CD56DIMCD16N,CD56DIMCD16P,CD56NCD16P,CD4NCD8N,CD4P,CD4PCD8P,CD8P,NKT,CD3PCD19P,CD19) %>% 
# adaptive panel
#df %>% select(CD4P,CD4PCD19P,CD19P,CD19hi,CD19lo,CD38PCD27P,PB,PC,IgDNCD27N,DN1,DN2,Naive,SM,USM) %>% 
  scale() %>%                  # scale to 0 mean and unit variance
  prcomp() ->pca               # store result as `pca`

# display the results from the PCA 
head(pca$x)
pca_data <- data.frame(pca$x, Group=df$Group)
head(pca_data)

percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100,2)
percentage <- paste( colnames(pca_data), "(", paste( as.character(percentage), "%", ")", sep="") )

theme<-theme(
  #panel.background = element_rect(fill="white"),
  #legend.position= "none",
  axis.line = element_line(size=1,colour="black"),
  axis.text.x = element_text(size=10,colour="black",face="bold"),
  axis.text.y = element_text(size=10,colour="black",face="bold"))

PCA<-ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
  scale_fill_manual(values = c("forestgreen","red2","cornflowerblue"))+
  geom_point(aes(colour=Group,fill=Group),shape=21,size = 2,color="black") +
  xlab(percentage[1]) + ylab(percentage[2])


PCA<-PCA+theme
PCA
#ggsave(PCA,  filename = "CovidWB/PCA/PCA_final_AD_parent.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")

# get the row corresponding to that specific point
point1 <- pca_data[16,] # deceased 1
point2 <- pca_data[17,] # deceased 2
point3 <- pca_data[23,] # deceased 3

PCAwithDead<-ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
  scale_fill_manual(values = c("forestgreen","red2","cornflowerblue"))+
  geom_point(aes(colour=Group,fill=Group),shape=21,size = 2,color="black") +
  xlab(percentage[1]) + ylab(percentage[2])+ 
  annotate("text", label='LL_026', x=point1$PC1+0.5, y=point1$PC2,size = 2,color="black")+
  annotate("text", label='LL_034', x=point2$PC1+0.5, y=point2$PC2,size = 2,color="black")+
  annotate("text", label='LL_097', x=point3$PC1+0.5, y=point3$PC2,size = 2,color="black")

PCAwithDead<-PCAwithDead+theme
PCAwithDead
#ggsave(PCAwithDead,  filename = "CovidWB/PCA/PCA_final_AD_parent_w_deceased.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")

#PCA for each subset-----------------------------------------------------------------------

# for innate panel
# In order to filter CD4P without CD4PCD8P change all the CD4PCD8P to 4P8P
colnames(df) <- gsub("CD4PCD8P","4P8P",colnames(df))
colnames(df) <- gsub("CD3PCD19P","3P19P",colnames(df))
list1<-list('Eosinophils','Neutrophil','Monocytes','CM','ITM','NCM','CD56BRTCD16BRT','CD56BRTCD16N','CD56DIMCD16N','CD56DIMCD16P','CD56NCD16P','CD4NCD8N','CD4P','4P8P','CD8P','NKT','3P19P','CD19')

# for adaptive panel
#colnames(df) <- gsub("CD4PCD19P","4P19P",colnames(df))
#list1<-list('CD4P','4P19P','CD19P','CD19hi','CD19lo','PB','PC','DN1','DN2','Naive','SM','USM')

pdf(file="CovidDec4/PCA/PCA_mod_subsets_extended.pdf")

for (i in list1) {
  print(i)
  
  #PCA for each parent
  df_sub<- df %>% dplyr:: select(starts_with(i))
  
  df_sub %>% 
    scale() %>%               # scale to 0 mean and unit variance
    prcomp() ->pca_sub        # store result as `pca`
  
  # display the results from the PCA analysis
  head(pca_sub$x)
  
  pca_data_sub <- data.frame(pca_sub$x, Group=df$Group)
  head(pca_data_sub)
  pca_sub
  
  percentage <- round(pca_sub$sdev^2 / sum(pca_sub$sdev^2) * 100,2)
  percentage <- paste( colnames(pca_data_sub), "(", paste( as.character(percentage), "%", ")", sep="") )
  theme<-theme(
    axis.line = element_line(size=1,colour="black"),
    axis.text.x = element_text(size=10,colour="black",face="bold"),
    axis.text.y = element_text(size=10,colour="black",face="bold"),
    legend.position="none"
    )
  
  p<-ggplot(pca_data_sub,aes(x=PC1,y=PC2,color=Group ))
  p<-p+theme+ xlab(percentage[1]) + ylab(percentage[2])+scale_fill_manual(values = c("forestgreen","red2","cornflowerblue"))+
    geom_point(aes(colour=Group,fill=Group),shape=21,size = 3,color="black") +
    ggtitle(paste0("group ",i))
  
  print(p)
}
dev.off()

