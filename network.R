library(dplyr)
library(ggplot2)
library(DescTools)
library(rstatix)
library(tibble)
library(broom)
library(ggpubr)
# data file
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/IN2 New NK cell gating(2 subsets) -December 20th 2020.csv")

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
zero<-df[, colSums(df[,-1]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df<-df %>% select(Group, names(nonzero))

# major immune markers

# for innate panel
#df_parent<-df%>%select(Eosinophils,Neutrophil,Monocytes,CM,ITM,NCM,CD56BRTCD16BRT,CD56BRTCD16N,CD56DIMCD16N,CD56DIMCD16P,CD56NCD16P,CD3P,CD14N,CD4NCD8N,CD4P,CD4PCD8P,CD8P,NKT,CD3PCD19P,CD19)
# for adaptive panel
df_parent<-df%>%select(CD3P,CD4P,CD4PCD19P,CD19P,CD19hi,CD19lo,CD38PCD27P,PB,PC,IgDNCD27N,DN1,DN2,Naive,SM,USM)

# input row index for each group
end_row_HD = 15 
start_row_ICU = 16
end_row_ICU =24
start_row_RD = 25
end_row_RD =32

df_parent_HD<-df_parent[(1:end_row_HD),]
df_parent_ICU<-df_parent[(start_row_ICU:end_row_ICU),]
df_parent_RD<-df_parent[(start_row_RD:end_row_RD),]

#network-------------------------------------------------
#---- change HD to ICU and RD to plot the network for ICU and RD ---------
library("qgraph")
network <- qgraph(cor(df_parent_HD,method="pearson")
                  ,label.cex=0.3
                  ,labels=colnames(df_parent_HD)
                  ,label.scale=F
                  ,layout = "circle"
                  #,details=T
                  ,edge.labels=F
                  ,alpha=0.05
                  ,minimum='sig'
                  ,sampleSize=nrow(df_parent_HD)
                  ,posCol = c("#0000FF","blue"),negCol = c("#BF0000","red")
                  ,negDashed =F
                  ,edge.label.cex=0.4
                  ,edge.label.color ="black"
                  ,edge.label.position=0.2
                  ,edge.label.bg=T
                  #,vTrans=100
                  ,filetype="png"
                  #,filetype="correlation_network_parents_r"
                  ,width=6
                  ,height=6
                  ,graph = 'cor')

# subgroup network _ not completed
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/Innate_study2/COVID -19 Paper 1 INT final data 9.14.2020.csv")
zero<-df[, colSums(df[,-1]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df<-df %>% select(Group, names(nonzero))

colnames(df) <- gsub("CD4PCD8P","4P8P",colnames(df))
colnames(df) <- gsub("CD3PCD19P","3P19P",colnames(df))
colnames(df) <- gsub(" ","",colnames(df))
colnames(df) <- gsub("-","",colnames(df))
colnames(df) <- gsub("(","",colnames(df),fixed=T)
colnames(df) <- gsub(")","",colnames(df),fixed=T)


grouplist<-list()

list1<-list('Eosinophils','Neutrophil','Monocytes','CM','ITM')
#'NCM','CD3P','CD56BRTCD16BRT','CD56BRTCD16N','CD56DIMCD16N','CD56DIMCD16P','CD56NCD16P','CD4NCD8N','CD4P','4P8P','CD8P','NKT','3P19P','CD19')
for (i in list1) {
  print(i)
  name <- paste(i,sep='')
  df_sub<- df %>% dplyr:: select(starts_with(i))
  
  print(colnames(df_sub))
  grouplist[[name]]<-as.factor(colnames(df_sub))
}


df<-df[,-1]
df_ICU<-df[c(12:31),c(1:148)]
zero<-df_ICU[, colSums(df_ICU) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df_ICU<-df_ICU %>% select(names(nonzero))

library("qgraph")
Q <- qgraph(cor(df_ICU,method="pearson"), minimum = 0.2,groups = grouplist, 
            legend = TRUE, borders = FALSE)

dev.off()


network <- qgraph(cor(df_ICU,method="pearson")
                  ,label.cex=0.5
                  #,labels=colnames(df_ICU)
                  ,label.scale=F
                  ,details=T
                  ,edge.labels=T
                  ,alpha=0.05
                  #,minimum='sig'
                  #,sampleSize=nrow(df_parent_ICU)
                  ,posCol = "blue", negCol = "red"
                  ,negDashed =F
                  ,edge.label.cex=0.5
                  ,edge.label.color ="black"
                  ,edge.label.position=0.2
                  ,edge.label.bg=F
                  ,groups = grouplist)
# for cytoscape
df<-df[,-1]
df_ICU<-df[c(12:31),c(1:103)]
zero<-df_ICU[, colSums(df_ICU) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df_ICU<-df_ICU %>% select(names(nonzero))

Corr=Hmisc::rcorr(as.matrix(df_ICU),type='pearson') # apply pearson correlation 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
C=flattenCorrMatrix(Corr$r, Corr$P) # extract correlation matrix and p-value matrix
C1=C[order(-abs(C$cor)),] # order data according to corretaion value (Descending)
C2=na.omit(C) # filter out NA values

C2$row <- as.character(C2$row)
C2$column <- as.character(C2$column)
library(dplyr)

list2<-list("Eosinophils","Neutrophil",'Monocytes')

data<-C2 %>%
  mutate(staterow = case_when(
      startsWith(row, "Eosinophils") ~ "Eosinophils",
      startsWith(row, "Neutrophil") ~ "Neutrophil",
      startsWith(row, "Monocytes") ~ "Monocytes",
      TRUE                      ~ NA_character_
    ))

data<-data %>%
  mutate(statecolumn = case_when(
    startsWith(column, "Eosinophils") ~ "Eosinophils",
    startsWith(column, "Neutrophil") ~ "Neutrophil",
    startsWith(column, "Monocytes") ~ "Monocytes",
    TRUE                      ~ NA_character_
  ))
write.csv(data,"for_cytospace.csv")

markers<-as.data.frame(colnames(df_ICU))
markers<-markers %>%
  mutate(group = case_when(
    startsWith(colnames(df_ICU), "Eosinophils") ~ "Eosinophils",
    startsWith(colnames(df_ICU), "Neutrophil") ~ "Neutrophil",
    startsWith(colnames(df_ICU), "Monocytes") ~ "Monocytes",
    TRUE                      ~ NA_character_
  ))
write.csv(markers,"for_cytospace_groups.csv")



  