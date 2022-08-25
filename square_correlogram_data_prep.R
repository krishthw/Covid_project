# correlogram 

library(dplyr)
library(corrplot)
# import data files
df1<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Cytokines Anna's data final- February 10th 2021.csv")
df2<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Isoplexis Codeplex final- February 10th 2021.csv")
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
df4<-df4 %>% select(c("Sample","Group","Age","Sex","Diabetes","State",'Eosinophils','Neutrophil','Monocytes','NCM','CM','ITM','CD3P','CD56BRTCD16N','CD56DIMCD16P','CD4NCD8N','CD4PCD8P','CD4P','CD8P','NKT','CD19','CD14N'))
# remove - from codeplex and anna data to identify if there is any same markers
colnames(df1) <- gsub("-","",colnames(df1)) 
colnames(df2) <- gsub("-","",colnames(df2))
colnames(df3) <- gsub("(pg/mL)","",colnames(df3))
colnames(df3) <- gsub("(","",colnames(df3),fixed=T)
colnames(df3) <- gsub(")","",colnames(df3),fixed=T)

(anti_join(df1, df2, by="Sample"))$Sample # 3 missing
(anti_join(df2, df1, by="Sample"))$Sample # 7 missing

# join Cytokine datasets
df<-inner_join(df1,df2, by="Sample") # innerjoin by 'Sample'
df1df2<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df2) <- gsub("\\.x$","",colnames(df1df2)) # remove .x from column names
# join Plasma dataset
dff<-inner_join(df1df2,df3, by="Sample")
df1df2df3<- dff %>% select(-ends_with(".y"))
colnames(df1df2df3) <- gsub("\\.x$","",colnames(df1df2df3))
# join innate dataset
dfff<-inner_join(df1df2df3,df4, by="Sample")
df1df2df3df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df3df4) <- gsub("\\.x$","",colnames(df1df2df3df4))

# cytokine innate combine
df<-inner_join(df1,df2, by="Sample") # innerjoin by 'Sample'
df1df2<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df2) <- gsub(".x","",colnames(df1df2)) # remove .x from column names
# join innate dataset
dfff<-inner_join(df1df2,df4, by="Sample")
df1df2df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df4) <- gsub(".x","",colnames(df1df2df4))
df1df2df4<- df1df2df4%>% select(-c(S)) # this contains cytokines and innate parents

data<-as.data.frame(df1df2df4) # convert to dataframe
CI_ICU<-data[data[, "Group"] == 'ICU',-c(1:5,41)]
CI_RD<-data[data[, "Group"] == 'RD',-c(1:5,41)]
CI_HD<-data[data[, "Group"] == 'HD',-c(1:5,41)]

subdata<-CI_RD # change this to CI_RD,CI_HD to get their correlograms

Corr=Hmisc::rcorr(as.matrix(subdata),type='pearson') # apply pearson correlation 
png(height=8, width=10, file="CovidFinalfeb/Mix_correlograms/Itest.png", type = "cairo",units="in",res=1800)
corrplot(Corr$r,tl.cex = 0.5,tl.col ="black",col=brewer.pal(n=10, name="PuOr"),
         method='color',
         diag=T)
corrRect(c(46,16),lwd=1)
corrRect(62,lwd=1)
dev.off()

library("qgraph")
network <- qgraph(cor(subdata,method="pearson")
                  ,label.cex=0.3
                  ,labels=colnames(subdata)
                  ,label.scale=F
                  ,layout = "circle"
                  , vsize = 4
                  #,details=T
                  ,edge.labels=F
                  ,alpha=0.05
                  ,minimum='sig'
                  ,sampleSize=nrow(subdata)
                  ,posCol = c("#4e258f"),negCol = c("#eb9431")
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



# plasma innate combine
dfff<-inner_join(df3,df4, by="Sample") # innerjoin by 'Sample'
df3df4<- dfff %>% select(-ends_with(".y")) # remove duplicates 
colnames(df3df4) <- gsub(".x","",colnames(df3df4)) # remove .x from column names

data<-as.data.frame(df3df4) # convert to dataframe
CI_ICU<-data[data[, "Group"] == 'ICU',-c(1:6)]
CI_RD<-data[data[, "Group"] == 'RD',-c(1:6)]
CI_HD<-data[data[, "Group"] == 'HD',-c(1:6)]

subdata<-CI_HD # change this to CI_RD,CI_HD to get their correlograms

Corr=Hmisc::rcorr(as.matrix(subdata),type='pearson') # apply pearson correlation 
png(height=8, width=10, file="CovidFinalfeb/Mix_correlograms/HD15_Innate_plasma.png", type = "cairo",units="in",res=1800)
corrplot(Corr$r,tl.cex = 0.5,tl.col ="black",col=brewer.pal(n=10, name="PuOr"),
         method='color',
         diag=T)
corrRect(c(8,16),lwd=1)
corrRect(24,lwd=1)
dev.off()


