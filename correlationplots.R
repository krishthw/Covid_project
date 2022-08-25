library(dplyr)
library(ggplot2)
library(DescTools)
library(rstatix)
library(tibble)
library(broom)
library(ggpubr)
library(corrplot)

# data file
df_data<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Cytokines Anna's data final- February 10th 2021.csv")
df_data<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Isoplexis Codeplex final- February 10th 2021.csv")
df_data<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
# select major immune markers for correlation analysis
df_data<-df4 %>% select(c("Sample","Group","Age","Sex","Diabetes","State",'Eosinophils','Neutrophil','Monocytes','NCM','CM','ITM','CD3P','CD56BRTCD16N','CD56DIMCD16P','CD4NCD8N','CD4PCD8P','CD4P','CD8P','NKT','CD19','CD14N'))

# __________________________________________________________________________________________
# Data Preperation 
dff<-df_data%>% select(c(2,7:length(colnames(df_data))))

#Look for duplicated columns
length(colnames(dff))-length(unique(colnames(dff))) # number of duplicates
duplicated(colnames(dff)) # if duplicated names present output "TRUE"

# remove unnecessary characters
colnames(dff) <- gsub("(ng/mL)","",colnames(dff))
colnames(dff) <- gsub("(pg/mL)","",colnames(dff))
colnames(dff) <- gsub("(u/mL)","",colnames(dff))
colnames(dff) <- gsub(" ","",colnames(dff))
colnames(dff) <- gsub("-","",colnames(dff))
colnames(dff) <- gsub("(","",colnames(dff),fixed=T)
colnames(dff) <- gsub(")","",colnames(dff),fixed=T)
colnames(dff) <- gsub("_","",colnames(dff))

# remove all zero columns
no_zero<-dff[, colSums(dff[,-1]) > 0] # get logical vector to show 
nonzero<-which(no_zero == TRUE)     # select non zero columns

names(which(colSums(is.na(dff))>0)) # names 
df<-dff %>% select(Group, names(nonzero)) # Group and non-zero columns

#major immune markers
data<-as.data.frame(df) # convert to dataframe
ICU<-data[data[, "Group"] == 'ICU',-1] 
RD<-data[data[, "Group"] == 'RD',-1]
HD<-data[data[, "Group"] == 'HD',-1]

# for adaptive panel
#df_parent<-df%>%select(CD3P,CD4P,CD19P,CD38PCD27P,Plasmablasts,Plasmacells,IgDNCD27N,DN1,DN2,Naive,SwitchedMemory,UnswitchedMemory)
CorrM=Hmisc::rcorr(as.matrix(RD),type='pearson') # apply pearson correlation 

png(height=10, width=10, file="CovidFinalfeb/Mix/all_correlations/RD/all_combined_RD.png", type = "cairo",units="in",res=1800)
correlogram<-corrplot(CorrM$r,type="upper",tl.cex = 0.7,tl.col ="black",
                           p.mat = CorrM$P, sig.level = 0.05, insig = "blank",
                           method='square', 
                           number.cex = 0.7,diag=F,
                          #addCoef.col = "black",
                           title= "RD", 
                           cex.main = 1, font.main= 2,
                           mar=c(0,0,3,0))
dev.off()
# function to extract corr value and p-value of the upper triangle ==================
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
#======================================================================================

C=flattenCorrMatrix(CorrM$r, CorrM$P) # extract correlation matrix and p-value matrix
C=C[order(-abs(C$cor)),] # order data according to corretaion value (Descending)
C=na.omit(C) # filter out NA values
#write.csv(C,'CovidFinalfeb/Innate/Correlations/innateparent_ICU_correlations.csv')

C_sig=C[abs(C$p)<0.05,]
write.csv(C_sig,'CovidFinalfeb/Mix/all_correlations/RD/all_combined_sig_correlations_RD.csv')

df1<-RD
df2<-C_sig
# ========================================================================================
for (i in seq(1:nrow(df2))){
  colname1 <-as.character(df2[i,1])
  colname2<-as.character(df2[i,2])
  #df3<-df1 %>% dplyr::select (all_of(colname1), all_of(colname2))
  cplot<-ggpubr::ggscatter(df1, x = colname2, y = colname1,
                           add = "reg.line",conf.int = TRUE,
                           xlab=colname2, ylab=colname1,
                           color="black",
                           title=paste0(colname1 ," vs " ,colname2))+
    ggpubr::stat_cor(method = "pearson", label.x.npc=0.25, label.y.npc= 1.0,r.accuracy=0.01,p.accuracy=0.001) +
    #scale_y_continuous(labels=function(x) paste0(x*100, "%"))+
    #scale_x_continuous(labels = scales::number_format(accuracy = 0.0001))+
    theme(axis.text.x=element_text(size=10),
          axis.title.y=element_text(size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=10))
  cplot<-cplot+geom_point( shape = 21,fill='cornflowerblue',color="black",lwd=3) # cornflowerblue for RD red2 for ICU forestgreen for HD
  print(cplot)
  ggsave(cplot, filename = paste0("CovidFinalfeb/Mix/all_correlations/RD/", colname1, "vs", colname2, ".png"),dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")
}
dev.off()

#---------- correaltions between different files--------------
lista<-colnames(data[,2:9])
listb<-colnames(data[,10:25])

for (i in seq(1:nrow(df2))){
  colname1 <-as.character(df2[i,1])
  colname2<-as.character(df2[i,2])
  #df3<-df1 %>% dplyr::select (all_of(colname1), all_of(colname2))
  if (((colname1 %in% lista) & (colname2 %in% listb)) || ((colname2 %in% lista) & (colname1 %in% listb))){
    print(colname1)
    print(colname2)
    
    cplot<-ggpubr::ggscatter(df1, x = colname2, y = colname1,
                             add = "reg.line",conf.int = TRUE,
                             xlab=colname2, ylab=colname1,
                             color="black",
                             title=paste0(colname1 ," vs " ,colname2))+
      ggpubr::stat_cor(method = "pearson", label.x.npc=0.25, label.y.npc= 1.0,r.accuracy=0.01,p.accuracy=0.001) +
      #scale_y_continuous(labels=function(x) paste0(x*100, "%"))+
      #scale_x_continuous(labels = scales::number_format(accuracy = 0.0001))+
      theme(axis.text.x=element_text(size=10),
            axis.title.y=element_text(size=12),
            axis.title.x=element_text(size=12),
            axis.text.y=element_text(size=10))
    cplot<-cplot+geom_point( shape = 21,fill='forestgreen',color="black",lwd=3) # cornflowerblue for RD red2 for ICU forestgreen for HD
    print(cplot)
    ggsave(cplot, filename = paste0("CovidFinalfeb/Mix/Correlations/PlasmaInnate/HD/", colname1, "vs", colname2, ".png"),dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")
  }
}
dev.off()

# multiple correlations

databubble<-ICU
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
# when databubble=ICU plotting for ICU
bubsub<-as.data.frame(databubble%>%select(CRP,sCD14,ITM,CD4P,CD8P,Eosinophils,Neutrophil,ITM,CD56DIMCD16P,CD19PCXCR3P,CD19PCXCR5P,TissueFactor,NKTPCCR7P,CD14N,ITMPPD1P))
# when databubble=RD Plotting for RD
#bubsub<-as.data.frame(databubble%>%select(Eosinophils,Neutrophil,FABP4,CD14N,TNFa,MonocytesPCXCR3P))

#bubsub<-as.data.frame(databubble%>%select(ITM,LBP,FABP4))

Corr=Hmisc::rcorr(as.matrix(bubsub),type='pearson') 
p<-Corr$P
r<-Corr$r

bubble<-ggplot(bubsub, aes(x=CD4P,y=CD8P,color = CD56DIMCD16P)) +
  geom_point(alpha=1, size=5)+
  scale_color_gradient(low="#fddde6", high="#DC143C")+
  geom_smooth(method=lm , color="black", se=FALSE,size=0.5) +
  #scale_size(range = c(.1, 10), name="Neutrophil")+
  theme_minimal()+
  theme(legend.position="right",
        axis.line= element_line(size = 0.5,color="darkgray"),
        legend.title=element_text(size=12,angle= 90),
        axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  #guides(color = guide_legend(title.position = 'left'))+
  #guides(color=guide_colourbar(title.vjust=1))+
  guides(color= guide_colourbar(title.position = "left"))+
  xlab("CD4P") +
  ylab("CD8P") 
 #scale_color_gradientn(colours = blues9(5))

bubble

ggsave(bubble, filename = "CovidFinalfeb/Mix/all_correlations/Multi/CD4P_CD8P_CD56DIMCD16P_multicorrICURED_fromall.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")




















dev.off()

library(reshape2) ## for melt()
dl  <- melt(data[,c(2,3,4)],id.var=1)

dcplot<-ggpubr::ggscatter(dl, y = "Neutrophil", x = "value",
                          add = "reg.line",conf.int = FALSE,
                          ylab="Neutrophil",xlab="",color = "variable")+
  stat_cor(aes(color = variable), label.x = 12)+
  
  theme(axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10))
ggsave(dcplot, filename = "Neutrophil_multi_corr_plots.png",dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")
dcplot


#bubbleplot 3 variables
dfall <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/twocytokine_innate_plasma_combined_Jan082021.csv")
dfin<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 12.27.2020_Krish.csv")

dfallICU<-dfall[c(13:31),]
dfinICU<-dfin[c(17:36),]

dfallRD<-dfall[c(32:55),]

data<-dfallICU
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
bubble1<-ggplot(data, aes(x=CD4P, y=CD8P,color = Neutrophil)) +
  geom_point(alpha=1, size=5)+
  scale_color_gradient(low="pink", high="red")+
  geom_smooth(method=lm , color="black", se=FALSE,size=0.5) +
  #scale_size(range = c(.1, 10), name="Neutrophil")+
  theme_minimal()+
  theme(legend.position="right",
        axis.line= element_line(size = 0.5,color="darkgray"),
        legend.title=element_text(size=12,angle= 90),
        axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  #guides(color = guide_legend(title.position = 'left'))+
  #guides(color=guide_colourbar(title.vjust=1))+
  guides(color= guide_colourbar(title.position = "left"))+
  xlab("FABP4") +
  ylab("MCP1") 
#scale_color_gradientn(colours = blues9(5))

bubble1

ggsave(bubble1, filename = "Eosinophils(RD)_multi_corr.png",dpi = 300, type = "cairo",  width = 5, height = 4, units = "in")

corrdf<-data%>% select(FABP4,MCP1,Eosinophils)
Corr<-Hmisc::rcorr(as.matrix(corrdf),type='pearson') # apply pearson correlation 
Corr$P
library(tidyverse)

# heat map
ggplot(data, aes(factor(CD4P), factor(CD8P), fill = Neutrophil)) + 
  geom_tile() +
  coord_fixed(ratio = 1)  + 
  scale_fill_gradientn(colours = rev(topo.colors(7))) # colourmap
ggsave(heat1, filename = "Neutrophil_multi_corr_plots3.png",dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")

# heat map2

data[,-1] %>% 
  slice(1:1000)%>% 
  ggplot(aes(x = CD4P, y = CD8P, color = Neutrophil))+
  geom_point(size = 5, alpha = 1) 

ggsave(heat2, filename = "Neutrophil_multi_corr_plots4.png",dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")


