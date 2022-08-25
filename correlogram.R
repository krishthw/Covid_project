library(dplyr)
library(ggplot2)
library(DescTools)
library(rstatix)
library(tibble)
library(broom)
library(ggpubr)

library(corrplot)

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

# HD correlogram
CorrP_HD=Hmisc::rcorr(as.matrix(df_parent_HD),type='pearson') # apply pearson correlation 
png(height=4, width=4, file="CovidWB/Figures/correlogramParent_HD_upper.png", type = "cairo",units="in",res=1800)
correlogramP_HD<-corrplot(CorrP_HD$r,type="upper",tl.cex = 0.4,tl.col ="black",
                          p.mat = CorrP_HD$P, sig.level = 0.05, insig = "blank",
                          method='square', 
                          number.cex = 0.3,diag=F,
                          addCoef.col = "black",
                          title= "Correlation for parent markers \n combined with \n the significance test(p<0.05) for HD", 
                          cex.main = 1, font.main= 2,
                          mar=c(0,0,3,0))
dev.off()

png(height=4, width=4, file="CovidWB/Figures/correlogramParent_HD_upper_pic2.png", type = "cairo",units="in",res=1800)
correlogramP_HD2<-corrplot(CorrP_HD$r,type="upper",tl.cex = 0.4,tl.col ="black",
                           p.mat = CorrP_HD$P ,insig = "label_sig",method='square',diag=F,
                           sig.level = c(.001, .01, .05),pch.cex = .5, pch.col = "white",
                           title="Correlation for parent markers \n combined with \n the significance test(p<0.05) for HD",
                           cex.main = 1, font.main= 2,mar=c(0,0,4,0))

dev.off()

# ICU correlogram
CorrP_ICU=Hmisc::rcorr(as.matrix(df_parent_ICU),type='pearson') # apply pearson correlation 
png(height=4, width=4, file="CovidWB/Figures/correlogramParent_ICU_upper.png", type = "cairo",units="in",res=1800)
correlogramP_ICU<-corrplot(CorrP_ICU$r,type="upper",tl.cex = 0.4,tl.col ="black",
                           p.mat = CorrP_ICU$P, sig.level = 0.05, insig = "blank",
                           method='square', 
                           number.cex = 0.3,diag=F,
                           addCoef.col = "black",
                           title= "Correlation for parent markers \n combined with \n the significance test(p<0.05) for ICU", 
                           cex.main = 1, font.main= 2,
                           mar=c(0,0,3,0))
dev.off()

png(height=4, width=4, file="CovidWB/Figures/correlogramParent_ICU_upper_pic2.png", type = "cairo",units="in",res=1800)
correlogramP_ICU2<-corrplot(CorrP_ICU$r,type="upper",tl.cex = 0.4,tl.col ="black",
                            p.mat = CorrP_ICU$P ,insig = "label_sig",method='square',diag=F,
                            sig.level = c(.001, .01, .05),pch.cex = .5, pch.col = "white",
                            title="Correlation for parent markers \n combined with \n the significance test(p<0.05) for ICU",
                            cex.main = 1, font.main= 2,mar=c(0,0,4,0))

dev.off()

# RD correlogram

CorrP_RD=Hmisc::rcorr(as.matrix(df_parent_RD),type='pearson') # apply pearson correlation 
png(height=4, width=4, file="CovidWB/Figures/correlogramParent_RD_upper.png", type = "cairo",units="in",res=1800)
correlogramP_ICU<-corrplot(CorrP_RD$r,type="upper",tl.cex = 0.4,tl.col ="black",
                           p.mat = CorrP_RD$P, sig.level = 0.05, insig = "blank",
                           method='square', 
                           number.cex = 0.3,diag=F,
                           addCoef.col = "black",
                           title= "Correlation for parent markers \n combined with \n the significance test(p<0.05) for RD", 
                           cex.main = 1, font.main= 2,
                           mar=c(0,0,3,0))
dev.off()

png(height=4, width=4, file="CovidWB/Figures/correlogramParent_RD_upper_pic2.png", type = "cairo",units="in",res=1800)
correlogramP_RD2<-corrplot(CorrP_RD$r,type="upper",tl.cex = 0.4,tl.col ="black",
                           p.mat = CorrP_RD$P ,insig = "label_sig",method='square',diag=F,
                           sig.level = c(.001, .01, .05),pch.cex = .5, pch.col = "white",
                           title="Correlation for parent markers \n combined with \n the significance test(p<0.05) for RD",
                           cex.main = 1, font.main= 2,mar=c(0,0,4,0))

dev.off()
# _________________________________________________________________________________________


