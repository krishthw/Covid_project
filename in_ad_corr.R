library(dplyr)
library(ggplot2)
library(DescTools)
library(rstatix)
library(tibble)
library(broom)
library(ggpubr)

library(corrplot)


# checking correlations between adaptive and innate panel for ICU group
# only samples in both were selected
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/adaptive_innate_trial2.csv")
zero<-df[, colSums(df[,-1]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
df<-df %>% select(Group, names(nonzero))

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
colnames(df)

Corr_in_ad=Hmisc::rcorr(as.matrix(df),type='pearson') # apply pearson correlation 
png(height=4, width=4, file="Correlogram_t2.png", type = "cairo",units="in",res=1800)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
correlogram<-corrplot(Corr_in_ad$r,type="upper",tl.cex = 0.4,tl.col ="black",
                      p.mat = Corr_in_ad$P, sig.level = 0.05, insig = "blank",
                      method='square', col = col(200),
                      number.cex = 0.2,diag=F,
                      addCoef.col = "black",
                      title= "Correlation for parent markers of innate and adaptive panels\n combined with \n the significance test(p<0.05) for ICU", 
                      cex.main = 0.5, font.main= 1,
                      mar=c(0,0,3,0))
dev.off()

png(height=4, width=4, file="correlogram_pic2.png", type = "cairo",units="in",res=1800)
correlogram2<-corrplot(Corr_in_ad$r,type="upper",tl.cex = 0.4,tl.col ="black",
                       p.mat = Corr_in_ad$P ,insig = "label_sig",method='square',diag=F,
                       sig.level = c(.001, .01, .05),pch.cex = .4, pch.col = "white",
                       title="Correlation for parent markers of innate and adaptive panels\n combined with \n the significance test(p<0.05) for ICU",
                       cex.main = 0.5, font.main= 1,mar=c(0,0,4,0))

dev.off()

#----------------------------------------------------------------------------------------
lista<-colnames(df[,1:20])
listb<-colnames(df[,21:35])

for (i in seq(1:nrow(df2))){
  colname1 <-as.character(df2[i,1])
  colname2<-as.character(df2[i,2])
  #df3<-df1 %>% dplyr::select (all_of(colname1), all_of(colname2))
  if (((colname1 %in% lista) & (colname2 %in% listb)) || ((colname2 %in% lista) & (colname1 %in% listb))){
    print(colname1)
    print(colname2)
    cplot<-ggpubr::ggscatter(df, x = colname2, y = colname1,
                             add = "reg.line",conf.int = TRUE,
                             xlab=colname2, ylab=colname1,
                             #color="red2",
                             title=paste0(colname1 ," vs " ,colname2))+
      ggpubr::stat_cor(method = "pearson", label.x.npc=0.25, label.y.npc= 1.0,r.accuracy=0.01,p.accuracy=0.001) +
      #scale_y_continuous(labels=function(x) paste0(x*100, "%"))+
      #scale_x_continuous(labels = scales::number_format(accuracy = 0.0001))+
      theme(axis.text.x=element_text(size=10),
            axis.title.y=element_text(size=12),
            axis.title.x=element_text(size=12),
            axis.text.y=element_text(size=10))
    cplot<-cplot+geom_point( shape = 21,fill='red2')
    print(cplot)
    ggsave(cplot, filename = paste0("CovidDec4/Figures/in_ad_correlations/", colname1, "vs", colname2, ".png"),dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")
  }
}