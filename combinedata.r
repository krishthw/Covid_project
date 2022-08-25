# This code was used to join Cytokine data, Plasma data and Innate Panel
library(dplyr)
# import data files
df1<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Cytokines Anna's data Updated 12.25.2020_Krish.csv")
df2<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Isoplexis Codeplex- 12.26.2020_Krish.csv")
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Plasma ELISA - 12.26.2020_Krish.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 12.27.2020_Krish.csv")

colnames(df1) <- gsub("-","",colnames(df1))
colnames(df2) <- gsub("-","",colnames(df2))

# make sample names the same in all data frames
df1$Sample <- gsub("_","",df1$Sample) 
df2$Sample <- gsub("_","",df2$Sample)
df3$Sample <- gsub("_","",df3$Sample)
df4$Sample <- gsub("_","",df4$Sample)

# join Cytokine datasets
df<-inner_join(df1,df2, by="Sample") # innerjoin by 'Sample'
df1df2<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df2) <- gsub(".x","",colnames(df1df2)) # remove .x from column names
# join Plasma dataset
dff<-inner_join(df1df2,df3, by="Sample")
df1df2df3<- dff %>% select(-ends_with(".y"))
colnames(df1df2df3) <- gsub(".x","",colnames(df1df2df3))
df1df2df3<- df1df2df3%>% select(-(S))
# join innate dataset
dfff<-inner_join(df1df2df3,df4, by="Sample")
df1df2df3df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df3df4) <- gsub(".x","",colnames(df1df2df3df4))
df1df2df3df4<- df1df2df3df4%>% select(-c(S))


write.csv(df1df2df3df4,'twocytokine_innate_plasma_combined_Krish.csv')

data<-df1df2df3df4

colnames(data) <- gsub("(ng/mL)","",colnames(data))
colnames(data) <- gsub("(pg/mL)","",colnames(data))
colnames(data) <- gsub("(u/mL)","",colnames(data))
colnames(data) <- gsub(" ","",colnames(data))
colnames(data) <- gsub("-","",colnames(data))
colnames(data) <- gsub("(","",colnames(data),fixed=T)
colnames(data) <- gsub(")","",colnames(data),fixed=T)
write.csv(data,'twocytokine_innate_plasma_combined_Krish1.csv')

dataICU<-data[c(13:31),] # ICU patient data
dataRD<-data[c(32:55),] # Recovered patient data

# following subsets were used to make correlograms
dataRDEosi<-dataRD[,c(64,9:63)] 
dataICU_Eosi_CM_ITM_NCM<-dataICU[,c(64,163,176,190,9:63)]
dataICU_Eosi_CM_ITM_NCM<-dataICU_Eosi_CM_ITM_NCM[,c(1:14,43:59)]

dataICU_NKT_NKTPCD69P<-dataICU[,c(325,328,9:63)]
dataRD_NKT_NKTPCD69P<-dataRD[,c(325,328,9:63)]




Corr=Hmisc::rcorr(as.matrix(dataRD_NKT_NKTPCD69P),type='pearson') # apply pearson correlation 
png(height=4, width=4, file="Eosi_CytoPlsma_RD_corr.png", type = "cairo",units="in",res=1800)
COL <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
correlogramP_HD<-corrplot(Corr$r,type="upper",tl.cex = 0.2,tl.col ="black",
                          p.mat = Corr$P, sig.level = 0.05, insig = "blank",
                          method='square', col = COL(200),
                          number.cex = 0.1,diag=F,
                          addCoef.col = "black",
                          #title= "Correlation for parent markers \n combined with \n the significance test(p<0.05) for HD", 
                          #cex.main = 1, font.main= 2,
                          mar=c(0,0,3,0))
dev.off()
df<-dataRD_NKT_NKTPCD69P


lista<-colnames(df[,1:2])
listb<-colnames(df[,3:57])

C=flattenCorrMatrix(Corr$r, Corr$P) # extract correlation matrix and p-value matrix
C=C[order(-abs(C$cor)),] # order data according to corretaion value (Descending)
C=na.omit(C) # filter out NA values
#C_sig=C[abs(C$cor)>0.5 & abs(C$p)<0.05,]
C_sig=C[abs(C$p)<0.05,]

df2<-C_sig


# 1/15 /2021 
# Dhanuja 27 RD patients me 24 RD patients
df3<-df3[c(1:67),]
# join CytokineAnna to Plasma
df<-inner_join(df1,df3, by="Sample") # innerjoin by 'Sample'
df1df3<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df3) <- gsub(".x","",colnames(df1df3)) # remove .x from column names

# join innate
dff<-inner_join(df1df3,df4, by="Sample")
df1df3df4<- dff %>% select(-ends_with(".y"))
colnames(df1df3df4) <- gsub(".x","",colnames(df1df3df4))
df1df3df4<- df1df3df4%>% select(-(S))
#df1df3df4<- df1df3df4%>% select(-c(V16,V17))

print(df1df3 %>% anti_join(dff)) 

dfff<-inner_join(df1df3df4,df2, by="Sample")

df1df2df3df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df3df4) <- gsub(".x","",colnames(df1df2df3df4))
df1df2df3df4<- df1df2df3df4%>% select(-c(S))


ICU58ia<-df1df3df4[c(13:31),c(9:48,57,77,137,156,169,183,200,201,213,230,231,251,297,278,318,350)]
RD58ia<-df1df3df4[c(32:58),c(9:48,57,77,137,156,169,183,200,201,213,230,231,251,297,278,318,350)]
HD58ia<-df1df3df4[c(1:12),c(9:48,57,77,137,156,169,183,200,201,213,230,231,251,297,278,318,350)]

ICU2<-df1df2df3df4[c(13:31),c(349,9:48,407:413,49:56)]
RD55<-df1df2df3df4[c(32:55),c(9:57,77,137,156,169,183,200,201,230,231,251,278,297,318,334,349,407:413)]
RD58<-df1df3df4[c(32:58),c(9:57,77,137,156,169,183,200,201,230,231,251,278,297,318,334,349)]

data<-HD58ia
colnames(data) <- gsub("(ng/mL)","",colnames(data))
colnames(data) <- gsub("(pg/mL)","",colnames(data))
colnames(data) <- gsub("(u/mL)","",colnames(data))
colnames(data) <- gsub(" ","",colnames(data))
colnames(data) <- gsub("-","",colnames(data))
colnames(data) <- gsub("(","",colnames(data),fixed=T)
colnames(data) <- gsub(")","",colnames(data),fixed=T)


heatmap(as.matrix(CorrICU$r))

lista<-c("Eosinophils")
listb<-colnames(data)
C=flattenCorrMatrix(Corr$r, Corr$P) # extract correlation matrix and p-value matrix
C=C[order(-abs(C$cor)),] # order data according to corretaion value (Descending)
C=na.omit(C) # filter out NA values
#C_sig=C[abs(C$cor)>0.5 & abs(C$p)<0.05,]
C_sig=C[abs(C$p)<0.05,]

df2<-C_sig


data<-data
Corr=Hmisc::rcorr(as.matrix(data),type='pearson') # apply pearson correlation 
png(height=8, width=10, file="HD15_Innate_Codeplex.png", type = "cairo",units="in",res=1800)

colour<- colorRampPalette(c("red", "white", "blue"))(20)
corrplot(Corr$r,tl.cex = 0.5,tl.col ="black",col=brewer.pal(n=10, name="PuOr"),
         method='color',
         diag=F)
corrRect(c(23,16),lwd=1)
corrRect(39,lwd=1)
dev.off()

# join innate to Plasma
df<-inner_join(df3,df4, by="Sample") # innerjoin by 'Sample'
df4df3<- df %>% select(-ends_with(".y")) # remove duplicates 
colnames(df4df3) <- gsub(".x","",colnames(df4df3)) # remove .x from column names

ICU65ip<-df4df3[c(17:36),c(8:16,36,96,115,128,142,159,160,172,189,190,210,256,237,277,308)]
RD65ip<-df4df3[c(37:66),c(8:16,36,96,115,128,142,159,160,172,189,190,210,256,237,277,308)]
HD65ip<-df4df3[c(1:15),c(8:16,36,96,115,128,142,159,160,172,189,190,210,256,237,277,308)]
data<-RD65ip

# join innate to codeplex
df<-inner_join(df2,df4, by="Sample") # innerjoin by 'Sample'
df4df2<- df %>% select(-ends_with(".y")) # remove duplicates 
colnames(df4df2) <- gsub(".x","",colnames(df4df2)) # remove .x from column names

ICU62ico<-df4df2[c(17:36),c(7:30,50,110,129,142,156,173,174,186,203,204,224,270,251,291,322)]
RD62ico<-df4df2[c(37:63),c(7:30,50,110,129,142,156,173,174,186,203,204,224,270,251,291,322)]
HD62ico<-df4df2[c(1:15),c(7:30,50,110,129,142,156,173,174,186,203,204,224,270,251,291,322)]
data<-HD62ico

