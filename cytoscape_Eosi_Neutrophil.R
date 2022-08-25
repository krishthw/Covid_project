

#Eosi ~ Eotoxin  correlations CI_ICU from density plot.r

ICU<-CI_ICU[,c(2,3,47:61)]
CorrICU<-Hmisc::rcorr(as.matrix(ICU),type='pearson') # apply pearson correlation 
C1<-flattenCorrMatrix(CorrICU$r, CorrICU$P) # extract correlation matrix and p-value matrix
C1 <-C1[order(-C1$cor),]
write.csv(C1,"CovidFinalfeb/Mix/EotoxinEosi/EotoxinEosi_ICUcorrelations.csv")

RD<-CI_RD[,c(2,3,47:61)]
CorrRD<-Hmisc::rcorr(as.matrix(RD),type='pearson') # apply pearson correlation 
C2<-flattenCorrMatrix(CorrRD$r, CorrRD$P) # extract correlation matrix and p-value matrix
C2 <-C2[order(-C2$cor),]
write.csv(C2,"CovidFinalfeb/Mix/EotoxinEosi/EotoxinEosi_RDcorrelations.csv")

HD<-CI_HD[,c(2,3,47:61)]
CorrHD<-Hmisc::rcorr(as.matrix(HD),type='pearson') # apply pearson correlation 
C3<-flattenCorrMatrix(CorrHD$r, CorrHD$P) # extract correlation matrix and p-value matrix
C3 <-C3[order(-C3$cor),]
write.csv(C3,"CovidFinalfeb/Mix/EotoxinEosi/EotoxinEosi_HDcorrelations.csv")

png(height=4, width=4, file="CovidFinalfeb/Mix/EotoxinEosi/correlogramEotoxin_eosi_HD.png", type = "cairo",units="in",res=1800)
correlogramHD<-corrplot(CorrHD$r,type="upper",tl.cex = 0.4,tl.col ="black",
                           p.mat = CorrHD$P, sig.level = 0.05, insig = "blank",
                           method='square', 
                           number.cex = 0.3,diag=F,
                           addCoef.col = "black",
                           title= "Eotoxin_eosi_HD", 
                           cex.main = 1, font.main= 2,
                           mar=c(0,0,3,0))
dev.off()
# correlation plots
C=C1[order(-abs(C1$cor)),] # order data according to corretaion value (Descending)
C=na.omit(C) # filter out NA values
#C_sig=C[abs(C$cor)>0.5 & abs(C$p)<0.05,]
C_sig=C[abs(C$p)<0.05,]

df1<-ICU
df2<-C_sig

#Eosi ~ cytokine correlations CI_ICU from density plot.r
ICU<-CI_ICU #[,1:35]
CorrICU<-Hmisc::rcorr(as.matrix(ICU),type='pearson') # apply pearson correlation 
C1<-flattenCorrMatrix(CorrICU$r, CorrICU$P) # extract correlation matrix and p-value matrix
C1 <-C1[order(-C1$cor),]
C1<-na.omit(C1) # filter out NA values
C1<-C1 %>% rename(metab1=row,matab2=column,pcor=cor,pval=p)

write.csv(C1,"for_cytospaceEosicytokineannaICU_test3.csv")

png(height=6, width=6, file="test3.png", type = "cairo",units="in",res=1800)
correlogram<-corrplot(CorrICU$r,type="upper",tl.cex = 0.4,tl.col ="black",
                        p.mat = CorrICU$P, sig.level = 0.05, insig = "blank",
                        method='square', 
                        number.cex = 0.3,diag=F,
                        #addCoef.col = "black",
                        title= "test", 
                        cex.main = 1, font.main= 2,
                        mar=c(0,0,3,0))
dev.off()

# correlations Eosi ~ with all df from PCAfinal.R for all data

data<-as.data.frame(df) # convert to dataframe
EICU<-data[data[, "Group"] == 'ICU',-c(1,58:71)] 
ERD<-data[data[, "Group"] == 'RD',-c(1,58:71)]
EHD<-data[data[, "Group"] == 'HD',-c(1,58:71)]

CorrICU<-Hmisc::rcorr(as.matrix(EICU),type='pearson') # apply pearson correlation 
C<-flattenCorrMatrix(CorrICU$r, CorrICU$P) # extract correlation matrix and p-value matrix
C <-C[order(-C$cor),]
C<-na.omit(C) # filter out NA values
C_sig=C[abs(C$p)<0.05,]

C_sig<-C_sig %>% rename(metab1=row,metab2=column,pcor=cor,pval=p)
Ctable<-C_sig[C_sig[,"metab1"] == 'Eosinophils' | C_sig[, "metab2"] == 'Eosinophils',] 
write.csv(Ctable,"for_cytospaceEosi_allother_ICU.csv")

data<-as.data.frame(df) # convert to dataframe
NICU<-data[data[, "Group"] == 'ICU',-c(1,73:87)] 
NRD<-data[data[, "Group"] == 'RD',-c(1,73:87)]
NHD<-data[data[, "Group"] == 'HD',-c(1,73:87)]

CorrICU<-Hmisc::rcorr(as.matrix(NHD),type='pearson') # apply pearson correlation 
C<-flattenCorrMatrix(CorrICU$r, CorrICU$P) # extract correlation matrix and p-value matrix
C <-C[order(-C$cor),]
C<-na.omit(C) # filter out NA values
C_sig=C[abs(C$p)<0.05,]

C_sig<-C_sig %>% rename(metab1=row,metab2=column,pcor=cor,pval=p)
Ctable<-C_sig[C_sig[,"metab1"] == 'Neutrophil' | C_sig[, "metab2"] == 'Neutrophil',] 
write.csv(Ctable,"for_cytospaceNeutrophil_allother_HD.csv")

# correlations NKT ~ with all df from PCAfinal.R for all data

NKTICU<-data[data[, "Group"] == 'ICU',-c(1,249:262)] 
NKTRD<-data[data[, "Group"] == 'RD',-c(1,249:262)]
NKTHD<-data[data[, "Group"] == 'HD',-c(1,249:262)]

CorrICU<-Hmisc::rcorr(as.matrix(NKTHD),type='pearson') # apply pearson correlation 
C<-flattenCorrMatrix(CorrICU$r, CorrICU$P) # extract correlation matrix and p-value matrix
C <-C[order(-C$cor),]
C<-na.omit(C) # filter out NA values
C_sig=C[abs(C$p)<0.05,]

C_sig<-C_sig %>% rename(metab1=row,metab2=column,pcor=cor,pval=p)
Ctable<-C_sig[C_sig[,"metab1"] == 'NKT' | C_sig[, "metab2"] == 'NKT',] 
write.csv(Ctable,"for_cytospaceNKT_allother_HD.csv")


