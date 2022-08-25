# this code prepare two CSV files; one conatining p and r values and one containing group 
# to feed to Cytoscappe to make correlation based network
library(dplyr)
df<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 12.27.2020_Krish.csv")
df<-df[,c(2,7:356)] # select group and variables

# check for nonzero columns
nozero<-df[, colSums(df[,-1]) > 0] # get logical vector to show nonzero columnss
nonzero<-which(nozero == TRUE) 

# replace n/a with group means
df<-df %>% 
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

# ICU Group
df<-as.data.frame(df)
df<-df[c(17:36),]

duplicated(colnames(df)) # if duplicated names present output "TRUE"

# remove unnecessary characters
colnames(df) <- gsub("(ng/mL)","",colnames(df))
colnames(df) <- gsub("(pg/mL)","",colnames(df))
colnames(df) <- gsub("(u/mL)","",colnames(df))
colnames(df) <- gsub(" ","",colnames(df))
colnames(df) <- gsub("-","",colnames(df))
colnames(df) <- gsub("(","",colnames(df),fixed=T)
colnames(df) <- gsub(")","",colnames(df),fixed=T)

# remove zero column #> df$CMPCD69PCD27P
df<-df[,-106]
#remove all CD3PCD19P
df<-df[,-c(278:292)]

# making  of dataframe 
df1<-data.frame(colnames(df[,-1]))
df1<-df1%>%rename(name=colnames.df....1..)

grouplist<-c('Eosinophils','Neutrophil','Monocytes','NCM','CM','ITM','CD3P','CD56BRTCD16BRT','CD56BRTCD16N','CD56DIMCD16N','CD56DIMCD16P','CD56NCD16P','CD4NCD8N','CD4PCD8P','CD4P','CD8P','NKT','CD19','CD14N')
library(stringr)
conditions <- purrr::map(grouplist, ~ quo(str_detect(name, fixed(!!.x, ignore_case = T))~!!.x))
df1<-df1 %>% mutate(cluster = case_when(!!!conditions) )

write.csv(df1,"for_cytospace_groups1.csv")


# Correlations
Corr=Hmisc::rcorr(as.matrix(df[,-1]),type='pearson') # apply pearson correlation 
# _____________________________________________
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
#  _____________________________________________

C=flattenCorrMatrix(Corr$r, Corr$P) # extract correlation matrix and p-value matrix
C1=C[order(-abs(C$cor)),] # order data according to corretaion value (Descending)
C2=na.omit(C1) # filter out NA values
C2<-C2 %>% rename(metab1=row,matab2=column,pcor=cor,pval=p)

write.csv(C2,"for_cytospace1.csv")




