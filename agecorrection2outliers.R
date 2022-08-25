library(dplyr)
library(ggplot2)

# Impoart the dataset (.csv)
df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")

df_innate<-as_tibble(df_innate) %>% 
  mutate(AgeM = Age - min(df_innate$Age))

dfh<-df_innate[df_innate[, "Group"] == 'HD',c(1:4,ncol(df_innate),7:(ncol(df_innate)-1))]

#dfh<-df_innate[c(1:15),c(1:4,ncol(df_innate),7:(ncol(df_innate)-1))]
#------------------------------------------------------
# function to identify outliers
outlier <- function(x) {
  return(x < quantile(x, 0.25) - 3 * IQR(x) | x > quantile(x, 0.75) + 3 * IQR(x))
}

#outlier eliminated data
dfh_nu<-as.data.frame(dfh[,-c(1:5)]) # numeric data
dfh_no_na<-dfh_nu[,] # create new dataframe just like dfh_nu
dfh_no_mean<-dfh_nu[,] 
for (i in seq(1:ncol(dfh_nu))) {
  dfh_no_na[,i] = ifelse(outlier(dfh_nu[,i]), as.numeric(NA), as.numeric(dfh_nu[,i]))
  
  dfh_no_mean[,i] = ifelse(outlier(dfh_nu[,i]), mean(dfh_no_na[,i], na.rm = TRUE), as.numeric(dfh_nu[,i]))
}

df<-cbind(dfh[,1:5],dfh_no_na)
df_m<-cbind(dfh[,1:5],dfh_no_mean) # outliers are replaced by mean (with out outlier mean)


dfall<-df_innate[,c(1:4,ncol(df_innate),7:(ncol(df_innate)-1))]

beta0<-c()
beta1<-c()
pval<-c()
fstat<-c()
model<-c()
marker<-c()
for (i in seq(1:(ncol(df)-5))){
  print(i)
  df2<-df %>% dplyr::select(c(3,4,5, i+5))
  colname1 <-as.character(colnames(df2)[4])
  print(colname1)
  
  m1<-lm(get(colname1) ~ AgeM,data=df)
  m1_p<-pf(summary(m1)$fstatistic[1L], summary(m1)$fstatistic[2L], summary(m1)$fstatistic[3L], lower.tail = FALSE)
  m1_f<-summary(m1)$fstatistic[1L]
  m1_beta0<-summary(m1)$coefficients[1L]
  m1_beta1<-summary(m1)$coefficients[2L]
  
  pval<-c(pval,m1_p)
  fstat<-c(fstat,m1_f)
  beta0<-c(beta0,m1_beta0)
  beta1<-c(beta1,m1_beta1)
  
  #a<-AIC(m1,m2,m3)$AIC
  #aic<-c(aic,a)
  
  modellist<-c(paste0(colname1,"_Age"))
  model<-c(model,modellist)
  markertolist<-c(paste0(colname1))
  marker<-c(marker,markertolist)
  
}
col0<-data.frame(marker)
col1<-data.frame(model)
col2<-data.frame(pval)
col3<-data.frame(fstat)
col4<-data.frame(beta0)
col5<-data.frame(beta1)
dataframe1<-cbind(col0,col1,col2,col3,col4,col5)
modelsig<-dataframe1 %>% filter(pval<0.01)
write.csv(modelsig,file="sig_models_Age_only_without_extreme_outliers.csv")

# for loop to update outliers with values from models
for (i in seq(1:ncol(dfh_nu))) {
  df2<-dfh %>% dplyr::select(c(1))
  colname1 <-as.character(colnames(df2)[1])
  if (colname1 %in% modelsig$marker) {
    dfh_no_mean[,i] = ifelse(outlier(dfh_nu[,i]), (fit$residuals + fit$coefficients["(Intercept)"]), as.numeric(dfh_nu[,i]))
  }
}


for (i in seq(1:(ncol(dfall)-5))){
  df2<-dfall %>% dplyr::select(c(3,4,5, i+5))
  colname1 <-as.character(colnames(df2)[4])
  if (colname1 %in% modelsig$marker) {
    if (modelsig[modelsig[, 1] == colname1 ,6 ]<0){
      print(paste0("- beta1_",colname1))
      dfall[[colname1]] <- dfall[[colname1]] -((modelsig[modelsig[, 1]==colname1,6 ])*dfall$AgeM)
    } else {
      
      print(paste0("+ beta1_",colname1))
      dfall[[colname1]] <- dfall[[colname1]] -((modelsig[modelsig[, 1]==colname1,6 ])*dfall$AgeM)
    }
  }
}
write.csv(dfall,file="corrected_for_age_bias_beta0_beta1_minAge_included_innate_extreme.csv")

age<-c(12 ,11 , 5 ,25  ,1 ,15 , 4 , 4 ,17 , 8 , 6 ,19, 32,19 ,20)
marker<-marker<-c(0.73,  1.66,  0.12 , 0.31 , 0.73 , 0.84 , 0.40 , 0.51 , 0.28 , 0.23 , 0.69 ,28.90 ,52.90,29.60   , NA)
dataf<-as.data.frame(cbind(age,marker))

m2<-lm( marker~ age,data=dataf)
summary(m2)


