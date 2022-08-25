library(dplyr)
library(ggplot2)

# Impoart the dataset (.csv)
df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February15th2021.csv")

df_innate<-as_tibble(df_innate) %>% 
  mutate(AgeM = Age - min(df_innate$Age))

df<-df_innate[df_innate[, "Group"] == 'HD',c(1:4,ncol(df_innate),7:(ncol(df_innate)-1))]

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
write.csv(modelsig,file="sig_models_.csv")

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
write.csv(dfall,file="corrected_for_age_bias_beta1_minAge_included_feb16.csv")

# with beta 0


for (i in seq(1:(ncol(dfall)-5))){
  df2<-dfall %>% dplyr::select(c(3,4,5, i+5))
  colname1 <-as.character(colnames(df2)[4])
  if (colname1 %in% modelsig$marker) {
    if (modelsig[modelsig[, 1] == colname1 ,6 ]<0){
      print(paste0("- beta1_",colname1))
      dfall[[colname1]] <- dfall[[colname1]] -((modelsig[modelsig[, 1]==colname1,5 ])+
                                                 (modelsig[modelsig[, 1]==colname1,6 ])*dfall$AgeM)
    } else {
  
      print(paste0("+ beta1_",colname1))
      dfall[[colname1]] <- dfall[[colname1]] -((modelsig[modelsig[, 1]==colname1,5 ])+
                                                 (modelsig[modelsig[, 1]==colname1,6 ])*dfall$AgeM)
    }
  }
}
write.csv(dfall,file="corrected_for_age_bias_beta1_minAge_included_forcorrectedinnatefeb15_beta0.csv")


Z <- c(0.9,1.2,3.0,4.5,0.8,0.4)
age <- c(30,22,45,60,33,20)

fit <-lm(Z~age)
summary(fit)
adjusted<- fit$residuals + fit$coefficients["(Intercept)"]
