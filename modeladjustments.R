# Impoart the dataset (.csv)
df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Manuja_made2/AD1 (WB) final - February 10th 2021.csv")
df_innate<-as_tibble(df_innate) %>% 
  mutate(AgeM = Age - min(df_innate$Age))

df<-df_innate[df_innate[, "Group"] == 'HD',c(1:4,ncol(df_innate),7:(ncol(df_innate)-1))]


beta=c()
pval<-c()
fstat<-c()
aic<-c()
model<-c()
marker<-c()
anova_p<-c()
model_comparison<-c()
for (i in seq(1:(ncol(df)-5))){
  print(i)
  df2<-df %>% dplyr::select(c(3,4,5, i+5))
  colname1 <-as.character(colnames(df2)[4])
  print(colname1)
  
  m1<-lm(get(colname1) ~ AgeM,data=df)
  m1_p<-pf(summary(m1)$fstatistic[1L], summary(m1)$fstatistic[2L], summary(m1)$fstatistic[3L], lower.tail = FALSE)
  m1_f<-summary(m1)$fstatistic[1L]
  m1_beta<-summary(m1)$coefficients[2L]

  m2<-lm(get(colname1) ~ Sex,data=df)
  m2_p<-pf(summary(m2)$fstatistic[1L], summary(m2)$fstatistic[2L], summary(m2)$fstatistic[3L], lower.tail = FALSE)
  m2_f<-summary(m2)$fstatistic[1L]
  m2_beta<-summary(m2)$coefficients[2L]
  
  m3<-lm(get(colname1) ~ AgeM+Sex,data=df)
  m3_p<-pf(summary(m3)$fstatistic[1L], summary(m3)$fstatistic[2L], summary(m3)$fstatistic[3L], lower.tail = FALSE)
  m3_f<-summary(m3)$fstatistic[1L]
  m3_agebeta_sexbeta<-paste0(summary(m3)$coefficients[2L],",",summary(m3)$coefficients[3L])

  pval<-c(pval,m1_p,m2_p,m3_p)
  fstat<-c(fstat,m1_f,m2_f,m3_f)
  beta<-c(beta,m1_beta,m2_beta,m3_agebeta_sexbeta)
  #a<-AIC(m1,m2,m3)$AIC
  #aic<-c(aic,a)
  
  modellist<-c(paste0(colname1,"_Age"),paste0(colname1,"_Sex"),paste0(colname1,"_Age+Sex"))
  model<-c(model,modellist)

  markertolist<-c(rep(paste0(colname1),3))
  marker<-c(marker,markertolist)
  cl<-c(paste0("Age_Sex"),paste0("Age_Age+Sex"),paste0("Sex_Age+Sex"))
  model_comparison<-c(model_comparison,cl)
  m1m2<-anova(m1,m2)[[6L]][2]
  m1m3<-anova(m1,m3)[[6L]][2]
  m2m3<-anova(m2,m3)[[6L]][2]
  anova_p<-c(anova_p,m1m2,m1m3,m2m3)
  
}
col0<-data.frame(marker)
col1<-data.frame(model)
col2<-data.frame(pval)
col3<-data.frame(fstat)
col4<-data.frame(beta)
col5<-data.frame(model_comparison)
col6<-data.frame(anova_p)
dataframe1<-cbind(col0,col1,col2,col3,col4)

dataframe2<-cbind(col0,col5,col6)
dataframe2<-dataframe2[!is.na(dataframe2$anova_p),]
dataframe2<-dataframe2%>%group_by(marker) %>% slice(which.min(anova_p))

modelsig<-dataframe1 %>% filter(pval<0.01)

write.csv(modelsig, file="wb_model_results.csv")

write.csv(dataframe2, file="wb_model_comparison.csv")







