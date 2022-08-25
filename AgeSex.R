
df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Manuja_made/IN2 Color Coded- February 2nd 2021.csv")
tags <- c("<30","[30-60)", "[60-100)") # define age brackets
cdf <- as_tibble(df_innate) %>% 
  mutate(AgeGroup = case_when(
    Age < 30 ~ tags[1],
    Age  >= 30 & Age  < 60 ~ tags[2],
    Age  >= 60 & Age  < 100 ~ tags[3]))


parent_cdf<-cdf%>% select(Eosinophils,Neutrophil,Monocytes,CM)#,ITM)#NCM,CD56BRTCD16N,CD56DIMCD16P,CD3P,CD14N,CD4NCD8N,CD4P,CD4PCD8P,CD8P,NKT,CD3PCD19P,CD19)

min<-as.data.frame(apply(parent_cdf,2, FUN = function(x) {min(x[x > 0])}))

# Sex log10
for (i in seq(1:(ncol(parent_cdf)))){
  df2<-parent_cdf %>% dplyr::select(c(i))
  a<-colnames(df2)
  colname <-as.character(a[1])
  print(min[i,])
  logdf<-parent_cdf%>%
   mutate(logpar = ifelse(get(colname) == 0, log10(0.001 * min[i,]), log10(get(colname))))
  
  plotdf<-cbind(cdf[,1:4],logdf)
  p<-ggplot(plotdf, aes(Group, logpar, colour = factor(Sex))) + geom_point() +
    stat_summary(aes(group = factor(Sex)), fun = mean, geom = "line")+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=0.7,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold"),
          legend.position = "bottom")+
    xlab("") + ylab(paste0("log10 ",colname))+
    labs(colour= "Sex")+
    scale_color_manual(values = c("#000080","#DE3163"))
  #ggsave(p,  filename = paste0("Innate_sex_02032021/log10/",colname, "_log10_sex.png"),
  #       dpi = 300, type = "cairo",  width = 6, height = 5, units = "in")
  new_col_name <- paste0("log", colname)
  parent_cdf<-parent_cdf%>%
    mutate(!!sym(new_col_name) := ifelse(get(colname) == 0, log10(0.001 * min[i,]), log10(get(colname))))
  }

logsex_mean<-parent_cdf_log[,c(2,4,6:22)] %>% 
  group_by(Group,Sex) %>% 
  summarise_all("mean")

# Sex raw
for (i in seq(1:(ncol(parent_cdf)))){
  df2<-parent_cdf %>% dplyr::select(c(i))
  a<-colnames(df2)
  colname <-as.character(a[1])
  plotdf_raw<-cbind(cdf[,c(1:4,326)],parent_cdf)
  p<-ggplot(plotdf_raw, aes(Group, get(colname), colour = factor(Sex))) + geom_point() +
    stat_summary(aes(group = factor(Sex)), fun = mean, geom = "line")+
    
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=0.7,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold"),
          legend.position = "bottom")+
    xlab("") + ylab(paste0(colname))+
    labs(colour= "Sex")+
    scale_color_manual(values = c("#000080","#DE3163"))
  
  #ggsave(p,  filename = paste0("Innate_sex_02032021/raw/",colname, "_raw_sex.png"),
   #      dpi = 300, type = "cairo",  width = 6, height = 5, units = "in")
}
# sex mean values to csv 
rawsex_mean<-plotdf_raw[,c(2,4,6:22)] %>% 
  group_by(Group,Sex) %>% 
  summarise_all("mean")

# Age log10
for (i in seq(1:(ncol(parent_cdf)))){
  df2<-parent_cdf %>% dplyr::select(c(i))
  a<-colnames(df2)
  colname <-as.character(a[1])
  logdf<-parent_cdf%>%
    mutate(logpar = ifelse(get(colname) == 0, log10(0.001 * min[i,]), log10(get(colname))))
  plotdf<-cbind(cdf[,c(1:4,326)],logdf)
  p<-ggplot(plotdf, aes(Group, logpar, colour = factor(AgeGroup))) + geom_point() +
    stat_summary(aes(group = factor(AgeGroup)), fun = mean, geom = "line")+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=0.7,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold"),
          legend.position = "bottom")+
    xlab("") + ylab(paste0("log10 ",colname))+
    labs(colour= "Age Group")+
    scale_color_manual(values = c("#750985","#FF985E","#326633"))
  #ggsave(p,  filename = paste0("Innate_age_02032021/log10/",colname, "_log10_age.png"),
   #      dpi = 300, type = "cairo",  width = 6, height = 5, units = "in")
  new_col_name <- paste0("log", colname)
  parent_cdf<-parent_cdf%>%
    mutate(!!sym(new_col_name) := ifelse(get(colname) == 0, log10(0.001 * min[i,]), log10(get(colname))))
  }

logage_mean<-parent_cdf_log[,c(2,5:22)] %>% 
  group_by(Group,AgeGroup) %>% 
  summarise_all("mean")

# Age raw
for (i in seq(1:(ncol(parent_cdf)))){
  df2<-parent_cdf %>% dplyr::select(c(i))
  a<-colnames(df2)
  colname <-as.character(a[1])
  plotdf_raw<-cbind(cdf[,c(1:4,326)],parent_cdf)
  p<-ggplot(plotdf_raw, aes(Group, get(colname), colour = factor(AgeGroup))) + geom_point() +
    stat_summary(aes(group = factor(AgeGroup)), fun = mean, geom = "line")+
    
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=0.7,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold"),
          legend.position = "bottom")+
    xlab("") + ylab(paste0(colname))+
    labs(colour= "AgeGroup")+
    scale_color_manual(values = c("#750985","#FF985E","#326633"))
  
 # ggsave(p,  filename = paste0("Innate_age_02032021/raw/",colname, "_raw_age.png"),
    #     dpi = 300, type = "cairo",  width = 6, height = 5, units = "in")
}


rawage_mean<-plotdf_raw[,c(2,5:22)] %>% 
  group_by(Group,AgeGroup) %>% 
  summarise_all("mean")

names(rawsex_mean)[c(1,2)] <- c("Group1","Group2")
names(rawage_mean)[c(1,2)] <- c("Group1","Group2")
names(logsex_mean)[c(1,2)] <- c("Group1","Group2")
names(logage_mean)[c(1,2)] <- c("Group1","Group2")
sex_age_raw_log_mean<-rbind(rawsex_mean,rawage_mean,logsex_mean,logage_mean)
write.csv(sex_age_raw_log_mean,"Innate_age_02032021/sex_age_raw_log_mean_02032021.csv")

log10data<-rbind(plotdf_raw,parent_cdf_log)
write.csv(log10data,"Innate_age_02032021/sex_age_raw_log_data_02032021.csv")
