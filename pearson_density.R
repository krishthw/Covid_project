library(ggplot2)
library(dplyr)
library(corrplot)
# import data files
df1<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Cytokines Anna's data final- February 10th 2021.csv")
df2<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Isoplexis Codeplex final- February 10th 2021.csv")
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
df4<-df4 %>% select(c("Sample","Group","Age","Sex","Diabetes","State",'Eosinophils','Neutrophil','Monocytes','NCM','CM','ITM','CD3P','CD56BRTCD16N','CD56DIMCD16P','CD4NCD8N','CD4PCD8P','CD4P','CD8P','NKT','CD19','CD14N'))

# combine cytokines
df<-inner_join(df1,df2, by="Sample") # innerjoin by 'Sample'
df1df2<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df2) <- gsub(".x","",colnames(df1df2)) # remove .x from column names
# join innate dataset
dfff<-inner_join(df1df2,df4, by="Sample")
df1df2df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df4) <- gsub(".x","",colnames(df1df2df4))
df1df2df4<- df1df2df4%>% select(-c(S)) # this contains cytokines and innate parents

data<-as.data.frame(df1df2df4) # convert to dataframe
CI_ICU<-data[data[, "Group"] == 'ICU',-c(1:5,41)] # 41 is the GranzymeB
CI_RD<-data[data[, "Group"] == 'RD',-c(1:5,41)]
CI_HD<-data[data[, "Group"] == 'HD',-c(1:5,41)]

ICUcors<- cor(CI_ICU,method="pearson")
ICUcorsdata<-as.data.frame(ICUcors[upper.tri(ICUcors, diag = FALSE)]) 
colnames(ICUcorsdata) <- "pearsonR"
ICUcorsdata['Group']= 'ICU'

RDcors<- cor(CI_RD,method="pearson")
RDcorsdata<-as.data.frame(RDcors[upper.tri(RDcors, diag = FALSE)] )
colnames(RDcorsdata) <- "pearsonR"
RDcorsdata['Group']= 'RD'

HDcors<- cor(CI_HD,method="pearson")
HDcorsdata<-as.data.frame(HDcors[upper.tri(HDcors, diag = FALSE)]  )
colnames(HDcorsdata) <- "pearsonR"
HDcorsdata['Group']= 'HD'

datadensity<-dplyr::bind_rows(ICUcorsdata, RDcorsdata,HDcorsdata)

d<-datadensity %>%
  ggplot( aes(x=pearsonR, color=Group)) +
  geom_density(alpha=0.6) +
  scale_color_manual(values=c("forestgreen","red2","cornflowerblue")) +
  theme(panel.background = element_rect(fill="white"),
      axis.line = element_line(size=0.5,colour="black"),
      axis.text.x = element_text(size=10,colour="black",face="bold"),
      axis.text.y = element_text(size=10,colour="black",face="bold")) +
  xlab("Pearson Correlation Coefficient") + ylab("Density")
print(d)
ggsave(d,  filename = paste0("CovidFinalfeb/Mix/density_plot_Cytokines+innate16.png"),
       dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")
