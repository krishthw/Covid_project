library(ggfortify)

#PCA final
# import data files
df1<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Cytokines Anna's data final- February 10th 2021.csv")
df2<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Isoplexis Codeplex final- February 10th 2021.csv")
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
df4<-df4 %>% select(c("Sample","Group","Age","Sex","Diabetes","State",'Eosinophils','Neutrophil','Monocytes','NCM','CM','ITM','CD3P','CD56BRTCD16N','CD56DIMCD16P','CD4NCD8N','CD4PCD8P','CD4P','CD8P','NKT','CD19','CD14N'))
# remove - from codeplex and anna data to identify if there is any same markers
colnames(df1) <- gsub("-","",colnames(df1)) 
colnames(df2) <- gsub("-","",colnames(df2))
colnames(df3) <- gsub("(pg/mL)","",colnames(df3))
colnames(df3) <- gsub("(","",colnames(df3),fixed=T)
colnames(df3) <- gsub(")","",colnames(df3),fixed=T)

(anti_join(df1, df2, by="Sample"))$Sample # 3 missing
(anti_join(df2, df1, by="Sample"))$Sample # 7 missing

# join Cytokine datasets
df<-inner_join(df1,df2, by="Sample") # innerjoin by 'Sample'
df1df2<- df %>% select(-ends_with(".y")) # remove duplicates (from Cytocodeplex data)
colnames(df1df2) <- gsub("\\.x$","",colnames(df1df2)) # remove .x from column names

# join Plasma dataset
dff<-inner_join(df1df2,df3, by="Sample")
df1df2df3<- dff %>% select(-ends_with(".y"))
colnames(df1df2df3) <- gsub("\\.x$","",colnames(df1df2df3))

# join innate dataset
dfff<-inner_join(df1df2df3,df4, by="Sample")
df1df2df3df4<- dfff %>% select(-ends_with(".y"))
colnames(df1df2df3df4) <- gsub("\\.x$","",colnames(df1df2df3df4))

data<-df1df2df3df4# change as you want

df<-data %>% 
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))%>% ungroup 

df<-df%>%
  mutate(Group= ifelse(State =="Deceased", paste0("ICU-Deceased"),Group))
  

df<- df%>% select(-c(1,3:6))

# remove all zero columns
#zero<-df[, colSums(df[,-c(1:6)]) > 0] # get logical vector to show 
#nonzero<-which(zero == TRUE)
#df<-df %>% select(Group, names(nonzero))

df %>% select(-Group)%>%
  scale() %>%                  # scale to 0 mean and unit variance
  prcomp() ->pca 

# display the results from the PCA 
head(pca$x)
pca_data <- data.frame(pca$x, Group=df$Group)
head(pca_data)

percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100,2)
percentage <- paste( colnames(pca_data), "(", paste( as.character(percentage), "%", ")", sep="") )

theme<-theme(
  #panel.background = element_rect(fill="white"),
  #legend.position= "none",
  axis.line = element_line(size=1,colour="black"),
  axis.text.x = element_text(size=10,colour="black",face="bold"),
  axis.text.y = element_text(size=10,colour="black",face="bold"))

PCA<-ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
  scale_fill_manual(limits=c("HD","ICU","ICU-Deceased","RD"),values = c("forestgreen","red2","black","cornflowerblue"))+
  geom_point(aes(colour=Group,fill=Group),shape=21,size = 3,color="black") +
  xlab(percentage[1]) + ylab(percentage[2])


PCA<-PCA+theme
PCA
ggsave(PCA,  filename = "CovidFinalfeb/PCA/allfiles_innate229_v2.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")

# clusters
library(ggbiplot)

data<-df3 # your choice
colnames(data) <- gsub(" ","",colnames(data))
colnames(data) <- gsub("(pg/mL)","",colnames(data))
colnames(data) <- gsub("(","",colnames(data),fixed=T)
colnames(data) <- gsub(")","",colnames(data),fixed=T)

df<-data %>% 
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))%>% ungroup 

df<- df%>% select(-c(1,3:6))
df %>% select(-Group)%>%
  scale() %>%                  # scale to 0 mean and unit variance
  prcomp() ->pca 

# factor Groups column
df$Group <- as.factor(df$Group)

cluster<-ggbiplot(pca, obs.scale=1, var.scale=1, groups=df$Group, ellipse=TRUE,var.axes=FALSE )+
  scale_colour_manual(name = "Group",
                      labels = c("HD","ICU","RD"),
                      values = c("forestgreen","red2","cornflowerblue")) +
  geom_point(aes(colour=df$Group), size = 2) +
  
  theme


ggsave(cluster,  filename = "CovidFinalfeb/PCA/plasma_clusters.png",dpi = 300, type = "cairo",  width = 6, height = 4, units = "in")



