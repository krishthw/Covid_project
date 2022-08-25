df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Manuja_made2/Plasma ELISA final - February 10th 2021.csv")
tags <- c("<30","[30-50)", "[50-70)","[70-100)") # define age brackets
tags1<- c("<50","[50-100)")
cdf <- as_tibble(df_innate) %>% 
  mutate(AgeGroup = case_when(
    Age < 30 ~ tags[1],
    Age  >= 30 & Age  < 50 ~ tags[2],
    Age  >= 50 & Age  < 70 ~ tags[3],
    Age  >= 70 & Age  < 100 ~ tags[4])) %>%
  mutate(AgeGroup1 = case_when(
    Age < 50 ~ tags1[1],
    Age  >= 50 & Age  < 100 ~ tags1[2]))


df<-cdf[,-c(1,3:6,(ncol(cdf)-1),ncol(cdf))] %>%
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)) %>%
  ungroup()


df <-df[,-1 ]%>% 
  scale() %>%                  # scale to 0 mean and unit variance
  prcomp() ->pca 

# display the results from the PCA 
head(pca$x)
pca_data <- data.frame(pca$x, Group=cdf$AgeGroup1)
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
  scale_fill_manual(values = c("pink","royalblue"))+
  
  #scale_fill_manual(values = c("#FFA500","#FF7F50","#FF4500","#FFD700"))+
  geom_point(aes(colour=Group,fill=Group),shape=21,size = 3,color="lightgray") +
  xlab(percentage[1]) + ylab(percentage[2])+
  theme
PCA
library("gg3D")

PCA3<-ggplot(pca_data, aes(x=PC1, y=PC2,z=PC3 ,color=Group)) + 
  theme_void() +
  axes_3D() +
  stat_3D()

PCA3
library(plotly)
plot_ly(x=pca_data$PC1, y=pca_data$PC2,z=pca_data$PC3, type="scatter3d", mode="markers", color=df_innate$Age)

plot_ly(x=pca_data$PC1, y=pca_data$PC2,z=pca_data$PC3, type="scatter3d", mode="markers", color=pca_data$Group)

