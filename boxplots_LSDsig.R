# 1.statistical analysis report of one-way ANOVA and Kruskal Wallis test 
#   to determine whether there are any statistically significant differences between the means 
#   of three or more independent (unrelated) groups. 
# 2.statistical analysis report of FisherLSD(unprotected) test follow up to ANOVA
# 3. Box-plots nd Violin plots with jitter
# -------------------------------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(broom)
library(tibble)
library(ggplot2)
library(DescTools)
library(rstatix)
library(ggpubr)

# Impoart the dataset (.csv)
df_innate<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 01.27.2021_Krish.csv")
df_anna <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Cytokines Anna's data Updated 12.25.2020_Krish.csv")
df_codeplex<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Isoplexis Codeplex- 01.27.2021_Krish.csv")
df_plasma <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/Plasma ELISA - 01.27.2021_Krish.csv")
df_pbmc <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/AD1 Color Coded(PBMC)- December 23rd  2020.csv")
df_wb <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/AD1 Color Coded(WB)- December 23rd  2020.csv")

# select one to analyse
df<-df_innate[,c(2,7:356)]  # innate panel
df<-df_anna[,c(2,10:49)]    # Cytokine anna
df<-df_codeplex[,c(2,7:29)] # Cytokine codeplex 
df<-df_plasma[,c(2,8:15)]   # Plasma
df<-df_pbmc[,c(2,10:95)]    # PBMC
df<-df_wb[,c(2,10:93)]      # Whole Blood

# __________________________________________________________________________________________
# Data Preperation 

#Look for duplicated columns
length(colnames(df))-length(unique(colnames(df))) # number of duplicates
duplicated(colnames(df)) # if duplicated names present output "TRUE"

# remove unnecessary characters
colnames(df) <- gsub("(ng/mL)","",colnames(df))
colnames(df) <- gsub("(pg/mL)","",colnames(df))
colnames(df) <- gsub("(u/mL)","",colnames(df))
colnames(df) <- gsub(" ","",colnames(df))
colnames(df) <- gsub("-","",colnames(df))
colnames(df) <- gsub("(","",colnames(df),fixed=T)
colnames(df) <- gsub(")","",colnames(df),fixed=T)
colnames(df) <- gsub("_","",colnames(df))

# remove all zero columns
no_zero<-df[, colSums(df[,-1]) > 0] # get logical vector to show 
nonzero<-which(no_zero == TRUE)     # select non zero columns
# important
# if columns with NA present, exsecute *** to get statistics and boxplots only. Not for correlations
names(which(colSums(is.na(df))>0))
df<-df %>% select(Group, names(nonzero)) # Group and non-zero columns

# ***
df<-df %>% select(Group, names(nonzero), names(which(colSums(is.na(df))>0))) # non-zero columns and columns with NA

# __________________________________________________________________________________________
# ANOVA and Kruskal -Wallis tests to identify  significant differences in means within 3 groups

anova<-c()   # create empty vector
kruskal<-c() # create empty vector
for (i in seq(1:(ncol(df)-1))){ # loop through the columns and perform the test
  df2<-df %>% dplyr::select(c(1, i+1))
  a<-colnames(df2)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  KW<-kruskal.test(get(colname2) ~ get(colname1), data = df2) 
  kruskal<-c(kruskal,KW$"p.value")
  AN<-broom::tidy(aov(get(colname2) ~ get(colname1), data = df2))
  anova<-c(anova,as.numeric(AN[1,6]))
}
names(kruskal)<-names(df[,2:ncol(df)]) # assign each p-value to marker type. this results a numeric vector
kruskaltest <-data.frame(kruskal) # convert to a data-frame
kruskaltest<-tibble::rownames_to_column(kruskaltest,var="Marker")
kruskalsig<-kruskaltest %>% filter(kruskal<0.05)
kruskaltest<-sort(kruskal)
write.csv(kruskaltest,'no17/PBMC/pbmc_kruskal.csv')

names(anova)<-names(df[,2:ncol(df)])
ANOVAtest <-data.frame(anova) # convert to a data-frame
ANOVAtest<-tibble::rownames_to_column(ANOVAtest,var="Marker")
ANOVAsig<-ANOVAtest %>% filter(anova<0.05)
ANOVAtest<-sort(anova)
write.csv(ANOVAtest,'no17/PBMC/pbmc_anova.csv')

joinsig<-inner_join(kruskalsig, ANOVAsig)
write.csv(joinsig,'no17/PBMC/pbmc_anovakruskal.csv')

kruskal_only<-kruskalsig %>% anti_join(ANOVAsig)
anova_only<-ANOVAsig %>% anti_join(kruskalsig)

# filter markers significant from both tests
dfsig<-df %>% select(one_of(dput(as.character(anova_only$Marker)))) # select columns match with rows of joinsig dataframe
dfsig<-cbind(df[,1],dfsig)
dfsig #  now can feed dfsig to loop to make boxplots for markers which have significance from ANOVA - Kruskal Wallis

# __________________________________________________________________________________________
# TukeyHSD and FisherLSD tests without any adjustment methods to identify the groups which have different means.

tukeyLSD_dfsig=NULL;
for (i in seq(1:(ncol(dfsig)-1))){
  dfsig1<-dfsig%>% dplyr::select(c(1, i+1))
  a<-colnames(dfsig1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  tukey<-tukey_hsd(aov(get(colname2)~get(colname1),dfsig1))
  fisher<-(PostHocTest(aov(get(colname2)~get(colname1),data=dfsig1), method = "lsd"))$`get(colname1)`
  pairsTL<-list("Tukey-HD_ICU","Tukey-HD_RD","Tukey-ICU_RD",
                "Fisher-HD_ICU","Fisher-HD_RD","Fisher-ICU_RD")
  valuesTL<-list(as.numeric(tukey[1,7]),as.numeric(tukey[2,7]),as.numeric(tukey[3,7]),
                 fisher[1,4],fisher[2,4],fisher[3,4])
  TL<-cbind(pairsTL,valuesTL) # combine pairs and values
  TLt<-data.frame(t(TL)) # make a dataframe after transposing
  tt<-data.frame(X1 = TLt[2,1], X2= TLt[2,2], X3=TLt[2,3],X4=TLt[2,4],X5=TLt[2,5],X6=TLt[2,6])
  tukeyLSD_dfsig<-rbind(tukeyLSD_dfsig,tt)
}

tukeyLSD_dfsig<-cbind(names(dfsig[,-1]),tukeyLSD_dfsig)
colnames(tukeyLSD_dfsig)<-c("marker","Tukey-HD_ICU","Tukey-HD_RD","Tukey-ICU_RD",
                            "Fisher-HD_ICU","Fisher-HD_RD","Fisher-ICU_RD")
write.csv(tukeyLSD_dfsig,'no17/PBMC/pbmc_TukeyFisher.csv')

# __________________________________________________________________________________________
#Boxplot with FisherLSD test significance

for (i in seq(1:(ncol(dfsig)-1))) {
  df1<-dfsig %>% dplyr::select(c(1, i+1))
  a <- colnames(df1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  tukey_p<-aov(as.formula(paste(a[2], a[1], sep="~")),df1) %>%
    tukey_hsd() %>%
    add_xy_position(x="Species",data=df1,fun = "max",
                    step.increase = 0.12,formula=as.formula(paste(a[2], a[1], sep="~")))
  lsd_test<-PostHocTest(aov(as.formula(paste(a[2], a[1], sep="~")),df1), method = "lsd")$`Group`
  lsd_test<-as.data.frame(lsd_test)
  lsd_p<-cbind(lsd_test,as.data.frame(tukey_p%>%select(2,3,9,10,11,12)))
  lsd_p<-lsd_p%>%mutate(p.adj.signif = ifelse(pval<= 0.0001,paste0("****"), 
                                              ifelse(pval <= 0.001,paste0("***"), 
                                                     ifelse(pval <= 0.01,paste0("**"),
                                                            ifelse(pval <= 0.05,paste0("*"),
                                                                   "ns")))))
  print(i)
  print(lsd_p)
  box_lsd<-ggplot(df1, aes(x = get(colname1), y = get(colname2 )))+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=1,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold")) +
    scale_x_discrete(limits=c("HD", "ICU","RD"),labels=c("HD" = "HD","ICU" = "ICU","RD"="RD"))+
    geom_boxplot(outlier.shape = NA)+
    stat_pvalue_manual(lsd_p, hide.ns = T,size=6)+
    geom_jitter(position=position_jitter(width=.30, height=0), 
                aes(colour=get(colname1),fill=get(colname1)),shape=21,size = 7,color="black")+
    scale_fill_manual(values = c("forestgreen","red2","cornflowerblue"))+
    xlab("") + ylab("")+  theme(legend.position = "none")
  print(box_lsd)
  ggsave(box_lsd,  filename = paste0("no17/PBMC/pbmc_boxplots/extras/", colname2, "_anova.png"),
         dpi = 300, type = "cairo",  width = 4, height = 6, units = "in")
}
dev.off()

# __________________________________________________________________________________________
# Violin plot with FisherLSD test significance

for (i in seq(1:(ncol(dfsig)-1))) {
  df1<-dfsig %>% dplyr::select(c(1, i+1))
  a <- colnames(df1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  tukey_p<-aov(as.formula(paste(a[2], a[1], sep="~")),df1) %>%
    tukey_hsd() %>%
    add_xy_position(x="Species",data=df1,fun = "max",
                    step.increase = 0.12,formula=as.formula(paste(a[2], a[1], sep="~")))
  lsd_test<-PostHocTest(aov(as.formula(paste(a[2], a[1], sep="~")),df1), method = "lsd")$`Group`
  lsd_test<-as.data.frame(lsd_test)
  lsd_p<-cbind(lsd_test,as.data.frame(tukey_p%>%select(2,3,9,10,11,12)))
  lsd_p<-lsd_p%>%mutate(p.adj.signif = ifelse(pval<= 0.0001,paste0("****"), 
                                              ifelse(pval <= 0.001,paste0("***"), 
                                                     ifelse(pval <= 0.01,paste0("**"),
                                                            ifelse(pval <= 0.05,paste0("*"),
                                                                   "ns")))))
  print(i)
  print(lsd_p)
  box_lsd<-ggplot(df1, aes(x = get(colname1), y = get(colname2 )))+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=1,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold")) +
    scale_x_discrete(limits=c("HD", "ICU","RD"),labels=c("HD" = "HD","ICU" = "ICU","RD"="RD"))+
    stat_pvalue_manual(lsd_p, hide.ns = T,size=6,bracket.nudge.y=0.8,step.increase=0.2)+
    geom_violin(aes(colour=get(colname1),fill=get(colname1)),trim=FALSE)+
    scale_color_manual(values = c("forestgreen","red2","cornflowerblue"))+
    
    geom_jitter(position=position_jitter(width=.1, height=0),shape=21,size=1,color="black")+
    scale_fill_manual(values = c("forestgreen","red2","cornflowerblue"))+
    
    xlab("") + ylab("")+  theme(legend.position = "none")
  b<-box_lsd+ stat_summary(fun.data=mean_sdl,geom="pointrange", color="black",size=0.1,mult=1)
  #b<-box_lsd+stat_summary(fun=mean, geom="point", shape="-",size=10,color='black')
  print(b)
  ggsave(b,  filename = paste0("no17/plasma/plasma_boxplots/", colname2, ".png"),
         dpi = 300, type = "cairo",  width = 4, height = 6, units = "in")
}

