library(dplyr)
library(data.table)
library(broom)
library(tibble)
library(ggplot2)
library(DescTools)
library(rstatix)
library(ggpubr)
# data file
# Impoart the dataset (.csv)
df1<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Cytokines Anna's data final- February 10th 2021.csv")
df2<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Isoplexis Codeplex final- February 10th 2021.csv")
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
df5<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/AD1 (PBMC) final- April 29th 2021.csv")

# select one to analyse
dff<-df4[df4$Diabetes != ""]    # innate panel
dff<-df1[df1$Diabetes != ""]    # Cytokine anna
dff<-df2[df2$Diabetes != ""]    # Cytokine codeplex 
dff<-df3[df3$Diabetes != ""]    # Plasma
dff<-df5[df5$Diabetes != ""]    # PBMC


# __________________________________________________________________________________________
# Data Preperation 

#Look for duplicated columns
length(colnames(dff))-length(unique(colnames(dff))) # number of duplicates
duplicated(colnames(dff)) # if duplicated names present output "TRUE"

# remove unnecessary characters
colnames(dff) <- gsub("(ng/mL)","",colnames(dff))
colnames(dff) <- gsub("(pg/mL)","",colnames(dff))
colnames(dff) <- gsub("(u/mL)","",colnames(dff))
colnames(dff) <- gsub(" ","",colnames(dff))
colnames(dff) <- gsub("-","",colnames(dff))
colnames(dff) <- gsub("(","",colnames(dff),fixed=T)
colnames(dff) <- gsub(")","",colnames(dff),fixed=T)
colnames(dff) <- gsub("_","",colnames(dff))

# remove all zero columns
no_zero<-dff[, colSums(dff[,-c(1:6)]) > 0] # get logical vector to show 
nonzero<-which(no_zero == TRUE)     # select non zero columns

# important
# if columns with NA present, exsecute *** to get statistics and boxplots only. Not for correlations
names(which(colSums(is.na(dff))>0))
df<-dff %>% select(Diabetes, names(nonzero)) # Group and non-zero columns

# ***
df<-dff %>% select(Diabetes, names(nonzero), names(which(colSums(is.na(dff))>0))) # non-zero columns and columns with NA

# __________________________________________________________________________________________
# ANOVA and Kruskal -Wallis tests to identify  significant differences in means within 3 groups

anova<-c()    # empty vector to store anova test 
kruskal<-c()  # empty vector to store kruskal wallis test
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
# get Kruskal test to a dataframe
names(kruskal)<-names(df[,2:ncol(df)]) # assign each p-value to marker type. this results a numeric vector
kruskaltest <-data.frame(kruskal) # convert to a data-frame
kruskaltest<-tibble::rownames_to_column(kruskaltest,var="Marker")
kruskalsig<-kruskaltest %>% filter(kruskal<0.05)
kruskaltest<-sort(kruskal)
write.csv(kruskaltest,'CovidFinalfeb/Diabetes/3Groups/Plasma/diabetes_plasma_kruskal_report.csv')
# get Anova test to a dataframe
names(anova)<-names(df[,2:ncol(df)]) # assign each p-value to marker type. this results a numeric vector
ANOVAtest <-data.frame(anova) # convert to a data-frame
ANOVAtest<-tibble::rownames_to_column(ANOVAtest,var="Marker")
ANOVAsig<-ANOVAtest %>% filter(anova<0.05)
ANOVAtest<-sort(anova)
write.csv(ANOVAtest,'CovidFinalfeb/Diabetes/3Groups/Plasma/diabetes_plasma_anova_report.csv')
#Filter out Markers significant from both tests
joinsig<-inner_join(kruskalsig, ANOVAsig)
write.csv(joinsig,'CovidFinalfeb/Diabetes/3Groups/Plasma/diabetes_plasma_anovakruskal_report.csv')

dfsig<-df %>% select(one_of(dput(as.character(joinsig$Marker)))) 
dfsig<-cbind(df[,1],dfsig)
dfsig 
#  now can feed dfsig to loop to make boxplots for markers which have significance from ANOVA - Kruskal Wallis

# __________________________________________________________________________________________
# TukeyHSD and FisherLSD tests without any adjustment methods to identify the groups which have different means.

#TukeyHSD and FisherLSD test results
tukeyLSD_dfsig=NULL;
for (i in seq(1:(ncol(dfsig)-1))){
  dfsig1<-dfsig%>% dplyr::select(c(1, i+1))
  a<-colnames(dfsig1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  tukey<-tukey_hsd(aov(get(colname2)~get(colname1),dfsig1))
  lsd<-(PostHocTest(aov(get(colname2)~get(colname1),data=dfsig1), method = "lsd"))$`get(colname1)`
  pairsTL<-list("Tukey-G1-G2","Tukey-G1-G3","Tukey-G2-G3","fisher-G1-G2","Fisher-G1-G3","Fisher-G2-G3")
  valuesTL<-list(as.numeric(tukey[1,7]),as.numeric(tukey[2,7]),as.numeric(tukey[3,7]),lsd[1,4],lsd[2,4],lsd[3,4])
  TL<-cbind(pairsTL,valuesTL) # combine pairs and values
  TLt<-data.frame(t(TL)) # make a dataframe after transposing
  tt<-data.frame(X1 = TLt[2,1], X2= TLt[2,2], X3=TLt[2,3],X4=TLt[2,4],X5=TLt[2,5],X6=TLt[2,6])
  tukeyLSD_dfsig<-rbind(tukeyLSD_dfsig,tt)
}

tukeyLSD_dfsig<-cbind(names(dfsig[,-1]),tukeyLSD_dfsig)
colnames(tukeyLSD_dfsig)<-c("marker","Tukey-G1-G2","Tukey-G1-G3","Tukey-G2-G3","fisher-G1-G2","Fisher-G1-G3","Fisher-G2-G3")
write.csv(tukeyLSD_dfsig,'CovidFinalfeb/Diabetes/3Groups/Plasma/diabetes_plasma_TukeyFisher_report.csv')

# __________________________________________________________________________________________
# Boxplots with FisherLSD test significance for diabetes
dead<-dff[dff$State == 'Deceased']
live<-dff[dff$State != 'Deceased']

for (i in seq(1:(ncol(dfsig)-1))) {
  df1<-dfsig %>% dplyr::select(c(1, i+1))
  a <- colnames(df1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  
  tukey_p<-aov(as.formula(paste(a[2], a[1], sep="~")),df1) %>%
    tukey_hsd() %>%
    add_xy_position(x="Species",data=df1,fun = "max",
                    step.increase = 0.12,formula=as.formula(paste(a[2], a[1], sep="~")))
  lsd_test<-PostHocTest(aov(as.formula(paste(a[2], a[1], sep="~")),df1), method = "lsd")$`Diabetes`
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
    scale_x_discrete(limits=c("CTRL", "No","Yes"),labels=c("CTRL" = "CTRL","No" = "Non Diabetes","Yes"="Diabetes"))+
    geom_boxplot(outlier.shape = NA)+
    stat_pvalue_manual(lsd_p, hide.ns = T,size=6)+
    geom_point(data=live, position=position_jitter(seed=1),
               aes(fill=get(colname1)),shape=21,size = 5)+
    geom_point(data=dead,position=position_jitter(seed=1), 
               shape=19,size = 3,color="black")+
    #geom_jitter(position=position_jitter(width=.30, height=0), aes(colour=get(colname1),fill=get(colname1)),shape=21,size = 7,color="black")+
    scale_fill_manual(values = c("forestgreen","orange","#8B0000"))+
    xlab("") + ylab("")+  theme(legend.position = "none")
  
  print(box_lsd)
  ggsave(box_lsd,  filename = paste0("CovidFinalfeb/Diabetes/3Groups/Plasma/", colname2, "_diabetes_3G.png"),dpi = 300, type = "cairo",  width = 4, height = 6, units = "in")
}

dev.off()


# _______________________________________________________
#               for two groups
#________________________________________________________
# select one to analyse
df<-df[df$Diabetes != "CTRL"]  
df$Diabetes[df$Diabetes == "No"]<- "Non Diabetes"
df$Diabetes[df$Diabetes == "Yes"] <- "Diabetes"

ttest<-c() 
wilcox<-c() # wilcox can not calculate EXACT p, when there is the same value 
for (i in seq(1:(ncol(df)-1))){ 
  df2<-df %>% dplyr::select(c(1, i+1))
  a<-colnames(df2)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  TT<-t.test(get(colname2) ~ get(colname1), data = df2)
  WIL<-wilcox.test(get(colname2) ~ get(colname1), data = df2)
  ttest<-c(ttest,TT$"p.value")
  wilcox<-c(wilcox,WIL$"p.value")
}
names(ttest)<-names(df[,2:ncol(df)])
ttestresults<-data.frame(ttest) # convert it to a data-frame
ttestresults<-tibble::rownames_to_column(ttestresults,var="Marker")
tsig<-ttestresults %>% filter(ttest<0.05)
ttestresults<-data.frame(sort(ttest)) # sort p-values in acsending order
write.csv(ttestresults,'CovidFinal/Diabetes/Plasma/2Groups/Plasmattest.csv')

names(wilcox)<-names(df[,2:ncol(df)])
wilcoxresults<-data.frame(wilcox) # convert it to a data-frame
wilcoxresults<-tibble::rownames_to_column(wilcoxresults,var="Marker")
wsig<-wilcoxresults %>% filter(wilcox<0.05)
wilcoxresults<-data.frame(sort(wilcox)) # sort p-values in acsending order
write.csv(wilcoxresults,'CovidFinal/Diabetes/Plasma/2Groups/Plasmawilcox.csv')

joinsig<-inner_join(tsig, wsig)
write.csv(joinsig,'CovidFinal/Diabetes/Plasma/2Groups/Plasmattestwilcox.csv')

dfsig<-df %>% select(one_of(dput(as.character(joinsig$Marker)))) 
dfsig<-cbind(df[,1],dfsig)
dfsig #  now can feed dfsig to loop to make boxplots for markers which have significance from ANOVA - Kruskal Wallis

for (i in seq(1:(ncol(dfsig)-1))) {
  df1<-dfsig %>% dplyr::select(c(1, i+1))
  a <- colnames(df1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  box_lsd<-ggboxplot(df1, x = colname1, y = colname2,outlier.shape = NA)+
    stat_compare_means(aes(label = ..p.signif..),label.x.npc=0.5,method = "t.test",size = 8)+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=1,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold")) +
    geom_jitter(position=position_jitter(width=.30, height=0), aes(colour=get(colname1),fill=get(colname1)),shape=21,size = 7,color="black")+
    scale_fill_manual(values = c("orange","#8B0000"))+
    xlab("") + ylab("")+  theme(legend.position = "none")
  
    print(box_lsd)
    ggsave(box_lsd,  filename = paste0("CovidFinal/Diabetes/Plasma/2Groups/", colname2, "_diabetes_2G.png"),dpi = 300, type = "cairo",  width = 4, height = 6, units = "in")
  
}

