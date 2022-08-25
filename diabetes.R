dat<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/IN2_Plasma_Cytokine_Diabetes.csv")

# Data Preperation------------------------- for two sets -----------------
#in2/plsma/cytokinecodeplex
dat<-dat[,c(2:383)]
#cytokineAnnas
dat<-dat[c(1:19),c(385:425)]

# na columns
dat<-dat[,c(2,155,157,174)]


# if NA needs to replace by group means
dff<-dat %>% 
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

# remove all zero columns
zero<-dat[, colSums(dat[,-1]) > 0] # get logical vector to show 
nonzero<-which(zero == TRUE)
dat<-dat %>% select(Group, names(nonzero))


dfin2<-dat[,c(1:347)]
dfplasma<-dat[,c(1,348:355)]
dfcyto<-dat[,c(1,356:378)]

df<-dff
#Look for duplicated columns
duplicated(colnames(df)) # if duplicated names present output "TRUE"


# remove all zero columns
#df %>% select_if(~ !is.numeric(.) || sum(.) != 0)

colnames(df) <- gsub("/","",colnames(df))
colnames(df) <- gsub(" ","",colnames(df))
colnames(df) <- gsub("-","",colnames(df))
colnames(df) <- gsub("(","",colnames(df),fixed=T)
colnames(df) <- gsub(")","",colnames(df),fixed=T)

ttest<-c() # create empty vector
wilcox<-c()
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
#write.csv(ttestresults,'ttest_cytoAnnas_diabetes.csv')

names(wilcox)<-names(df[,2:ncol(df)])
wilcoxresults<-data.frame(wilcox) # convert it to a data-frame
wilcoxresults<-tibble::rownames_to_column(wilcoxresults,var="Marker")
wsig<-wilcoxresults %>% filter(wilcox<0.05)

wilcoxresults<-data.frame(sort(wilcox)) # sort p-values in acsending order
#write.csv(wilcoxresults,'wilcox_cytoAnnas_diabetes.csv')

joinsig<-inner_join(tsig, wsig)
#write.csv(joinsig,'ttestwilcox_cytoAnnas_diabetes.csv')

dfsig<-df %>% select(one_of(dput(as.character(tsig$Marker)))) # select columns match with rows of anovasig dataframe
dfsig<-cbind(df[,1],dfsig)
dfsig #  now can feed dfsig to loop to make boxplots for markers which have significance from ANOVA - Kruskal Wallis

for (i in seq(1:(ncol(dfsig)-1))) {
  df1<-dfsig %>% dplyr::select(c(1, i+1))
  a <- colnames(df1)
  colname1 <-as.character(a[1])
  colname2 <-as.character(a[2])
  box_lsd<-ggboxplot(df1, x = colname1, y = colname2,outlier.shape = NA)+
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(size=1,colour="black"),
          axis.text.x = element_text(size=10,colour="black",face="bold"),
          axis.text.y = element_text(size=10,colour="black",face="bold")) +
    scale_x_discrete(limits=c("ND", "D"),labels=c("ND" = "Non Diabetes","D" = "Diabetes"))+
    geom_jitter(position=position_jitter(width=.30, height=0), aes(colour=get(colname1),fill=get(colname1)),shape=21,size = 7,color="black")+
    scale_fill_manual(values = c("orange","#8B0000"))+
    xlab("") + ylab("")+  theme(legend.position = "none")
  p<-box_lsd+stat_compare_means(aes(label = paste0("", ..p.signif..)),label.x.npc=0.5,method = "t.test")
  print(p)
  ggsave(p,  filename = paste0(colname2, ".png"),dpi = 300, type = "cairo",  width = 4, height = 6, units = "in")
  
}
