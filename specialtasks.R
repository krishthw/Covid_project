# special tasks

#replace na with group means
dff<-df %>% 
  group_by(Group) %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

#get group means
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Isoplexis Codeplex- December 23rd  2020.csv")
df2<-aggregate(df[, 2:24], list(df$`Donor Group`), mean)
df_num <- as.matrix(df2[,2:24])
rownames(df_num) <- df2$`Group.1`

  
  
  
  
  

  