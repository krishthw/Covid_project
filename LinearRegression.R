library(plot3D)
library(ggplot2)

#Linear regression model 
df3<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/Plasma ELISA final - February 18th 2021.csv")
df3<-df3 %>% select(c("Sample","Group","LBP","FABP4"))
df4<-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/final_corrected/IN2_final_corrected_dhanuja_February17th2021.csv")
df4<-df4 %>% select(c("Sample" ,"Group",'ITM'))

(anti_join(df3, df4, by="Sample"))$Sample 
(anti_join(df4, df3, by="Sample"))$Sample 

# plasma innate combine
dfff<-inner_join(df3,df4, by="Sample") # innerjoin by 'Sample'
df3df4<- dfff %>% select(-ends_with(".y")) # remove duplicates 
colnames(df3df4) <- gsub(".x","",colnames(df3df4)) # remove .x from column names

data<-as.data.frame(df3df4) # convert to dataframe
ICUdata<-data[data[, "Group"] == 'ICU', -c(1:2)] # select ICU set


Corr=Hmisc::rcorr(as.matrix(ICUdata),type='pearson') 
correlations<-Corr$r
P<-Corr$P
write.csv(correlations,"ITM_LBP_FABP4corr.csv")
write.csv(P,"ITM_LBP_FABP4_pvalues.csv")


#ICUdata_s<-as.data.frame(scale(ICUdata[,-c(1:2)])) #scaling is not necessary for linear regression

lrmodel<-lm(ITM ~ LBP + FABP4, data = ICUdata)
sink(file = "ITM_LBP_FABP4_regression_output.txt")
summary(lrmodel)
sink(file = NULL)

plot(lrmodel)

# x, y, z variables for regression model
x <- ICUdata$LBP
y <- ICUdata$FABP4
z <- ICUdata$ITM
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
png("regression_ITM_LBP_FABP4_plot1.png", width = 1600, height = 1200, units = "px")


scatter3D(x, y, z, 
          pch = 19,  cex = 3,colvar = NULL, col="red", # point shape, size and color
          theta = 20, phi = 10,  # oriantation of the grid
          bty="b", # background grid
          xlab = "LBP", ylab = "FABP4", zlab = "ITM", 
          #colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
          cex.lab=3, # axis label font size
          cex.axis=1.5, # axis tick marks font size
          ticktype = "detailed", 
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      col=ramp.col (col = c("dodgerblue3","seagreen2"), n = 300, alpha=0.7),
                      facets = TRUE, fit = fitpoints, border="grey42"))
#scatter3D(x = 12806.631, y = 22991.485, z = 5.12, add = TRUE, colkey =FALSE, 
         # pch = 19, cex = 1, col = "black")


dev.off()

png("regression_ITM_LBP_FABP4_plot2.png", width = 1600, height = 1200, units = "px")

scatter3D(x, y, z, pch = 19, cex = 3,col = RColorBrewer::brewer.pal(4, "RdBu"),
          bty ="b",
          #colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
          cex.lab=3,cex.axis=2,
          theta = 20, phi = 10, 
          #ticktype = "detailed",
          xlab = "LBP", ylab = "FABP4", zlab = "ITM", 
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = TRUE, fit = fitpoints,alpha=0.6,border="grey42"))
dev.off()



png("regression_ITM_LBP_FABP4_plot3.png", width = 1600, height = 1200, units = "px")

scatter3D(x, y, z, pch = 19, cex = 3,col = RColorBrewer::brewer.pal(9, "YlGnBu"),
          bty ="b",
          colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
          cex.lab=3,cex.axis=2,
          theta = 20, phi = 10, 
          ticktype = "detailed",
          xlab = "LBP", ylab = "FABP4", zlab = "ITM", 
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA, fit = fitpoints))
dev.off()

png("regression_ITM_LBP_FABP4_plot4.png", width = 1600, height = 1200, units = "px")

plot3D::scatter3D(x, y, z, pch = 19, cex = 3,col =topo.colors(100),#RColorBrewer::brewer.pal(9, "YlGnBu"),
                  bty ="b2", alpha=1,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 2),
                  cex.lab=2,cex.axis=2,
                  theta = 40, phi = 30, ticktype = "detailed",
                  xlab = "LBP", ylab = "FABP4", zlab = "ITM", 
                  surf = list(x = x.pred, y = y.pred, z = z.pred, 
                              col =topo.colors(100),alpha=1,
                              facets = NA, fit = fitpoints))
dev.off()

png("regression_ITM_LBP_FABP4_plot5.png", width = 1600, height = 1200, units = "px")

scatter3D(x, y, z, pch = 19, cex = 3,col = RColorBrewer::brewer.pal(4, "RdBu"),
          bty ="b",
          #colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
          cex.lab=3,cex.axis=2,
          theta = 20, phi = 10, 
          #ticktype = "detailed",
          xlab = "LBP", ylab = "FABP4", zlab = "ITM", 
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA, fit = fitpoints,border="grey42"))
dev.off()
          