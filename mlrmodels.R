# Multiple liner regression models

# 1. Eosinophils,NKTPCCR7P,CD56DIMCD16PNKp44PMFI
library(plot3D)
df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 12.27.2020_Krish.csv")

dfICU<-df[c(17:36),]

data<-dfICU %>% select(Sample,Eosinophils,NKTPCCR7P,CD56DIMCD16PNKp44PMFI) # 20 observations
data<-dfICU %>% select(Sample,Neutrophil,CD8P,CD4P)

data<-data[complete.cases(data), ] # remove 2 rows with NA 

datas<-as.data.frame(scale(data[,-1]))

dataprint<-cbind(data,datas)
write.csv(dataprint,"MLR1.csv")

Corr=Hmisc::rcorr(as.matrix(datas),type='pearson') 
correlatons<-Corr$r
write.csv(correlatons,"MLR1corr.csv")
# normality
ggpubr::ggdensity(data$NKTPCCR7P ,main= "Density plot",xlab = "NKTPCCR7P")
ggpubr::ggdensity(data$CD56DIMCD16PNKp44PMFI ,main= "Density plot",xlab = "CD56DIMCD16PNKp44PMFI")

ggpubr::ggqqplot(data$NKTPCCR7P)
ggpubr::ggqqplot(data$CD56DIMCD16PNKp44PMFI )

model1 <- lm(Eosinophils ~ NKTPCCR7P + CD56DIMCD16PNKp44PMFI, data = datas)
sink(file = "lm_output.txt")
summary(model1)
sink(file = NULL)

par(mfrow = c(2, 2))
plot(model1)
dev.off()

res = resid(model1)
plot(datas$CD56DIMCD16PNKp44PMFI, res,ylab="Residuals", xlab="CD56DIMCD16PNKp44PMFI") 
abline(0, 0) 
plot(datas$NKTPCCR7P, res,ylab="Residuals", xlab="NKTPCCR7P") 
abline(0, 0) 

## qq plot
stdres = rstandard(model1)
ggpubr::ggqqplot(stdres)


model2 <- lm(NKTPCCR7P ~ Eosinophils + CD56DIMCD16PNKp44PMFI, data = data) # not good
summary(model2)

model3 <- lm(CD56DIMCD16PNKp44PMFI ~ Eosinophils + NKTPCCR7P, data = data) # not good
summary(model3)


model4<-lm(Eosinophils ~ CD56DIMCD16PNKp44PMFI+ NKTPCCR7P+I(CD56DIMCD16PNKp44PMFI*NKTPCCR7P),data = datas)
summary(model4)
plot(model4)
# x, y, z variables for model1
x <- datas$NKTPCCR7P
y <- datas$CD56DIMCD16PNKp44PMFI
z <- datas$Eosinophils
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
scatter3D(x, y, z, pch = 18, cex = 1,col = RColorBrewer::brewer.pal(9, "RdYlGn"),
          bty ="g",
          colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
          cex.lab=0.7,cex.axis=0.3,
          theta = 40, phi = 30, ticktype = "detailed",
          xlab = "NKTPCCR7P", ylab = "CD56DIMCD16PNKp44PMFI", zlab = "Eosinophils")#,  
          #surf = list(x = x.pred, y = y.pred, z = z.pred, 
                      #facets = NA, fit = fitpoints))
# 2. eosinophils, FABP4 and MCP.1
data <-data.table::fread("/Users/krishanthi/twocytokine_innate_plasma_combined_Krish.csv")
colnames(data) <- gsub("(ng/mL)","",colnames(data))
colnames(data) <- gsub("(pg/mL)","",colnames(data))
colnames(data) <- gsub("(u/mL)","",colnames(data))
colnames(data) <- gsub(" ","",colnames(data))
colnames(data) <- gsub("-","",colnames(data))
colnames(data) <- gsub("(","",colnames(data),fixed=T)
colnames(data) <- gsub(")","",colnames(data),fixed=T)

df2<-data

df2RD<-df2[c(32:55),] # Recovered patient data

data2<-df2RD%>% select(Sample,Eosinophils,MCP1,FABP4)
data2<-data2[complete.cases(data2), ] # remove 2 rows with NA 

data2s<-as.data.frame(scale(data2[,-1]))

dataprint2<-cbind(data2,data2s)
write.csv(dataprint2,"MLR2.csv")

model1 <- lm(Eosinophils ~ MCP1 + FABP4, data = data2s)
sink(file = "lm2_output.txt")
summary(model1)
sink(file = NULL)

par(mfrow = c(2, 2))
plot(model1)
dev.off()

model2 <- lm(MCP1 ~ Eosinophils + FABP4, data = data2s) # not good
summary(model2)
plot(model1)

model3 <- lm(FABP4 ~ Eosinophils + MCP1, data = data2s) # not good
summary(model3)

x <- data2s$MCP1
y <- data2s$FABP4
z <- data2s$Eosinophils

# 3. CD4P, CD8P Neutrophil

df <-data.table::fread("/Users/krishanthi/Documents/Liyanage_lab/Liyanage_data/Covid-3/IN2-final/Krish_made/IN2 Color Coded- 12.27.2020_Krish.csv")

dfICU<-df[c(17:36),]

data<-dfICU %>% select(Sample,Neutrophil,CD8P,CD4P,CD56DIMCD16P)
data<-data[complete.cases(data), ] # remove 2 rows with NA 
datas<-as.data.frame(scale(data[,-1]))
model1 <- lm(Neutrophil ~ CD4P+CD8P, data = data)
summary(model1)
model2 <- lm(CD56DIMCD16P ~ CD4P+CD8P, data = data)
summary(model2)

par(mfrow = c(2, 2))
plot(model1)
dev.off()
x <- data$CD4P
y <- data$CD8P
z <- data$Neutrophil
z<- data$CD56DIMCD16P

png("regression_CD56DIMCD16P.png", width = 1600, height = 1200, units = "px")

plot3D::scatter3D(x, y, z, pch = 19, cex = 3,col =topo.colors(100),#RColorBrewer::brewer.pal(9, "YlGnBu"),
                  bty ="b2", alpha=1,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 2),
                  cex.lab=2,cex.axis=2,
                  theta = 40, phi = 30, ticktype = "detailed",
                  xlab = "CD4P", ylab = "CD8P", zlab = "CD56DIMCD16P",
                  surf = list(x = x.pred, y = y.pred, z = z.pred, 
                             col =topo.colors(100),alpha=1,
                             facets = NA, fit = fitpoints))
dev.off()

library(reshape2) ## for melt()
dl  <- melt(data[,c(2,3,4)],id.var=1)

dcplot<-ggpubr::ggscatter(dl, y = "Neutrophil", x = "value",
                  add = "reg.line",conf.int = FALSE,
                  ylab="Neutrophil",xlab="",color = "variable")+
  stat_cor(aes(color = variable), label.x = 12)+
  
  theme(axis.text.x=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=10))
ggsave(dcplot, filename = "Neutrophil_multi_corr_plots.png",dpi = 300, type = "cairo",  width = 4, height = 4, units = "in")
dcplot

