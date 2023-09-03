setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 5")

#Load packages
library(tidyverse)
library(car)
library(statmod)

#Load data
cholla<-read.csv("cholla_demography_corrected.csv") 

#creates reproduction variable that is binary
cholla$Reproduced <- ifelse(cholla$Goodbuds_t>0, 1, 0)

#plotting reproduction vs. size
plot(x=cholla$logvol_t, y=cholla$Reproduced, 
     xlab = "Volume in year t", ylab = "Reproduced in t+1")


rep_bin <- cholla %>% 
  mutate(rep_bin = cut_interval(logvol_t, n=10)) %>% #Makes 8 size bins based on log(vol_t)
  group_by(rep_bin) %>% 
  summarise(mean_size = mean(logvol_t,na.rm=T), #Mean size in each bin
            mean_rep = mean(Reproduced,na.rm=T)) #Mean reproduction in each bin

plot(x=cholla$logvol_t, y=cholla$Reproduced, 
     xlab = "Volume in year t", ylab = "Reproduced in t+1", col = "grey")
points(rep_bin$mean_size, rep_bin$mean_rep, pch=16, cex=2)

#more plots

invlogit <- function(x){exp(x)/(1+exp(x))}

x_dummy <- seq(min(cholla$logvol_t,na.rm=T),max(cholla$logvol_t,na.rm=T),0.1)
#add predicted line using dummy x
lines(x_dummy,invlogit(coef(reprod_model)[1]+coef(reprod_model)[2]*x_dummy),col="blue",lwd=3) 



#glm 
reprod_model <- glm(Reproduced ~ logvol_t, family = "binomial", data=cholla)

plot(reprod_model)

qres<-qresid(reprod_model)#calculates quantilized residuals with some introduced randomization to overcome discreteness
qqnorm(qres)#make qqplot
abline(0,1)


resid.ssq <- sum(residuals(reprod_model,type="pearson")^2)  
resid.df <- nrow(subset(cholla,!is.na(logvol_t) & !is.na(Reproduced)))-length(coef(reprod_model)) 
resid.ssq/resid.df #overdispersion!!!s


summary(reprod_model)
scatter.smooth(qres~fitted(reprod_model))
