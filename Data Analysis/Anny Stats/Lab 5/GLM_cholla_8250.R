#Cholla demography as GLM example
#Adapted from TEX Miller 2018
#AC 210131

#Housekeeping----------------
#Make sure clean environment

#Set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 5")

#Load packages
library(tidyverse)
library(car)
library(statmod)

#Load data
cholla<-read.csv("cholla_demography_corrected.csv")

#1.1 Growth: prelim visualization--------------------
#Note we will work directly with log(volume) instead of raw data

#Plot log(vol) change from year to year
plot(x=cholla$logvol_t, y=cholla$logvol_t1, 
     xlab = "Volume in year t", ylab = "Volume in year t+1")
#Looks linear

#Do they grow or shrink? Add 1:1 line
abline(a=0, b=1) #a=intercept, b=slope
legend("topleft",bty="n",lty=1,legend=c("y=x"))
#Both, but it seems like there's more growth when plants are small

#1.2 Growth: Fit model---------------------------
#State model
growth_model <- lm(logvol_t1 ~ logvol_t, data = cholla)

#Check assumptions
qqPlot(growth_model$residuals) #looks bad?
hist(growth_model$residuals) #actually not terrible
plot(growth_model,which=1) #good

#Fitted model results
summary(growth_model)

#1.3 Growth: Plot fitted model---------------------
plot(x=cholla$logvol_t, y=cholla$logvol_t1, 
     xlab = "Volume in year t", ylab = "Volume in year t+1")
abline(coef(growth_model)[1],coef(growth_model)[2],col="red",lwd=3) #Fitted line
abline(a=0, b=1) 
legend("topleft",bty="n",lty=1,col=c("black","red"),lwd=c(1,3),legend=c("y=x","fitted model"))

#2.1 Survival: Prelim visualization--------------------
#Try simple x-y plot
plot(cholla$Survival_t1 ~ cholla$logvol_t, xlab = "Volume in year t", ylab = "Survival in year t+1")
#Seems like larger cholla survive more, but hard to tell what's happening

#Try binning by size and plotting mean survival for each size bin for easier visualization
surv_bin <- cholla %>% 
  mutate(size_bin = cut_interval(logvol_t, n=8)) %>% #Makes 8 size bins based on log(vol_t)
  group_by(size_bin) %>% 
  summarise(mean_size = mean(logvol_t,na.rm=T), #Mean size in each bin
            mean_surv = mean(Survival_t1,na.rm=T)) #Mean survival in each bin
#Visualize
plot(cholla$Survival_t1 ~ cholla$logvol_t, xlab = "Volume in year t", ylab = "Survival in year t+1", col="gray")        
points(surv_bin$mean_size, surv_bin$mean_surv, pch=16, cex=2) #Calculated means based on bins
#Looks like survival does indeed increase as size increases

#2.2 Survival: fit model----------------------------
#use glm with family=binomial (because bernoulli is special case of binomial)
#check out glm help and list of families, also binomial documentation
survival_model <- glm(Survival_t1 ~ logvol_t, family = "binomial", data=cholla)

#Check model assumptions
#Try our usual method
plot(survival_model)
#None of this makes sense because these are designed for lm() models!
#So what are the assumptions and how do we check?

#Check for outliers/large residuals
qres<-qresid(survival_model)#calculates quantilized residuals with some introduced randomization to overcome discreteness
qqnorm(qres)#make qqplot
abline(0,1)#add expected line

#Check for constant dispersion by eye
scatter.smooth(qres~fitted(survival_model))
#Check for overdispersion using stats
resid.ssq <- sum(residuals(survival_model,type="pearson")^2)  
resid.df <- nrow(subset(cholla,!is.na(logvol_t) & !is.na(Survival_t1)))-length(coef(survival_model)) 
resid.ssq/resid.df 
#As a rule of thumb, this ratio of ss:df should be close to 1 
#if it's much over 1, that would indicate overdispersion

#Look at results
summary(survival_model)

#2.3 Survival: plot fitted model----------------------
#Our model predicts logit(survival), so survival=inverse logit(predictors)
#write inverse logit function
invlogit <- function(x){exp(x)/(1+exp(x))}

#plot data as before
plot(cholla$Survival_t1 ~ cholla$logvol_t, xlab = "Volume in year t", ylab = "Survival in year t+1", col="gray")        
points(surv_bin$mean_size, surv_bin$mean_surv, pch=16, cex=2)

#add predicted model using our function
#First make dummy x variable across range of sizes
x_dummy <- seq(min(cholla$logvol_t,na.rm=T),max(cholla$logvol_t,na.rm=T),0.1)
#add predicted line using dummy x
lines(x_dummy,invlogit(coef(survival_model)[1]+coef(survival_model)[2]*x_dummy),col="red",lwd=3) 
#or you can do it by calling $family$linkinv from the model
lines(x_dummy,survival_model$family$linkinv(coef(survival_model)[1]+coef(survival_model)[2]*x_dummy),col="blue",lty=2,lwd=3)

#This plot is also how we evaluate whether we used appropriate family/link...is it a good fit?

#3 Reproduction (try on your own)------------------------
#Hint, before you start, convert the "Goodbuds" variable to a binary outcome to form the response variable for reproduction

#4.1 Fertility: preliminary visualization------------------
plot(cholla$Goodbuds_t ~ cholla$logvol_t, xlab = "Volume in year t", ylab = "Fertility of flowering plants in year t", col="gray") 
#Fertility is strongly dependent on size

#4.2 Fertility: fit model----------------------------
fert_model <- glm(Goodbuds_t ~ logvol_t, family = "poisson", data=cholla)

#It's even harder to judge assumptions for poisson families using diagnositic plots
qqPlot(resid(fert_model))#but we expect the S shape with poisson resids
#Check for constant dispersion by eye
scatter.smooth(rstandard(fert_model)~sqrt(fitted(fert_model)))#standardized residuals against square root of predicted
#Check for overdispersion using stats
resid.ssq <- sum(residuals(fert_model,type="pearson")^2)  
resid.df <- nrow(subset(cholla,!is.na(logvol_t) & !is.na(Goodbuds_t)))-length(coef(fert_model)) 
resid.ssq/resid.df 
#As a rule of thumb, this ratio of ss:df should be close to 1 
#so we are way overdispersed....
#Let's just keep going for now

#Look at results
summary(fert_model)

#4.3 Fertility: plot fitted model---------------------------
plot(cholla$Goodbuds_t ~ cholla$logvol_t, xlab = "Volume in year t", ylab = "Fertility of flowering plants in year t", col="gray") 
lines(x_dummy,exp(coef(fert_model)[1]+coef(fert_model)[2]*x_dummy),col="red",lwd=3)
#The back-transforming function for the log link function in a poisson family is exp()
#It's easy to see where the overdispersion comes in, and it's hard to avoid with this data
#Small plants just don't reproduce
#Reasonable fit, but take the exact fitted "slope" coefficient with a grain of salt given overdispersion
