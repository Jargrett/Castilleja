#Potato virus example
#Data from Clafkin et al. 2016
#AC 210210

#Housekeeping----------------
#Make sure clean environment

#Set working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Anny Stats/Lab 6")

#Load packages
library(tidyverse)
library(car)
library(statmod)
library(emmeans)

#Load data
PVY<-read.csv("potato virus data.csv")

#Take a peek at data structure
str(PVY)

#Change Year to factor
PVY$Year<-as.factor(PVY$Year)

#Preliminary visualizations------------
#Example question: Does PVY prevalence depend on year and %agricultural land cover at 1500m radius?

#We are interested in the proportion of infected plants, since total plants sampled differ

#Can calculate prevalence as the ratio of infected plants
PVY$Prev<-PVY$FinalInfected/PVY$NPlants

PVY$APT <- PVY$Naphids/PVY$NTraps

#Take a look at distribution of prevalence
hist(PVY$APT)
#Distribution of data is very right-skewed; lots of zeros

#What percentage of data is zero?
100*sum(PVY$APT == 0)/nrow(PVY)
#37.5%; seems high

#Take a look at potential relationships with predictors
plot(x=PVY$Ag1500, y=PVY$Prev)
#This definitely cannot be modeled with a linear model, but we knew that
boxplot(PVY$Prev~PVY$Year) #Skewed distribution, hard to tell if there's a difference

#Go ahead with a binomial glm, and check if residuals are overdispersed

#Binomial glm---------------------------------
#For binomial glms, there are 3 ways to specify the response
#See "Details" under help(family)
#For these data, we can use options 2 or 3

#Specify model using option 2
m1<-glm(Prev~Year*Ag1500, weights=NPlants, data=PVY, family="binomial")#note the inclusion of "weights"


model <- glm(Prev~Year*Ag1500, data=PVY, family="binomial")

#Check for outliers/large residuals
qres1<-qresid(m1)
qqnorm(qres1)#make qqplot
abline(0,1)#add expected line
#Doesn't look great

#Check for constant dispersion by eye
scatter.smooth(qres1~fitted(m1))
#doesn't actually look terrible, but larger positive residuals than negative residuals

#Check for overdispersion using stats
resid.ssq1 <- sum(residuals(m1,type="pearson")^2)  
resid.df1 <- nrow(subset(PVY,!is.na(Ag1500) & !is.na(Prev)))-length(coef(m1)) 
resid.ssq1/resid.df1 
#6.2, much above 1
#Overdispersed

#Try hurdle model instead

#Step 1: Infected yes/no-------------------------------
#Make new response variable
PVY$Infection<-ifelse(PVY$Prev>0,1,0)

#Model bionmial response of infection
m2<-glm(Infection~Year*Ag1500, data=PVY, family="binomial")#bernoulli response so no need for weights

#Check for outliers/large residuals
qres2<-qresid(m2)
qqnorm(qres2)#make qqplot
abline(0,1)#add expected line
#Looks much better

#Check for constant dispersion by eye
scatter.smooth(qres2~fitted(m2))
#more even distribution of negative and bositive residuals

#Check for overdispersion using stats
resid.ssq2 <- sum(residuals(m2,type="pearson")^2)  
resid.df2 <- nrow(subset(PVY,!is.na(Ag1500) & !is.na(Infection)))-length(coef(m2)) 
resid.ssq2/resid.df2 
#1.08, this is much more reasonable

#Look at model results
summary(m2) #none of the coefficients are significantly different from zero

#Can you do analysis-of-variance type tests on glm?
#Technically yes, using analysis-of-deviance as an analog
#see help(anova.glm) for details
#It is essentially comparing each effect to the null
anova(m2,test="Chisq") #need to specify test type for glm models

#Pseudo R-sq
1-(49.577/52.925)
#0.06325933

#Conclusion: Year and Ag1500 have no significant effects on the presence/absence of infection

#Step 2: If infected, how prevalent?----------------------------
#Subset data to only include prevalence>0 for this analysis
posPVY<-subset(PVY, Prev>0)

#Specify model for prevalence given infection
m3<-glm(Prev~Year*Ag1500, weights=NPlants, data=posPVY, family="binomial")#note the inclusion of "weights"

#Check for outliers/large residuals
qres3<-qresid(m3)
qqnorm(qres3)#make qqplot
abline(0,1)#add expected line
#Not perfect but better than before

#Check for constant dispersion by eye
scatter.smooth(qres3~fitted(m3))
#more even distribution of negative and bositive residuals

#Check for overdispersion using stats
resid.ssq3 <- sum(residuals(m3,type="pearson")^2)  
resid.df3 <- nrow(subset(posPVY,!is.na(Ag1500) & !is.na(Infection)))-length(coef(m3)) 
resid.ssq3/resid.df3 
#4.8, still pretty high, so there is overdispersion beyond zero-inflation
#Let's keep going with it and assess model fit by plotting predictions later

#Check resutls
summary(m3)
#Looks like there may be a significant slope for Ag1500 in Year1, and the slope in Y2 is not significantly different
anova(m3, test="Chisq") #Checks out. Significant slope for Ag1500, no interaction or year effects
#Pseudo R-sq
1-(88.756/113.728)
#0.2195765


#Plot results of hurdle model-----------------------------------
#We used base R plotting last week; let's try ggplot this time

#Plot infection presence/absence
#Calculate the fitted values based on our data
fitted2 <- predict(m2, type = "response", se = TRUE)
fitted2 #this includes 3 things: $fit, $se.fit, and $residual.scale
PVY$fitted2 <- fitted2$fit
PVY$fitted2.se <- fitted2$se.fit
#Calculate fitted +/- SE values
PVY$fitted.2.upper<-PVY$fitted2+PVY$fitted2.se
PVY$fitted.2.lower<-PVY$fitted2-PVY$fitted2.se
#Make plot
p1<-ggplot(data = PVY, aes(y = fitted2, x = Ag1500, col=Year) ) +
  geom_line() +
  geom_point(aes(y = Infection, x=Ag1500), alpha=0.5) + 
  geom_ribbon(aes(ymin = fitted.2.lower, ymax = fitted.2.upper),fill = "grey", alpha = 0.3)  + 
  theme_classic()+
  labs(x="%Ag at 1500m radius", y="Probability of infection presence")+
  annotate("text", label = "all effects P>0.05", x = 60, y = 0.9) 

#Plot prevalence given infection
#Calculate the fitted values based on our data
fitted3 <- predict(m3, type = "response", se = TRUE)
posPVY$fitted3 <- fitted3$fit
posPVY$fitted3.se <- fitted3$se.fit
#Calculate fitted +/- SE value
posPVY$fitted.3.upper<-posPVY$fitted3+posPVY$fitted3.se
posPVY$fitted.3.lower<-posPVY$fitted3-posPVY$fitted3.se
#Make plot
p2<-ggplot(data = posPVY, aes(y = fitted3, x = Ag1500, col=Year) ) +
  geom_line() +
  geom_point(aes(y = Prev, x=Ag1500), alpha=0.5) + 
  geom_ribbon(aes(ymin = fitted.3.lower, ymax = fitted.3.upper),fill = "grey", alpha = 0.3)  + 
  theme_classic()+
  labs(x="%Ag at 1500m radius", y="Prevalence of PVY infection")+
  annotate("text", label = c("Year P>0.05","%Ag P<0.0001", "Year:%Ag P>0.05"), x = c(60,60,60), y = c(0.75,0.7, 0.65)) 

#Make a two-paneled plot with labels
library(ggpubr)
ggarrange(p1, p2, labels = c("A", "B"),
          ncol = 2, nrow = 1)

