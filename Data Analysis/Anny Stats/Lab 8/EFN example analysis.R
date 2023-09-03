
#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 8")
#load libraries
library(car)
library(lme4)#this is new
library(lmerTest)
library(performance)#this is new
library(see)#this is new
library(ggpubr)
library(emmeans)

#Import data
EFN<-read.csv("EFN data.csv")

str(EFN)

#Make Plant a factor
EFN$Plant<-as.factor(EFN$Plant)
table(EFN$Plant,EFN$EFNtype)#Each plant measured different numbers of times

#Do foliar and floral nectaries differ in nectar volume?-----------------
m1<-lmer(Sugar~EFNtype + (1|Plant), data=EFN)
check_model(m1)
summary(m1)

Anova(m1)
anova(m1)
#Do foliar and floral nectaries differ in nectar calories?-----------------
m2<-lmer(Sugar~EFNtype + (1|Plant), data=EFN)

