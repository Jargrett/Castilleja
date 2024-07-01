# Parasitic plant status
# Initiated: 5/13/2024

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Parasite Conservation")

library(tidyverse) # for data working
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)
library(lme4)

gap<- read.csv("Georgia's Parasitic Plants Copy.csv")

#-----------------Checking data structure-----------------#
str(gap)

#change all to factors
gap[sapply(gap, is.character)] <- lapply(gap[sapply(gap, is.character)], as.factor)
#-----------------Basic plotting-----------------#
ggplot(gap, aes(x=Georgia.Status..DNR., y = INat.Observations)) + 
  geom_point() +



