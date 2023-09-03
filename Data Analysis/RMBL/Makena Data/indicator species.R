#Intitated: 8/28/2023
#For the purpose of analyzing indicator species by treatment

#Set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

#Loading necessary packages
library(indicspecies) # indicator species
library(tidyverse) # data working
library(vegan) # diversity metrix
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # regression analysis

cali <- read.csv("Combined deer creek data - Individuals.csv")
case <- read.csv("Combined septentrionalis plant data - Cover.csv")

cali.ind <- cali[ -c(1:4,15)]
case.ind <- case[ -c(1:4, 7, 9)]

cali.cover = cali.ind[,2:ncol(cali.ind)]
cali.treatment = cali.ind$Treatment

case.cover = case.ind[,3:ncol(case.ind)]
case.Site = case.ind$Treatment

#indicator species code
inv = multipatt(cali.cover, cali.treatment, func = "r.g", control = how(nperm=9999))
summary(inv)

inv2 = multipatt(case.cover, case.treatment, func = "r.g", control = how(nperm=9999))
summary(inv2)
