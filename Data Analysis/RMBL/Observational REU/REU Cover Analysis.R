#Castilleja Diveristy Analysis

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(glmm)
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(indicspecies)
library(statmod)

# #Load in 2023 + 2024 datasets (For cover analysis)
# CALI.24 <- read.csv("Cali 2024 - Cover.csv")
# CALI.23 <- read.csv("Cali 2023 - Cover.csv")
# CALI.Johnson <- read.csv("Cali Johnson Hill 2024 - Cover.csv")
# 
# #shifting structure to long format rather than matrix format
# cali.24.long<- pivot_longer(CALI.24, cols = Bare.ground:Wyethia.amplexicaulis,
#                          names_to = "species",
#                         values_to = "cover")
# 
# cali.23.long<- pivot_longer(CALI.23, cols = Bare.ground:Wyethia.amplexicaulis,
#                             names_to = "species",
#                             values_to = "cover")
# 
# cali.johnson.long<- pivot_longer(CALI.Johnson, cols = Bare.ground:Viola.praemorsa,
#                             names_to = "species",
#                             values_to = "cover")
# #Removing 0s from dateset
# cali.24.cover <- filter(cali.24.long, cover > 0)
# cali.23.cover <- filter(cali.23.long, cover > 0)
# cali.johnson.cover <- filter(cali.johnson.long, cover > 0)
# 
# cali.cover.comb <- rbind(cali.23.cover,cali.24.cover,cali.johnson.cover)
# cali.cover.zero <- rbind(cali.23.l,cali.24.long,cali.johnson.long)
# 
# write.csv(cali.cover.comb, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Species Cover.csv", row.names=FALSE)

#-----------------------------COVER ANALYSIS BEGINS----------------------------#

cali.cover <- read.csv("CALI Species Cover.csv")

str(cali.cover)
cali.cover$site <- as.factor(cali.cover$site)
cali.cover$year <- as.factor(cali.cover$year)
cali.cover$castilleja <- as.factor(cali.cover$castilleja)
cali.cover$species <- as.factor(cali.cover$species)

hist(cali.cover$cover)
boxplot(cali.cover$cover~cali.cover$castilleja)

#----------------Individual Species analysis--------------------#
#Does the presence of Castilleja alter the cover of species?
#Species list:
#want to understand whether the cover of a species changes, if it exists when Castilleja is present 
#should only consider a pair of plots where that species exists in both plots

#We need to get the data into that format, 
#i.e comparison of cover values between paired plots with a give species present

#First we subset the data
indv.species = subset(cali.cover, select = -c(2,5,7,8,11:13))
#The we pivot the columns wider creating two cover columns based on a given species
species.pair <- indv.species %>%
  pivot_wider(names_from = castilleja,
              values_from = c(cover)) %>% 
  drop_na()

#subset data for pair 1
#pair <- indv.species %>% dplyr::filter(pair %in% c("1"))
# pair3 <- pair2 %>%
#   group_by(species) %>% 
#   fill(Castilleja, Control, .direction = 'up') %>% 
#   fill(Castilleja, Control) %>% 
#   distinct()

#final <- pair3 %>% dplyr::filter(plot %in% c("1"))

#Birds eye view
plant.lm <- lm(plant.cover ~ castilleja*site, data = cali.cover)
summary(plant.lm)
Anova(plant.lm)
emmeans(plant.lm, pairwise ~ castilleja|site)
emmip(plant.lm, castilleja ~ site)

