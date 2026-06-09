setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")

#----------Packages----------#
library(tidyverse)
library(magrittr)
library(plyr)
library(conflicted)
library(car)
library(lme4)
library(MuMIn)
library(emmeans)
library(indicspecies)
library(pwr)

#----------------------------------------------------------#
#-------------------RARITY RESPONSE RATIO------------------#
#----------------------------------------------------------#

rr <- readRDS("Processed Data/Robinhood Summary.rds")
rr %<>% filter(!is.na(spatial_rarity), !is.na(temporal_rarity), !is.na(response_ratio))

#spatial rarity
sr.lm <- lm(response_ratio ~ spatial_rarity, data = rr)
summary(sr.lm)
qqnorm(resid(sr.lm))
qqline(resid(sr.lm))
plot(sr.lm)
Anova(sr.lm)#spatial_rarity F = 3.286, df = 1,44, p = 0.077

#temporal rarity
tr.lm <- lm(response_ratio ~ temporal_rarity, data = rr)
summary(tr.lm)
qqnorm(resid(tr.lm))
qqline(resid(tr.lm))
plot(tr.lm)
Anova(tr.lm)#temporal_rarity F = 11.098, df = 1,44, p = 0.002

#----------------------------------------------------------#
#--------------------DELTA OCCUPANCY-----------------------#
#----------------------------------------------------------#

#spatial rarity
sr.occ <- lm(occupancy_shift ~ spatial_rarity, data = occ)
summary(sr.occ)
qqnorm(resid(sr.occ))#Check passed
qqline(resid(sr.occ))#Check passed
plot(sr.occ)#Check passed
Anova(sr.occ)#spatial_rarity F = 0.669, df = 1,44, p = 0.418

#temporal rarity
tr.occ <- lm(occupancy_shift ~ temporal_rarity, data = occ)
summary(tr.occ)
qqnorm(resid(tr.occ))#Check passed
qqline(resid(tr.occ))#Check passed
plot(tr.occ)#Check passed
Anova(tr.occ)#temporal_rarity F = 0.475, df = 1,44, p = 0.494



