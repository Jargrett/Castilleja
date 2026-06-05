setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")
#----------Data importing, cleaning, and resctructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)
library(statmod)
library(lme4)
library(emmeans) # for comparison of means
library(betareg)

#----------Biomass Analysis----------#
biomass <- readRDS("Processed Data/Plant Biomass.rds")
#Log total plant community biomass (Castilleja Excluded)
total.bio.lmm <- lmer(log(total_no_cas) ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(total.bio.lmm)
Anova(total.bio.lmm)# No significance removal Chisq = 1.5329, 1, p = 0.216
qqnorm(resid(total.bio.lmm))
qqline(resid(total.bio.lmm))
plot(total.bio.lmm)#
emmeans(total.bio.lmm, pairwise ~ litter|removal)
