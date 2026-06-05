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



#-------Bare Ground Analysis------#
envi <- readRDS("Processed Data/Environmental Cover.rds")
#Total Environment
total.lmm <- lmer(total_cover ~ year*removal*litter + (1|year) + (1|block) + (1|pair), data = envi)
summary(total.lmm)
Anova(total.lmm)#removal chisq = 12.788, df = 1, p < 0.001 
emmeans(total.lmm, pairwise ~ year|removal)
emmip(total.lmm, year ~ removal)