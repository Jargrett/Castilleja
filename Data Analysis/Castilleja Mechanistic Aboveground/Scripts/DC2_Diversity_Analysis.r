setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(rstatix)
library(magrittr)#for data wrangling and restructuring

conflicted::conflicts_prefer(lme4::lmer)

#----------------------------------------------------------#
#---------------------TEMPORAL DIVERSITY-------------------#
#----------------------------------------------------------#

diversity <- readRDS("Processed Data/Plant Diversity Full.rds")
diversity %<>% filter(year != "0")

#diversity
div.lmm <- lmer(div ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(div.lmm)
qqnorm(resid(div.lmm))#Check passed
qqline(resid(div.lmm))#Check passed
plot(div.lmm)#Check passed
Anova(div.lmm)#removal:year chisq = 10.546, df = 3, p = 0.005
emmeans(div.lmm, pairwise ~ removal|year)
emmip(div.lmm, removal ~ year)

#richness
rich.lmm <- lmer(rich ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(rich.lmm)
qqnorm(resid(rich.lmm))#Check passed
qqline(resid(rich.lmm))#Check passed
plot(rich.lmm)#Check passed
Anova(rich.lmm)#removal:year chisq = 14.927, df = 2, p < 0.001
emmeans(rich.lmm, pairwise ~ removal|year)
emmip(rich.lmm, removal ~ year)

#eveness
even.lmm <- lmer(even ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(even.lmm)
qqnorm(resid(even.lmm))#Check passed
qqline(resid(even.lmm))#Check passed
plot(even.lmm)#Check passed
Anova(even.lmm)#n.s.  
emmeans(even.lmm, pairwise ~ removal|year)
emmip(even.lmm, removal ~ year)

#----------------------------------------------------------#
#---------------------CHANGE IN DIVERSITY------------------#
#----------------------------------------------------------#

#---------delta diversity Analysis---------#
delta_diversity <- readRDS("Processed Data/Delta Diversity.rds")

#Change in diversity from year 1 to year 3
delta.rich <- lmer(delta_rich ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.rich)
Anova(delta.rich)#removal chisq = 16.9246, df = 1, p < 0.001
emmip(delta.rich, litter~removal)
emmeans(delta.rich, pairwise ~ litter|removal)

#Change in richness from year 1 to year 3
delta.div <- lmer(delta_div ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.div)
Anova(delta.div)#removal chisq = 17.0773, df = 1, p < 0.001
emmip(delta.div, litter~removal)
emmeans(delta.div, pairwise ~ litter|removal)

#Change in evenness from year 1 to year 3
delta.even <- lmer(delta_even ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.even)
Anova(delta.even)#litter:removal chisq = 6.5534, df = 3, p = 0.088
emmip(delta.even, litter~removal)
emmeans(delta.even, pairwise ~ litter|removal)


