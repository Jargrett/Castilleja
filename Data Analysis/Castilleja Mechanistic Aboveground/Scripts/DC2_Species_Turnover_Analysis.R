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
library(conflicted)
library(codyn)
library(labdsv)
library(vegan)


#----------------------------------------------------------#
#----------------COMPOSITIONAL CHANGE----------------------#
#----------------------------------------------------------#
comp_change <- readRDS("Processed Data/Compositional Change.rds")
comp.lmm <- lmer(distance ~ removal*litter + (1|block) + (1|pair), data = comp_change)
summary(comp.lmm)
Anova(comp.lmm)#removal: chisq = 3.5907, df= 1, p = 0.058*

#----------------------------------------------------------#
#--------------------SPECIES TURNOVER----------------------#
#----------------------------------------------------------#
total_turn <- readRDS("Processed Data/Species Turnover.rds")

turn.lmm <- lmer(total ~ removal + (1|block) + (1|pair), data = total_turn)
summary(turn.lmm)
Anova(turn.lmm)#removal: Chisq = 0.3679, df= 1, p = 0.544

gain.lmm <- lmer(appearance ~ removal + (1|block) + (1|pair), data = total_turn)
summary(gain.lmm)
Anova(gain.lmm)#removal: Chisq = 13.098 , df= 1, p < 0.001

loss.lmm <- lmer(disappearance ~ removal + (1|block) + (1|pair), data = total_turn)
summary(loss.lmm)
Anova(loss.lmm)#removal: Chisq = 12.6935, df= 1, p < 0.001
