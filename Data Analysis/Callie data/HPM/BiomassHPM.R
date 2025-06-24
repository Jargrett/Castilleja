#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/HPM")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis

biomass <- read.csv("HPM Biomass - Biomass.csv")
str(biomass)
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)


#Total Belowground biomass analysis

#agalinis
agalinis.biomass <- filter(biomass, species == "AGPU")

#belowground
ag.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.bg.lm)
Anova(ag.bg.lm)
emmeans(ag.bg.lm, pairwise ~ type|treatment)
emmip(ag.bg.lm, ~ type ~ treatment)

#heterotheca
hetero.biomass <- filter(biomass, species == "HESU")

#belowground
he.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.bg.lm)
Anova(he.bg.lm)
emmeans(he.bg.lm, pairwise ~ type|treatment)
emmip(he.bg.lm, ~ type ~ treatment)


#Standard error calculations
ag.biomass <- agalinis.biomass %>% drop_na(below_total)
AGPU.below <- ag.biomass %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(below_total),
                   se = sd(below_total)/sqrt(n()))

he.biomass <- hetero.biomass %>% drop_na(below_total)
HESU.below <- he.biomass %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(below_total),
                   se = sd(below_total)/sqrt(n()))

