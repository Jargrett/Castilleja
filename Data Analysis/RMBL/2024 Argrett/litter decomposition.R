#Litter Decomposition
#Initiated 8/23/24

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

remotes::install_github("cornwell-lab-unsw/litterfitter")

library(litterfitter)#for k-curve fitting
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(lme4)#for modeling linear mixed effect models
library(nlme)#alternative for modeling linear mixed effect models
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis

#load-in the data
winter <- read.csv("Litter Decomposition - 2023-2024 Overwinter.csv")
#check structure
str(winter)
winter$litter <- as.factor(winter$litter)
winter$removal <- as.factor(winter$removal)
winter$deployment_duration <- as.factor(winter$deployment_duration)
#remove missing bags
winter <- winter %>% 
  filter(missing != "Yes") %>% 
  select(-c((missing)))

#first lets see if decomp is different between treatment
over.lm <- lm(final_dry_weight ~ litter, data = winter)
summary(over.lm)
Anova(over.lm)
emmip(over.lm, removal ~ litter)
emmeans(over.lmm, pairwise ~ litter)

