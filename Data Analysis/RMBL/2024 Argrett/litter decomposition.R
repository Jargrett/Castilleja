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
within <- read.csv("Litter Decomposition - 2024 Within Season.csv")
#yearlong <- read.csv("Litter Decomposition - 2023-2024 Yearlong.csv")
#bind
decomp <- rbind.fill(winter, within)
#check structure
str(decomp)
decomp <- as.data.frame(unclass(decomp),stringsAsFactors=TRUE)
within <- as.data.frame(unclass(within),stringsAsFactors=TRUE)
winter <- as.data.frame(unclass(winter),stringsAsFactors=TRUE)
#yearlong <- as.data.frame(unclass(yearlong),stringsAsFactors=TRUE)
#remove missing bags
decomp <- decomp %>% 
  filter(missing != "Yes") %>% 
  select(-c((missing)))

within <- within %>% 
  filter(missing != "Yes") %>% 
  select(-c((missing)))

winter <- winter %>% 
  filter(missing != "Yes") %>% 
  select(-c((missing)))

#first lets see if decomp is different between treatment
over.lm <- lm(final_dry_weight ~ litter*removal, data = within)
summary(over.lm)
Anova(over.lm)
emmip(over.lm, litter ~ removal)
emmeans(over.lm, pairwise ~ litter)

