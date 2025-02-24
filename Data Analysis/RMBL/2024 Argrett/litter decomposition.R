#Litter Decomposition
#Initiated 8/23/24

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

#remotes::install_github("cornwell-lab-unsw/litterfitter")

library(litterfitter)#for k-curve fitting
library(tidyverse)#for data wrangling and restructuring
library(dplyr)#for data wrangling and restructuring
library(lme4)#for modeling linear mixed effect models
library(nlme)#alternative for modeling linear mixed effect models
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(dplyr)
library(plyr)

#load-in the data
winter <- read.csv("Litter Decomposition - 2023-2024 Overwinter.csv")
within <- read.csv("Litter Decomposition - 2024 Within Season.csv")
yearlong <- read.csv("Litter Decomposition - 2023-2024 Full Year.csv")
#yearlong <- read.csv("Litter Decomposition - 2023-2024 Yearlong.csv")
#bind
decomp <- rbind.fill(winter, within, yearlong)
#check structure
str(decomp)
decomp <- as.data.frame(unclass(decomp),stringsAsFactors=TRUE)
within <- as.data.frame(unclass(within),stringsAsFactors=TRUE)
winter <- as.data.frame(unclass(winter),stringsAsFactors=TRUE)
yearlong <- as.data.frame(unclass(yearlong),stringsAsFactors=TRUE)
#yearlong <- as.data.frame(unclass(yearlong),stringsAsFactors=TRUE)
#remove missing bags
decomp <- decomp %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

within <- within %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

winter <- winter %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

yearlong <- yearlong %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

decomp <- decomp %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
within <- within %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
winter <- winter %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
yearlong <- yearlong %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)

decomp <- decomp %>% mutate(time = deployment_duration/365)
within <- within %>% mutate(time = deployment_duration/365)
winter <- winter %>% mutate(time = deployment_duration/365)
yearlong <- yearlong %>% mutate(time = deployment_duration/365)


#first lets see if decomp is different between treatment

#Removing Outliers for testing
decomp <- decomp[-c(57,69,157),]

over.lm <- lm(mass_remaining ~ removal*litter + deployment_period, data = decomp)
summary(over.lm)
Anova(over.lm)
emmip(over.lm, litter ~ deployment_period)
emmeans(over.lm, pairwise ~ litter|deployment_period)

litter <- fit_litter(time = decomp$time,
                     mass.remaining = decomp$mass_remaining,
                     model = "weibull",
                     iters=1000)

plot_multiple_fits(time = decomp$time,
                   mass.remaining = decomp$mass_remaining,
                   model=c("neg.exp","weibull"),
                   iters=500)

summary(litter)

class(litter)
plot(litter)

boxplot(decomp$mass_remaining ~ decomp$litter)

ggplot(decomp, aes(x = litter, y = mass_remaining)) +
  geom_jitter(aes(color = (litter))) +
  facet_wrap(~deployment_period) +
  labs(x = "Litter Treatment", y = "Mass remaining")

