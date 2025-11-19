setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")
#load in relevant packages
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(rstatix)
library(sjPlot)

#Diversity and Cover data file (combined in excel)
castilleja.cover <- read.csv("Processed Data/castilleja cover complete.csv")
castilleja.cover$castilleja[castilleja.cover$castilleja == "Control"] <- "Absent"
castilleja.cover$castilleja[castilleja.cover$castilleja == "Castilleja"] <- "Present"
castilleja.cover <- as.data.frame(unclass(castilleja.cover),stringsAsFactors=TRUE)

#Diversity Analysis
div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(div)
Anova(div)
emmip(div, ~ castilleja ~ year)
emmeans(div, pairwise ~ castilleja|year)

rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(rich)
Anova(rich) 
emmip(rich, castilleja ~ year)
emmeans(rich, pairwise ~ species|year)

even <- lmer(even ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(even)
Anova(even) 
emmip(even, castilleja ~ year)
emmeans(even, pairwise ~ castilleja|year)

#Summary Output
tab_model(div, p.val = "kr", show.df = TRUE)
tab_model(rich, p.val = "kr", show.df = TRUE)
tab_model(even, p.val = "kr", show.df = TRUE)
