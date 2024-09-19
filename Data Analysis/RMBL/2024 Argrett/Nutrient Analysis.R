#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models

soil <- read.csv("Soil Nutrients - Full.csv")

str(soil)
soil$pair <- as.factor(soil$pair)
soil$plot <- as.factor(soil$plot)
soil$litter <- as.factor(soil$litter)
soil$removal <- as.factor(soil$removal)
soil$block <- as.factor(soil$block)

nitrate <- ggplot(soil, aes(x = litter, y = P)) +
  geom_point(aes(color = (removal))) +
  facet_wrap(~burial) +
  labs(x = "Litter Treatment", y = "Nitrate")

nitrate

nitrate.lmm <- lmer(P ~ litter*removal + burial + (1|block) + (1|pair), data = soil)
summary(nitrate.lmm)
Anova(nitrate.lmm)
emmip(nitrate.lmm, litter ~ removal)
emmeans(nitrate.lmm, pairwise ~  removal|litter)
        