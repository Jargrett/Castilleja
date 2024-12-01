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
soil <- as.data.frame(unclass(soil),stringsAsFactors=TRUE)


nitrate <- ggplot(soil, aes(x = litter, y = K)) +
  geom_point(aes(color = (removal))) +
  facet_wrap(~burial) +
  labs(x = "Litter Treatment", y = "Nitrate")

nitrate

K.lmm <- lmer(K ~ litter*removal + burial + (1|block) + (1|pair), data = soil)
summary(K.lmm)
Anova(K.lmm)
emmip(K.lmm, litter ~ removal)
emmeans(K.lmm, pairwise ~  removal|litter)
        