#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2023 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models

soil <- read.csv("Soil Nutrients - Overwinter 2023.csv")
str(soil)
soil$Pair <- as.factor(soil$Pair)
soil$Plot <- as.factor(soil$Plot)
soil$Litter <- as.factor(soil$Litter)
soil$Removal <- as.factor(soil$Removal)
soil$Block <- as.factor(soil$Block)

nitrate <- ggplot(soil, aes(x = Litter, y = P)) +
  geom_point(aes(color = (Removal))) +
  labs(x = "Litter Treatment", y = "P")
  
nitrate


k.lmm <- lmer(K ~ Litter*Removal + (1|Block) + (1|Pair), data = soil)
summary(k.lmm)
Anova(k.lmm)
emmip(k.lmm, Litter ~ Removal)
emmeans(k.lmm, pairwise ~  Removal|Litter)

ggplot(data = soil, aes(x = Litter, y = K, fill = Removal)) +
  geom_bar(stat ="identity", position=position_dodge(1)) +
  labs(x = "Litter Treatment", y = "Potassium") +
  scale_fill_manual(values=c('black','lightgray')) +
  theme_classic()


Mn.lmm <- lmer(Mn ~ Removal + (1|Block) + (1|Pair), data = soil)
summary(Mn.lmm)
Anova(Mn.lmm)
emmip(Mn.lmm, Removal ~ Litter)
emmeans(Mn.lmm, pairwise ~ Removal)
