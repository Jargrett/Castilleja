#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models

soil.overwinter <- read.csv("Soil Nutrients - Overwinter 2023.csv")
soil.within <- read.csv("Soil Nutrients - Within 2024.csv")
str(soil.within)
soil.within$Pair <- as.factor(soil.within$Pair)
soil.within$Plot <- as.factor(soil.within$Plot)
soil.within$Litter <- as.factor(soil.within$Litter)
soil.within$Removal <- as.factor(soil.within$Removal)
soil.within$Block <- as.factor(soil.within$Block)

nitrate <- ggplot(soil.within, aes(x = Litter, y = P)) +
  geom_point(aes(color = (Removal))) +
  labs(x = "Litter Treatment", y = "P")
  
nitrate

nitrate.lmm <- lmer(Nitrate ~ Litter*Removal + (1|Block) + (1|Pair), data = soil.within)
summary(nitrate.lmm)
Anova(nitrate.lmm)
emmip(nitrate.lmm, Litter ~ Removal)
emmeans(nitrate.lmm, pairwise ~  Removal|Litter)

k.lmm <- lmer(P ~ Litter*Removal + (1|Block) + (1|Pair), data = soil.within)
summary(k.lmm)
Anova(k.lmm)
emmip(k.lmm, Litter ~ Removal)
emmeans(k.lmm, pairwise ~  Litter)

ggplot(data = soil.within, aes(x = Litter, y = Nitrate, fill = Removal)) +
  geom_bar(stat ="identity", position=position_dodge(1)) +
  labs(x = "Litter Treatment", y = "Potassium") +
  scale_fill_manual(values=c('black','lightgray')) +
  theme_classic()


Mn.lmm <- lmer(Mn ~ Removal + (1|Block) + (1|Pair), data = soil)
summary(Mn.lmm)
Anova(Mn.lmm)
emmip(Mn.lmm, Removal ~ Litter)
emmeans(Mn.lmm, pairwise ~ Removal)
