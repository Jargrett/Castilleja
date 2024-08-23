# EL Soil Moisture Calulations
# Initiated: 8/23/2023
# Complted:
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

# Load-in packages
library(tidyverse) # for data working
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # Regression Analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models

soil.moisture.23 <- read.csv("Emerald Lake Soil Moisture Data - 2023.csv")
soil.moisture.24 <- read.csv("Emerald Lake Soil Moisture Data - 2024.csv")

#checking data structure
str(soil.moisture.23)

soil.moisture.23$date <- as.factor(soil.moisture.23$date)
soil.moisture.24$date <- as.factor(soil.moisture.24$date)
soil.moisture.23$plot <- as.factor(soil.moisture.23$plot)
soil.moisture.24$plot <- as.factor(soil.moisture.24$plot)
moisture.23 <- na.omit(soil.moisture.23)
moisture.24 <- na.omit(soil.moisture.24)

#Create a new DF that groups by averages per plot
ASM.24 <- moisture.24 %>%
  group_by(plot) %>%
  mutate(asm24=mean(average_soil_moisture)) %>%
  select(plot,asm24) %>%
  distinct(plot, asm24, .keep_all = TRUE)

ASM.23 <- moisture.23 %>%
  group_by(plot) %>%
  mutate(asm23=mean(average_soil_moisture)) %>%
  select(plot, asm23) %>%
  distinct(plot, asm23, .keep_all = TRUE)

#Data Visualizations
el.plot <- read.csv("Emerald Lake Plot Data - Info.csv", header=TRUE)
el.plot$plot <- as.factor(el.plot$plot)
el.data.23 <- merge(el.plot,ASM.23, by = "plot")
el.data.full <- merge(el.data.23,ASM.24, by = "plot")


el.moisture <- el.data.full %>% select(plot,pair,block,removal,asm23,asm24)
str(el.moisture)
el.moisture$removal <- as.factor(el.moisture$removal)
el.moisture$block <- as.factor(el.moisture$block)
el.moisture$pair <- as.factor(el.moisture$pair)
#linear regression and ANOVA to to see if Block and Elevation effect SM
sm.lm <- lmer(asm24 ~ removal*block + (1|pair), data = el.moisture)
summary(sm.lm)
Anova(sm.lm)
emmeans(sm.lm, pairwise ~ removal|block)

