setwd("~/Desktop/Castilleja/Data Analysis/Conner Data")

#Packages
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis

#load data
biomass <- read.csv("HPL Data - Biomass New.csv")
str(biomass)
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)
growth <- read.csv("HPL Data - Growth.csv")
str(growth)
growth <- as.data.frame(unclass(growth),stringsAsFactors=TRUE)
#-------------------------BIOMASS----------------------------------#
#Questions: How does parasite attachment influence host biomass?
#How does host identity influence parasite biomass?

#filter biomass between species
agalinis.biomass <- filter(biomass, focal_code == "AGPU")
little.biomass <- filter(biomass, focal_code == "SCSC")
aster.biomass <- filter(biomass, focal_code == "SYRA")


#Agalinis
#aboveground
ag.above.lm <- lmer(above ~ species*light + (1|replicate), data = agalinis.biomass)
summary(ag.above.lm)
Anova(ag.above.lm) #species Chisq = 7.5035, df = 3, p = 0.02348*
emmeans(ag.above.lm, pairwise ~ species|light)
emmip(ag.above.lm, ~ species ~ light)

#belowground
ag.below.lm <- lmer(below ~ species*light + (1|replicate), data = agalinis.biomass)
summary(ag.below.lm)
Anova(ag.below.lm) #species Chisq = 18.8086, df = 2, p < 0.001
emmeans(ag.below.lm, pairwise ~ species|light)
emmip(ag.below.lm, ~ species ~ light)

#total
ag.total.lm <- lmer(total ~ species*light + (1|replicate), data = agalinis.biomass)
summary(ag.total.lm)
Anova(ag.total.lm) #
emmeans(ag.total.lm, pairwise ~ species|light)
emmip(ag.total.lm, ~ species ~ light)



sd.little.bio <- little.biomass %>% drop_na(total)
sd.little.bio <- sd.little.bio %>% 
  group_by(species, light) %>% 
  dplyr::summarise(mean= mean(total),
                   se = sd(total)/sqrt(n()))


little.plot <- ggplot(data = sd.little.bio, aes(x = species, y = mean, fill = light)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  labs(x = "Pot Treatment", y = "Total Biomass (g)") +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 15))+
  ylim(0,5)

little.plot


#Agalinis
#aboveground
little.lm <- lmer(above ~ species*light + (1|replicate), data = little.biomass)
summary(little.lm)
Anova(little.lm) #species Chisq = 7.5035, df = 3, p = 0.02348*
emmeans(little.lm, pairwise ~ species|light)
emmip(little.lm, ~ species ~ light)

sd.aster.bio <- aster.biomass %>% drop_na(total)
sd.aster.bio <- sd.aster.bio %>% 
  group_by(species, light) %>% 
  dplyr::summarise(mean= mean(total),
                   se = sd(total)/sqrt(n()))

aster.plot <- ggplot(data = sd.aster.bio, aes(x = species, y = mean, fill = light)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  labs(x = "Pot Treatment", y = "Total Biomass (g)") +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 15))+
  ylim(0,5)

aster.plot

#aboveground
aster.lm <- lm(above ~ species*light, data = aster.biomass)
summary(aster.lm)
Anova(aster.lm)
emmeans(aster.lm, pairwise ~ species|light)
emmip(aster.lm, ~ species ~ light)


agalinis.biomass <- filter(agalinis.biomass, species != "P1")
sd.agalinis.bio <- agalinis.biomass %>% drop_na(total)
sd.agalinis.bio <- sd.agalinis.bio %>% 
  group_by(species, light) %>% 
  dplyr::summarise(mean= mean(total),
                   se = sd(total)/sqrt(n()))

ag.plot <- ggplot(data = sd.agalinis.bio, aes(x = species, y = mean, fill = light)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  labs(x = "Pot Treatment", y = "Total Biomass (g)") +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 15))+
  ylim(0,5)

ag.plot

