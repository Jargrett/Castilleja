#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/HPM")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis

biomass <- read.csv("HPM Biomass - Biomass.csv")
str(biomass)
biomass$treatment[biomass$treatment == "innoculated"] <- "AMF"
biomass$treatment[biomass$treatment == "sterilized"] <- "Control"
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)


#Total Belowground biomass analysis

#agalinis
agalinis.biomass <- filter(biomass, species == "AGPU")


#aboveground
ag.ag.lm <- lmer(above_total ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.ag.lm)
Anova(ag.ag.lm)
emmeans(ag.ag.lm, pairwise ~ type|treatment)
emmip(ag.ag.lm, ~ type ~ treatment)

#belowground
ag.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.bg.lm)
Anova(ag.bg.lm)
emmeans(ag.bg.lm, pairwise ~ type|treatment)
emmip(ag.bg.lm, ~ type ~ treatment)

#Total
ag.tb.lm <- lmer(total_biomass ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.tb.lm)
Anova(ag.tb.lm)
emmeans(ag.tb.lm, pairwise ~ type|treatment)
emmip(ag.tb.lm, ~ type ~ treatment)

#Root:shoot
ag.rs.lm <- lmer(root_shoot ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.rs.lm)
Anova(ag.rs.lm)
emmeans(ag.rs.lm, pairwise ~ type|treatment)
emmip(ag.rs.lm, ~ type ~ treatment)

#heterotheca
hetero.biomass <- filter(biomass, species == "HESU")

#aboveground
he.ag.lm <- lmer(above_total ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.ag.lm)
Anova(he.ag.lm)
emmeans(he.ag.lm, pairwise ~ type|treatment)
emmip(he.ag.lm, ~ type ~ treatment)

#belowground
he.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.bg.lm)
Anova(he.bg.lm)
emmeans(he.bg.lm, pairwise ~ type|treatment)
emmip(he.bg.lm, ~ type ~ treatment)

#total
he.tb.lm <- lmer(total_biomass ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.tb.lm)
Anova(he.tb.lm)
emmeans(he.tb.lm, pairwise ~ type|treatment)
emmip(he.tb.lm, ~ type ~ treatment)

#Root:shoot
he.rs.lm <- lmer(root_shoot ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.rs.lm)
Anova(he.rs.lm)
emmeans(he.rs.lm, pairwise ~ type|treatment)
emmip(he.rs.lm, ~ type ~ treatment)


#Standard error calculations
sd.above <- biomass %>% drop_na(above_total)
sd.above <- sd.above %>% 
  group_by(species,treatment, type) %>% 
  dplyr::summarise(mean= mean(above_total),
                   se = sd(above_total)/sqrt(n()))

sd.below <- biomass %>% drop_na(below_total)
sd.below <- sd.below %>% 
  group_by(species,treatment, type) %>% 
  dplyr::summarise(mean= mean(below_total),
                   se = sd(below_total)/sqrt(n()))

sd.total <- biomass %>% drop_na(total_biomass)
sd.total <- sd.total %>% 
  group_by(species, treatment, type) %>% 
  dplyr::summarise(mean= mean(total_biomass),
                   se = sd(total_biomass)/sqrt(n()))

sd.rs <- biomass %>% drop_na(root_shoot)
sd.rs <- sd.rs %>% 
  group_by(species, treatment, type) %>% 
  dplyr::summarise(mean= mean(root_shoot),
                   se = sd(root_shoot)/sqrt(n()))

hetero.above <- filter(sd.above, species == "HESU")
hetero.above.plot <- ggplot(data = hetero.above, aes(x = type, y = mean, fill = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Hemiparasite presence", y = "Aboveground Biomass") +
  facet_wrap(~treatment) +
  ylim(0,1.5)
hetero.above.plot 

agalinis.above <- filter(sd.above, species == "AGPU")
agpu.above.plot <- ggplot(data = agalinis.above, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Hemiparasite presence", y = "Aboveground Biomass") +
  facet_wrap(~treatment) +
  ylim(0,1.5)
agpu.above.plot 

hetero.below <- filter(sd.below, species == "HESU")
hetero.below.plot <- ggplot(data = hetero.below, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Hemiparasite presence", y = "Belowground Biomass") +
  facet_wrap(~treatment) +
  ylim(0,3)
hetero.below.plot 

agalinis.below <- filter(sd.below, species == "AGPU")
agpu.below.plot <- ggplot(data = agalinis.below, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Hemiparasite presence", y = "Belowground Biomass") +
  facet_wrap(~treatment) +
  ylim(0,3)
agpu.below.plot 


biomass.plots <- ggarrange(agpu.above.plot,hetero.above.plot, agpu.below.plot, hetero.below.plot,
                          labels = c("A", "B","C","D"), 
                          nrow = 2, ncol = 2, legend = FALSE)
biomass.plots

agalinis.total <- filter(sd.total, species == "AGPU")
agpu.total.plot <- ggplot(data = agalinis.total, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Hemiparasite presence", y = "Total Biomass") +
  facet_wrap(~treatment) +
  ylim(0,5)
agpu.total.plot 

hetero.total <- filter(sd.total, species == "HESU")
hetero.total.plot <- ggplot(data = hetero.total, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Hemiparasite presence", y = "Total Biomass") +
  facet_wrap(~treatment) +
  ylim(0,5)
hetero.total.plot

agalinis.rs<- filter(sd.rs, species == "AGPU")
agpu.rs.plot <- ggplot(data = agalinis.rs, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Hemiparasite presence", y = "Root/Shoot Ratio") +
  facet_wrap(~treatment) +
  ylim(0,3.5)
agpu.rs.plot 

hetero.rs<- filter(sd.rs, species == "HESU")
hesu.rs.plot <- ggplot(data = hetero.rs, aes(x = type, y = mean, color = type)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Hemiparasite presence", y = "Root/Shoot Ratio") +
  facet_wrap(~treatment) +
  ylim(0,3.5)
hesu.rs.plot 

rootshoot.plots <- ggarrange(agpu.total.plot, hetero.total.plot, agpu.rs.plot, hesu.rs.plot,
                           labels = c("A", "B","C","D"), 
                           nrow = 2, ncol = 2, legend = FALSE)
rootshoot.plots 
