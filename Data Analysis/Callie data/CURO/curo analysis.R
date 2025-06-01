#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/CURO")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis


data <- read.csv("Callie Data - Master.csv")
data$parasite[data$parasite == "no"] <- "absent"
data$parasite[data$parasite == "yes"] <- "present"
data <- as.data.frame(unclass(data),stringsAsFactors=TRUE)
data <- filter(data, haustoria != "no")
data$solidago_rs <- data$solidago_below/data$solidago_above
data$agalinis_rs <- data$agalinis_below/data$agalinis_above

#Above
soal.agb.lm <- lmer(solidago_above ~ treatment*parasite + (1|soil_id), data = data)
summary(soal.agb.lm)
Anova(soal.agb.lm) #treatment P = 0.01074 (significant difference of 1.1 gram in parasite presence plots)
emmeans(soal.agb.lm, pairwise ~ treatment|parasite)

agpu.agb.lm <- lmer(agalinis_above ~ treatment + (1|soil_id), data = data)
summary(agpu.agb.lm)
Anova(agpu.agb.lm)
emmeans(agpu.agb.lm, pairwise ~ treatment|parasite)

#Below
soal.bgb.lm <- lmer(solidago_below ~ treatment*parasite + (1|soil_id), data = data)
summary(soal.bgb.lm)
Anova(soal.bgb.lm) #treatment P = 0.007758 (difference of 1.3 gram in parasite absence plots)
emmeans(soal.bgb.lm, pairwise ~ treatment|parasite)

#Root:shoot
soal.rs.lm <- lmer(solidago_rs ~ treatment*parasite + (1|soil_id), data = data)
summary(soal.rs.lm)
Anova(soal.rs.lm) #treatment P = 0.007758 (difference of 1.3 gram in parasite absence plots)
emmeans(soal.rs.lm, pairwise ~ treatment|parasite)

agpu.rs.lm <- lmer(agalinis_rs ~ treatment + (1|soil_id), data = data)
summary(agpu.rs.lm)
Anova(agpu.rs.lm) #treatment P = 0.007758 (difference of 1.3 gram in parasite absence plots)
emmeans(agpu.rs.lm, pairwise ~ treatment|parasite)

#Total
soal.total.lm <- lmer(solidago_total ~ treatment*parasite + (1|soil_id), data = data)
summary(soal.total.lm)
Anova(soal.total.lm) #treatment P = 0.016.04 (difference of 1.3 gram in parasite absence plots)
emmeans(soal.total.lm, pairwise ~ treatment|parasite)

agpu.total.lm <- lmer(agalinis_total ~ treatment + (1|soil_id), data = data)
summary(agpu.total.lm)
Anova(agpu.total.lm)
emmeans(agpu.total.lm, pairwise ~ treatment)

#Growth
soal.growth.lm <- lmer(solidago_final_height ~ treatment*parasite + (1|soil_id), data = data)
summary(soal.growth.lm)
Anova(soal.growth.lm)
emmeans(soal.growth.lm, pairwise ~ treatment|parasite)

agpu.growth.lm <- lmer(agalinis_final_height ~ treatment + (1|soil_id), data = data)
summary(agpu.growth.lm)
Anova(agpu.growth.lm)
emmeans(agpu.growth.lm, pairwise ~ treatment)


#Standard error calculations
AGB.s <- data %>% drop_na(solidago_above)
SOAL.AGB <- AGB.s %>% 
  group_by(treatment, parasite) %>% 
  dplyr::summarise(mean= mean(solidago_above),
                   se = sd(solidago_above)/sqrt(n()))

BGB.s <- data %>% drop_na(solidago_below)
SOAL.BGB <- BGB.s %>% 
  group_by(treatment, parasite) %>% 
  dplyr::summarise(mean= mean(solidago_below),
                   se = sd(solidago_below)/sqrt(n()))

total.s <- data %>% drop_na(solidago_total)
SOAL.total <- total.s %>% 
  group_by(treatment, parasite) %>% 
  dplyr::summarise(mean= mean(solidago_total),
                   se = sd(solidago_total)/sqrt(n()))

AGB.a <- data %>% drop_na(agalinis_above)
AGPU.AGB <- AGB.a %>% 
  group_by(treatment) %>% 
  dplyr::summarise(mean= mean(agalinis_above),
                   se = sd(agalinis_above)/sqrt(n()))

total.a <- data %>% drop_na(agalinis_total)
AGPU.total <- AGB.a %>% 
  group_by(treatment) %>% 
  dplyr::summarise(mean= mean(agalinis_total),
                   se = sd(agalinis_total)/sqrt(n()))

growth.s <- data %>% drop_na(solidago_final_height)
SOAL.growth <- growth.s %>% 
  group_by(treatment, parasite) %>% 
  dplyr::summarise(mean= mean(solidago_final_height),
                   se = sd(solidago_final_height)/sqrt(n()))

growth.a <- data %>% drop_na(agalinis_final_height)
AGPU.growth <- growth.a %>% 
  group_by(treatment) %>% 
  dplyr::summarise(mean= mean(agalinis_final_height),
                   se = sd(agalinis_final_height)/sqrt(n()))

rs.s <- data %>% drop_na(solidago_rs)
SOAL.rootshoot <- rs.s %>% 
  group_by(treatment, parasite) %>% 
  dplyr::summarise(mean= mean(solidago_rs),
                   se = sd(solidago_rs)/sqrt(n()))

rs.a <- data %>% drop_na(agalinis_rs)
AGPU.rootshoot<- rs.a %>% 
  group_by(treatment) %>% 
  dplyr::summarise(mean= mean(agalinis_rs),
                   se = sd(agalinis_rs)/sqrt(n()))

#Graphs
solidago.above.plot <- ggplot(data = SOAL.AGB, aes(x = reorder(treatment, mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~parasite) +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Solidago Aboveground Biomass") +
  theme(legend.position="none") +
  ylim(2,6)

solidago.above.plot


solidago.below.plot <- ggplot(data = SOAL.BGB, aes(x = reorder(treatment, mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~parasite) +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Solidago Belowground Biomass") +
  theme(legend.position="none") +
  ylim(2,6)

solidago.below.plot

solidago.total.plot <- ggplot(data = SOAL.total, aes(x = reorder(treatment, mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~parasite) +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Solidago Belowground Biomass") +
  theme(legend.position="none") +
  ylim(2,12)

solidago.total.plot


solidago.plots <- ggarrange(solidago.above.plot, solidago.below.plot,
                             labels = c("A", "B"), 
                             nrow = 2)
solidago.plots

agalinis.above.plot <- ggplot(data = AGPU.AGB, aes(x = reorder(treatment, -mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Agalinis Aboveground Biomass") +
  theme(legend.position="none") +
  ylim(0,0.075)

agalinis.above.plot

#Growth plots
solidago.growth.plot <- ggplot(data = SOAL.growth, aes(x = reorder(treatment, mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~parasite) +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Solidago Final Height") +
  theme(legend.position="none") +
  ylim(25,50)

solidago.growth.plot

agalinis.growth.plot <- ggplot(data = AGPU.growth, aes(x = reorder(treatment, -mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Agalinis Final Height") +
  theme(legend.position="none") +
  ylim(0,15)

agalinis.growth.plot

solidago.plots <- ggarrange(solidago.growth.plot, solidago.above.plot, solidago.below.plot,
                            labels = c("A", "B", "C"), 
                            nrow = 3)
solidago.plots


agalinis.plots <- ggarrange(agalinis.growth.plot, agalinis.above.plot,
                            labels = c("A", "B"), 
                            nrow = 2)
agalinis.plots


solidago.rs.plot <- ggplot(data = SOAL.rootshoot, aes(x = reorder(treatment, mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~parasite) +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Solidago Root/Shoot") +
  theme(legend.position="none") +
  ylim(0,1.5)
solidago.rs.plot

agalinis.rs.plot <- ggplot(data = AGPU.rootshoot, aes(x = reorder(treatment, -mean), y = mean, color = treatment)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  scale_color_manual( values=c("#40b0a6", "#E1BE6A")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Innoculumn", y = "Agalinis Root/Shoot") +
  theme(legend.position="none") +
  ylim(0,0.2)

agalinis.rs.plot
