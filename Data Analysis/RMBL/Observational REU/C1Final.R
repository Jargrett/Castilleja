setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")
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
library(ggpubr)

castilleja.cover <- read.csv("castilleja cover complete.csv")
castilleja.cover$castilleja[castilleja.cover$castilleja == "Control"] <- "Absent"
castilleja.cover$castilleja[castilleja.cover$castilleja == "Castilleja"] <- "Present"

#Diversity Analysis
div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(div)
Anova(div)
div.em <- emmip(div, ~ castilleja, plotit = FALSE)
emmeans(div, pairwise ~ castilleja)

rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(rich)
Anova(rich) 
emmip(rich, castilleja ~ year)
emmeans(rich, pairwise ~ castilleja|year)

even <- lmer(even ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(even)
Anova(even) 
emmip(even, castilleja ~ year)
emmeans(even, pairwise ~ castilleja|year)

#Standard error calculations
castilleja.div <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(div),
                   se = sd(div)/sqrt(n()))
castilleja.rich <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(rich),
                   se = sd(rich)/sqrt(n()))
castilleja.even <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(even),
                   se = sd(even)/sqrt(n()))

#Graphs
ggplot(data = castilleja.div, aes(x = reorder(castilleja, -mean), y = mean, fill = castilleja)) +
geom_point(aes (fill = castilleja), shape=18, size = 4) +
scale_fill_manual(values=c("#d8b365", "#5ab4ac")) 
  
  
div.plot <- ggplot(data = castilleja.div, aes(x = reorder(castilleja, -mean), y = mean)) +
  #geom_jitter(data = castilleja.cover, aes(x= castilleja, y = div), width = .15, alpha = .2) +
  stat_summary(fun=mean, colour="grey82", geom="line", aes(group = 1)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_fill_manual( values=c("#d8b365", "#5ab4ac")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
       color = "gray", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Shannon diversity of co-occuring species") +
  ylim(1,2.2)

div.plot

rich.plot <- ggplot(data = castilleja.rich, aes(x = reorder(castilleja, -mean), y = mean)) +
  geom_jitter(data = castilleja.cover, aes(x= castilleja,y = rich), width = .15, alpha = .2) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~year) +
  labs(x = "Castilleja", y = "Species Richness") +
  ylim(0,20) +
  theme_pubr()

rich.plot

even.plot <- ggplot(data = castilleja.even, aes(x = reorder(castilleja, -mean), y = mean)) +
  geom_jitter(data = castilleja.cover, aes(x= castilleja, y = even), width = .15, alpha = .2) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07,) +
  theme_pubr() +
  facet_wrap(~year) +
  labs(x = "Castilleja", y = "Species Evenness") +
  ylim(0.25,1) +
  theme_pubr()

even.plot
