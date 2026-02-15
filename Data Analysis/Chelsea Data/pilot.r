setwd("~/Desktop/Castilleja/Data Analysis/Chelsea Data")

#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
library(car)
library(lme4)
library(emmeans) # for comparison of means
library(ggpubr)
library(ggplot2)

#important cover data (raw)
pilot <- read.csv("pilot.csv")
cain.pilot <- filter(pilot, species == "parasite")
scsc.pilot <- filter(pilot, species == "host")
cain.pilot <- na.omit(cain.pilot)
scsc.pilot <- na.omit(scsc.pilot)
#SCSC Final Height
scsc.pilot.lm <- lmer(height ~ light*pot.type + (1|replicate), data = scsc.pilot)
summary(scsc.pilot.lm)
Anova(scsc.pilot.lm) #p < 0.001
emmeans(scsc.pilot.lm, pairwise ~ light*pot.type)
emmip(scsc.pilot.lm, ~ light ~ pot.type)

#CAIN Final Height
cain.pilot.lm <- lmer(height ~ light*pot.type + (1|replicate), data = cain.pilot)
summary(cain.pilot.lm)
Anova(cain.pilot.lm) #p < 0.001
emmeans(cain.pilot.lm, pairwise ~ light*pot.type)
emmip(cain.pilot.lm, ~ light ~ pot.type)


scsc.mean <- scsc.pilot %>% 
  group_by(light,pot.type) %>% 
  dplyr::summarise(mean = mean(height),
                   se = sd(height)/sqrt(n()))

cain.mean <- cain.pilot %>% 
  group_by(light,pot.type) %>% 
  dplyr::summarise(mean = mean(height),
                   se = sd(height)/sqrt(n()))

#graphs

ggplot(cain.mean, aes(x=pot.type, y=mean, fill=func)) + 
  geom_bar(stat="identity") +
  labs(x = "Castilleja", y = "Proportion of Total Biomass") +
  scale_fill_manual(values=c("#582f0e", "#936639","#b6ad90","#656d4a","#414833" )) +
  theme_pubr()

ggplot(cain.mean, aes(x=pot.type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Height (cm)") +
  scale_fill_manual(values=c("#936639","#656d4a")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()

ggplot(scsc.mean, aes(x=pot.type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Height (cm)") +
  scale_fill_manual(values=c("#414833","#b6ad90")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()

