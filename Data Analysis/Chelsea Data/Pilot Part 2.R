#This sets our working directory so that we have access to our files
setwd("~/Desktop/Castilleja/Data Analysis/Chelsea Data")

#The goal of this is to assess how light influences biomas growth for host and parasite

#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
library(car)
library(lme4)
library(emmeans)# for comparison of means
library(ggpubr)
library(ggplot2)

#conflicts
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

#import data
full_data <- read.csv("Pilot - Total.csv")
full_data <- as.data.frame(unclass(full_data),stringsAsFactors=TRUE)
full_data %<>% drop_na(biomass)

#Filter datasets by species so that we can conduct analysis and graph
cain_data <- filter(full_data, plant == "CAIN")#parasite data
cain_data %<>% mutate(pot_type = recode(pot_type,
                          "HP" = "with host",
                          "P" = "alone"))
scsc_data <- filter(full_data, plant == "SCSC")#host data
scsc_data %<>% mutate(pot_type = recode(pot_type,
                                        "HP" = "with parasite",
                                        "H" = "alone"))
cain_above <- filter(cain_data, bio_type == "AGB")#parasite data
cain_below <- filter(cain_data, bio_type == "BGB")#parasite data


#remove NA values from rows and columns

#----------Graphing----------#
cain_above <- filter(cain_data, bio_type == "AGB")#parasite data
cain_below <- filter(cain_data, bio_type == "BGB")#parasite data

scsc_above <- filter(scsc_data, bio_type == "AGB")#parasite data
scsc_below <- filter(scsc_data, bio_type == "BGB")#parasite data

cain_above_mean <- cain_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

cain_below_mean <- cain_below %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

scsc_above_mean <- scsc_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

scsc_below_mean <- scsc_below %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

#We will graph Paraite biomass by our two treatments (pot_type & light)
cain.above <- ggplot(cain_above_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Aboveground Biomass") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()
cain.above

cain.below <- ggplot(cain_below_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Belowground Biomass") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()
cain.below

scsc.above <- ggplot(scsc_above_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Aboveground Biomass") +
  scale_fill_manual(values=c("#45463e","#b3b1be")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()
scsc.above

scsc.below <- ggplot(scsc_below_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Belowground Biomass") +
  scale_fill_manual(values=c("#45463e","#b3b1be")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()
scsc.below

biomass.plots <- ggarrange(cain.above, scsc.above, cain.below, scsc.below,
                           labels = c("A", "B","C", "D"), 
                           nrow = 2, ncol = 2)
biomass.plots
