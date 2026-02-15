setwd("~/Desktop/Castilleja/Data Analysis/TK data")

library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis
library(gasanalyzer)

microscopy <- read.csv("HPM Microscopy - AMF.csv")
microscopy <- as.data.frame(unclass(microscopy),stringsAsFactors = TRUE)
str(microscopy)

#Filter the whole dataset into sub dataframes for slide type
host_micro <- microscopy %>% filter(type == "host")
parasite_micro <- microscopy %>% filter(type == "parasite")
host_parasite_micro <- microscopy %>% filter(type == "host-parasite")

#How does sterilization influence fungal presence
parasite_mean <- parasite_micro %>% 
  filter(prop_total != 0) %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(prop_total), se = sd(prop_total)/sqrt(n()))

host_mean <- host_parasite_micro %>% 
  filter(prop_total != 0) %>% 
  group_by(treatment, root_type) %>% 
  summarise(mean = mean(prop_total), se = sd(prop_total)/sqrt(n()))

parasite_innoc <- ggplot(host_mean, aes(x = root_type, y = mean, fill = treatment)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Prop colonized") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90", "#45463e")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr()
parasite_innoc
