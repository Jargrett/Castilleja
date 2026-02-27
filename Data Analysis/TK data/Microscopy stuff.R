setwd("~/Desktop/Castilleja/Data Analysis/TK data")

library(tidyverse)#for data wrangling and restructuring
library(statmod) 
library(lme4)
library(emmeans) 
library(car) 
library(ggpubr)
library(ggplot2)
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages

microscopy <- read.csv("HPM Microscopy - AMF.csv")
microscopy <- as.data.frame(unclass(microscopy),stringsAsFactors = TRUE)
str(microscopy)


#Filter the whole dataset into sub dataframes for slide type
parasite_comp <- microscopy %>% 
  filter(root_type == "parasite") %>% 
  drop_na(views)

host_comp <- microscopy %>% 
  filter(root_type == "host") %>% 
  drop_na(views)


#How does sterilization influence fungal presence
parasite_mean <- parasite_comp %>% 
  #filter(prop_total != 0) %>% 
  group_by(treatment, type) %>% 
  summarise(mean = mean(prop_hyphae), se = sd(prop_hyphae)/sqrt(n()))

host_mean <- host_comp %>% 
  #filter(prop_total != 0) %>% 
  group_by(treatment, type) %>% 
  summarise(mean = mean(prop_hyphae), se = sd(prop_hyphae)/sqrt(n()))

parasite_innoc <- ggplot(parasite_mean, aes(x = type, y = mean, fill = treatment)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Prop colonized") +
  scale_fill_manual(values = c("#4b3b40","#b6ad90")) +
  theme(panel.background = element_rect(fill = 'transparent'), #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA), #transparent plot bg
        legend.background = element_rect(fill = 'transparent'), #transparent legend bg
        legend.box.background = element_rect(fill = 'transparent')) + #transparent legend pane
  theme_pubr()
parasite_innoc

host_innoc <- ggplot(host_mean, aes(x = type, y = mean, fill = treatment)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Prop colonized") +
  scale_fill_manual(values = c("#45463e","#c6ac8f")) +
  theme(panel.background = element_rect(fill = 'transparent'), #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA), #transparent plot bg
        legend.background = element_rect(fill = 'transparent'), #transparent legend bg
        legend.box.background = element_rect(fill = 'transparent')) + #transparent legend pane
  theme_pubr()
host_innoc

#combines plots together
colonization.plots <- ggarrange(parasite_innoc, host_innoc,
                           labels = c("A", "B"), 
                           nrow = 1, ncol = 2)
colonization.plots



