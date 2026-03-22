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
library(rstatix)
library(sjPlot)
library(ggpmisc)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarise)

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


#Biomass by colonization plot
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/HPM")
biomass <- read.csv("HPM Biomass - Biomass.csv")
str(biomass)
biomass$treatment[biomass$treatment == "innoculated"] <- "AMF"
biomass$treatment[biomass$treatment == "sterilized"] <- "Control"
agalinis.biomass <- filter(biomass, species == "AGPU")
agalinis.biomass$type[agalinis.biomass$type == "host-parasite"] <- "With Host"
agalinis.biomass$type[agalinis.biomass$type == "parasite"] <- "Alone"
hetero.biomass <- filter(biomass, species == "HESU")
hetero.biomass$type[hetero.biomass$type == "host-parasite"] <- "With Parasite"
hetero.biomass$type[hetero.biomass$type == "host"] <- "Alone"
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)
hetero.biomass <- as.data.frame(unclass(hetero.biomass),stringsAsFactors=TRUE)
agalinis.biomass <- as.data.frame(unclass(agalinis.biomass),stringsAsFactors=TRUE)

parasite_bio_myc <- full_join(agalinis.biomass,parasite_comp, by="pot_id")
host_bio_myc <- full_join(hetero.biomass,host_comp, by="pot_id")

parasite_scatter <- ggplot(parasite_bio_myc, aes(x = prop_hyphae, y = total_biomass, color = type.x)) +
  geom_point() +
  geom_smooth(method=lm , color="#472d30", fill = "cornsilk3", se = TRUE) + 
  scale_color_manual(values = c("#4b3b40","#b6ad90")) +
  labs(x = "Parasite Root Colonization (%)", y = "Parasite Total Biomass (g)") +
  stat_poly_eq(aes(label = after_stat(rr.label), group = type.x),
               formula = y ~ x,
               parse   = TRUE,
               size    = 3) +
  theme_pubr()
parasite_scatter


host_scatter <- ggplot(host_bio_myc, aes(x = prop_hyphae, y = total_biomass, color = type.x)) +
  geom_point() +
  geom_smooth(method=lm , color="#472d30", fill = "cornsilk3", se = TRUE) + 
  scale_color_manual(values = c("#4b3b40","#b6ad90")) +
  labs(x = "Host Root Colonization (%)", y = "Host Total Biomass (g)") +
  stat_poly_eq(aes(label = after_stat(rr.label), group = type.x),
               formula = y ~ x,
               parse   = TRUE,
               size    = 3) +
  theme_pubr()
host_scatter


#combines plots together
biomyc.plots <- ggarrange(parasite_scatter, host_scatter,
                                labels = c("A", "B"), 
                                nrow = 1, ncol = 2)
biomyc.plots

  
