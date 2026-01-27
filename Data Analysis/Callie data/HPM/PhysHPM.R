#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/HPM")

#Question (1): How is the physiology of Agalinis effected by host attachment and innoculum?
#Question (2): How is the physiology of Heterotheca effected by host attachment and innoculum?
#Question (3): Does gas exchange differ between Host and parasites, both alone and within the same pot?

#load in relevant packages
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis
library(gasanalyzer)

physiology <- read.csv("HPM Physiology - Phys.csv")
str(physiology)
physiology$treatment[physiology$treatment == "innoculated"] <- "AMF"
physiology$treatment[physiology$treatment == "sterilized"] <- "Control"
physiology <- as.data.frame(unclass(physiology),stringsAsFactors=TRUE)


#Standard error calculations
phys <- physiology %>% drop_na()
ass.phys <- phys %>% 
  group_by(species, treatment, plant_id) %>% 
  dplyr::summarise(mean= mean(A),
                   se = sd(A)/sqrt(n()))

trans.phys <- phys %>% 
  group_by(species, treatment, plant_id) %>% 
  dplyr::summarise(mean= mean(E),
                   se = sd(E)/sqrt(n()))

sto.phys <- phys %>% 
  group_by(species, treatment, plant_id) %>% 
  dplyr::summarise(mean= mean(gsw),
                   se = sd(gsw)/sqrt(n()))

#graphs
ass.graph <- ggplot(ass.phys, aes(x = plant_id, y = mean, color = plant_id, group = plant_id)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 1, width = 0.09) +
  geom_point(aes(colour=plant_id),shape = 18, size = 5) +
  labs(x = "Plant Identity", y = "Carbon Assimilation (µmol m-2 s-1)") +
  scale_color_manual( values=c("#e07a5f", "#D6A839","#71A4A0", "#3d405b")) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(legend.position="none") +
  facet_wrap(~treatment) +
  ylim(0,10)

ass.graph

trans.graph <- ggplot(trans.phys, aes(x = plant_id, y = mean, color = plant_id)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 1, width = 0.09) +
  geom_point(aes(colour=plant_id),shape = 18, size = 5) +
  labs(x = "Plant Identity", y = "Transpiration rate (mol m-2 s-1)") +
  scale_color_manual( values=c("#e07a5f", "#D6A839","#71A4A0", "#3d405b")) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(legend.position="none") +
  facet_wrap(~treatment) +
  ylim(0,0.03)

trans.graph

sto.graph <- ggplot(sto.phys, aes(x = plant_id, y = mean, color = plant_id)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 1, width = 0.09) +
  geom_point(aes(colour=plant_id),shape = 18, size = 5) +
  labs(x = "Plant Identity", y = "Stomatal Conductance (mol m-2 s-1)") +
  scale_color_manual( values=c("#e07a5f", "#D6A839","#71A4A0", "#3d405b")) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(legend.position="none") +
  facet_wrap(~treatment) +
  ylim(0,1)


sto.graph

phys.plots <- ggarrange(ass.graph, trans.graph, sto.graph,
                             labels = c("A", "B","C"), 
                             nrow = 1)
phys.plots  
