setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and resctructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)
library(statmod)
library(lme4)
library(emmeans) # for comparison of means
library(betareg)
library(ggplot2)
library(ggpubr)#extended functions for plottinglibrary(remotes)
library(ggpattern)


#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::filter)

#import and restructure  biomass data (raw) 
biomass <- read.csv("Raw Data/EL Biomass - Biomass.csv")
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)
biomass$pair <- as.factor(biomass$pair)
biomass$plot <- as.factor(biomass$plot)
biomass$block <- as.factor(biomass$block)
str(biomass)

#----------Biomass Analysis----------#
#total plant community biomass (Castilleja Excluded)
total.bio.lmm <- lmer(total_no_cas ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(total.bio.lmm)
Anova(total.bio.lmm)# No significance 
emmeans(total.bio.lmm, pairwise ~ litter|removal)

#Graminoid Biomass 
gram.lmm <- lmer(gram_final ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(gram.lmm)
Anova(gram.lmm)# No significance 
emmeans(gram.lmm, pairwise ~ litter|removal)

#Graminoid Biomass Proportion
gram.prop.lmm <- lmer(gram_prop ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(gram.prop.lmm)
Anova(gram.prop.lmm)# No significance (removal: p = 0.067)
emmeans(gram.prop.lmm, pairwise ~ litter|removal)

#Forb Biomass 
forb.lmm <- lmer(forb_final ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(forb.lmm)
Anova(forb.lmm)# No significance 
emmeans(forb.lmm, pairwise ~ litter|removal)

#Forb Biomass Proportion
forb.prop.lmm <- lmer(forb_prop ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(forb.prop.lmm)
Anova(forb.prop.lmm)# No significance
emmeans(forb.prop.lmm, pairwise ~ litter|removal)

#Legume Biomass 
leg.lmm <- lmer(leg_final ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(leg.lmm)
Anova(leg.lmm)# No significance (litter p = 0.082, removal p = 0.098)
emmeans(leg.lmm, pairwise ~ litter|removal)

#Legume Biomass Proportion
leg.prop.lmm <- lmer(leg_prop ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(leg.prop.lmm)
Anova(leg.prop.lmm)# No significance
emmeans(leg.prop.lmm, pairwise ~ litter|removal)

#Shrub Biomass 
shurb.lmm <- lmer(shrub_final ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(shurb.lmm)
Anova(shurb.lmm)# Litter p = 0.027
emmeans(shurb.lmm, pairwise ~ removal|litter)

#Shrub Biomass Proportion
shrub.prop.lmm <- lmer(shrub_prop ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(shrub.prop.lmm)
Anova(shrub.prop.lmm)# Litter p = 0.020
emmeans(shrub.prop.lmm, pairwise ~ litter|removal)


biomass.long <- biomass %>%
  pivot_longer(
    cols = ends_with("_prop"),      # columns to pivot
    names_to = "func",  # new column for former col names
    values_to = "biomass_prop"        # new column for the numeric values
  )

biomass.long <- biomass.long %>% 
  select(-c(9:22))

#total biomass graphs
func.biomass <- biomass.long %>% 
  group_by(removal,func) %>% 
  dplyr::summarise(mean = mean(biomass_prop),
                   se = sd(biomass_prop)/sqrt(n()))

total.biomass <- biomass %>% 
  group_by(removal) %>% 
  dplyr::summarise(mean = mean(total_no_cas),
                   se = sd(total_no_cas)/sqrt(n()))

total.biomass %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "R", "stripe", "none")
  )

ggplot(total.biomass, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.8, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),color = "black", width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Total Biomass (No Castilleja)")

ggplot(biomass, aes(x = removal, y = total_no_cas,fill = removal, color = removal)) +
  geom_boxplot( lwd = 0.7, outlier.shape = NA,
              position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                      dodge.width = 0.10),
                                      alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c("#333d29", "#4A3D21")) +
  scale_fill_manual(values = c("#c5c6af", "#D3BC8D")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  ylim(0,425) +
  labs(x = "Parasite", y = "Total Biomass (No Castilleja)")
