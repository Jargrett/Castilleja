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
emmeans(shurb.lmm, pairwise ~ litter|removal)

#Shrub Biomass Proportion
shrub.prop.lmm <- lmer(shrub_prop ~ litter*removal + (1|block) + (1|pair), data = biomass)
summary(shrub.prop.lmm)
Anova(shrub.prop.lmm)# No significance
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

ggplot(func.biomass, aes(x=removal, y=mean, fill=func)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Castilleja", y = "Proportion of Total Biomass") +
  scale_fill_manual(values=c("#582f0e", "#936639","#b6ad90","#656d4a","#414833" )) +
  theme_pubr()
  
  

