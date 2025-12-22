#Set the working directory
setwd("~/Desktop/Castilleja/Data Analysis/Host Parasite Microbe Project")

#----------Packages----------#
#data cleaning and structuring
library(tidyverse)
library(conflicted)
library(magrittr)
library(hablar)
#Statistical Analysis
library(statmod)
library(lme4)
library(emmeans)
library(car)
library(ggstatsplot)
#Figures and Tables
library(ggpubr)
library(ggplot2)


#----------Data Working----------#
#Reading in the data from a downloaded csv (text) file
microscopy <- read.csv("Raw Data/HPM Microscopy - AMF.csv")
#Removing all slides that have not been measured and removing parasite alone slides
micro <- microscopy %>%  
  drop_na(views) %>% 
  filter(type != "parasite")
#checking data structure
str(micro)
#Converting strings to factors for analaysis
micro <- as.data.frame(unclass(micro),stringsAsFactors=TRUE)
micro %<>% convert(hablar::num(col_total)) %>% 
  convert(hablar::num(col_hyphae)) %>% 
  convert(hablar::num(col_extra)) %>% 
  convert(hablar::num(col_vesicle)) %>% 
  convert(hablar::num(col_arb)) %>% 
  convert(hablar::num(col_spore)) %>% 
  convert(hablar::num(col_dse)) %>% 
  convert(hablar::num(col_microsclerotia)) 


#----------Analysis----------#
#Effect of innoculuation on colonization
col.lmm <- lmer(col_total ~ treatment*type + (1|replicate_id), data = micro)
summary(col.lmm)
Anova(col.lmm) #Treatment p < 0.001
emmeans(col.lmm, pairwise ~ type|treatment)
emmip(col.lmm, ~ type ~ treatment)

col.plot <- ggbetweenstats(data = micro, x = treatment, y = col_total) 
col.plot

hyphae.lmm <- lmer(col_hyphae ~ treatment*type + (1|replicate_id), data = micro)
summary(hyphae.lmm)
Anova(hyphae.lmm)
emmeans(hyphae.lmm, pairwise ~ treatment|type)
emmip(hyphae.lmm, ~ type ~ treatment)

hyphae.plot <- ggbetweenstats(data = micro, x = treatment, y = col_hyphae) 
hyphae.plot

extra.lmm <- lmer(col_extra ~ treatment*type + (1|replicate_id), data = micro)
summary(extra.lmm)
Anova(extra.lmm)
emmeans(extra.lmm, pairwise ~ treatment|type)
emmip(extra.lmm, ~ type ~ treatment)

extra.plot <- ggbetweenstats(data = micro, x = treatment, y = col_extra) 
extra.plot

vesicle.lmm <- lmer(col_vesicle ~ treatment*type + (1|replicate_id), data = micro)
summary(vesicle.lmm)
Anova(vesicle.lmm)
emmeans(vesicle.lmm, pairwise ~ treatment|type)
emmeans(vesicle.lmm, pairwise ~ type|treatment)
emmip(vesicle.lmm, ~ type ~ treatment)

vesicle.plot <- ggbetweenstats(data = micro, x = treatment, y = col_vesicle) 
vesicle.plot

arb.lmm <- lmer(col_arb ~ treatment*type + (1|replicate_id), data = micro)
summary(arb.lmm)
Anova(arb.lmm)
emmeans(arb.lmm, pairwise ~ treatment|type)
emmip(arb.lmm, ~ type ~ treatment)

arb.plot <- ggbetweenstats(data = micro, x = treatment, y = col_arb) 
arb.plot

spore.lmm <- lmer(col_spore ~ treatment*type + (1|replicate_id), data = micro)
summary(spore.lmm)
Anova(spore.lmm)
emmeans(spore.lmm, pairwise ~ treatment|type)
emmip(spore.lmm, ~ type ~ treatment)

spore.plot <- ggbetweenstats(data = micro, x = type, y = col_spore) 
spore.plot

dse.lmm <- lmer(col_dse ~ treatment*type + (1|replicate_id), data = micro)
summary(dse.lmm)
Anova(dse.lmm)
emmeans(dse.lmm, pairwise ~ treatment|type)
emmip(dse.lmm, ~ type ~ treatment)

dse.plot <- ggbetweenstats(data = micro, x = treatment, y = col_dse) 
dse.plot


#Standard error calculations
total.col <- micro %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(col_total),
                   se = sd(col_total)/sqrt(n()))

#graphs
colonization.plot <- ggplot(total.col, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Total Colonization") +
  scale_color_manual( values=c("#994F00", "#7C99B3")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,0.5)
colonization.plot
