#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/HPM")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(ggpubr)#post-hoc analysis
library(ggplot2)#post-hoc analysis
library(patchwork)

nutrients <- read.csv("HPM Nutrient - Tissue.csv")
nutrients$treatment[nutrients$treatment == "innoculated"] <- "AMF"
nutrients$treatment[nutrients$treatment == "sterilized"] <- "Control"
nutrients <- as.data.frame(unclass(nutrients),stringsAsFactors=TRUE)
summary(nutrients)

nutrients %<>% 
  drop_na(N) %>% 
  mutate(CN = C / N) %>% 
  mutate(replicate_id = as.factor(str_extract(pot_id, "^\\d+")))

hesu_above <- nutrients %>% filter(species == "HESU", biomass_type == "Above")
hesu_above <- nutrients %>% filter(species == "HESU", biomass_type == "Below")

#HESU Analysis
#Above
hesu.cn.above <- lmer(CN ~ type*treatment + (1|replicate_id), data = hesu_above)
summary(hesu.cn.above)
Anova(hesu.cn.above)# pot_type p < 0.0001
emmeans(hesu.cn.above, pairwise ~ type|treatment)
emmip(hesu.cn.above, ~ type ~ treatment) 
# Growing with a parasite increases above ground CN ratio

hesu.c.above <- lmer(C ~ type*treatment + (1|replicate_id), data = hesu_above)
summary(hesu.c.above)
Anova(hesu.c.above) # AMF p = 0.0074
emmeans(hesu.c.above, pairwise ~ type|treatment)
emmip(hesu.c.above, ~ type ~ treatment)
# AMF presence decreases %C 

hesu.n.above <- lmer(N ~ type*treatment + (1|replicate_id), data = hesu_above)
summary(hesu.n.above)
Anova(hesu.n.above) # pot_type p < 0.0001
emmeans(hesu.n.above, pairwise ~ type|treatment)
emmip(hesu.n.above, ~ type ~ treatment)
# Growing with Agalinis reduces host foliar N

hesu.p.above <- lmer(P ~ type*treatment + (1|replicate_id), data = hesu_above)
summary(hesu.p.above)
Anova(hesu.p.above) # pot_type p < 0.0001, AMF p < 0.0001
emmeans(hesu.p.above, pairwise ~ type|treatment)
emmip(hesu.p.above, ~ type ~ treatment)
# AMF increases above ground P and Host reduces total P (No interaction)

hesu.k.above <- lmer(K ~ type*treatment + (1|replicate_id), data = hesu_above)
summary(hesu.k.above)
Anova(hesu.k.above) # pot_type p < 0.0001
emmeans(hesu.k.above, pairwise ~ type|treatment)
emmip(hesu.k.above, ~ type ~ treatment)
#Agalinis reduces host %K aboveground

#Below
hesu.cn.below <- lmer(CN ~ type*treatment + (1|replicate_id), data = hesu_below)
summary(hesu.cn.below)
Anova(hesu.cn.below) # pot_type:AMF p = 0.0134
emmeans(hesu.cn.below, pairwise ~ type|treatment)
emmip(hesu.cn.below, ~ type ~ treatment)
#AMF presence leads to differentiation between root CN with Agalinis attachment leading to lower CN ratio

hesu.c.below <- lmer(C ~ type*treatment + (1|replicate_id), data = hesu_below)
summary(hesu.c.below)
Anova(hesu.c.below) # pot_type:AMF p = 0.0018
emmeans(hesu.c.below, pairwise ~ type|treatment)
emmip(hesu.c.below, ~ type ~ treatment)
#AMF presence leads to differentiation between root C with Agalinis attachment leading to lower C

hesu.n.below <- lmer(N ~ type*treatment + (1|replicate_id), data = hesu_below)
summary(hesu.n.below)
Anova(hesu.n.below) # pot_type:AMF p = 0.048
emmeans(hesu.n.below, pairwise ~ type|treatment)
emmip(hesu.n.below, ~ type ~ treatment)
#AMF presence leads to differentiation between root N with Agalinis attachment leading to increased N

hesu.p.below <- lmer(P ~ type*treatment + (1|replicate_id), data = hesu_below)
summary(hesu.p.below)
Anova(hesu.p.below) # pot_type p = 0.0054, AMF p < 0.0001
emmeans(hesu.p.below, pairwise ~ type|treatment)
emmip(hesu.p.below, ~ type ~ treatment)
#AMF presence increases root P. Agalinis attachement leads to a decrease in root P

hesu.k.below <- lmer(K ~ type*treatment + (1|replicate_id), data = hesu_below)
summary(hesu.k.below)
Anova(hesu.k.below) # pot_type p = 0.0098,
emmeans(hesu.k.below, pairwise ~ type|treatment)
emmip(hesu.k.below, ~ type ~ treatment)
#Agalinis attachment leads to increased K in host roots

#AGPU Analysis

agpu_above <- nutrients %>% filter(species == "AGPU", biomass_type == "Above")
agpu_below <- nutrients %>% filter(species == "AGPU", biomass_type == "Below")

#Above
agpu.cn.above <- lmer(CN ~ type*treatment + (1|replicate_id), data = agpu_above)
summary(agpu.cn.above)
Anova(agpu.cn.above) # pot_type p < 0.0001, AMF p = 0.0398
emmeans(agpu.cn.above, pairwise ~ type|treatment)
emmip(agpu.cn.above, ~ type ~ treatment)
#Agalinis has a higher Shoot CN ratio when attached to host plants. AMF increases Agalinis Shoot CN ratio.

agpu.c.above <- lmer(C ~ type*treatment + (1|replicate_id), data = agpu_above)
summary(agpu.c.above)
Anova(agpu.c.above) # pot_type:AMF p = 0.0083
emmeans(agpu.c.above, pairwise ~ type|treatment)
emmip(agpu.c.above, ~ type ~ treatment)
#host attachment lead to lower Shoot C. This effect is exacerbated when AMF is present

agpu.n.above <- lmer(N ~ type*treatment + (1|replicate_id), data = agpu_above)
summary(agpu.n.above)
Anova(agpu.n.above) # pot_type p < 0.0001, AMF p = 0.0165
emmeans(agpu.n.above, pairwise ~ type|treatment)
emmip(agpu.n.above, ~ type ~ treatment)
# Host attachment leads to lower Shoot N. AMF presence reduces Agalinis Shoot N


#Below
agpu.cn.below <- lmer(CN ~ type*treatment + (1|replicate_id), data = agpu_below)
summary(agpu.cn.below)
Anova(agpu.cn.below) # pot_type p < 0.0001
emmeans(agpu.cn.below, pairwise ~ type|treatment)
emmip(agpu.cn.below, ~ type ~ treatment)
#Agalinis has higher Root CN when grown with a host.

agpu.c.below <- lmer(C ~ type*treatment + (1|replicate_id), data = agpu_below)
summary(agpu.c.below)
Anova(agpu.c.below) # pot_type:AMF p = 0.021
emmeans(agpu.c.below, pairwise ~ type|treatment)
emmip(agpu.c.below, ~ type ~ treatment)
# The effect of AMF on Agalinis Root C is dependent upon host presence. WHen grown alone AMF reduces Agalinis Root C, however when grown with a host AMF presence increases Root C.

agpu.n.below <- lmer(N ~ type*treatment + (1|replicate_id), data = agpu_below)
summary(agpu.n.below)
Anova(agpu.n.below) # pot_type:AMF p = 0.000013
emmeans(agpu.n.below, pairwise ~ type|treatment)
emmip(agpu.n.below, ~ type ~ treatment)
#Agalinis root N AMF reduces the magnitude of difference (Increased N with Host, decreased N alone) in root N between Agalinis grown alone and with a host.

agpu_above_P <- nutrients %>% 
  filter(species == "AGPU", biomass_type == "Above") %>% 
  drop_na(P)

agpu_below_P <- nutrients %>% 
  filter(species == "AGPU", biomass_type == "Below") %>% 
  drop_na(P)

agpu_above_K <- nutrients %>% 
  filter(species == "AGPU", biomass_type == "Above") %>% 
  drop_na(K)

agpu_below_K <- nutrients %>% 
  filter(species == "AGPU", biomass_type == "Below") %>% 
  drop_na(K)

#graphing
CN_summary <- nutrients %>%
  drop_na(CN) %>%
  group_by(species, treatment, type, biomass_type) %>%
  summarise(mean = mean(CN), se = sd(CN)/sqrt(n()), .groups = "drop") %>%
  mutate(nutrient = "C:N")

C_summary <- nutrients %>%
  drop_na(C) %>%
  group_by(species, treatment, type, biomass_type) %>%
  summarise(mean = mean(C), se = sd(C)/sqrt(n()), .groups = "drop") %>%
  mutate(nutrient = "C")

N_summary <- nutrients %>%
  drop_na(N) %>%
  group_by(species, treatment, type, biomass_type) %>%
  summarise(mean = mean(N), se = sd(N)/sqrt(n()), .groups = "drop") %>%
  mutate(nutrient = "N")

P_summary <- nutrients %>%
  drop_na(P) %>%
  group_by(species, treatment, type, biomass_type) %>%
  summarise(mean = mean(P), se = sd(P)/sqrt(n()), .groups = "drop") %>%
  mutate(nutrient = "P")

K_summary <- nutrients %>%
  drop_na(K) %>%
  group_by(species, treatment, type, biomass_type) %>%
  summarise(mean = mean(K), se = sd(K)/sqrt(n()), .groups = "drop") %>%
  mutate(nutrient = "K")

plotting_mean <- bind_rows(CN_summary, C_summary, N_summary, P_summary, K_summary) %>%
  mutate(
    nutrient    = factor(nutrient,    levels = c("C:N", "C", "N", "P", "K")),
    biomass_type = factor(biomass_type, levels = c("Above", "Below")),
    species     = factor(species,     levels = c("HESU", "AGPU")),
    lower = mean - se,
    upper = mean + se)

