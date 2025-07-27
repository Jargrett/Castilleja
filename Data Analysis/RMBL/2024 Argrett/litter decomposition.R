#Litter Decomposition
#Initiated 8/23/24

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

#remotes::install_github("cornwell-lab-unsw/litterfitter")

library(litterfitter)#for k-curve fitting
library(tidyverse)#for data wrangling and restructuring
library(dplyr)#for data wrangling and restructuring
library(lme4)#for modeling linear mixed effect models
library(nlme)#alternative for modeling linear mixed effect models
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(dplyr)
library(plyr)

#load-in the data
winter1 <- read.csv("Litter Decomposition - 2024 Overwinter.csv")
winter2 <- read.csv("Litter Decomposition - 2025 Overwinter.csv")
within <- read.csv("Litter Decomposition - 2024 Within Season.csv")
yearlong <- read.csv("Litter Decomposition - 2024 Full Year.csv")


#bind

winter1 <- winter1 %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing))) %>%
  dplyr::select(-c((redo_coin_litter_dry_weight))) %>% 
  dplyr::select(-c((diff))) 

winter2 <- winter2 %>% 
  filter(ups_missing != "Yes") %>% 
  dplyr::select(-c((ups_missing))) %>% 
  dplyr::select(-c((missing)))

winter <- rbind.fill(winter1, winter2)
winter <- winter %>% drop_na(mass_loss)

decomp <- rbind.fill(within, winter, yearlong)
#check structure
str(decomp)
decomp <- as.data.frame(unclass(decomp),stringsAsFactors=TRUE)
within <- as.data.frame(unclass(within),stringsAsFactors=TRUE)
winter <- as.data.frame(unclass(winter),stringsAsFactors=TRUE)
yearlong <- as.data.frame(unclass(yearlong),stringsAsFactors=TRUE)
#yearlong <- as.data.frame(unclass(yearlong),stringsAsFactors=TRUE)
#remove missing bags
decomp <- decomp %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

within <- within %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

yearlong <- yearlong %>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

decomp <- decomp %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
within <- within %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
winter <- winter %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)
yearlong <- yearlong %>% mutate(mass_remaining = final_dry_weight/initial_dry_weight)

decomp <- decomp %>% mutate(time = deployment_duration/365)
within <- within %>% mutate(time = deployment_duration/365)
winter <- winter %>% mutate(time = deployment_duration/365)
yearlong <- yearlong %>% mutate(time = deployment_duration/365)


#first lets see if decomp is different between treatment

#Removing Outliers for testing
decomp <- decomp[-c(57,69,157),]

over.lm <- lmer(mass_remaining ~ removal*litter + year + (1|block) + (1|location), data = winter)
summary(over.lm)
Anova(over.lm)
emmip(over.lm, litter ~ removal)
emmeans(over.lm, pairwise ~ litter|removal)

litter <- fit_litter(time = decomp$time,
                     mass.remaining = decomp$mass_remaining,
                     model = "weibull",
                     iters=1000)

plot_multiple_fits(time = decomp$time,
                   mass.remaining = decomp$mass_remaining,
                   model=c("neg.exp","weibull"),
                   iters=500)

summary(litter)

class(litter)
plot(litter)

boxplot(decomp$mass_remaining ~ decomp$litter)

#Standard error calculations
decomp.mean <- decomp %>% 
  group_by(deployment_period, litter) %>% 
  dplyr::summarise(mean= mean(mass_remaining),
                   se = sd(mass_remaining)/sqrt(n()))


decomp.plot <- ggplot(data = decomp.mean, aes(x = reorder (deployment_period, -mean), y = mean, color = litter)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~litter) + 
  scale_color_manual( values=c("#004D40", "#C52812", "#FFC107")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Deployment Period", y = "Proportion of Mass Remaining") +
  ylim(0,1)

decomp.plot 

full.decomp.plot <- ggplot(data = decomp.mean, aes(x = reorder (deployment_period, -mean), y = mean, color = litter)) +
  geom_point(shape=18, size = 4,position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  scale_color_manual( values=c("#004D40", "#C52812", "#FFC107")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Deployment Period", y = "Proportion of Mass Remaining") +
  ylim(0.25 ,0.75)

full.decomp.plot 

