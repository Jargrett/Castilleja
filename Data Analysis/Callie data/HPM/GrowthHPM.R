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

height <- read.csv("HPM Growth - Height.csv")
str(height)
height <- as.data.frame(unclass(height),stringsAsFactors=TRUE)

leaf <- read.csv("HPM Growth - Leaf_Number.csv")
str(leaf)
leaf <- as.data.frame(unclass(leaf),stringsAsFactors=TRUE)

#Height analysis
agalinis.height <- filter(height, species == "AGPU")
ag.height.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = agalinis.height)
summary(ag.height.lm)
Anova(ag.height.lm)
emmeans(ag.height.lm, pairwise ~ type|treatment)
emmip(ag.height.lm, ~ type ~ treatment)

hetero.height <- filter(height, species == "HESU")
he.height.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = hetero.height)
summary(he.height.lm)
Anova(he.height.lm)
emmeans(he.height.lm, pairwise ~ type|treatment)
emmip(he.height.lm, ~ type ~ treatment)

#leaf number analysis
agalinis.leaf <- filter(leaf, species == "AGPU")
ag.leaf.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = agalinis.leaf)
summary(ag.leaf.lm)
Anova(ag.leaf.lm)
emmeans(ag.leaf.lm, pairwise ~ type|treatment)
emmip(ag.leaf.lm, ~ type ~ treatment)

hetero.leaf <- filter(leaf, species == "HESU")
he.leaf.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = hetero.leaf)
summary(he.leaf.lm)
Anova(he.leaf.lm)
emmeans(he.leaf.lm, pairwise ~ type|treatment)
emmip(he.leaf.lm, ~ type ~ treatment)

#height graphs
height.long <- height %>% pivot_longer(cols=8:15, names_to = "week",values_to="height")
agalinis.height.long <- filter(height.long, species == "AGPU")
agalinis.height.long$type <- as.character(agalinis.height.long$type)
agalinis.height.long$type[agalinis.height.long$type == "host-parasite"] <- "With Host"
agalinis.height.long$type[agalinis.height.long$type == "parasite"] <- "Alone"
agalinis.height.long$type <- as.factor(agalinis.height.long$type)
str(agalinis.height.long)

#Standard error calculations
a.h.l <- agalinis.height.long %>% drop_na()
AGPU.Height <- a.h.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(height),
                   se = sd(height)/sqrt(n()))

ag.height <- ggplot(AGPU.Height, aes(x = week, y = mean, color = type, group = type)) +
  geom_point(aes(shape = type), position =  position_dodge(width = 0.5)) +
  geom_line(aes(linetype = type), position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
              position =  position_dodge(width = 0.5), width = 0.07) +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Time (Week)", y = "Agalinis Height (cm)") +
  theme_pubr() +
  facet_wrap(~treatment)
ag.height

#height graphs
height.long <- height %>% pivot_longer(cols=8:15, names_to = "week",values_to="height")
hetero.height.long <- filter(height.long, species == "HESU")
hetero.height.long$type <- as.character(hetero.height.long$type)
hetero.height.long$type[hetero.height.long$type == "host-parasite"] <- "With Parasite"
hetero.height.long$type[hetero.height.long$type == "host"] <- "Alone"
hetero.height.long$type <- as.factor(hetero.height.long$type)
str(hetero.height.long)

#Standard error calculations
h.h.l <- hetero.height.long %>% drop_na()
HESU.Height <- h.h.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(height),
                   se = sd(height)/sqrt(n()))

hes.height <- ggplot(HESU.Height, aes(x = week, y = mean, color = type, group = type)) +
  geom_point(aes(shape = type), position =  position_dodge(width = 0.5)) +
  geom_line(aes(linetype = type), position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Time (Week)", y = "Heterotheca Height (cm)") +
  theme_pubr() +
  facet_wrap(~treatment)
hes.height
#leaf graphs
leaf.long <- leaf %>% pivot_longer(cols=8:15, names_to = "week",values_to="leaves")
hetero.leaf.long <- filter(leaf.long, species == "HESU")
hetero.leaf.long$type <- as.character(hetero.leaf.long$type)
hetero.leaf.long$type[hetero.leaf.long$type == "host-parasite"] <- "With Parasite"
hetero.leaf.long$type[hetero.leaf.long$type == "host"] <- "Alone"
hetero.leaf.long$type <- as.factor(hetero.leaf.long$type)
str(hetero.leaf.long)

#Standard error calculations
h.l.l <- hetero.leaf.long %>% drop_na()
HESU.leaf <- h.l.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(leaves),
                   se = sd(leaves)/sqrt(n()))

hes.leaf <- ggplot(HESU.leaf, aes(x = week, y = mean, color = type, group = type)) +
  geom_point(aes(shape = type), position =  position_dodge(width = 0.5)) +
  geom_line(aes(linetype = type), position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  labs(x = "Time (Week)", y = "Heterotheca leaf number") +
  theme_pubr() +
  facet_wrap(~treatment)
hes.leaf

#leaf graphs
leaf.long <- leaf %>% pivot_longer(cols=8:15, names_to = "week",values_to="leaves")
agalinis.leaf.long <- filter(leaf.long, species == "AGPU")
agalinis.leaf.long$type <- as.character(agalinis.leaf.long$type)
agalinis.leaf.long$type[agalinis.leaf.long$type == "host-parasite"] <- "With Host"
agalinis.leaf.long$type[agalinis.leaf.long$type == "parasite"] <- "Alone"
agalinis.leaf.long$type <- as.factor(agalinis.leaf.long$type)
str(hetero.leaf.long)

#Standard error calculations
a.l.l <- agalinis.leaf.long %>% drop_na()
AGPU.leaf <- a.l.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(leaves),
                   se = sd(leaves)/sqrt(n()))

ag.leaf <- ggplot(AGPU.leaf, aes(x = week, y = mean, color = type, group = type)) +
  geom_point(aes(shape = type), position =  position_dodge(width = 0.5)) +
  geom_line(aes(linetype = type), position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  labs(x = "Time (Week)", y = "Agalinis leaf number") +
  theme_pubr() +
  facet_wrap(~treatment)
ag.leaf

growth.plots <- ggarrange(ag.height, ag.leaf, hes.height, hes.leaf,
                            labels = c("A", "B","C","D"), 
                            nrow = 2, ncol = 2)
growth.plots
