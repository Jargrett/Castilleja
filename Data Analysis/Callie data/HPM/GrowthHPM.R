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
height$treatment[height$treatment == "innoculated"] <- "AMF"
height$treatment[height$treatment == "sterilized"] <- "Control"
agalinis.height <- filter(height, species == "AGPU")
hetero.height <- filter(height, species == "HESU")
height <- as.data.frame(unclass(height),stringsAsFactors=TRUE)
hetero.height <- as.data.frame(unclass(hetero.height),stringsAsFactors=TRUE)
agalinis.height <- as.data.frame(unclass(agalinis.height),stringsAsFactors=TRUE)

leaf <- read.csv("HPM Growth - Leaf_Number.csv")
str(leaf)
leaf$treatment[leaf$treatment == "innoculated"] <- "AMF"
leaf$treatment[leaf$treatment == "sterilized"] <- "Control"
leaf <- as.data.frame(unclass(leaf),stringsAsFactors=TRUE)


flower <- read.csv("HPM Growth - Flowering.csv")
str(flower)
flower$treatment[flower$treatment == "innoculated"] <- "AMF"
flower$treatment[flower$treatment == "sterilized"] <- "Control"
flower <- as.data.frame(unclass(flower),stringsAsFactors=TRUE)

SLA<- read.csv("HPM SLA - SLA.csv")
str(SLA)
SLA$treatment[SLA$treatment == "innoculated"] <- "AMF"
SLA$treatment[SLA$treatment == "sterilized"] <- "Control"
SLA <- as.data.frame(unclass(SLA),stringsAsFactors=TRUE)

#Height analysis
agalinis.height <- filter(height, species == "AGPU")
ag.height.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = agalinis.height)
summary(ag.height.lm)
Anova(ag.height.lm) #treatment:type p = 0.029
emmeans(ag.height.lm, pairwise ~ type*treatment)
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

#final time point graphs
tiff("effect.tiff", units = "in", width = 10, length = 8, res = 300)
dev.off()

a.h <- agalinis.height %>% drop_na()
agpu.Height <- a.h %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
a.l <- agalinis.leaf %>% drop_na()
agpu.Leaf <- a.l %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
h.h <- hetero.height %>% drop_na()
agpu.Height <- h.h %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
h.l <- hetero.leaf %>% drop_na()
agpu.Leaf <- h.l %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))

t8.ag.height.graph <- ggplot(agpu.Height, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final Height (cm)") +
  scale_color_manual(values=c("#D6A839", "#71A4A0")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,50)

t8.ag.height.graph

t8.ag.leaf.graph <- ggplot(agpu., aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final Leaf Number") +
  scale_color_manual(values=c("#D6A839", "#71A4A0")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,300)

t8.ag.height.graph

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


#bud analysis
ag.flower.lm <- lmer(total_bud ~ treatment*type + (1|replicate_id), data = flower)
summary(ag.flower.lm)
Anova(ag.flower.lm) #Treatment by type interaction: Chisq: 5.806, p = 0.0159
emmeans(ag.flower.lm, pairwise ~ type*treatment)
emmip(ag.flower.lm, ~treatment ~ type)


#
flower$type <- as.character(flower$type)
flower$type[flower$type == "host-parasite"] <- "With Host"
flower$type[flower$type == "parasite"] <- "Alone"
flower$type <- as.factor(flower$type)

flower$treatment <- as.character(flower$treatment)
flower$treatment[flower$treatment == "innoculated"] <- "AMF"
flower$treatment[flower$treatment == "sterilized"] <- "Control"
flower$treatment <- as.factor(flower$treatment)

ag.bud <- flower %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(total_bud),
                   se = sd(total_bud)/sqrt(n()))

bud.graph <- ggplot(ag.bud, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Number of flower buds") +
  scale_color_manual(values=c("#D6A839", "#71A4A0")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,15)

bud.graph

#SLA analysis
agalinis.sla <- filter(SLA, species == "AGPU")
ag.sla.lm <- lmer(sla ~ treatment*type + (1|replicate_id), data = agalinis.sla)
summary(ag.sla.lm)
Anova(ag.sla.lm)
emmeans(ag.sla.lm, pairwise ~ type|treatment)
emmip(ag.sla.lm, ~ type ~ treatment)

hetero.sla <- filter(SLA, species == "HESU")
he.sla.lm <- lmer(sla ~ treatment*type + (1|replicate_id), data = hetero.sla)
summary(he.sla.lm)
Anova(he.sla.lm)
emmeans(he.sla.lm, pairwise ~ type|treatment)
emmip(he.sla.lm, ~ type ~ treatment)

#Standard error calculations
AGPU.sla <- agalinis.sla %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(sla),
                   se = sd(sla)/sqrt(n()))

HESU.sla <- hetero.sla %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(sla),
                   se = sd(sla)/sqrt(n()))

ag.SLA.graph <- ggplot(AGPU.sla, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Specific leaf Area (cm^2/g)") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(100,200)

ag.SLA.graph


he.SLA.graph <- ggplot(HESU.sla, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Specific leaf Area (cm^2/g)") +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(200,600)

he.SLA.graph

sla.plots <- ggarrange(ag.SLA.graph, he.SLA.graph,
                          labels = c("A", "B"), 
                          nrow = 1, ncol = 2)
sla.plots
