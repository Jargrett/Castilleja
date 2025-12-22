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

#-------------------------GROWTH----------------------------------#
#------Growth Data Importing and restructuring------#

#Height: import, clean and rename height data
height <- read.csv("HPM Growth - Height.csv")
str(height)
height$treatment[height$treatment == "innoculated"] <- "AMF"
height$treatment[height$treatment == "sterilized"] <- "Control"
agalinis.height <- filter(height, species == "AGPU")
agalinis.height$type[agalinis.height$type == "host-parasite"] <- "With Host"
agalinis.height$type[agalinis.height$type == "parasite"] <- "Alone"
hetero.height <- filter(height, species == "HESU")
hetero.height$type[hetero.height$type == "host-parasite"] <- "With Parasite"
hetero.height$type[hetero.height$type == "host"] <- "Alone"
height <- as.data.frame(unclass(height),stringsAsFactors=TRUE)
hetero.height <- as.data.frame(unclass(hetero.height),stringsAsFactors=TRUE)
agalinis.height <- as.data.frame(unclass(agalinis.height),stringsAsFactors=TRUE)

#Height over time data conversions

agalinis.height.long <- agalinis.height %>% pivot_longer(cols=8:15, names_to = "week",values_to = "height")
hetero.height.long <- hetero.height %>% pivot_longer(cols=8:15, names_to = "week",values_to = "height")
height.long <- rbind(agalinis.height.long, hetero.height.long)

#Height over time MLE logistic fitting
height.long$week <- as.numeric(gsub("t","",height.long$week))

fit_logistic_model <- function(height.long) {
  time <- height.long$week
  height <- height.long$height
  
  logistic_growth <- function(time,r,maxH) {
    N0 <- height[1]
    N <- N0 * maxH / (N0 + (maxH - N0) * exp(-r * time))
    return(N)
  }
  
  #MLE 
  model <- optim(par = c(r = 0.05, maxH = 40),
                 fn = function(par) sum((height - logistic_growth(time,par[1],par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  maxH_hat <- model$par[2]
  
  #predictions
  predicted_height <- logistic_growth(time,r_hat,maxH_hat)
  
  return(data.frame(species = unique(height.long$species),
                    type = unique(height.long$type),
                    treatment = unique(height.long$treatment),
                    week = time,
                    observed_height = height,
                    predicted_height = predicted_height,
                    r = r_hat,
                    maxH = maxH_hat,
                    N0 = height[1]))
}

results_height <- height.long %>%
  filter(!is.na(height)) %>%
  group_by(treatment,species,type) %>%
  do(fit_logistic_model(.))

results_height %>%
  mutate(unique_ID = paste(species,treatment,type,sep="_")) %>%
  ggplot(aes(x = week, y = predicted_height)) +
  geom_line() +
  facet_wrap(~unique_ID,scales="free") +
  geom_point(aes(x=week,y=observed_height),size=2,alpha=0.3) +
  labs(x="Week",y="Height") +
  theme_classic()


results_height2 <- height.long %>%
  filter(!is.na(height)) %>%
  group_by(treatment,species,type,replicate_id) %>%
  do(fit_logistic_model(.))

results_height2 %>%
  mutate(unique_ID = paste(species,treatment,type,replicate_id,sep="_")) %>%
  ggplot(aes(x = week, y = predicted_height)) +
  geom_line() +
  facet_wrap(~unique_ID,scales="free") +
  geom_point(aes(x=week,y=observed_height),size=2,alpha=0.3) +
  labs(x="Week",y="Height") +
  theme_classic()

#Leaf Number: import, clean and rename height data
leaf <- read.csv("HPM Growth - Leaf_Number.csv")
str(leaf)
leaf$treatment[leaf$treatment == "innoculated"] <- "AMF"
leaf$treatment[leaf$treatment == "sterilized"] <- "Control"
agalinis.leaf <- filter(leaf, species == "AGPU")
agalinis.leaf$type[agalinis.leaf$type == "host-parasite"] <- "With Host"
agalinis.leaf$type[agalinis.leaf$type == "parasite"] <- "Alone"
hetero.leaf <- filter(leaf, species == "HESU")
hetero.leaf$type[hetero.leaf$type == "host-parasite"] <- "With Parasite"
hetero.leaf$type[hetero.leaf$type == "host"] <- "Alone"
leaf <- as.data.frame(unclass(leaf),stringsAsFactors=TRUE)
hetero.leaf <- as.data.frame(unclass(hetero.leaf),stringsAsFactors=TRUE)
agalinis.leaf <- as.data.frame(unclass(agalinis.leaf),stringsAsFactors=TRUE)

#leaf number over time data conversions
leaf.long <- leaf %>% pivot_longer(cols=8:15, names_to = "week",values_to="height")
agalinis.leaf.long <- agalinis.leaf %>% pivot_longer(cols=8:15, names_to = "week",values_to="leaf")
hetero.leaf.long <- hetero.leaf %>% pivot_longer(cols=8:15, names_to = "week",values_to="leaf")

#Flower/Bud data: import, clean and rename bud data (Agalinis only)
flower <- read.csv("HPM Growth - Flowering.csv")
str(flower)
flower$treatment[flower$treatment == "innoculated"] <- "AMF"
flower$treatment[flower$treatment == "sterilized"] <- "Control"
flower$type[flower$type == "host-parasite"] <- "With Host"
flower$type[flower$type == "parasite"] <- "Alone"
flower <- as.data.frame(unclass(flower),stringsAsFactors=TRUE)

#Specific Leaf Area: import, clean and rename 
SLA<- read.csv("HPM SLA - SLA.csv")
str(SLA)
SLA$treatment[SLA$treatment == "innoculated"] <- "AMF"
SLA$treatment[SLA$treatment == "sterilized"] <- "Control"
agalinis.sla <- filter(SLA, species == "AGPU")
agalinis.sla$type[agalinis.sla$type == "host-parasite"] <- "With Host"
agalinis.sla$type[agalinis.sla$type == "parasite"] <- "Alone"
hetero.sla <- filter(SLA, species == "HESU")
hetero.sla$type[hetero.sla$type == "host-parasite"] <- "With Parasite"
hetero.sla$type[hetero.sla$type == "host"] <- "Alone"
SLA <- as.data.frame(unclass(SLA),stringsAsFactors=TRUE)
hetero.sla <- as.data.frame(unclass(hetero.sla),stringsAsFactors=TRUE)
agalinis.sla <- as.data.frame(unclass(agalinis.sla),stringsAsFactors=TRUE)

#------Growth Data Analysis------#

#Agalinis Final Height, Leaf number, Budding and SLA
#Eventually plan to add over time data analysis and visualization here as well

#AGPU Final Height
ag.height.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = agalinis.height)
summary(ag.height.lm)
Anova(ag.height.lm) #treatment:type  Chisq = 4.719, p = 0.029
emmeans(ag.height.lm, pairwise ~ type*treatment)
emmip(ag.height.lm, ~ type ~ treatment)

#AGPU Final Leaf Number
ag.leaf.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = agalinis.leaf)
summary(ag.leaf.lm)
Anova(ag.leaf.lm) #treatment:type  Chisq = 3.847, p = 0.049
emmeans(ag.leaf.lm, pairwise ~ type|treatment)
emmip(ag.leaf.lm, ~ type ~ treatment)

#AGPU Final Bud Count
ag.flower.lm <- lmer(total_bud ~ treatment*type + (1|replicate_id), data = flower)
summary(ag.flower.lm)
Anova(ag.flower.lm) #Treatment:type Chisq = 5.806, p = 0.0159
emmeans(ag.flower.lm, pairwise ~ type*treatment)
emmip(ag.flower.lm, ~treatment ~ type)

#AGPU SLA
ag.sla.lm <- lmer(sla ~ treatment*type + (1|replicate_id), data = agalinis.sla)
summary(ag.sla.lm)
Anova(ag.sla.lm)# treatment: Chisq = 3.097, p = 0.0785
emmeans(ag.sla.lm, pairwise ~ type|treatment)
emmip(ag.sla.lm, ~ type ~ treatment)

#Heterotheca Final Height, Leaf number and SLA
#Eventually plan to add over time data analysis and visualization here as well

#HESU Final Height
he.height.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = hetero.height)
summary(he.height.lm)
Anova(he.height.lm)
emmeans(he.height.lm, pairwise ~ type|treatment)
emmip(he.height.lm, ~ type ~ treatment)

#HESU Final Leaf Number
he.leaf.lm <- lmer(t8 ~ treatment*type + (1|replicate_id), data = hetero.leaf)
summary(he.leaf.lm)
Anova(he.leaf.lm) #Treatment:type Chisq = 3.847, p = 0.049
emmeans(he.leaf.lm, pairwise ~ type|treatment)
emmip(he.leaf.lm, ~ type ~ treatment)

#HESU SLA
he.sla.lm <- lmer(sla ~ treatment*type + (1|replicate_id), data = hetero.sla)
summary(he.sla.lm)
Anova(he.sla.lm) #Treatment:type Chisq = 3.916, p = 0.048
emmeans(he.sla.lm, pairwise ~ type|treatment)
emmip(he.sla.lm, ~ type ~ treatment)

#------Growth Data Standard error Calculations------#

#Agalinis Final time point calculations
#Height
a.h <- agalinis.height %>% drop_na()
agpu.height <- a.h %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
#height/time
a.h.l <- agalinis.height.long %>% drop_na()
agpu.height.time <- h.h.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(height),
                   se = sd(height)/sqrt(n()))
#Leaf Number
a.l <- agalinis.leaf %>% drop_na()
agpu.leaf <- a.l %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
#Flower
ag.bud <- flower %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(total_bud),
                   se = sd(total_bud)/sqrt(n()))
#SLA
ag.sla <- agalinis.sla %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(sla),
                   se = sd(sla)/sqrt(n()))

#Heterotheca Final time point calculations
#Height
h.h <- hetero.height %>% drop_na()
hesu.height <- h.h %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
#height/time
h.h.l <- hetero.height.long %>% drop_na()
hesu.height.time <- h.h.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(height),
                   se = sd(height)/sqrt(n()))
#Leaf Number
h.l <- hetero.leaf %>% drop_na()
hesu.Leaf <- h.l %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(t8),
                   se = sd(t8)/sqrt(n()))
#leaf number/time
h.l.l <- hetero.leaf.long %>% drop_na()
hesu.leaf.time <- h.h.l %>% 
  group_by(week, treatment, type) %>% 
  dplyr::summarise(mean= mean(leaf),
                   se = sd(leaf)/sqrt(n()))
#SLA
he.sla <- hetero.sla %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(sla),
                   se = sd(sla)/sqrt(n()))

#------Growth Graphs------#
#Unlike before this section is split by response variable
#For parasite: Alone graphs are coded in light blue and with host graphs in yellow
#For Host: Alone graphs are coded in dark blue and with host graphs are coded in orange

#Height Graphs
t8.ag.height.graph <- ggplot(agpu.height, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 6, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final Height (cm)") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,50)

t8.ag.height.graph

t8.he.height.graph <- ggplot(hesu.height, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 6, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final Height (cm)") +
  scale_color_manual( values=c("#3d405b", "#e07a5f")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,3)

t8.he.height.graph

#Leaf Number Graphs
t8.ag.leaf.graph <- ggplot(agpu.leaf, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 6, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final Leaf Number") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,300)

t8.ag.leaf.graph

t8.he.leaf.graph <- ggplot(hesu.leaf, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 6, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Final leaf Number") +
  scale_color_manual( values=c("#3d405b", "#e07a5f")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,50)

t8.he.leaf.graph

#SLA graphs
ag.SLA.graph <- ggplot(ag.sla, aes(x = treatment, y = mean, color = type)) +
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

he.SLA.graph <- ggplot(he.sla, aes(x = treatment, y = mean, color = type)) +
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

#bud count graphs
ag.bud.graph <- ggplot(ag.bud, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Bud count") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,20)

ag.bud.graph

ag.growth.plots <- ggarrange(t8.ag.height.graph, t8.ag.leaf.graph,
                          labels = c("A","B"), 
                          nrow = 1, ncol = 2, common.legend = TRUE)
ag.growth.plots

he.growth.plots <- ggarrange(t8.he.height.graph, t8.he.leaf.graph,
                             labels = c("C","D"), 
                             nrow = 1, ncol = 2, common.legend = TRUE)
he.growth.plots

growth.plots <- ggarrange(ag.growth.plots,he.growth.plots,
                          nrow = 2, ncol = 1)

growth.plots

ggsave(plot = growth.plots, filename = 'growth.png',
       width = 10 ,height = 8, units = "in", dpi = 600, 
       bg = "transparent")



#-------------------------BIOMASS----------------------------------#
#------Growth Data Importing and restructuring------#
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

#------Biomass Data Analysis------#

#Agalinis
#aboveground
ag.ag.lm <- lmer(above_total ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.ag.lm)
Anova(ag.ag.lm) #Treatment:type Chisq = 31.258, p < 0.0001
emmeans(ag.ag.lm, pairwise ~ type|treatment)
emmip(ag.ag.lm, ~ type ~ treatment)

#belowground
ag.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.bg.lm)
Anova(ag.bg.lm) #Treatment:type Chisq = 31.309, p < 0.0001
emmeans(ag.bg.lm, pairwise ~ type|treatment)
emmip(ag.bg.lm, ~ type ~ treatment)

#Total
ag.tb.lm <- lmer(total_biomass ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.tb.lm)
Anova(ag.tb.lm)
emmeans(ag.tb.lm, pairwise ~ type|treatment)
emmip(ag.tb.lm, ~ type ~ treatment)

#Root:shoot
ag.rs.lm <- lmer(root_shoot ~ treatment*type + (1|replicate_id), data = agalinis.biomass)
summary(ag.rs.lm)
Anova(ag.rs.lm)
emmeans(ag.rs.lm, pairwise ~ type|treatment)
emmip(ag.rs.lm, ~ type ~ treatment)

#Heterotheca
#aboveground
he.ag.lm <- lmer(above_total ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.ag.lm)
Anova(he.ag.lm)
emmeans(he.ag.lm, pairwise ~ type|treatment)
emmip(he.ag.lm, ~ type ~ treatment)

#belowground
he.bg.lm <- lmer(below_total ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.bg.lm)
Anova(he.bg.lm)
emmeans(he.bg.lm, pairwise ~ type|treatment)
emmip(he.bg.lm, ~ type ~ treatment)

#total
he.tb.lm <- lmer(total_biomass ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.tb.lm)
Anova(he.tb.lm)
emmeans(he.tb.lm, pairwise ~ type|treatment)
emmip(he.tb.lm, ~ type ~ treatment)

#Root:shoot
he.rs.lm <- lmer(root_shoot ~ treatment*type + (1|replicate_id), data = hetero.biomass)
summary(he.rs.lm)
Anova(he.rs.lm)
emmeans(he.rs.lm, pairwise ~ type|treatment)
emmip(he.rs.lm, ~ type ~ treatment)


#------Biomass Data Standard error Calculations------#

sd.ag.above <- agalinis.biomass %>% drop_na(above_total)
sd.ag.above <- sd.ag.above %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(above_total),
                   se = sd(above_total)/sqrt(n()))

sd.ag.below <- agalinis.biomass %>% drop_na(below_total)
sd.ag.below <- sd.ag.below %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(below_total),
                   se = sd(below_total)/sqrt(n()))

sd.ag.total <- agalinis.biomass %>% drop_na(total_biomass)
sd.ag.total <- sd.ag.total %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(total_biomass),
                   se = sd(total_biomass)/sqrt(n()))

sd.ag.rootshoot <- agalinis.biomass %>% drop_na(root_shoot)
sd.ag.rootshoot <- sd.ag.rootshoot %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(root_shoot),
                   se = sd(root_shoot)/sqrt(n()))

sd.he.above <- hetero.biomass %>% drop_na(above_total)
sd.he.above <- sd.he.above %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(above_total),
                   se = sd(above_total)/sqrt(n()))

sd.he.below <- hetero.biomass %>% drop_na(below_total)
sd.he.below <- sd.he.below %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(below_total),
                   se = sd(below_total)/sqrt(n()))

sd.he.total <- hetero.biomass %>% drop_na(total_biomass)
sd.he.total <- sd.he.total %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(total_biomass),
                   se = sd(total_biomass)/sqrt(n()))

sd.he.rootshoot <- hetero.biomass %>% drop_na(root_shoot)
sd.he.rootshoot <- sd.he.rootshoot %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean= mean(root_shoot),
                   se = sd(root_shoot)/sqrt(n()))

#------Biomass Graphs------#

hetero.above.plot <- ggplot(data = sd.he.above, aes(x = treatment, y = mean, fill = type)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  labs(x = "Fungi", y = "Aboveground Biomass (g)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,2.5)
hetero.above.plot 

hetero.below.plot <- ggplot(data = sd.he.below, aes(x = treatment, y = mean, fill = type)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  
  labs(x = "Fungi", y = "Belowground Biomass (g)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,2.5)
hetero.below.plot 

hetero.total.plot <- ggplot(data = sd.he.total, aes(x = treatment, y = mean, fill = type)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#3d405b", "#e07a5f")) +
  
  labs(x = "Fungi", y = "Belowground Biomass (g)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,10)
hetero.below.plot 

agalinis.above.plot <- ggplot(data = sd.ag.above, aes(x = treatment, y = mean, fill = type)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#71A4A0", "#D6A839")) +
  labs(x = "Fungi", y = "Aboveground Biomass (g)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,2.5)
agalinis.above.plot 

agalinis.below.plot <- ggplot(data = sd.ag.below, aes(x = treatment, y = mean, fill = type)) +
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  theme_pubr() +
  scale_fill_manual( values=c("#71A4A0", "#D6A839")) +
  labs(x = "Fungi", y = "Belowground Biomass (g)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0,2.5)
agalinis.below.plot 



ag.biomass.plots <- ggarrange(agalinis.above.plot, agalinis.below.plot,
                             labels = c("A","B"), 
                             nrow = 1, ncol = 2, common.legend = TRUE)
ag.biomass.plots

he.biomass.plots <- ggarrange(hetero.above.plot, hetero.below.plot,
                             labels = c("C","D"), 
                             nrow = 1, ncol = 2, common.legend = TRUE)
he.biomass.plots

biomass.plots <- ggarrange(ag.biomass.plots, he.biomass.plots,
                          nrow = 2, ncol = 1)

biomass.plots
ggsave(plot = biomass.plots, filename = 'biomass.png',
       width = 10 ,height = 8, units = "in", dpi = 600, 
       bg = "transparent")

#-------------------------PHYSIOLOGY----------------------------------#
library(gasanalyzer)
physiology <- read.csv("HPM Physiology - Phys.csv")
str(physiology)
physiology$treatment[physiology$treatment == "innoculated"] <- "AMF"
physiology$treatment[physiology$treatment == "sterilized"] <- "Control"
agalinis.physiology <- filter(physiology, species == "AGPU")
agalinis.physiology$type[agalinis.physiology$type == "host-parasite"] <- "With Host"
agalinis.physiology$type[agalinis.physiology$type == "parasite"] <- "Alone"
hetero.physiology <- filter(physiology, species == "HESU")
hetero.physiology$type[hetero.physiology$type == "host-parasite"] <- "With Parasite"
hetero.physiology$type[hetero.physiology$type == "host"] <- "Alone"
physiology <- as.data.frame(unclass(physiology),stringsAsFactors=TRUE)
hetero.physiology <- as.data.frame(unclass(hetero.physiology),stringsAsFactors=TRUE)
agalinis.physiology <- as.data.frame(unclass(agalinis.physiology),stringsAsFactors=TRUE)

#------Phys Data Analysis------#
#Agalinis
#Carbon assimilation (photosynthesis)
ag.carb.lm <- lmer(A ~ treatment*type + (1|replicate_id), data = agalinis.physiology)
summary(ag.carb.lm)
Anova(ag.carb.lm) 
emmeans(ag.carb.lm, pairwise ~ type|treatment)
emmip(ag.carb.lm, ~ type ~ treatment)

#Transpiration
ag.trans.lm <- lmer(E ~ treatment*type + (1|replicate_id), data = agalinis.physiology)
summary(ag.trans.lm)
Anova(ag.trans.lm) 
emmeans(ag.trans.lm, pairwise ~ type|treatment)
emmip(ag.trans.lm, ~ type ~ treatment)

#Stomotal conductance
ag.sto.lm <- lmer(gsw ~ treatment*type + (1|replicate_id), data = agalinis.physiology)
summary(ag.sto.lm)
Anova(ag.sto.lm) 
emmeans(ag.sto.lm, pairwise ~ type|treatment)
emmip(ag.sto.lm, ~ type ~ treatment)

#Heterotheca
#Carbon assimilation (photosynthesis)
he.carb.lm <- lmer(A ~ treatment*type + (1|replicate_id), data = hetero.physiology)
summary(he.carb.lm)
Anova(he.carb.lm) # type: chisq = 1.978, p = 0.0022**
emmeans(he.carb.lm, pairwise ~ type|treatment)
emmip(he.carb.lm, ~ type ~ treatment)

#Transpiration
he.trans.lm <- lmer(E ~ treatment*type + (1|replicate_id), data = hetero.physiology)
summary(he.trans.lm)
Anova(he.trans.lm) #type: chisq = 3.578, p = 0.059. treatment:type: chisq = 2.829, p = 0.0926.
emmeans(he.trans.lm, pairwise ~ type|treatment)
emmip(he.trans.lm, ~ type ~ treatment)

#Stomotal conductance
he.sto.lm <- lmer(gsw ~ treatment*type + (1|replicate_id), data = hetero.physiology)
summary(he.sto.lm)
Anova(he.sto.lm) #no significance
emmeans(he.sto.lm, pairwise ~ type|treatment)
emmip(he.sto.lm, ~ type ~ treatment)

#Comparative Work
#Carbon assimilation (photosynthesis)
comp.carb.lm <- lmer(A ~ species*type + species*treatment + (1|replicate_id), data = physiology)
summary(comp.carb.lm)
Anova(comp.carb.lm) # species: chisq: 5.982, p = 0.0146
emmeans(comp.carb.lm, pairwise ~ treatment|species)
emmip(comp.carb.lm, ~ species ~ treatment)

#Transpiration
comp.trans.lm <- lmer(E ~ species*type + species*treatment + (1|replicate_id), data = physiology)
summary(comp.trans.lm)
Anova(comp.trans.lm) 
emmeans(comp.trans.lm, pairwise ~ type|treatment)
emmip(comp.trans.lm, ~ type ~ treatment)

#Stomotal conductance
comp.sto.lm <- lmer(gsw ~ species*type + species*treatment + (1|replicate_id), data = physiology)
summary(comp.sto.lm)
Anova(comp.sto.lm) 
emmeans(comp.sto.lm, pairwise ~ type|treatment)
emmip(comp.sto.lm, ~ type ~ treatment)

#Standard error calculations
phys <- physiology %>% drop_na()
he.phys <- hetero.physiology %>% drop_na()
ag.phys <- agalinis.physiology %>% drop_na()

ag.carb.phys <- ag.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(A),
                   se = sd(A)/sqrt(n()))

he.carb.phys <- he.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(A),
                   se = sd(A)/sqrt(n()))

ag.trans.phys <- ag.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(E),
                   se = sd(E)/sqrt(n()))

he.trans.phys <- he.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(E),
                   se = sd(E)/sqrt(n()))

ag.sto.phys <- ag.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(gsw),
                   se = sd(gsw)/sqrt(n()))

he.sto.phys <- he.phys %>% 
  group_by(treatment, type) %>% 
  dplyr::summarise(mean = mean(gsw),
                   se = sd(gsw)/sqrt(n()))

#------Physiology graphs Graphs------#

ag.carb.plot <- ggplot(ag.carb.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Carbon Assimilation (µmol m-2 s-1)") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,10)

ag.carb.plot

ag.trans.plot <- ggplot(ag.trans.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Transpiration (mol m-2 s-1)") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,0.04)

ag.trans.plot

ag.sto.plot <- ggplot(ag.sto.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Stomatal Conductance (mol m-2 s-1)") +
  scale_color_manual(values=c("#71A4A0", "#D6A839")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,1)

ag.sto.plot

ag.phys.plots <- ggarrange(ag.carb.plot, ag.trans.plot, ag.sto.plot,
                           labels = c("A", "B","C"), 
                           nrow = 1,
                           common.legend = TRUE)
ag.phys.plots  

he.carb.plot <- ggplot(he.carb.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Carbon Assimilation (µmol m-2 s-1)") +
  scale_color_manual( values=c("#3d405b", "#e07a5f")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,10)

he.carb.plot

he.trans.plot <- ggplot(he.trans.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Transpiration (mol m-2 s-1)") +
  scale_color_manual( values=c("#3d405b", "#e07a5f")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,0.04)

he.trans.plot

he.sto.plot <- ggplot(he.sto.phys, aes(x = treatment, y = mean, color = type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Stomatal Conductance (mol m-2 s-1)") +
  scale_color_manual( values=c("#3d405b", "#e07a5f")) +
  scale_shape_manual(values=c(18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,1)

he.sto.plot

he.phys.plots <- ggarrange(he.carb.plot, he.trans.plot, he.sto.plot,
                        labels = c("A", "B","C"), 
                        nrow = 1,
                        common.legend = TRUE)
he.phys.plots  
#---------------------HAUSTORIA----------------------#

#Haustoria: import, clean and rename
haustoria <- read.csv("HPM Haustoria - Haustoria.csv")
str(haustoria)
haustoria$treatment[haustoria$treatment == "innoculated"] <- "AMF"
haustoria$treatment[haustoria$treatment == "sterilized"] <- "Control"
haustoria <- as.data.frame(unclass(haustoria),stringsAsFactors=TRUE)
#Pallette: Host: #E1BE6A, Parasite: #40B0A6, HP: #7D0050

#---------Haustoria data analysis--------#

#Proportion of Haustoria Attached
haustoria.lm <- lmer(prop_attached ~ root_type*treatment + (1|replicate_id), data = haustoria)
summary(haustoria.lm)
Anova(haustoria.lm) # root_type:treatment p = 0.002
emmeans(haustoria.lm, pairwise ~ treatment|root_type)
emmip(haustoria.lm, ~ root_type ~ treatment)

#Standard error calculations
haus <- haustoria %>% drop_na()

attached <- haus %>% 
  group_by(treatment, root_type) %>% 
  dplyr::summarise(mean = mean(prop_attached),
                   se = sd(prop_attached)/sqrt(n()))
#------ graphs ------#
haustoria.plot <- ggplot(attached, aes(x = treatment, y = mean, color = root_type)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), linewidth = 1, width = 0.09) +
  geom_point(aes(shape = root_type), size = 7, position =  position_dodge(width = 0.5)) +
  labs(x = "Fungi", y = "Proportion of Attached Haustoria") +
  scale_color_manual(values=c("#71A4A0", "#E1BE6A","#E26688")) +
  scale_shape_manual(values=c(18,18,18)) +
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  ylim(0,1)

haustoria.plot

