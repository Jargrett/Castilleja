#cleaned and simplified anlaysis for ease of use

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(ggpubr)

castilleja.cover <- read.csv("castilleja cover complete.csv")
castilleja.count <- read.csv("castilleja count complete.csv")

#Diversity Analysis
div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(div)
Anova(div)
emmip(div, castilleja ~ year)
emmeans(div, pairwise ~ castilleja|year)

rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(rich)
Anova(rich) 
emmip(rich, castilleja ~ species)
emmeans(rich, pairwise ~ castilleja|species)

even <- lmer(even ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(even)
Anova(even) 
emmip(even, castilleja ~ species)
emmeans(even, pairwise ~ castilleja|species)

diversity.plot <- ggplot(castilleja.cover, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Shannon Diversity") +
  ylim(0,3)

diversity.plot

richness.plot <- ggplot(castilleja.cover, aes(x = castilleja, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Species Richness") +
  ylim(0,20)

richness.plot

#Productivity Analysis

#taking the log (plant cover)
castilleja.cover$log_no_plant <- log(castilleja.cover$no_plant)

bare <- lmer(bare ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(bare)
Anova(bare) #p = 0.0001613
emmeans(bare, pairwise ~ castilleja|year) #higher in control by 7.7%
check_model(bare)

plant <- lmer(log_no_plant ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(plant)
Anova(plant)
emmip(plant, castilleja)
emmeans(plant, pairwise ~ castilleja)

total <- lmer(no_total ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(total)
Anova(total)
emmeans(total, pairwise ~ castilleja)


castilleja.bare <- castilleja.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
castilleja.plant <- castilleja.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(log_no_plant),
            se = sd(log_no_plant)/sqrt(n()))
castilleja.total <- castilleja.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_total),
            se = sd(no_total)/sqrt(n()))

bare.plot <- ggplot(data = castilleja.bare, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Percent Bareground") +
  geom_bracket(data = castilleja.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.65,
               label = "***") +
  ylim(0,1)
bare.plot 

plant.plot <- ggplot(data = castilleja.plant, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Percent Plant Cover") +
  geom_bracket(data = castilleja.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.85,
               label = "ns") +
  ylim(0,1)
plant.plot 


#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compositional differences between our sites,castilleja,etc.
library(ggrepel)
library(vegan)

#Import Datasets: we will continue using the castilleja complete = castilleja

#we are working towards Matrix format so we can take our castilleja matrix as our starting point
species.matrix <- castilleja.cover[ -c(1:16)]
species.env <- subset(castilleja.cover, select=c(1:8))
set.seed(20)

#First calculate distance matrix
dist <-vegdist(species.matrix, method="bray")


#Run NMDS on distance matrix
nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                k=2, #2 dimensions
                try=500) #for publication I recommend 500)

nmds#stress value 0.14 which is below .2 so we need to investigate

ordiplot(nmds, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds))
NMDS <- cbind(species.env,nmds.scores) #final dataset

adonis2(dist~species*castilleja + castilleja*year + castilleja*site, data = NMDS, permutations=9999)

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape=species)) +
  coord_equal() +
  theme_bw()

