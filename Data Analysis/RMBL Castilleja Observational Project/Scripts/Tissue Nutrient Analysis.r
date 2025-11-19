setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(rstatix)
library(sjPlot)
library(ggcharts)
library(ggthemes)
library(ggplot2)
library(ggpubr)

#Diversity and Cover data file (combined in excel)
litter <- read.csv("Raw Data/Foliar Analysis - litter.csv")
senescence <- read.csv("Raw Data/Foliar Analysis - senescence.csv")
tissue <- read.csv("Raw Data/Foliar Analysis - tissue.csv")

#------------Plant Litter Bag Nutrient Analysis -----------#
#Total Nitrogen (N)
n.litter.lm <- lm(n ~ litter_bag, data = litter)
summary(n.litter.lm)
Anova(n.litter.lm)

n.senescence.lm <- lm(n ~ treatment, data = senescence)
summary(n.senescence.lm)
Anova(n.senescence.lm)

n.tissue.lm <- lm(n ~ organ, data = tissue)
summary(n.tissue.lm)
Anova(n.tissue.lm)

#Total Carbon (C)
c.litter.lm <- lm(c ~ litter_bag, data = litter)
summary(c.litter.lm)
Anova(c.litter.lm)

c.senescence.lm <- lm(c ~ treatment, data = senescence)
summary(c.senescence.lm)
Anova(c.senescence.lm)

c.tissue.lm <- lm(total_c ~ organ, data = tissue)
summary(c.tissue.lm)
Anova(c.tissue.lm)

#Carbon to Nitrogen Ratio (C/N)
cn.litter.lm <- lm(cn ~ litter_bag, data = litter)
summary(cn.litter.lm)
Anova(cn.litter.lm)

cn.senescence.lm <- lm(cn ~ treatment, data = senescence)
summary(cn.senescence.lm)
Anova(cn.senescence.lm)

cn.tissue.lm <- lm(cn ~ organ, data = tissue)
summary(cn.tissue.lm)
Anova(cn.tissue.lm)

#Phosphorus (P)
p.litter.lm <- lm(p ~ litter_bag, data = litter)
summary(p.litter.lm)
Anova(p.litter.lm)

p.senescence.lm <- lm(p ~ treatment, data = senescence)
summary(p.senescence.lm)
Anova(p.senescence.lm)

p.tissue.lm <- lm(p ~ organ, data = tissue)
summary(p.tissue.lm)
Anova(p.tissue.lm)

#Potassium (K)
k.litter.lm <- lm(k ~ litter_bag, data = litter)
summary(k.litter.lm)
Anova(k.litter.lm)

k.senescence.lm <- lm(k ~ treatment, data = senescence)
summary(k.senescence.lm)
Anova(k.senescence.lm)

k.tissue.lm <- lm(k ~ organ, data = tissue)
summary(k.tissue.lm)
Anova(k.tissue.lm)

#statplotting
install.packages('ggstatsplot')
library(ggstatsplot)

litter.plot <- ggbetweenstats(data = litter, x = litter_bag, y = p) 
litter.plot

tissue.plot <- ggbetweenstats(data = tissue, x = organ, y = p)
tissue.plot
senescence.plot <- ggbetweenstats(data = senescence, x = treatment, y = cn_ratio)
senescence.plot


#Multivariate
#-------------------------Multivariate analysis-------------------------#
library(ggrepel)
library(vegan)
library(ggordiplots)
library(labdsv)#enables restructuring for ecological analysis
conflicts_prefer(plyr::mutate)
litter.matrix <- litter %>%
  mutate(lab_id = as.factor(lab_id)) %>% 
  subset(select = c('lab_id','n', 'c', 'p', 'k'))

#calculate distance matrix
dist.litter <-vegdist(litter.matrix, method="bray")
beta.23 <- betadisper(dist.litter, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
nmds.23 <- metaMDS(dist.23, distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=500) #for publication I recommend 500)
nmds.23#stress value 0.14 which is below .2 so we need to investigate


ordiplot(nmds.23, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds.23))

NMDS <- cbind(plot,nmds.scores) #final dataset

perm <- adonis2(dist ~ removal*litter, data = NMDS, permutations=9999)
perm

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=removal , shape = litter), size = 2.2, alpha = 0.8) +
  scale_color_manual(values=c("#dda15e", "#606c38")) +
  coord_equal() +
  theme_bw()

