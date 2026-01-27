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

n.senescence.lm <- lm(n ~ senescence, data = senescence)
summary(n.senescence.lm)
Anova(n.senescence.lm)

n.tissue.lm <- lm(n ~ organ, data = tissue)
summary(n.tissue.lm)
Anova(n.tissue.lm)

#Total Carbon (C)
c.litter.lm <- lm(c ~ litter_bag, data = litter)
summary(c.litter.lm)
Anova(c.litter.lm)

c.senescence.lm <- lm(c ~ senescence, data = senescence)
summary(c.senescence.lm)
Anova(c.senescence.lm)

c.tissue.lm <- lm(c ~ organ, data = tissue)
summary(c.tissue.lm)
Anova(c.tissue.lm)

#Carbon to Nitrogen Ratio (C/N)
cn.litter.lm <- lm(cn ~ litter_bag, data = litter)
summary(cn.litter.lm)
Anova(cn.litter.lm)

cn.senescence.lm <- lm(cn ~ senescence, data = senescence)
summary(cn.senescence.lm)
Anova(cn.senescence.lm)

cn.tissue.lm <- lm(cn ~ organ, data = tissue)
summary(cn.tissue.lm)
Anova(cn.tissue.lm)

#Phosphorus (P)
p.litter.lm <- lm(p ~ litter_bag, data = litter)
summary(p.litter.lm)
Anova(p.litter.lm)

p.senescence.lm <- lm(p ~ senescence, data = senescence)
summary(p.senescence.lm)
Anova(p.senescence.lm)

p.tissue.lm <- lm(p ~ organ, data = tissue)
summary(p.tissue.lm)
Anova(p.tissue.lm)

#Potassium (K)
k.litter.lm <- lm(k ~ litter_bag, data = litter)
summary(k.litter.lm)
Anova(k.litter.lm)

k.senescence.lm <- lm(k ~ senescence, data = senescence)
summary(k.senescence.lm)
Anova(k.senescence.lm)

k.tissue.lm <- lm(k ~ organ, data = tissue)
summary(k.tissue.lm)
Anova(k.tissue.lm)

#statplotting
library(ggstatsplot)

litter.plot <- ggbetweenstats(data = litter, x = litter_bag, y = k) 
litter.plot

tissue.plot <- ggbetweenstats(data = tissue, x = organ, y = k)
tissue.plot

senescence.plot <- ggbetweenstats(data = senescence, x = senescence, y = k)
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
  subset(select = c('n', 'c', 'p', 'k'))

#calculate distance matrix
dist.litter <-vegdist(litter.matrix, method="bray")
beta.litter <- betadisper(dist.litter, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
litter.nmds <- metaMDS(dist.litter, distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=1000) #for publication I recommend 500)
litter.nmds#stress value 0.14 which is below .2 so we need to investigate


ordiplot(litter.nmds, type="text", display="sites")

litter.scores <- as.data.frame(vegan::scores(litter.nmds))

litter.NMDS <- cbind(litter,litter.scores) #final dataset

perm <- adonis2(dist.litter ~ litter, data = NMDS, permutations=9999)
perm

ggplot(litter.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=litter_bag , shape = litter_bag), size = 2.2, alpha = 1) +
  scale_color_manual(values=c("#6b9080", "#c6ac8f","#22333b")) +
  stat_ellipse(geom = "polygon", segments = 20, linetype = 2, alpha = 0.1, aes(group = litter_bag)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = litter_bag)) +
  coord_equal() +
  theme_bw()

#-------Tissue Type------#
tissue.matrix <- tissue %>%
  mutate(lab_id = as.factor(lab_id)) %>% 
  subset(select = c('n', 'c', 'p', 'k'))

#calculate distance matrix
dist.tissue <-vegdist(tissue.matrix, method="bray")
beta.litter <- betadisper(dist.litter, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
tissue.nmds <- metaMDS(dist.tissue, distance="bray", #use bray-curtis distance
                       k=3, #2 dimensions
                       try=1000) #for publication I recommend 500)
tissue.nmds#stress value 0.14 which is below .2 so we need to investigate


ordiplot(tissue.nmds, type="text", display="sites")

tissue.scores <- as.data.frame(vegan::scores(tissue.nmds))

tissue.NMDS <- cbind(tissue,tissue.scores) #final dataset

perm <- adonis2(dist.tissue ~ litter, data = tissue.NMDS, permutations=9999)
perm

ggplot(tissue.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=organ , shape = organ), size = 2.2, alpha = 1) +
  scale_color_manual(values=c("#6b9080", "#c6ac8f","#22333b")) +
  stat_ellipse(geom = "polygon", segments = 20, linetype = 2, alpha = 0.1, aes(group = organ)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = organ)) +
  coord_equal() +
  theme_bw()
#-------Senescence------#
senesc.matrix <- senescence %>%
mutate(lab_id = as.factor(lab_id)) %>% 
subset(select = c('n', 'c', 'p', 'k'))

#calculate distance matrix
dist.senesc <-vegdist(senesc.matrix, method="bray")
beta.litter <- betadisper(dist.litter, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
senesc.nmds <- metaMDS(dist.senesc, distance="bray", #use bray-curtis distance
                       k=3, #2 dimensions
                       try=1000) #for publication I recommend 500)
senesc.nmds#stress value 0.14 which is below .2 so we need to investigate


ordiplot(senesc.nmds, type="text", display="sites")

senesc.scores <- as.data.frame(vegan::scores(senesc.nmds))

senesc.NMDS <- cbind(senescence,senesc.scores) #final dataset

perm <- adonis2(dist.tissue ~ litter, data = senesc.NMDS, permutations=9999)
perm

ggplot(senesc.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = senescence , shape = senescence), size = 2.2, alpha = 1) +
  scale_color_manual(values=c("#c6ac8f","#22333b")) +
  stat_ellipse(geom = "polygon", segments = 20, linetype = 2, alpha = 0.1, aes(group = senescence)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = senescence)) +
  coord_equal() +
  theme_bw()
