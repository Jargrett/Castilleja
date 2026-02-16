setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")

#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compositional differences between our sites,castilleja,etc.
#we are working towards Matrix format so we can take our castilleja matrix as our starting point
library(ggrepel)
library(vegan)
library(ggordiplots)
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
library(rstatix)
library(permute)

castilleja.cover <- readRDS("Processed Data/Total Castilleja Cover.rds")
species.matrix <- castilleja.cover[ -c(1:16)]
species.env <- subset(castilleja.cover, select=c(1:8))

#First calculate distance matrix
dist <-vegdist(species.matrix, method="bray")

set.seed(20)
#Run NMDS on distance matrix
nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                k=2, #2 dimensions
                try=500) #for publication I recommend 500)
nmds#stress value 0.14 which is below .2 so we need to investigate


ordiplot(nmds, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds))

NMDS <- cbind(species.env,nmds.scores) #final dataset
saveRDS(NMDS,"NMDS.rds")

perm <- adonis2(dist ~ castilleja*species + castilleja*year + castilleja*site, 
                data = NMDS, permutations = 9999)
perm


cap.mod <- capscale(dist ~ castilleja*species + castilleja*year + castilleja*site + Condition(pair), 
                    data = NMDS)

perm <- how(blocks = NMDS$pair, nperm = 9999)

anova(cap.mod, permutations = perm, by = "terms")

