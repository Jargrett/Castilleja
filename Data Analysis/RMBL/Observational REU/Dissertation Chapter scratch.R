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

castilleja <- read.csv("Castilleja.csv")


#Diversity Analysis
div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(div)
Anova(div)
emmip(div, castilleja ~ year)
emmeans(div, pairwise ~ castilleja|year)

rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(rich)
Anova(rich) 
emmip(rich, castilleja ~ species)
emmeans(rich, pairwise ~ castilleja|species)

even <- lmer(even ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(even)
Anova(even) 
emmip(even, castilleja ~ species)
emmeans(even, pairwise ~ castilleja|species)

diversity.plot <- ggplot(castilleja, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Shannon Diversity") +
  ylim(0,3)

diversity.plot

richness.plot <- ggplot(castilleja, aes(x = castilleja, y = rich)) +
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
bare <- lmer(bare ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(bare)
Anova(bare) #p = 0.0001613
emmeans(bare, pairwise ~ castilleja|year) #higher in control by 7.7%
check_model(bare)

plant <- lmer(no_plant ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(plant)
Anova(plant)
emmip(plant, castilleja)
emmeans(plant, pairwise ~ castilleja)

total <- lmer(no_total ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja)
summary(total)
Anova(total)
emmeans(total, pairwise ~ castilleja)


castilleja.bare <- castilleja %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
castilleja.plant <- castilleja %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_plant),
            se = sd(no_plant)/sqrt(n()))
castilleja.total <- castilleja %>% 
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

  
  #--------------------------Multivariate analysis-------------------------------#
  #we will now run a (Multivariate analysis)
  #This allows us to look at the compositional differences between our sites,castilleja,etc.
  library(ggrepel)
  
  #Import Datasets: We will be using the 23_24 Combined dataset for our analysis
  case <- read.csv("Case Cover 23_24 Combined - Cover.csv")
  cali <- read.csv("Cali Cover 23_24 Combined - Cover.csv")
  cach <- read.csv("Cacr Cover 24 - Cover.csv")
  
  case.matrix <- case[ -c(1:10)]
  cali.matrix <- cali[ -c(1:10)]
  cach.matrix <- cach[ -c(1:11)]
  
  case.matrix <- case.matrix %>% 
    mutate(species = "C. septentrionalis")
  cali.matrix <- cali.matrix %>% 
    mutate(species = "C. linariifolia")
  cach.matrix <- cach.matrix %>% 
    mutate(species = "C. chromosa")
  
  species.matrix <- plyr::rbind.fill(case.matrix, cali.matrix, cach.matrix)
  species.matrix[is.na(species.matrix)] <- 0
  
  str(species.matrix)
  write.csv(species.matrix, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/castilleja species matrix.csv", row.names=FALSE)
  #Isolating the species matrix, we also will remove Castilleja from the analysis to assess the background community
  species <- read.csv("castilleja species matrix.csv")
  species.order <- species %>% 
    dplyr::select(order(names(species))) %>% 
    select(species, everything())
  nospecies.matrix <- species.order[ -c(36,37,38)]
  write.csv(nospecies.matrix, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/no castilleja species matrix.csv", row.names=FALSE)
  
  complete <- read.csv(".csv")
  #we are working towards Matrix format so we can take our castilleja matrix as our starting point
  set.seed(20)
  
  #First calculate distance matrix
  dist <-vegdist(cover.matrix, method="bray")

  
  #Run NMDS on distance matrix
  nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                          k=2, #2 dimensions
                          try=500) #for publication I recommend 500)

  nmds#stress value 0.27 which is above .2 so we need to investigate
  
  ordiplot(nmds, type="text", display="sites")
  
  nmds.scores <- as.data.frame(vegan::scores(nmds))
  NMDS <- cbind(emerald.env,emerald.nmds.scores) #final dataset
  
  adonis2(emerald.dist~species*castilleja*year*elevation, data = emerald.NMDS, permutations=999)
  
  ggplot(NMDS, aes(NMDS1, NMDS2)) +
    geom_point(aes(color=castilleja , shape=species)) +
    coord_equal() +
    theme_bw()

