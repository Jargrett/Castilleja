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

cali.count <- read.csv("combined linariifolia counts.csv")
cali.count <- cali.count %>% 
  mutate(species = "C. linariifolia")
case.count <- read.csv("combined septentrionalis counts.csv")
case.count <- case.count %>% 
  mutate(species = "C. septentrionalis")

cali.count[is.na(cali.count)] <- 0
case.count[is.na(case.count)] <- 0

comba.count <- rbind.fill(case.count, cali.count)
comba.count[is.na(comba.count)] <- 0

#castilleja.count <- read.csv("castilleja count complete.csv")

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

evenness.plot <- ggplot(castilleja.cover, aes(x = castilleja, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Species Evenness") +
  ylim(0,1)

evenness.plot

castilleja.div <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(div),
                   se = sd(div)/sqrt(n()))
castilleja.rich <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(rich),
                   se = sd(rich)/sqrt(n()))
castilleja.even <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(even),
                   se = sd(even)/sqrt(n()))

div.plot <- ggplot(data = castilleja.div, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  facet_wrap(~year) +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  labs(x = "Hemiparasite presence", y = "Species Diversity") +
  ylim(1,2)
div.plot

#raincloudplot attempts
install.packages("ggthemes")
library(tidyquant)
library(ggdist)
library(ggthemes)
install.packages("ggrain")
library(ggrain)

ggplot(castilleja.cover, aes(1, div, fill = castilleja, color = castilleja)) +
  geom_rain(alpha = .5, rain.side = 'l',
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .1, width = 0.12), width = 0.09
            )) +
  theme_classic() +
  facet_wrap(~year) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2') +
  guides(fill = 'none', color = 'none')



#Productivity Analysis
install.packages("fitdistrplus")
install.packages("glmmTMB")
library(glmmTMB)
library(fitdistrplus)
#taking the log (plant cover)
castilleja.cover$log_no_plant <- log(castilleja.cover$no_plant)

bare <- lmer(bare ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(bare)
Anova(bare) #p = 0.0001613
emmeans(bare, pairwise ~ castilleja|year) #higher in control by 7.7%
check_model(bare)

gl.bare <- glmer(no_plant ~ castilleja*species + castilleja*year + (1|pair) + (1|site), family = "beta", data = castilleja.cover)
summary(gl.bare)
Anova(gl.bare)

plant <- lmer(no_plant ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
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
  dplyr::summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
castilleja.plant <- castilleja.cover %>% 
  group_by(castilleja) %>% 
  dplyr::summarise(mean= mean(no_plant),
            se = sd(no_plant)/sqrt(n()))
castilleja.total <- castilleja.cover %>% 
  group_by(castilleja) %>% 
  dplyr::summarise(mean= mean(no_total),
            se = sd(no_total)/sqrt(n()))

bare.plot <- ggplot(data = castilleja.bare, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  labs(x = "Hemiparasite presence", y = "Percent Bareground") +
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
  labs(x = "Hemiparasite presence", y = "Percent Plant Cover") +
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


case.cover <- castilleja.cover %>%filter((species == "C. septentrionalis"))
case.matrix <- case.cover[ -c(1:16)]
case.env <- subset(case.cover, select=c(1:8))

cali.cover <- castilleja.cover %>%filter((species == "C. linariifolia"))
cali.matrix <- cali.cover[ -c(1:16)]
cali.env <- subset(cali.cover, select=c(1:8))

cacr.cover <- castilleja.cover %>%filter((species == "C. chromosa"))
cacr.matrix <- cacr.cover[ -c(1:16)]
cacr.env <- subset(cacr.cover, select=c(1:8))

#First calculate distance matrix
dist <-vegdist(species.matrix, method="bray")
case.dist <-vegdist(case.matrix, method="bray")
cali.dist <-vegdist(cali.matrix, method="bray")
cacr.dist <-vegdist(cacr.matrix, method="bray")

set.seed(20)
#Run NMDS on distance matrix
nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                k=2, #2 dimensions
                try=500) #for publication I recommend 500)
case.nmds <- metaMDS(case.dist, distance="bray", #use bray-curtis distance
                k=2, #2 dimensions
                try=500) #for publication I recommend 500)
cali.nmds <- metaMDS(cali.dist, distance="bray", #use bray-curtis distance
                     k=2, #2 dimensions
                     try=500) #for publication I recommend 500)
cacr.nmds <- metaMDS(cacr.dist, distance="bray", #use bray-curtis distance
                     k=2, #2 dimensions
                     try=500) #for publication I recommend 500)

cacr.nmds#stress value 0.14 which is below .2 so we need to investigate

ordiplot(nmds, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds))
NMDS <- cbind(species.env,nmds.scores) #final dataset

adonis2(dist~species*castilleja + castilleja*year + castilleja*site, data = NMDS, permutations=9999)

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape=species)) +
  coord_equal() +
  theme_bw()


# calculate the difference for each pair
castilleja_pres <- castilleja.cover %>% dplyr::select(pair,site,div,rich,castilleja) %>% 
  filter(castilleja == "Castilleja") 

castilleja_abs <- castilleja.cover %>% dplyr::select(pair,site,div,rich,castilleja) %>% 
  filter(castilleja == "Control") 

castilleja_biplot <- cbind(castilleja_pres,castilleja_abs)

castilleja_biplot$diff_rich <- castilleja_biplot[4] - castilleja_biplot[9]
colnames(castilleja_biplot)[6:10] <- c("pair_2","site_2","div_2","rich_2","castilleja_2","diff_rich")

ggplot(castilleja_biplot,aes(x=div,y=div_2)) +
  geom_point() +
  labs(x="richness w/ Castilleja",y="richness w/o Castilleja") +
  theme_classic()



