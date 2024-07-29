#Castilleja Diveristy Analysis

#What I need:
#list of species (5NN/Host, and 5 Dominant species)

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(glmm)
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(indicspecies)
library(statmod)

# #Load in 2023 + 2024 datasets (For cover analysis)
# CALI.24 <- read.csv("Cali 2024 Cover.csv")
# CALI.23 <- read.csv("Cali 2023 - Cover.csv")
# 
# #shifting structure to long format rather than matrix format
# cali.24.long<- pivot_longer(CALI.24, cols = Bare.ground:Wyethia.amplexicaulis,
#                          names_to = "species",
#                         values_to = "cover")
# 
# cali.23.long<- pivot_longer(CALI.23, cols = Bare.ground:Wyethia.amplexicaulis,
#                             names_to = "species",
#                             values_to = "cover")
# 
# #Removing 0s from dateset
# cali.24.cover <- filter(cali.24.long, cover > 0)
# cali.23.cover <- filter(cali.23.long, cover > 0)
# cali.johnson.cover <- filter(cali.johnson.long, cover > 0)
# 
# cali.cover.comb <- rbind(cali.23.cover,cali.24.cover,cali.johnson.cover)
# cali.cover.zero <- rbind(cali.23.l,cali.24.long,cali.johnson.long)
# 
# write.csv(cali.cover.comb, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Species Cover.csv", row.names=FALSE)

#Load in 2023 + 2024 datasets for septentrionalis(For cover analysis)
CASE.24 <- read.csv("Case 2024 Cover.csv")
CASE.23 <- read.csv("Case 2023 Cover.csv")

#shifting structure to long format rather than matrix format
case.24.long<- pivot_longer(CASE.24, cols = Bare.ground:Viola.adunca,
                         names_to = "species",
                        values_to = "cover")

case.23.long<- pivot_longer(CASE.23, cols = Bare.ground:Viola.adunca,
                            names_to = "species",
                            values_to = "cover")

#Removing 0s from dateset
case.24.cover <- filter(case.24.long, cover > 0)
case.23.cover <- filter(case.23.long, cover > 0)


case.cover.comb <- rbind(case.23.cover,case.24.cover)
case.cover.zero <- rbind(case.23.long,case.24.long)

write.csv(case.cover.comb,"/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Species Cover.csv", row.names=FALSE)
install.packages("openxlsx")
library(openxlsx)
write.xlsx(case.cover.comb,"/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Species Cover.xlsx")
#-----------------------------linariifolia COVER ANALYSIS BEGINS----------------------------#

cali.cover <- read.csv("CALI Species Cover.csv")

str(cali.cover)
cali.cover$site <- as.factor(cali.cover$site)
cali.cover$year <- as.factor(cali.cover$year)
cali.cover$castilleja <- as.factor(cali.cover$castilleja)
cali.cover$species <- as.factor(cali.cover$species)

hist(cali.cover$cover)
boxplot(cali.cover$cover~cali.cover$castilleja)

#----------------Individual Species analysis--------------------#
#Does the presence of Castilleja alter the cover of species?
#Species list:
#want to understand whether the cover of a species changes, if it exists when Castilleja is present 
#should only consider a pair of plots where that species exists in both plots

#We need to get the data into that format, 
#i.e comparison of cover values between paired plots with a give species present

#First we subset the data
indv.species = subset(cali.cover, select = -c(2,5,7,8))
#The we pivot the columns wider creating two cover columns based on a given species
species.pair <- indv.species %>%
  pivot_wider(names_from = castilleja,
              values_from = c(cover)) %>% drop_na() %>%
  filter(if_any(species, ~ !(.x %in% c("Bare.ground"))))

species.pair$cover.difference <- species.pair$Castilleja-species.pair$Control

#analysis
summary(species.pair)
ARTR.pair <- species.pair%>% filter (species == "Artemisia.tridentada")

t.test(ARTR.pair$Castilleja, ARTR.pair$Control,
       paired = TRUE,   
       conf.level = 0.95)

grass.pair <- species.pair%>% filter (functional.group == "grass")

t.test(grass.pair$Castilleja, grass.pair$Control,
       paired = TRUE,   
       conf.level = 0.95)

ERCO.pair <- species.pair%>% filter (species == "Eremogone.congesta")

t.test(ERCO.pair$Castilleja, ERCO.pair$Control,
       paired = TRUE,   
       conf.level = 0.95)


#subset data for pair 1
#pair <- indv.species %>% dplyr::filter(pair %in% c("1"))
# pair3 <- pair2 %>%
#   group_by(species) %>% 
#   fill(Castilleja, Control, .direction = 'up') %>% 
#   fill(Castilleja, Control) %>% 
#   distinct()

#final <- pair3 %>% dplyr::filter(plot %in% c("1"))

func.diff <- ggplot(species.pair, aes(x = functional.group, y = cover.difference)) +
  geom_point(stat='identity',position = position_jitter(width = 0.1, height = 0.1)) +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "red") +
  facet_wrap(~site) +
  ylim(-0.5,0.5)
func.diff

presence.bar <- ggplot(species.pair, aes(x = functional.group, y = Castilleja)) +
  geom_point(stat='identity',position = position_jitter(width = 0.1, height = 0.1)) +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "red") +
  theme_pubr() +
  facet_wrap(~site) +
  ylim(0,1)
presence.bar

absence.bar <- ggplot(species.pair, aes(x = functional.group, y = Control)) +
  geom_point(stat='identity',position = position_jitter(width = 0.1, height = 0.1)) +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "red") +
  theme_pubr() +
  facet_wrap(~site) +
  ylim(0,1)
absence.bar

#Birds eye view
plant.lm <- lm(plant.cover ~ castilleja*site, data = cali.cover)
summary(plant.lm)
Anova(plant.lm)
emmeans(plant.lm, pairwise ~ castilleja|site)
emmip(plant.lm, castilleja ~ site)


grass.cover <- cali.cover %>% filter (functional.group == "grass")
grass.lm <- lm(cover ~ castilleja*site, data = grass.cover)
summary(grass.lm)
Anova(grass.lm)
emmeans(grass.lm, pairwise ~ castilleja|site)
emmip(grass.lm, castilleja ~ site)

ARTR.cover <- cali.cover%>% filter (species == "Artemisia.tridentada")
ARTR.lm <- lm(cover ~ castilleja*site, data = ARTR.cover)
summary(ARTR.lm)
Anova(ARTR.lm)
emmeans(ARTR.lm, pairwise ~ castilleja|site)
emmip(ARTR.lm, castilleja ~ site)


#-----------------------------linariifolia COVER ANALYSIS BEGINS----------------------------#
case.cover <- read.csv("Case Species Cover.csv")

str(case.cover)

case.cover$site <- as.factor(case.cover$site)
case.cover$year <- as.factor(case.cover$year)
case.cover$castilleja <- as.factor(case.cover$castilleja)
case.cover$species <- as.factor(case.cover$species)