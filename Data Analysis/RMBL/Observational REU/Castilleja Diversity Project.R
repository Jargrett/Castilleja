#Castilleja Diversity Project
#Initiated: 8/27/24

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)

#Import Datasets: We will be using the 23_24 Combined dataset for our analysis
case.cover <- read.csv("Case Cover 23_24 Combined - Cover.csv")
cali.cover <- read.csv("Cali Cover 23_24 Combined - Cover.csv")
cacr.cover <- read.csv("Cacr Cover 24 - Cover.csv")

#sperating the species matrix from the environmental data
#also will isolate bareground with the environmental data
cali.env <- subset(cali.cover, select=c(1:3,5:10))
case.env <- subset(case.cover, select=c(1:3,5:10))
cacr.env <- subset(cacr.cover, select=c(1:3,6:11))

#Isolating the species matrix, we also will remove Castilleja from the analysis to assess the background community
nocase.cover.matrix <- case.cover[ -c(1:10,26)]
nocali.cover.matrix <- cali.cover[ -c(1:10,25)]
nocacr.cover.matrix <- cacr.cover[ -c(1:11,30)]
#Isolating the species matric with Castilleja included we will use this to run the secondary analysis
case.cover.matrix <- case.cover[ -c(1:10)]
cali.cover.matrix <- cali.cover[ -c(1:10)]
cacr.cover.matrix <- cacr.cover[ -c(1:11)]
#--------------------Calculating Diversity Values--------------------#
library(vegan)#for diversity analysis
#Using the vegan package we will calcuate Shannon Diversity, Species Richness, and Pielou's evenness
#We are using our cover data for this analysis, we will use the conservative no Castilleja matrix 

#cali
cali.cover.div <- diversity(nocali.cover.matrix, index = "shannon")
cali.cover.rich <- specnumber(nocali.cover.matrix)
cali.cover.even <- diversity(nocali.cover.matrix, index = "shannon") / log(specnumber(nocali.cover.matrix)) 
#case
case.cover.div <- diversity(nocase.cover.matrix, index = "shannon")
case.cover.rich <- specnumber(nocase.cover.matrix)
case.cover.even <- diversity(nocase.cover.matrix, index = "shannon") / log(specnumber(nocase.cover.matrix))
#cacr
cacr.cover.div <- diversity(nocacr.cover.matrix, index = "shannon")
cacr.cover.rich <- specnumber(nocacr.cover.matrix)
cacr.cover.even <- diversity(nocacr.cover.matrix, index = "shannon") / log(specnumber(nocacr.cover.matrix))

#Combine calculated values with our environmental data
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)
cacr.cover.diversity <- cbind(cacr.env,cacr.cover.div,cacr.cover.rich,cacr.cover.even)
#renaming columns
case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich"))

cacr.cover.diversity <- plyr::rename(cacr.cover.diversity, c("cacr.cover.div" = "div",
                                                             "cacr.cover.even" = "even",
                                                             "cacr.cover.rich" = "rich"))
#creating a new csv for later use and sharing
write.csv(cali.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
cali.cover.diversity <- read.csv("Combined linariifolia Cover Diversity.csv")
write.csv(case.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
case.cover.diversity <- read.csv("Combined septentrionalis Cover Diversity.csv")
write.csv(cacr.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/chromosa Cover Diversity.csv", row.names=FALSE)
cacr.cover.diversity <- read.csv("chromosa Cover Diversity.csv")
#----------------------------Diversity Analysis--------------------------------#
#We will do analysis by species starting with septentrionalis first
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(ggpubr)
#Lets first check the data structure
str(case.cover.diversity)
#We will want to change some of our integers and characters to factors
#this will allow us to group and compare by these factor groupings
case.cover.diversity <- as.data.frame(unclass(case.cover.diversity),stringsAsFactors=TRUE)
#Now we can begin by running a linear mixed effect model
#Fixed effects = castilleja and site
#Random effects = year and pair

#Shannon Diversity
case.div <- lmer(div ~ castilleja*site + year + (1|pair), data = case.cover.diversity)
summary(case.div)
Anova(case.div) #Castilleja p = 0.02415, Site p =  < 2e-16
emmip(case.div, castilleja ~ site)
emmeans(case.div, pairwise ~ castilleja|site) #Seems that difference is driven by EL
check_model(case.div)

#Species Richness
case.rich <- lmer(rich ~ castilleja*site + year + (1|pair), data = case.cover.diversity)
summary(case.rich)
Anova(case.rich) #Castilleja p = 0.002825, Site p = < 2.2e-16
emmip(case.rich, castilleja ~ site)
emmeans(case.rich, pairwise ~ castilleja|site) #Seems that difference is again driven by EL
check_model(case.rich)

#Pielou's evenness
case.even <- lmer(even ~ castilleja*site + year + (1|pair), data = case.cover.diversity)
summary(case.even)
Anova(case.even) #Castilleja p = 0.0001496, Site p = < 2.2e-16, insignificant interaction p = 0.0878680
emmip(case.even, castilleja ~ site)
emmeans(case.even, pairwise ~ castilleja|site) #Seems that difference is again driven by Avery control is more even
check_model(case.even)

#visuals
case.diversity.plot <- ggplot(case.cover.diversity, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Shannon Diversity") +
  ylim(1,2.7)

case.diversity.plot

case.richness.plot <- ggplot(case.cover.diversity, aes(x = castilleja, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Richness") +
  ylim(0,20)

case.richness.plot

case.evenness.plot <- ggplot(case.cover.diversity, aes(x = castilleja, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Evenness") +
  ylim(0.2,1)

case.evenness.plot

septentrionalis.plot <- ggarrange(case.diversity.plot, case.richness.plot, case.evenness.plot,
                                  labels = c("A", "B","C"), 
                                  nrow = 1, common.legend = TRUE, legend = "bottom")
septentrionalis.plot

#Now we will do the same for linariifolia, code should be more or less the same
cali.cover.diversity <- as.data.frame(unclass(cali.cover.diversity),stringsAsFactors=TRUE)
#Shannon Diverisity
cali.div <- lmer(div ~ castilleja*site + year + (1|pair), data = cali.cover.diversity)
summary(cali.div)
Anova(cali.div) #Castilleja p = 0.006585, Castilleja by site interaction p = 0.018733, insignificant site p = 0.055909
emmip(cali.div, castilleja ~ site) # looks like magnitudes are different between sites
emmeans(cali.div, pairwise ~ castilleja|site) #Seems that difference is driven by DC2

#Species Richness
cali.rich <- lmer(rich ~ castilleja*site + year + (1|pair), data = cali.cover.diversity)
summary(cali.rich)
Anova(cali.rich)#Castilleja p = 0.006421, Castilleja by site interaction p = 0.007023
emmip(cali.rich, castilleja ~ site)# magnitude differences as well as Johnson inversing the relationship
emmeans(cali.rich, pairwise ~ castilleja|site) #Seems that difference is again driven by DC2

#Pielou's evenness
cali.even <- lmer(even ~ castilleja*site + year + (1|pair), data = cali.cover.diversity)
summary(cali.even)
Anova(cali.even) #Site p = 0.003408
emmip(cali.even, castilleja ~ site)# lot going on here, looks like no consistant trend between sites
emmeans(cali.even, pairwise ~ castilleja|site)#DC1 difference is sig p = 0.0223 but the retalionship is higher in control similar to Case data

#Visuals
cali.diversity.plot <- ggplot(cali.cover.diversity, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Shannon Diversity") +
  ylim(0,3)

cali.diversity.plot

cali.richness.plot <- ggplot(cali.cover.diversity, aes(x = castilleja, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Richness") +
  ylim(3,18)

cali.richness.plot

cali.evenness.plot <- ggplot(cali.cover.diversity, aes(x = castilleja, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Evenness") +
  ylim(.4,1)

cali.evenness.plot

linariifolia.plot <- ggarrange(cali.diversity.plot, cali.richness.plot, cali.evenness.plot,
                               labels = c("A", "B","C"), 
                               nrow = 1, common.legend = TRUE, legend = "bottom")
linariifolia.plot

#lastly chromosa
cacr.cover.diversity <- as.data.frame(unclass(cacr.cover.diversity),stringsAsFactors=TRUE)

cacr.div <- lmer(div ~ castilleja + (1|pair), data = cacr.cover.diversity)
summary(cacr.div)
Anova(cacr.div) #Not significant higher in castilleja plots
emmeans(cacr.div, pairwise ~ castilleja) #

cacr.rich <- lmer(rich ~ castilleja + (1|pair), data = cacr.cover.diversity)
summary(cacr.rich)
Anova(cacr.rich) #Not significant marginally (p = 0.067) 
emmeans(cacr.rich, pairwise ~ castilleja) 

cacr.even <- lmer(even ~ castilleja + (1|pair), data = cacr.cover.diversity)
summary(cacr.even)
Anova(cacr.even) #Not significant However higher in control plots
emmeans(cacr.even, pairwise ~ castilleja) #

#Visuals
cacr.diversity.plot <- ggplot(cacr.cover.diversity, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c( "darkred", "burlywood4")) +
  labs(x = "Castilleja chromosa", y = "Shannon Diversity") +
  ylim(0,3)

cacr.diversity.plot

cacr.richness.plot <- ggplot(cacr.cover.diversity, aes(x = castilleja, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c( "darkred", "burlywood4")) +
  labs(x = "Castilleja chromosa", y = "Species Richness") +
  ylim(0,20)

cacr.richness.plot

cacr.evenness.plot <- ggplot(cacr.cover.diversity, aes(x = castilleja, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c( "darkred", "burlywood4")) +
  labs(x = "Castilleja chromosa", y = "Species Evenness") +
  ylim(.4,1)

cacr.evenness.plot

chromosa.plot <- ggarrange(cacr.diversity.plot, cacr.richness.plot, cacr.evenness.plot,
                               labels = c("A", "B","C"), 
                               nrow = 1, common.legend = TRUE, legend = "bottom")
chromosa.plot

#combined Castilleja diverisity
case.cover.diversity <- case.cover.diversity %>% 
  mutate(species = "C. septentrionalis")
cali.cover.diversity <- cali.cover.diversity %>% 
  mutate(species = "C. linariifolia")
cacr.cover.diversity <- cacr.cover.diversity %>% 
  mutate(species = "C. chromosa")
castilleja.diversity <- rbind(case.cover.diversity, cali.cover.diversity,cacr.cover.diversity)
castilleja.diversity <- as.data.frame(unclass(castilleja.diversity),stringsAsFactors=TRUE)

castilleja.div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.diversity)
summary(castilleja.div)
Anova(castilleja.div)
emmip(castilleja.div, castilleja ~ year)
emmeans(castilleja.div, pairwise ~ castilleja|year)

castilleja.rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.diversity)
summary(castilleja.rich)
Anova(castilleja.rich) 
emmip(castilleja.rich, castilleja ~ species)
emmeans(castilleja.rich, pairwise ~ castilleja|species)

castilleja.even <- lmer(even ~ castilleja*species + year + (1|pair) + (1|site), data = castilleja.diversity)
summary(castilleja.even)
Anova(castilleja.even) 
emmip(castilleja.even, castilleja ~ species)
emmeans(castilleja.even, pairwise ~ castilleja|species)

#Visuals
castilleja.diversity.plot <- ggplot(castilleja.diversity, aes(x = castilleja, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~species) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Shannon Diversity") +
  ylim(0,3)

castilleja.diversity.plot

castilleja.richness.plot <- ggplot(castilleja.diversity, aes(x = castilleja, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~species) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Species Richness") +
  ylim(0,20)

castilleja.richness.plot


castilleja.evenness.plot <- ggplot(castilleja.diversity, aes(x = castilleja, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~species) +
  scale_color_manual(values=c( "darkseagreen4", "burlywood4")) +
  labs(x = "Castilleja", y = "Pelou's Evenness") +
  ylim(0,1)

castilleja.evenness.plot

castilleja.plot <- ggarrange(castilleja.diversity.plot, castilleja.richness.plot, castilleja.evenness.plot,
                           labels = c("A", "B","C"), 
                           nrow = 1, common.legend = TRUE, legend = "bottom")
castilleja.plot
#-------------------------Indicator Species Analysis---------------------------#
#In this analysis we will be assessing species that are found more often in one treatment group compared to another.
#This package takes into account both the relative abundance in a given plot as well as presence absence so we will use cover data
#we will assess this by species and sites but across years

#Loading in necessary packages
library(indicspecies)

#We will need a species matrix (castilleja removed) and vector that contains our presence absence plot info

#Avery
case.avery <- filter(case.cover, site == "Avery")#filtering for a specific site
case.avery.matrix <- case.avery %>%#this segment selects just our species matrix and removes castilleja
  select(11:82) %>% 
  select (-c(Castilleja.septentrionalis))

case.avery.cast = case.avery$castilleja

case.avery.inv = multipatt(case.avery.matrix, case.avery.cast, func = "r.g", control = how(nperm=9999))

summary(case.avery.inv)#Nothing significant for Avery

#Emerald Lake
case.emerald <- filter(case.cover, site == "Emerald Lake")
case.emerald.matrix <- case.emerald %>%
  select(11:82) %>% 
  select (-c(Castilleja.septentrionalis))

case.emerald.cast = case.emerald$castilleja

case.emerald.inv = multipatt(case.emerald.matrix, case.emerald.cast, func = "r.g", control = how(nperm=9999))

summary(case.emerald.inv)#Castilleja Group: Fragaria.virginiana p = 0.0419

#Copper Creek
case.copper <- filter(case.cover, site == "Copper Creek")
case.copper.matrix <- case.copper %>%
  select(11:82) %>% 
  select (-c(Castilleja.septentrionalis))

case.copper.cast = case.copper$castilleja

case.copper.inv = multipatt(case.copper.matrix, case.copper.cast, func = "r.g", control = how(nperm=9999))

summary(case.copper.inv)#Nothing significant for Copper Creek

#Now we do linariifolia
#Deer Creek 1
cali.dc1 <- filter(cali.cover, site == "Deer Creek 1")
cali.dc1.matrix <- cali.dc1 %>%
  select(11:76) %>% 
  select (-c(Castilleja.linariifolia))

cali.dc1.cast = cali.dc1$castilleja

cali.dc1.inv = multipatt(cali.dc1.matrix, cali.dc1.cast, func = "r.g", control = how(nperm=9999))

summary(cali.dc1.inv)#Nothing significant for Deer Creek 1

#Deer Creek 2
cali.dc2 <- filter(cali.cover, site == "Deer Creek 2")
cali.dc2.matrix <- cali.dc2 %>%
  select(11:76) %>% 
  select (-c(Castilleja.linariifolia))

cali.dc2.cast = cali.dc2$castilleja

cali.dc2.inv = multipatt(cali.dc2.matrix, cali.dc2.cast, func = "r.g", control = how(nperm=9999))

summary(cali.dc2.inv)#Castilleja Group: Delphinum.nuttalliianum p = 0.0048, Koeleria.macrantha = 0.0388

#Johnson Hill
cali.johnson <- filter(cali.cover, site == "Johnson Hill")
cali.johnson.matrix <- cali.johnson %>%
  select(11:76) %>% 
  select (-c(Castilleja.linariifolia))

cali.johnson.cast = cali.johnson$castilleja

cali.johnson.inv = multipatt(cali.johnson.matrix, cali.johnson.cast, func = "r.g", control = how(nperm=9999))

summary(cali.johnson.inv)#Nothing significant for Johnson Hill

#Almont
cacr.cover.matrix <- cacr.cover %>%
  select(12:57) %>% 
  select (-c(Castilleja.chromosa))


cacr.cover.matrix[is.na(cacr.cover.matrix)] <- 0

cacr.cover.cast = cacr.cover$castilleja

cacr.cover.inv = multipatt(cacr.cover.matrix, cacr.cover.cast, func = "r.g", control = how(nperm=9999))

summary(cacr.cover.inv)#Nothing significant for Almont

#--------------------------Species Composition Analysis-----------------------------#
#Now we want to look at differences between our paired plots in terms of cover(and sometimes count)
#We can start with our combined dataset
#We can fist look into differences in Bare ground, total cover, and plant cover following a similar model structure to before
#we will aslo first remove castilleja from both datasets

#septentrionalis
nocase.cover <- case.cover%>% select (-c(Castilleja.septentrionalis))
nocase.cover$no_case_plant <- rowSums(nocase.cover[11:81])
nocase.cover$no_case_total <- nocase.cover$no_case_plant +nocase.cover$bare


case.bare <- lmer(bare ~ castilleja*site + year + (1|pair), data = nocase.cover)
summary(case.bare)
Anova(case.bare) #castilleja p = 1.868e-07, site p = 1.481e-14
emmip(case.bare, castilleja ~ site)#looks like sites are different, but control plots have higher bareground consistently
emmeans(case.bare, pairwise ~ castilleja|site)#higher in control by 7-9%

case.plant<- lmer(no_case_plant ~ castilleja*site + year + (1|pair), data = nocase.cover)
summary(case.plant)
Anova(case.plant) #castilleja p = 0.002466, site p = 1.168e-14
emmip(case.plant, castilleja ~ site)#Castilleja plots have significantly higher plant cover 
emmeans(case.plant, pairwise ~ castilleja|site)#looks like between 3-6% higher

case.total<- lmer(no_case_total ~ castilleja*site + year + (1|pair), data = nocase.cover)
summary(case.total)
Anova(case.total) #castilleja p = 0.04872, site p = 0.02171
emmip(case.total, castilleja ~ site)#Castilleja plots have significantly higher total cover 
emmeans(case.total, pairwise ~ castilleja|site) #looks like between 1-6% higher

#Standard error calc
mean_nocase.bare <- nocase.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
mean_nocase.plant <- nocase.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_case_plant),
            se = sd(no_case_plant)/sqrt(n()))
mean_nocase.total <- nocase.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_case_total),
            se = sd(no_case_total)/sqrt(n()))
#Visials

case.bareplot <- ggplot(data = mean_nocase.bare, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual(values=c("cornsilk2", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Percent Bareground") +
  geom_bracket(data = mean_nocali.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.50,
               label = "***") +
  ylim(0,1)
  
case.bareplot

case.plantplot <- ggplot(data = mean_nocase.plant, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja septentrionalis", y = "Percent Plant Cover") +
  scale_fill_manual(values=c("cornsilk2", "burlywood4")) +
  geom_bracket(data = mean_nocali.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.9,
               label = "**") +
  ylim(0,1)

case.plantplot

case.totalplot <- ggplot(data = mean_nocase.total, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja septentrionalis", y = "Percent Total Cover") +
  scale_fill_manual(values=c("cornsilk2", "burlywood4")) +
  ylim(0,1.25)

case.totalplot

case.plants.plot <- ggarrange(case.plantplot, case.bareplot,
                              labels = c("A", "B"), 
                              nrow = 1, common.legend = TRUE, legend = "bottom")
case.plants.plot

#linariifolia
nocali.cover <- cali.cover%>% select (-c(Castilleja.linariifolia))
nocali.cover$no_cali_plant <- rowSums(nocali.cover[11:75])
nocali.cover$no_calie_total <- nocali.cover$no_cali_plant + nocali.cover$bare

cali.bare<- lmer(bare ~ castilleja*site + year + (1|pair), data = nocali.cover)
summary(cali.bare)
Anova(cali.bare) #castilleja p = 0.000199, site p = 0.012024
emmip(cali.bare, castilleja ~ site)#sites are different, DC1 and DC2 seem to be driving
emmeans(cali.bare, pairwise ~ castilleja|site)# higher in control plots by 2-14%
check_model(cali.bare)

cali.plant<- lmer(no_cali_plant ~ castilleja*site + year + (1|pair), data = nocali.cover)
summary(cali.plant)
Anova(cali.plant) #site p = 0.003679, year p = 0.006288
emmip(cali.plant, castilleja ~ site)#no significance in plant cover
emmeans(cali.plant, pairwise ~ castilleja|site)#looks like between 4-10% higher

cali.total<- lmer(no_cali_total ~ castilleja*site + year + (1|pair), data = nocali.cover)
summary(cali.total)
Anova(cali.total) #site p = 0.004993, Castilleja p = 4.108e-06 
emmip(cali.total, castilleja ~ site)
emmeans(cali.total, pairwise ~ castilleja|site)

#Standard error calc
mean_nocali.bare <- nocali.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
mean_nocali.plant <- nocali.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_cali_plant),
            se = sd(no_cali_plant)/sqrt(n()))
mean_nocali.total <- nocali.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_cali_total),
            se = sd(no_cali_total)/sqrt(n()))

#Visials
cali.bareplot <- ggplot(data = mean_nocali.bare, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Percent Bareground") +
  geom_bracket(data = mean_nocali.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.65,
               label = "***") +
  ylim(0,1)
cali.bareplot

cali.plantplot <- ggplot(data = mean_nocali.plant, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.2) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja linariifolia", y = "Percent Plant Cover") +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  geom_bracket(data = mean_nocali.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.75,
               label = "ns") +
  ylim(0,1)

cali.plantplot

cali.totalplot <- ggplot(data = mean_nocali.total, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.2) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja linariifolia", y = "Percent Plant Cover") +
  scale_fill_manual(values=c("indianred3", "burlywood4")) +
  ylim(0,1.25)

cali.totalplot

cali.plants.plot <- ggarrange(cali.plantplot, cali.bareplot,
                              labels = c("A", "B"), 
                              nrow = 1, common.legend = TRUE, legend = "bottom")
cali.plants.plot

#chromosa
nocacr.cover <- cacr.cover%>% select (-c(Castilleja.chromosa))
nocacr.cover$no_cacr_plant <- rowSums(nocacr.cover[12:56])
nocacr.cover$no_cacr_total <- nocacr.cover$no_cacr_plant + nocacr.cover$bare

cacr.bare<- lmer(bare ~ castilleja + (1|pair), data = nocacr.cover)
summary(cacr.bare)
Anova(cacr.bare) #no significance
emmeans(cacr.bare, pairwise ~ castilleja)
check_model(cacr.bare)

cacr.plant<- lmer(no_cacr_plant ~ castilleja + (1|pair), data = nocacr.cover)
summary(cacr.plant)
Anova(cacr.plant) #no significance
emmip(cacr.plant, castilleja)
emmeans(cacr.plant, pairwise ~ castilleja)

cacr.total<- lmer(no_cacr_total ~ castilleja + (1|pair), data = nocacr.cover)
summary(cacr.total)
Anova(cacr.total) #castilleja p = 0.03736
emmeans(cacr.total, pairwise ~ castilleja) #2% higher in castilleja plots

#Standard error calc
mean_nocacr.bare <- nocacr.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(bare),
            se = sd(bare)/sqrt(n()))
mean_nocacr.plant <- nocacr.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_cacr_plant),
            se = sd(no_cacr_plant)/sqrt(n()))
mean_nocacr.total <- nocacr.cover %>% 
  group_by(castilleja) %>% 
  summarise(mean= mean(no_cacr_total),
            se = sd(no_cacr_total)/sqrt(n()))

cacr.bareplot <- ggplot(data = mean_nocacr.bare, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.2) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja chromosa", y = "Percent Bareground") +
  scale_fill_manual(values=c("darkred", "burlywood4")) +
  geom_bracket(data = mean_nocacr.bare,
               xmin = "Castilleja", xmax = "Control", y.position = 0.75,
               label = "ns") +
  ylim(0,1)

cacr.bareplot

cacr.plantplot <- ggplot(data = mean_nocacr.plant, aes(x = castilleja, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = .7, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.2) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.15) +
  theme_pubr() +
  labs(x = "Castilleja chromosa", y = "Percent Plant Cover") +
  scale_fill_manual(values=c("darkred", "burlywood4")) +
  geom_bracket(data = mean_nocacr.plant,
               xmin = "Castilleja", xmax = "Control", y.position = .75,
               label = "ns") +
  ylim(0,1)

cacr.plantplot

cacr.plants.plot <- ggarrange(cacr.plantplot, cacr.bareplot,
                              labels = c("A", "B"), 
                              nrow = 1, common.legend = TRUE, legend = "bottom")
cacr.plants.plot

#Now we will look at compararisons between individual species
#We will focus on the three species identified by our indicator species analysis
#CALI = Delphinum.nuttalliianum, Koeleria.macrantha
#CASE = Fragaria.virginiana

#load in our pivot_longer data
case.long <- read.csv("CASE Species Cover.csv")
cali.long <- read.csv("CALI Species Cover.csv")

str(cali.long)
str(case.long)
case.long <- as.data.frame(unclass(case.long),stringsAsFactors=TRUE)
cali.long <- as.data.frame(unclass(cali.long),stringsAsFactors=TRUE)

cali.species = subset(cali.long, select = -c(2,5,7,8,12,13))
case.species = subset(case.long, select = -c(2,5,7,8,12,13))
#The we pivot the columns wider creating two cover columns based on a given species
case.pair <- case.species %>%
  pivot_wider(names_from = castilleja,
              values_from = c(cover)) %>% drop_na() %>%
  filter(if_any(species, ~ !(.x %in% c("Bare.ground"))))

cali.pair <- cali.species %>%
  pivot_wider(names_from = castilleja,
              values_from = c(cover)) %>% drop_na() %>%
  filter(if_any(species, ~ !(.x %in% c("Bare.ground"))))

case.pair$cover.difference <- case.pair$Castilleja-case.pair$Control
cali.pair$cover.difference <- cali.pair$Castilleja-cali.pair$Control

FRVI.pair <- case.pair%>% filter (species == "Fragaria.virginiana")

t.test(FRVI.pair$Castilleja, FRVI.pair$Control,#higher in Castilleja plots by 2%, p = 0.006285
       paired = TRUE,   
       conf.level = 0.95)

DENU.pair <- cali.pair%>% filter (species == "Delphinum.nuttalliianum")

t.test(DENU.pair$Castilleja, DENU.pair$Control,#not signficant
       paired = TRUE,   
       conf.level = 0.95)

KOMA.pair <- cali.pair%>% filter (species == "Koeleria.macrantha")#not enough paired plots

#Now lets investigate functional group differences

#Grasses
case.grass <- case.long %>% 
  filter (functional.group == "grass") %>%
  pivot_wider(names_from = species,
              values_from = c(cover),
              values_fill = 0) %>% 
  select(,-c(total.cover, plant.cover, functional.group, life.history, family))

case.grass$grass.cover <- rowSums(case.grass[7:20])

case.grass.lm <- lmer(grass.cover ~ castilleja*site + year + (1|pair), data = case.grass)
summary(case.grass.lm)
Anova(case.grass.lm)
emmip(case.grass.lm, castilleja ~ site)
emmeans(case.grass.lm, pairwise ~ castilleja|site)

#Legumes
case.legume <- case.long %>% 
  filter (functional.group == "legume") %>%
  filter (site != "Emerald") %>%
  pivot_wider(names_from = species,
              values_from = c(cover),
              values_fill = 0) %>% 
  select(,-c(total.cover, plant.cover, functional.group, life.history, family))

case.legume$legume.cover <- rowSums(case.legume[7:10])

case.legume.lm <- lmer(legume.cover ~ castilleja*site + year + (1|pair), data = case.legume)
summary(case.legume.lm)
Anova(case.legume.lm)
emmip(case.legume.lm, castilleja ~ site)
emmeans(case.legume.lm, pairwise ~ castilleja|site)

#Forbs
case.forb <- case.long %>% 
  filter (functional.group == "forb") %>%
  pivot_wider(names_from = species,
              values_from = c(cover),
              values_fill = 0) %>% 
  select(,-c(total.cover, plant.cover, functional.group, life.history, family))

case.forb$forb.cover <- rowSums(case.forb[7:50])

case.forb.lm <- lmer(forb.cover ~ castilleja*site + year + (1|pair), data = case.forb)
summary(case.forb.lm)
Anova(case.forb.lm)
emmip(case.forb.lm, castilleja ~ site)
emmeans(case.forb.lm, pairwise ~ castilleja|site)

#Sedge
case.sedge <- case.long %>% 
  filter (functional.group == "sedge") %>%
  filter (site != "Avery") %>%
  pivot_wider(names_from = species,
              values_from = c(cover),
              values_fill = 0) %>% 
  select(,-c(total.cover, plant.cover, functional.group, life.history, family))

case.sedge$sedge.cover <- rowSums(case.sedge[7:10])

case.sedge.lm <- lmer(sedge.cover ~ castilleja*site + year + (1|pair), data = case.sedge)
summary(case.sedge.lm)
Anova(case.sedge.lm)
emmip(case.sedge.lm, castilleja ~ site)
emmeans(case.sedge.lm, pairwise ~ castilleja|site)

#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compoisitinoal differences between our sites,castilleja,etc.
library(ggrepel)

#First we will subset data by sites
emerald.cover <- case.cover %>%filter((site == "Emerald Lake"))
emerald.cover.matrix <- emerald.cover[ -c(1:10,26)]
emerald.env <- subset(emerald.cover, select=c(1:3,5:10))

avery.cover <- case.cover %>%filter((site == "Avery"))
avery.cover.matrix <- avery.cover[ -c(1:10,26)]
avery.env <- subset(avery.cover, select=c(1:3,5:10))

copper.cover <- case.cover %>%filter((site == "Copper Creek"))
copper.cover.matrix <- copper.cover[ -c(1:10,26)]
copper.env <- subset(copper.cover, select=c(1:3,5:10))
#we are working towards Matrix format so we can take our castilleja matrix as our starting point
set.seed(20)

#First calculate distance matrix
emerald.dist <-vegdist(emerald.cover.matrix, method="bray")
avery.dist <-vegdist(avery.cover.matrix, method="bray")
copper.dist <-vegdist(copper.cover.matrix, method="bray")

cali.dist <-vegdist(nocali.cover.matrix, method="bray")

#Run NMDS on distance matrix
emerald.nmds <- metaMDS(emerald.dist, distance="bray", #use bray-curtis distance
                      k=3, #2 dimensions
                      try=500) #for publication I recommend 500)

avery.nmds <- metaMDS(avery.dist, distance="bray", #use bray-curtis distance
                        k=2, #2 dimensions
                        try=500) #for publication I recommend 500)

copper.nmds <- metaMDS(copper.dist, distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=500) #for publication I recommend 500)
emerald.nmds#stress value 0.15 which is below .2 so we are good!
avery.nmds
copper.nmds
cali.nmds#stress value 0.27 which is above .2 so we need to investigate

ordiplot(emerald.nmds, type="text", display="sites")
ordiplot(avery.nmds, type="text", display="sites")

emerald.nmds.scores <- as.data.frame(vegan::scores(emerald.nmds))
emerald.NMDS <- cbind(emerald.env,emerald.nmds.scores) #final dataset

adonis2(emerald.dist~castilleja*year*pair, data = emerald.NMDS, permutations=999)
adonis2(avery.dist~castilleja*year*pair, data = emerald.NMDS, permutations=999)
adonis2(copper.dist~castilleja*year*pair, data = emerald.NMDS, permutations=999)

ggplot(emerald.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape=site)) +
  coord_equal() +
  theme_bw()

cali.nmds.scores <- as.data.frame(vegan::scores(cali.nmds))
cali.NMDS <- cbind(cali.env,cali.nmds.scores) #final dataset

adonis2(cali.dist~castilleja*site*pair, data = cali.NMDS, permutations=999)

ggplot(cali.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape=site)) +
  coord_equal() +
  theme_bw()

#---------------------------species count data analysis------------------------#
# Our key question here is looking at the differences between functional groups and key speices
case.23.count <- read.csv("Case 2023 - Individuals.csv")
case.24.count <- read.csv("Case 2024 - Individuals.csv")
cali.23.count <- read.csv("Cali 2023 - Individuals.csv")
cali.24.count <- read.csv("Cali 2024 - Individuals.csv")

#combining datesets by castilleja species
case.count <- rbind.fill(case.23.count,case.24.count)
cali.count <- rbind.fill(cali.23.count,cali.24.count)

#remove set merged NA values to 0
cali.count[is.na(cali.count)] <- 0
case.count[is.na(case.count)] <- 0

cali.count.long<- pivot_longer(cali.count, cols = Agastache.urticifolia:Viola.adunca,
                                                      names_to = "species",
                                                     values_to = "counts")

case.count.long<- pivot_longer(case.count, cols = Achnatherum.sp.:Vicia.americana,
                               names_to = "species",
                               values_to = "counts")

#Removing 0s from dateset
case.count.long <- filter(case.count.long, counts > 0)
cali.count.long <- filter(cali.count.long, counts > 0)

#write and reload data with functional groups added in manually
write.csv(case.count.long, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/combined septentrionalis counts.csv", row.names=FALSE)
write.csv(cali.count.long, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/combined linariifolia counts.csv", row.names=FALSE)
case.counts <- read.csv("combined septentrionalis counts - sheet1.csv")
cali.counts <- read.csv("combined linariifolia counts - sheet1.csv")


cali.counts$site <- as.factor(cali.counts$site)
cali.counts$year <- as.factor(cali.counts$year)
cali.counts$castilleja <- as.factor(cali.counts$castilleja)
cali.counts$species <- as.factor(cali.counts$species)
cali.counts$functional_group <- as.factor(cali.counts$functional_group)
case.counts$site <- as.factor(case.counts$site)
case.counts$year <- as.factor(case.counts$year)
case.counts$castilleja <- as.factor(case.counts$castilleja)
case.counts$species <- as.factor(case.counts$species)
case.counts$functional_group <- as.factor(case.counts$functional_group)

case.species.counts = subset(case.counts, select = -c(2,4,6,10))
cali.species.counts = subset(cali.counts, select = -c(2,4,6,10))


case.pair.counts <- case.species.counts %>%
  pivot_wider(names_from = castilleja,
              values_from = c(counts)) %>% drop_na()


cali.pair.counts <- cali.species.counts %>%
  pivot_wider(names_from = castilleja,
              values_from = c(counts)) %>% drop_na()

summary(cali.species.counts)

LALA.pair <- cali.pair.counts %>% filter (species == "")

t.test(LALA.pair$Castilleja, LALA.pair$Control,#higher in Castilleja plots by 2%, p = 0.006285
       paired = TRUE,   
       conf.level = 0.95)

#--------------------------Nearest Neighbor Analysis---------------------------#
#
#
library(ggplot2)
library(ggpubr)
library(hrbrthemes)

case.nn <- read.csv("NN - Case.csv")
case.nn$rel_abund_cover <- case.nn$rel_abund_cover / 100
case.nn$NN_freq <- case.nn$NN_freq / 100

cali.nn <- read.csv("NN - Cali.csv")
cali.nn$rel_abund_cover <- cali.nn$rel_abund_cover / 100
cali.nn$NN_freq <- cali.nn$NN_freq / 100

cacr.nn <- read.csv("NN - Cacr.csv")
cacr.nn$rel_abund_cover <- cacr.nn$rel_abund_cover / 100
cacr.nn$NN_freq <- cacr.nn$NN_freq / 100

#cacr.nn <- read.csv("NN - Cacr.csv")
#cacr.nn$rel_abund_cover <- cacr.nn$rel_abund_cover / 100
#cacr.nn$NN_freq <- cacr.nn$NN_freq / 100

emerald.nn <- filter(case.nn, site == "Emerald Lake")
emerald.lm <- lm(NN_freq ~ rel_abund_cover, data = emerald.nn)
summary(emerald.lm)

emerald.predict <- as.data.frame(predict(emerald.lm, newdata = emerald.nn, interval = "prediction", level = 0.95))

EL.nn <- cbind(emerald.nn, emerald.predict)


avery.nn <- filter(case.nn, site == "Avery")
avery.lm <- lm(NN_freq ~ rel_abund_cover, data = avery.nn)
summary(avery.lm)

avery.predict <- as.data.frame(predict(avery.lm, newdata = avery.nn, interval = "prediction", level = 0.95))

AO.nn <- cbind(avery.nn, avery.predict)

copper.nn <- filter(case.nn, site == "Copper Creek")
copper.lm <- lm(NN_freq ~ rel_abund_cover, data = copper.nn)
summary(copper.lm)

copper.predict <- as.data.frame(predict(copper.lm, newdata = copper.nn, interval = "prediction", level = 0.95))

CC.nn <- cbind(copper.nn, copper.predict)

emerald.nearest <- ggplot(EL.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  geom_text(aes(0.1071574642, 0.24050633), label = "Deschampsia cespitosa", color = "grey22", nudge_y = 0.015) +
  geom_text(aes(0.1251533742, 0.18987342), label = "Fragaria virginiana", color = "grey22", nudge_y = -0.01) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal()
emerald.nearest

avery.nearest <- ggplot(AO.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  #geom_point(aes(y = fit), col = "steelblue", size = 2.5) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  geom_text(aes(0.1531226486, 0.29411765), label = "Poa pratensis", color = "grey22", nudge_y = 0.015) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal()
avery.nearest

copper.nearest <- ggplot(CC.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  #geom_point(aes(y = fit), col = "steelblue", size = 2.5) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0778150748, 0.13333333), label = "Thalictrum fendleri", color = "grey22", nudge_y = 0.01) +
  geom_text(aes(0.0234589564, 0.06666667), label = "Mertensia brevistyla", color = "grey22", nudge_y = 0.01) +
  theme_minimal()
copper.nearest

case.nearestplots <- ggarrange(emerald.nearest, avery.nearest, copper.nearest,
                              labels = c("A", "B","C"), 
                              nrow = 1, common.legend = FALSE, legend = "bottom")

case.nearestplots 

johnson.nn <- filter(cali.nn, site == "Johnson Hill")
johnson.nn <- johnson.nn[-c(7), ]
johnson.lm <- lm(NN_freq ~ rel_abund_cover, data = johnson.nn)
summary(johnson.lm)

johnson.predict <- as.data.frame(predict(johnson.lm, newdata = johnson.nn, interval = "prediction", level = 0.95))

JH.nn <- cbind(johnson.nn, johnson.predict)

deer.nn <- filter(cali.nn, site == "Deer Creek")
deer.lm <- lm(NN_freq ~ rel_abund_cover, data = deer.nn)
summary(deer.lm)

deer.predict <- as.data.frame(predict(deer.lm, newdata = deer.nn, interval = "prediction", level = 0.95))

DC.nn <- cbind(deer.nn, deer.predict)

deercreek.nearest <- ggplot(DC.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  #geom_point(aes(y = fit), col = "steelblue", size = 2.5) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0746228926, 0.155038760), label = "Eremogone congesta", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0143414410, 0.121212121), label = "Bromus inermis", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0402707664, 0.111111111), label = "Carex sp.", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0383203304, 0.090909091), label = "Achnatherum sp.", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0463176575, 0.085271318), label = "Carex sp.", color = "grey22", nudge_y = -0.005) +
  theme_minimal()
deercreek.nearest

johnson.nearest <- ggplot(JH.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  geom_text(aes(0.0738952297, 0.14492754), label = "Lathyrus lanszwertii", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0496049166, 0.10144928), label = "Viola praemorsa", color = "grey22", nudge_y = 0.007) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal()
johnson.nearest

cali.nearestplots <- ggarrange(deercreek.nearest, johnson.nearest,
                               labels = c("A", "B"), 
                               nrow = 1, common.legend = FALSE, legend = "bottom")

cali.nearestplots 

almont.lm <- lm(NN_freq ~ rel_abund_cover, data = cacr.nn)
summary(almont.lm)

almont.predict <- as.data.frame(predict(almont.lm, newdata = cacr.nn, interval = "prediction", level = 0.95))

AL.nn <- cbind(cacr.nn, almont.predict)

almont.nearest <- ggplot(AL.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point() +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "cadetblue") +
  geom_line(aes(y = upr), linetype = "dashed", col = "cadetblue") +
  geom_text(aes(0.0574450029, 0.16417910), label = "Stipeae", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0168418304, 0.08955224), label = "Crepis sp.", color = "grey22", nudge_y = 0.007) +
  geom_text(aes(0.0577061166, 0.01492537), label = "Artemisia arbuscula", color = "grey22", nudge_y = -0.006) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal()
almont.nearest
