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

#sperating the species matrix from the environmental data
#also will isolate bareground with the environmental data
cali.env <- subset(cali.cover, select=c(1:3,5:10))
case.env <- subset(case.cover, select=c(1:3,5:10))

#Isolating the species matrix, we also will remove Castilleja from the analysis to assess the background community
nocase.cover.matrix <- case.cover[ -c(1:10,26)]
nocali.cover.matrix <- cali.cover[ -c(1:10,25)]
#Isolating the species matric with Castilleja included we will use this to run the secondary analysis
case.cover.matrix <- case.cover[ -c(1:10)]
cali.cover.matrix <- cali.cover[ -c(1:10)]

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

#Combine calculated values with our environmental data
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)

#renaming columns
case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich"))

#creating a new csv for later use and sharing
write.csv(cali.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
cali.cover.diversity <- read.csv("Combined linariifolia Cover Diversity.csv")
write.csv(case.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
case.cover.diversity <- read.csv("Combined septentrionalis Cover Diversity.csv")

#----------------------------Diversity Analysis--------------------------------#
#We will do analysis by species starting with septentrionalis first
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
#Lets first check the data structure
str(case.cover.diversity)
#We will want to change some of our integers and characters to factors
#this will allow us to group and compare by these factor groupings

case.cover.diversity$site <- as.factor(case.cover.diversity$site)
case.cover.diversity$year <- as.factor(case.cover.diversity$year)
case.cover.diversity$castilleja <- as.factor(case.cover.diversity$castilleja)
case.cover.diversity$pair <- as.factor(case.cover.diversity$pair)
case.cover.diversity$plot <- as.factor(case.cover.diversity$plot)

#Now we can begin by running a linear mixed effect model
#Fixed effects = castilleja and site
#Random effects = year and pair

#Shannon Diverisity
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

#Now we will do the same for linariifolia, code should be more or less the same
cali.cover.diversity$site <- as.factor(cali.cover.diversity$site)
cali.cover.diversity$year <- as.factor(cali.cover.diversity$year)
cali.cover.diversity$castilleja <- as.factor(cali.cover.diversity$castilleja)
cali.cover.diversity$pair <- as.factor(cali.cover.diversity$pair)
cali.cover.diversity$plot <- as.factor(cali.cover.diversity$plot)

#Now we can begin by running a linear mixed effect model
#Fixed effects = castilleja and site
#Random effects = year and pair

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

#--------------------------Species Composition Analysis-----------------------------#
#Now we want to look at differences between our paired plots in terms of cover(and sometimes count)
#We can start with our combined dataset
#We can fist look into differences in Bare ground, total cover, and plant cover following a similar model structure to before
#we will aslo first remove castilleja from both datasets

#septentrionalis
nocase.cover <- case.cover%>% select (-c(Castilleja.septentrionalis))
nocase.cover$no_case_plant <- rowSums(nocase.cover[11:81])

case.bare<- lmer(bare ~ castilleja*site + year + (1|pair), data = nocase.cover)
summary(case.bare)
Anova(case.bare) #castilleja p = 1.868e-07, site p = 1.481e-14
emmip(case.bare, castilleja ~ site)#looks like sites are different, but control plots have higher bareground consistently
emmeans(case.bare, pairwise ~ castilleja|site)#higher in control by 7-9%

case.plant<- lmer(no_case_plant ~ castilleja*site + year + (1|pair), data = nocase.cover)
summary(case.plant)
Anova(case.plant) #castilleja p = 0.002466, site p = 1.168e-14
emmip(case.plant, castilleja ~ site)#Castilleja plots have significantly higher plant cover 
emmeans(case.plant, pairwise ~ castilleja|site)#looks like between 3-6% higher

#linariifolia
nocali.cover <- cali.cover%>% select (-c(Castilleja.linariifolia))
nocali.cover$no_cali_plant <- rowSums(nocase.cover[11:75])

cali.bare<- lmer(bare ~ castilleja*site + year + (1|pair), data = nocali.cover)
summary(cali.bare)
Anova(cali.bare) #castilleja p = 0.000199, site p = 0.012024
emmip(cali.bare, castilleja ~ site)#sites are different, DC1 and DC2 seem to be driving
emmeans(cali.bare, pairwise ~ castilleja|site)# higher in control plots by 2-14%
check_model(cali.bare)

cali.plant<- lmer(no_cali_plant ~ castilleja*site + year + (1|pair), data = nocali.cover)
summary(cali.plant)
Anova(cali.plant) #castilleja p = 0.005498, site p = 2.644e-05
emmip(cali.plant, castilleja ~ site)#Castilleja plots have significantly higher plant cover, Johnson hill in driver
emmeans(cali.plant, pairwise ~ castilleja|site)#looks like between 4-10% higher

#Now we will look at compararisons between individual species
#We will focus on the three species identified by our indicator species analysis
#CALI = Delphinum.nuttalliianum, Koeleria.macrantha
#CASE = Fragaria.virginiana

#load in our pivot_longer data
case.long <- read.csv("CASE Species Cover.csv")
cali.long <- read.csv("CALI Species Cover.csv")

str(cali.long)
str(case.long)
cali.long$site <- as.factor(cali.long$site)
cali.long$year <- as.factor(cali.long$year)
cali.long$castilleja <- as.factor(cali.long$castilleja)
cali.long$species <- as.factor(cali.long$species)
cali.long$functional.group <- as.factor(cali.long$functional.group)
case.long$site <- as.factor(case.long$site)
case.long$year <- as.factor(case.long$year)
case.long$castilleja <- as.factor(case.long$castilleja)
case.long$species <- as.factor(case.long$species)
case.long$functional.group <- as.factor(case.long$functional.group)


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

#we are working towards Matrix format so we can take our castilleja matrix as our starting point
set.seed(20)

#First calculate distance matrix
case.dist <-vegdist(nocase.cover.matrix, method="bray")
cali.dist <-vegdist(nocali.cover.matrix, method="bray")

#Run NMDS on distance matrix
case.nmds <- metaMDS(case.dist, distance="bray", #use bray-curtis distance
                      k=2, #2 dimensions
                      try=500) #for publication I recommend 500)
cali.nmds <- metaMDS(cali.dist, distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=500) #for publication I recommend 500)
case.nmds#stress value 0.15 which is below .2 so we are good!
cali.nmds#stress value 0.27 which is above .2 so we need to investigate

ordiplot(case.nmds, type="text", display="sites")
ordiplot(cali.nmds, type="text", display="sites")

case.nmds.scores <- as.data.frame(vegan::scores(case.nmds))
case.NMDS <- cbind(case.env,case.nmds.scores) #final dataset

adonis2(case.dist~castilleja*site, data = case.NMDS, permutations=999)

ggplot(case.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape=site)) +
  coord_equal() +
  theme_bw()

cali.nmds.scores <- as.data.frame(vegan::scores(cali.nmds))
cali.NMDS <- cbind(cali.env,cali.nmds.scores) #final dataset

adonis2(cali.dist~castilleja*site, data = cali.NMDS, permutations=999)

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
case.nn <- read.csv("NN - Case.csv")
cali.nn <- read.csv("NN - Cali.csv")


