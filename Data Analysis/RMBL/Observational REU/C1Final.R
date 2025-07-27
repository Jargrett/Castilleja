setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")
#load in relevant packages
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


castilleja.cover <- read.csv("castilleja cover complete.csv")
castilleja.cover$castilleja[castilleja.cover$castilleja == "Control"] <- "Absent"
castilleja.cover$castilleja[castilleja.cover$castilleja == "Castilleja"] <- "Present"
castilleja.cover <- as.data.frame(unclass(castilleja.cover),stringsAsFactors=TRUE)
cover.overview <- read.csv("average cover.csv")
cover.overview <- as.data.frame(unclass(cover.overview),stringsAsFactors=TRUE)
castilleja.cover$year = as.factor(castilleja.cover$year)

#Diversity Analysis
div <- lmer(div ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(div)
Anova(div)
emmip(div, ~ castilleja ~ year)
emmeans(div, pairwise ~ castilleja|year)

rich <- lmer(rich ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(rich)
Anova(rich) 
emmip(rich, castilleja ~ year)
emmeans(rich, pairwise ~ species|year)

even <- lmer(even ~ castilleja*species + castilleja*year + (1|pair) + (1|site), data = castilleja.cover)
summary(even)
Anova(even) 
emmip(even, castilleja ~ year)
emmeans(even, pairwise ~ castilleja|year)

#Standard error calculations
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

castilleja.div$colcast = castilleja.div$castilleja
castilleja.rich$colcast = castilleja.rich$castilleja
castilleja.even$colcast = castilleja.even$castilleja

#Diversity Graphs
  
div.plot <- ggplot(data = castilleja.div, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.6, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.6, colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Shannon Diversity of co-occuring species") +
  theme(legend.position="none") +
  ylim(1.4,2)

div.plot

rich.plot <- ggplot(data = castilleja.rich, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.7, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.6, colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species Richness of co-occuring species") +
  theme(legend.position="none") +
  ylim(6,12)

rich.plot

even.plot <- ggplot(data = castilleja.even, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                   position =  position_dodge(width = 0.5), size = 0.7, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.2,colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species evenness of co-occuring species") +
  theme(legend.position="none") +
  ylim(0.6,1)

even.plot

diversity.plots <- ggarrange(div.plot, rich.plot, even.plot,
                            labels = c("a", "b","c"), 
                            nrow = 1)
diversity.plots  

#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compositional differences between our sites,castilleja,etc.
#we are working towards Matrix format so we can take our castilleja matrix as our starting point
library(ggrepel)
library(vegan)
library(ggordiplots)

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

perm <- adonis2(dist ~ castilleja*species + castilleja*year + castilleja*site, data = NMDS, permutations=9999)
perm

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape = species), size = 2.2, alpha = 0.8) +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  stat_ellipse(geom = "polygon", segments = 20, linetype = 2, alpha = 0.1, aes(group = site)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = site)) +
  geom_text(aes(-1.53,0.72), label = "EL", color = "grey22", size = 3.5) +
  geom_text(aes(-0.9,0.1), label = "AP", color = "grey22", size = 3.5) +
  geom_text(aes(-0.3,0.58), label = "CC", color = "grey22", size = 3.5) +
  geom_text(aes(0.85,1.85), label = "AL", color = "grey22", size = 3.5) +
  geom_text(aes(0.85,0.9), label = "DC1", color = "grey22", size = 3.5) +
  geom_text(aes(0.23,0), label = "DC2", color = "grey22", size = 3.5) +
  geom_text(aes(-0.14,-0.5), label = "JH", color = "grey22", size = 3.5) +
  coord_equal() +
  theme_bw()


#-------------------------Indicator Species Analysis---------------------------#
#In this analysis we will be assessing species that are found more often in one treatment group compared to another.
#This package takes into account both the relative abundance in a given plot as well as presence absence so we will use cover data
#we will assess this by species and sites but across years

#Loading in necessary packages
library(indicspecies)

#We will need a species matrix (castilleja removed) and vector that contains our presence absence plot info

#Avery
case.avery <- filter(castilleja.cover, site == "Avery")#filtering for a specific site
case.avery.matrix <- case.avery %>%#this segment selects just our species matrix and removes castilleja
  select(17:158)

case.avery.cast = case.avery$castilleja

case.avery.inv = multipatt(case.avery.matrix, case.avery.cast, func = "r.g",control = how(nperm=9999))

summary(case.avery.inv, alpha=1)#Nothing significant for Avery

#Emerald Lake
case.emerald <- filter(castilleja.cover, site == "Emerald Lake")
case.emerald.matrix <- case.emerald %>%
  select(17:158)

case.emerald.cast = case.emerald$castilleja

case.emerald.inv = multipatt(case.emerald.matrix, case.emerald.cast, func = "r.g", control = how(nperm=9999))

summary(case.emerald.inv, alpha=1)#Castilleja Group: Fragaria.virginiana p = 0.0422, 
#Present: Erigeron.coulteri p = 0.0553

#Copper Creek
case.copper <- filter(castilleja.cover, site == "Copper Creek")
case.copper.matrix <- case.copper %>%
  select(17:158)

case.copper.cast = case.copper$castilleja

case.copper.inv = multipatt(case.copper.matrix, case.copper.cast, func = "r.g", control = how(nperm=9999))

summary(case.copper.inv, alpha=1)#Nothing significant for Copper Creek
#Absent: Cymopterus.lemmonii p = 0.0893 Present: Thalictrum.fendleri p = 0.0682, Campanula.petiolata p = 0.0849

#Now we do linariifolia
#Deer Creek 1
cali.dc1 <- filter(castilleja.cover, site == "Deer Creek 1")
cali.dc1.matrix <- cali.dc1 %>%
  select(17:158)

cali.dc1.cast = cali.dc1$castilleja

cali.dc1.inv = multipatt(cali.dc1.matrix, cali.dc1.cast, func = "r.g", control = how(nperm=9999))

summary(cali.dc1.inv, alpha=1)#Nothing significant for Deer Creek 1
#Absent: Ipomopsis.aggregatta p = 0.0551, Elymus.elymoides p = 0.0864
#Present: Heterotheca.sp. p = 0.0836

#Deer Creek 2
cali.dc2 <- filter(castilleja.cover, site == "Deer Creek 2")
cali.dc2.matrix <- cali.dc2 %>%
  select(17:158)

cali.dc2.cast = cali.dc2$castilleja

cali.dc2.inv = multipatt(cali.dc2.matrix, cali.dc2.cast, func = "r.g", control = how(nperm=9999))

summary(cali.dc2.inv, alpha=1)#Castilleja Group: Delphinum.nuttalliianum p = 0.0050, Koeleria.macrantha = 0.0437
#Absent: Gayophytum.sp. p = 0.0512, Erigeron.speciosus p = 0.0661

#Johnson Hill
cali.johnson <- filter(castilleja.cover, site == "Johnson Hill")
cali.johnson.matrix <- cali.johnson %>%
  select(17:158)

cali.johnson.cast = cali.johnson$castilleja

cali.johnson.inv = multipatt(cali.johnson.matrix, cali.johnson.cast, func = "r.g", control = how(nperm=9999))

summary(cali.johnson.inv, alpha=1)#Nothing significant for Johnson Hill

#Almont
cacr.Almont <- filter(castilleja.cover, site == "Almont")
cacr.cover.matrix <- cacr.Almont %>%
  select(17:158)

cacr.cover.cast = cacr.Almont$castilleja

cacr.cover.inv = multipatt(cacr.cover.matrix, cacr.cover.cast, func = "r.g", control = how(nperm=9999))

summary(cacr.cover.inv, alpha=1)#Nothing significant for Almont

#--------------------------Nearest Neighbor Analysis---------------------------#
#
#
library(ggplot2)
library(ggpubr)
library(hrbrthemes)

case.nn <- read.csv("NN - Case.csv")
case.nn$rel_abund_cover <- case.nn$rel_abund_cover / 100
case.nn$NN_freq <- case.nn$NN_freq / 100
case.nn <- as.data.frame(unclass(case.nn),stringsAsFactors=TRUE)
case.nn$year <- factor(case.nn$year)

cali.nn <- read.csv("NN - Cali.csv")
cali.nn$rel_abund_cover <- cali.nn$rel_abund_cover / 100
cali.nn$NN_freq <- cali.nn$NN_freq / 100
cali.nn <- as.data.frame(unclass(cali.nn),stringsAsFactors=TRUE)
cali.nn$year <- factor(cali.nn$year)

cacr.nn <- read.csv("NN - Cacr.csv")
cacr.nn$rel_abund_cover <- cacr.nn$rel_abund_cover / 100
cacr.nn$NN_freq <- cacr.nn$NN_freq / 100
cacr.nn <- as.data.frame(unclass(cacr.nn),stringsAsFactors=TRUE)
cacr.nn$year <- factor(cacr.nn$year)

emerald.nn <- filter(case.nn, site == "Emerald Lake")
emerald.lm <- lm(NN_freq ~ rel_abund_cover, data = emerald.nn)


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
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.1071574642, 0.24050633), label = "Deschampsia cespitosa", color = "grey22", nudge_y = 0.015, size = 3.5) +
  geom_text(aes(0.1251533742, 0.18987342), label = "Fragaria virginiana", color = "grey22", nudge_y = -0.01, size = 3.5) +
  geom_text(aes(0.1480333652, 0.20430108), label = "Fragaria virginiana", color = "grey22", nudge_y = 0.01, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
emerald.nearest

avery.nearest <- ggplot(AO.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.1531226486, 0.29411765), label = "Poa pratensis", color = "grey22", nudge_y = 0.015, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
avery.nearest

copper.nearest <- ggplot(CC.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0778150748, 0.13333333), label = "Thalictrum fendleri", color = "grey22", nudge_y = 0.01, size = 3.5) +
  geom_text(aes(0.0234589564, 0.06666667), label = "Mertensia brevistyla", color = "grey22", nudge_y = 0.01, size = 3.5) +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
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
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0746228926, 0.155038760), label = "Eremogone congesta", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0143414410, 0.121212121), label = "Bromus inermis", color = "grey22",nudge_x = 0.007, nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0402707664, 0.111111111), label = "Carex sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0383203304, 0.090909091), label = "Achnatherum sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0463176575, 0.085271318), label = "Carex sp.", color = "grey22", nudge_y = -0.005, size = 3.5) +
  theme_minimal() +
  ggtitle("Castilleja linariifolia") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
deercreek.nearest

johnson.nearest <- ggplot(JH.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.0738952297, 0.14492754), label = "Lathyrus lanszwertii", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0496049166, 0.10144928), label = "Viola praemorsa", color = "grey22", nudge_y = 0.007, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja linariifolia") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
johnson.nearest

cali.cacr.nearestplots <- ggarrange(deercreek.nearest, johnson.nearest, almont.nearest,
                               labels = c("C", "D", "E"), 
                               nrow = 1, common.legend = FALSE, legend = "bottom")

cali.nearestplots 

almont.lm <- lm(NN_freq ~ rel_abund_cover, data = cacr.nn)
summary(almont.lm)

almont.predict <- as.data.frame(predict(almont.lm, newdata = cacr.nn, interval = "prediction", level = 0.95))

AL.nn <- cbind(cacr.nn, almont.predict)

almont.nearest <- ggplot(AL.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.0574450029, 0.16417910), label = "Stipeae", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0168418304, 0.08955224), label = "Crepis sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0577061166, 0.01492537), label = "Artemisia arbuscula", color = "grey22", nudge_y = -0.006, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja chromosa") +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
almont.nearest

nearestplots <- ggarrange(emerald.nearest, avery.nearest, copper.nearest, deercreek.nearest, johnson.nearest, almont.nearest,
                               labels = c("a", "b", "c", "d", "e", "f"), 
                               nrow = 2, ncol = 3, common.legend = TRUE, legend = "right")
nearestplots

