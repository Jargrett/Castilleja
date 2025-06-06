#Castilleja Diversity Analysis

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(lme4)#for modeling linear mixed effect models
library(nlme)#alternative for modeling linear mixed effect models
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(vegan)#for diversity analysis
library(emmeans)#post-hoc analysis

#Load in 2023 + 2024 datasets (For diversity analysis we will use species cover)
CALI.24.cover <- read.csv("Cali 2024 Cover.csv")
CALI.23.cover <- read.csv("Cali 2023 Cover.csv")
CASE.23.cover <- read.csv("Case 2023 Cover.csv")
CASE.24.cover <- read.csv("Case 2024 Cover.csv")

#combining datesets by castilleja species
cali.cover <- rbind.fill(CALI.23.cover,CALI.24.cover)
case.cover <- rbind.fill(CASE.23.cover,CASE.24.cover)

#changing NA to 0
cali.cover[is.na(cali.cover)] <- 0
case.cover[is.na(case.cover)] <- 0

#sperating the species matrix from the environmental data
cali.env <- subset(cali.cover, select=c(1:3,5:7)) #gathering enviornmental data
case.env <- subset(case.cover, select=c(1:3,5:7))

#Isolating the species matric with castilleja included 
case.cover.matrix <- case.cover[ -c(1:10)] 
cali.cover.matrix <- cali.cover[ -c(1:10)]

#Isolating the species matrix without castilleja included in the analysis
nocase.cover.matrix <- case.cover[ -c(1:10,12)]
nocali.cover.matrix <- cali.cover[ -c(1:10,20)]
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating Shannon diversity, richness, and evenness for 2023 plots (Cover data)
#linariifolia
cali.cover.div <- diversity(nocali.cover.matrix, index = "shannon")
cali.cover.rich <- specnumber(nocali.cover.matrix)
cali.cover.even <- diversity(nocali.cover.matrix, index = "shannon") / log(specnumber(cali.cover.matrix)) 
#septentrionalis
case.cover.div <- diversity(nocase.cover.matrix, index = "shannon")
case.cover.rich <- specnumber(nocase.cover.matrix)
case.cover.even <- diversity(nocase.cover.matrix, index = "shannon") / log(specnumber(case.cover.matrix)) 

#combined data set with environmental and calculated values
#linariifolia
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
#septentrionalis
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)

#renaming columns
case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich",
                                                             "Treatment" = "Castilleja"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich",
                                                             "Treatment" = "Castilleja"))

#creating a new csv for later use and sharing
write.csv(cali.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
cali.cover.diversity <- read.csv("Combined linariifolia Cover Diversity.csv")
write.csv(case.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
case.cover.diversity <- read.csv("Combined septentrionalis Cover Diversity.csv")

#-----------------linariifolia-----------------#
#-----------------Checking data structure-----------------#
cali.cover.diversity$Site <- as.factor(cali.cover.diversity$Site)
cali.cover.diversity$Year <- as.factor(cali.cover.diversity$Year)
cali.cover.diversity$Castilleja <- as.factor(cali.cover.diversity$Castilleja)

#--------Models------#

cali.cover.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.cover.diversity)
summary(cali.cover.div)
Anova(cali.cover.div)
emmip(cali.cover.div, Castilleja ~ Site)
emmeans(cali.cover.div, pairwise ~ Castilleja|Site)

cali.cover.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.cover.diversity)
summary(cali.cover.rich)
Anova(cali.cover.rich)
emmip(cali.cover.rich, Castilleja ~ Site)
emmeans(cali.cover.rich, pairwise ~ Castilleja|Site)

cali.cover.even <- lmer(even ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.cover.diversity)
summary(cali.cover.even)
Anova(cali.cover.even)
emmip(cali.cover.even, Castilleja ~ Site)
emmeans(cali.cover.even, pairwise ~ Castilleja|Site)

#-------plotting------#
cali.diversity.plot <- ggplot(cali.cover.diversity, aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Shannon Diversity") +
  ylim(0,3)

cali.diversity.plot

cali.richness.plot <- ggplot(cali.cover.diversity, aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Richness") +
  ylim(3,18)

cali.richness.plot

cali.evenness.plot <- ggplot(cali.cover.diversity, aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Evenness") +
  ylim(.4,1)

cali.evenness.plot

linariifolia.plot <- ggarrange(cali.diversity.plot, cali.richness.plot, cali.evenness.plot,
                         labels = c("A", "B","C"), 
                         nrow = 1, common.legend = TRUE, legend = "bottom")
linariifolia.plot

#-----------------septentrionalis-----------------#
#-----------------Checking data structure-----------------#
case.cover.diversity$Site <- as.factor(case.cover.diversity$Site)
case.cover.diversity$Year <- as.factor(case.cover.diversity$Year)
case.cover.diversity$Castilleja <- as.factor(case.cover.diversity$Castilleja)

#--------Models------#
case.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.cover.diversity)
summary(case.div)
Anova(case.div)
emmip(case.div, Castilleja ~ Site)
emmeans(case.div, pairwise ~ Castilleja|Site)

case.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.cover.diversity)
summary(case.rich)
Anova(case.rich)
emmip(case.rich, Castilleja ~ Site)
emmeans(case.rich, pairwise ~ Castilleja|Site)

case.even <- lmer(even ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.cover.diversity)
summary(case.even)
Anova(case.even)
emmip(case.even, Castilleja ~ Site)
emmeans(case.even, pairwise ~ Castilleja|Site)

#-------plotting------#
case.diversity.plot <- ggplot(case.cover.diversity, aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Shannon Diversity") +
  ylim(1,2.7)

case.diversity.plot

case.richness.plot <- ggplot(case.cover.diversity, aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Richness") +
  ylim(0,20)

case.richness.plot

case.evenness.plot <- ggplot(case.cover.diversity, aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Evenness") +
  ylim(0.2,1)

case.evenness.plot

septentrionalis.plot <- ggarrange(case.diversity.plot, case.richness.plot, case.evenness.plot,
                               labels = c("A", "B","C"), 
                               nrow = 1, common.legend = TRUE, legend = "bottom")
septentrionalis.plot


#--------Composition Analysis (Multivariate analysis) -------#
# non-metric multidimensional scaling (NMDS)

