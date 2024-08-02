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

#Load in 2023 + 2024 datasets (For diversity analysis we will use species counts)
CALI.24.count <- read.csv("Cali 2024 Count.csv")
CALI.23.count<- read.csv("Cali 2023 Count.csv")
CASE.23.count <- read.csv("Case 2023 Count.csv")
CASE.24.count <- read.csv("Case 2024 Count.csv")
#Load in 2023 + 2024 datasets (For diversity analysis we will use species counts)
CALI.24.cover <- read.csv("Cali 2024 Cover.csv")
CALI.23.cover <- read.csv("Cali 2023 Cover.csv")
CASE.23.cover <- read.csv("Case 2023 Cover.csv")
CASE.24.cover <- read.csv("Case 2024 Cover.csv")

cali.count <- rbind.fill(CALI.23.count,CALI.24.count)
case.count <- rbind.fill(CASE.23.count,CASE.24.count)

cali.cover <- rbind.fill(CALI.23.cover,CALI.24.cover)
case.cover <- rbind.fill(CASE.23.cover,CASE.24.cover)

#changing NA to 0
cali.count[is.na(cali.count)] <- 0
case.count[is.na(case.count)] <- 0

cali.cover[is.na(cali.cover)] <- 0
case.cover[is.na(case.cover)] <- 0
#sperating the species matrix from the environmental data
cali.env <- subset(cali.count, select=c(1:3,5:7)) #gathering enviornmental data
case.env <- subset(case.count, select=c(1:3,5:7))

#Isolating the species matric with castilleja included 
case.count.matrix <- case.count[ -c(1:7)] 
cali.count.matrix <- cali.count[ -c(1:7)]

#for cover
case.cover.matrix <- case.cover[ -c(1:10)] 
cali.cover.matrix <- cali.cover[ -c(1:10)]
#Isolating the species matrix without castilleja included in the analysis
nocase.count.matrix <- case.count[ -c(1:7,9)]
nocali.count.matrix <- cali.count[ -c(1:7,17)]

nocase.cover.matrix <- case.cover[ -c(1:10,12)]
nocali.cover.matrix <- cali.cover[ -c(1:10,20)]
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating Shannon diversity, richness, and evenness for 2023 plots (Cover data)
#linariifolia
cali.cover.div <- diversity(nocali.cover.matrix, index = "shannon")
cali.cover.rich <- specnumber(nocali.cover.matrix)
cali.cover.even <- diversity(nocali.cover.matrix, index = "shannon") / log(specnumber(nocali.cover.matrix)) 
#septentrionalis
case.cover.div <- diversity(nocase.cover.matrix, index = "shannon")
case.cover.rich <- specnumber(nocase.cover.matrix)
case.cover.even <- diversity(nocase.cover.matrix, index = "shannon") / log(specnumber(nocase.cover.matrix)) 
# Calculating Shannon diversity, richness, and evenness for 2023 plots (Count data)
#linariifolia
cali.count.div <- diversity(nocali.count.matrix, index = "shannon")
cali.count.rich <- specnumber(nocali.count.matrix)
cali.count.even <- diversity(nocali.count.matrix, index = "shannon") / log(specnumber(nocali.count.matrix)) 
#septentrionalis
case.count.div <- diversity(nocase.count.matrix, index = "shannon")
case.count.rich <- specnumber(nocase.count.matrix)
case.count.even <- diversity(nocase.count.matrix, index = "shannon") / log(specnumber(nocase.count.matrix)) 

#combined data set with environmental and calculated values
#linariifolia
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
cali.count.diversity <- cbind(cali.env,cali.count.div,cali.count.rich,cali.count.even)
#septentrionalis
case.count.diversity <- cbind(case.env,case.count.div,case.count.rich,case.count.even)
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)

#renaming columns
case.count.diversity <- plyr::rename(case.count.diversity, c("case.count.div" = "div",
                                                 "case.count.even" = "even",
                                                 "case.count.rich" = "rich",
                                                 "Treatment" = "Castilleja"))

cali.count.diversity <- plyr::rename(cali.count.diversity, c("cali.count.div" = "div",
                                                 "cali.count.even" = "even",
                                                 "cali.count.rich" = "rich",
                                                 "Treatment" = "Castilleja"))

case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich",
                                                             "Treatment" = "Castilleja"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich",
                                                             "Treatment" = "Castilleja"))

#creating a new csv for later use and sharing
write.csv(cali.count.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Count Diversity.csv", row.names=FALSE)
cali.count.diversity <- read.csv("Combined linariifolia Count Diversity.csv")
write.csv(case.count.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Count Diversity.csv", row.names=FALSE)
case.count.diversity <- read.csv("Combined septentrionalis Count Diversity.csv")

write.csv(cali.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
cali.cover.diversity <- read.csv("Combined linariifolia Cover Diversity.csv")
write.csv(case.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
case.cover.diversity <- read.csv("Combined septentrionalis Cover Diversity.csv")

#-----------------linariifolia-----------------#
#-----------------Checking data structure-----------------#
cali.count.diversity$Site <- as.factor(cali.count.diversity$Site)
cali.count.diversity$Year <- as.factor(cali.count.diversity$Year)
cali.count.diversity$Castilleja <- as.factor(cali.count.diversity$Castilleja)
cali.cover.diversity$Site <- as.factor(cali.cover.diversity$Site)
cali.cover.diversity$Year <- as.factor(cali.cover.diversity$Year)
cali.cover.diversity$Castilleja <- as.factor(cali.cover.diversity$Castilleja)

#--------Models------#
cali.count.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.count.diversity)
summary(cali.count.div)
Anova(cali.count.div)
emmip(cali.count.div, Castilleja ~ Site)
emmeans(cali.count.div, pairwise ~ Castilleja|Site)

cali.cover.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.cover.diversity)
summary(cali.cover.div)
Anova(cali.cover.div)
emmip(cali.cover.div, Castilleja ~ Site)
emmeans(cali.cover.div, pairwise ~ Castilleja|Site)

cali.count.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.count.diversity)
summary(cali.count.rich)
Anova(cali.count.rich)
emmip(cali.count.rich, Castilleja ~ Site)
emmeans(cali.count.rich, pairwise ~ Castilleja|Site)

cali.cover.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.cover.diversity)
summary(cali.cover.rich)
Anova(cali.cover.rich)
emmip(cali.cover.rich, Castilleja ~ Site)
emmeans(cali.cover.rich, pairwise ~ Castilleja|Site)

cali.count.even <- lmer(even ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.count.diversity)
summary(cali.count.even)
Anova(cali.count.even)
emmip(cali.count.even, Castilleja ~ Site)
emmeans(cali.count.even, pairwise ~ Castilleja|Site)

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

cali.richness.plot <- ggplot(cali.comb.div, aes(x = Castilleja, y = rich)) +
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

cali.evenness.plot <- ggplot(cali.comb.div, aes(x = Castilleja, y = even)) +
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
str(case.comb.div)
case.comb.div$Site <- as.factor(case.comb.div$Site)
case.comb.div$Year <- as.factor(case.comb.div$Year)
case.comb.div$Castilleja <- as.factor(case.comb.div$Castilleja)

#--------Models------#
case.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.comb.div)
summary(case.div)
Anova(case.div)
emmip(case.div, Castilleja ~ Site)
emmeans(case.div, pairwise ~ Castilleja|Site)

case.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.comb.div)
summary(case.rich)
Anova(case.rich)
emmip(case.rich, Castilleja ~ Site)
emmeans(case.rich, pairwise ~ Castilleja|Site)

case.even <- lmer(even ~ Castilleja*Site + (1|Year) + (1|Pair), data = case.comb.div)
summary(case.even)
Anova(case.even)
emmip(case.even, Castilleja ~ Site)
emmeans(case.even, pairwise ~ Castilleja|Site)

#-------plotting------#
case.diversity.plot <- ggplot(case.comb.div, aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Shannon Diversity") +
  ylim(1,2.7)

case.diversity.plot

case.richness.plot <- ggplot(case.comb.div, aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Richness") +
  ylim(3,16)

case.richness.plot

case.evenness.plot <- ggplot(case.comb.div, aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Evenness") +
  ylim(.5,1)

case.evenness.plot

septentrionalis.plot <- ggarrange(case.diversity.plot, case.richness.plot, case.evenness.plot,
                               labels = c("A", "B","C"), 
                               nrow = 1, common.legend = TRUE, legend = "bottom")
septentrionalis.plot

