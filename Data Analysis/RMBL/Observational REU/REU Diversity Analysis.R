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
CALI.24 <- read.csv("Cali 2024 Count.csv")
CALI.23 <- read.csv("Cali 2023 Count.csv")
CASE.23 <- read.csv("Case 2023 Count.csv")
CASE.24 <- read.csv("Case 2024 Count.csv")

#merging dataframes by pivot longer then rbind
cali <- rbind.fill(CALI.23,CALI.24)
case <- rbind.fill(CASE.23,CASE.24)

#changing NA to 0
cali[is.na(cali)] <- 0
case[is.na(case)] <- 0

#sperating the species matrix from the environmental data
cali.env <- subset(cali, select=c(1:3,5:7)) #gathering enviornmental data
case.env <- subset(case, select=c(1:3,5:7))

#Isolating the species matric with castilleja included 
case.matrix <- case[ -c(1:7)] 
cali.matrix <- cali[ -c(1:7)]

#Isolating the species matrix without castilleja included in the analysis
nocase.matrix <- case[ -c(1:7,9)]
nocali.matrix <- cali[ -c(1:7,17)]
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating Shannon diversity, richness, and evenness for 2023 plots
#linariifolia
cali.div <- diversity(cali.matrix, index = "shannon")
cali.rich <- specnumber(cali.matrix)
cali.even <- diversity(cali.matrix, index = "shannon") / log(specnumber(cali.matrix)) 
#septentrionalis
case.div <- diversity(case.matrix, index = "shannon")
case.rich <- specnumber(case.matrix)
case.even <- diversity(case.matrix, index = "shannon") / log(specnumber(case.matrix)) 

#combined data set with environmental and calculated values
#linariifolia
cali.diversity <- cbind(cali.env,cali.div,cali.rich,cali.even)

#septentrionalis
case.diversity <- cbind(case.env,case.div,case.rich,case.even)


#renaming columns
case.diversity <- plyr::rename(case.diversity, c("case.div" = "div",
                                                 "case.even" = "even",
                                                 "case.rich" = "rich",
                                                 "Treatment" = "Castilleja"))

cali.diversity <- plyr::rename(cali.diversity, c("cali.div" = "div",
                                                 "cali.even" = "even",
                                                 "cali.rich" = "rich",
                                                 "Treatment" = "Castilleja"))

#creating a new csv for later use and sharing
write.csv(cali.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Diversity.csv", row.names=FALSE)
cali.comb.div <- read.csv("Combined linariifolia Diversity.csv")
write.csv(case.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Diversity.csv", row.names=FALSE)
case.comb.div <- read.csv("Combined septentrionalis Diversity.csv")

#-----------------linariifolia-----------------#
#-----------------Checking data structure-----------------#
str(cali.comb.div)
cali.comb.div$Site <- as.factor(cali.comb.div$Site)
cali.comb.div$Year <- as.factor(cali.comb.div$Year)
cali.comb.div$Castilleja <- as.factor(cali.comb.div$Castilleja)

#--------Models------#
cali.div <- lmer(div ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.comb.div)
summary(cali.div)
Anova(cali.div)
emmip(cali.div, Castilleja ~ Site)
emmeans(cali.div, pairwise ~ Castilleja|Site)

cali.rich <- lmer(rich ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.comb.div)
summary(cali.rich)
Anova(cali.rich)
emmip(cali.rich, Castilleja ~ Site)
emmeans(cali.rich, pairwise ~ Castilleja|Site)

cali.even <- lmer(even ~ Castilleja*Site + (1|Year) + (1|Pair), data = cali.comb.div)
summary(cali.even)
Anova(cali.even)
emmip(cali.even, Castilleja ~ Site)
emmeans(cali.even, pairwise ~ Castilleja|Site)

#-------plotting------#
cali.diversity.plot <- ggplot(cali.comb.div, aes(x = Castilleja, y = div)) +
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
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
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
