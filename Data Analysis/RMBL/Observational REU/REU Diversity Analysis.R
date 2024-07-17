#Castilleja Diversity Project

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in packages
library(tidyverse)#for data wrangling and restructuring
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
CASE.24 <- read.csv("Partial Case 2024 Count.csv")
#Adding a Column for year (2023)
CALI.23 <- CALI.23 %>%
  mutate(Year = 2023)
CASE.23 <- CASE.23 %>%
  mutate(Year = 2023)
#Adding a Column for year (2024)
CALI.24 <- CALI.24 %>%
  mutate(Year = 2024)
CASE.24 <- CASE.24 %>%
  mutate(Year = 2024)
#sperating the species matrix from the environmental data
CALI.23.env <- subset(CALI.23, select=c(1,2,4:6,50)) #gathering enviornmental data
CALI.24.env <- subset(CALI.24, select=c(1,2,4:6,70))
CASE.23.env <- subset(CASE.23, select=c(1,2,4:6,41))
CASE.24.env <- subset(CASE.24, select=c(1,2,4:6,36))
#Isolating the species matrix without castilleja included in the analysis
#CALI.23.species <- CALI.23[ -c(1:6,16, 50)] 
#CALI.24.species <- CALI.24[ -c(1:6,20,70)]

#Isolating the species matric with castilleja included 
CALI.23.species <- CALI.23[ -c(1:6, 50)] 
CALI.24.species <- CALI.24[ -c(1:6,70)]
CASE.23.species <- CASE.23[ -c(1:6, 41)] 
CASE.24.species <- CASE.24[ -c(1:6,36)]
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating Shannon diversity,richness, and evenness for 2023 plots
#linariifolia
cali.div.23 <- diversity(CALI.23.species, index = "shannon")
cali.rich.23 <- specnumber(CALI.23.species)
cali.even.23 <- diversity(CALI.23.species, index = "shannon") / log(specnumber(CASE.23.species)) 
#septentrionalis
case.div.23 <- diversity(CASE.23.species, index = "shannon")
case.rich.23 <- specnumber(CASE.23.species)
case.even.23 <- diversity(CASE.23.species, index = "shannon") / log(specnumber(CASE.23.species)) 

# Calculating Shannon diversity,richness, and evenness for 2024 plots
#linariifolia
cali.div.24 <- diversity(CALI.24.species, index = "shannon")
cali.rich.24 <- specnumber(CALI.24.species)
cali.even.24 <- diversity(CALI.24.species, index = "shannon") / log(specnumber(CALI.24.species)) 
#septentrionalis
case.div.24 <- diversity(CASE.24.species, index = "shannon")
case.rich.24 <- specnumber(CASE.24.species)
case.even.24 <- diversity(CASE.24.species, index = "shannon") / log(specnumber(CASE.24.species)) 

#combined data set with environmental and calculated values
#linariifolia
CALI.23.div <- cbind(CALI.23.env,cali.div.23,cali.rich.23,cali.even.23)
CALI.24.div <- cbind(CALI.24.env,cali.div.24,cali.rich.24,cali.even.24)
#septentrionalis
CASE.23.div <- cbind(CASE.23.env,case.div.23,case.rich.23,case.even.23)
CASE.24.div <- cbind(CASE.24.env,case.div.24,case.rich.24,case.even.24)

CALI.23.div <- CALI.23.div %>% #renaming columns to prepare for merging
  rename("div" = "cali.div.23",
         "even" = "cali.even.23",
         "rich" = "cali.rich.23")

CASE.23.div <- CASE.23.div %>% #renaming columns to prepare for merging
  rename("div" = "case.div.23",
         "even" = "case.even.23",
         "rich" = "case.rich.23")

CALI.24.div <- CALI.24.div %>% #renaming columns to prepare for merging
  rename("div" = "cali.div.24",
         "even" = "cali.even.24",
         "rich" = "cali.rich.24")

CASE.24.div <- CASE.24.div %>% #renaming columns to prepare for merging
  rename("div" = "case.div.24",
         "even" = "case.even.24",
         "rich" = "case.rich.24")

#merging dataframes
cali <- rbind(CALI.23.div,CALI.24.div)
cali <- cali %>%
  relocate(Year)
cali <- cali %>% 
  rename("Castilleja" = "Treatment")

case <- rbind(CASE.23.div,CASE.24.div)
case <- case %>%
  relocate(Year)
case <- case %>% 
  rename("Castilleja" = "Treatment")

castilleja <- rbind(case,cali)

castilleja <- castilleja %>%
  mutate(species = case_when(
    (Site == "Deer Creek 1") ~ "Castilleja linarifolia",
    (Site == "Deer Creek 2") ~ "Castilleja linarifolia",
    (Site == "Avery") ~ "Castilleja septentrionalis",
    (Site == "Emerald Lake") ~ "Castilleja septentrionalis",
  ))
#creating a new csv for later use and sharing
write.csv(cali, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Diversity.csv", row.names=FALSE)
cali.comb.div <- read.csv("Combined linariifolia Diversity.csv")
write.csv(cali, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Diversity.csv", row.names=FALSE)
case.comb.div <- read.csv("Combined septentrionalis Diversity.csv")
write.csv(cali, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined Castilleja Diversity.csv", row.names=FALSE)
case.comb.div <- read.csv("Combined Castilleja Diversity.csv")
#-----------------Checking data structure-----------------#
str(cali.comb.div)
cali.comb.div$Site <- as.factor(cali.comb.div$Site)
cali.comb.div$Year <- as.factor(cali.comb.div$Year)
cali.comb.div$Castilleja <- as.factor(cali.comb.div$Castilleja)

#--------Models------#
cali.div <- lmer(div ~ Castilleja*Site*Year + (1|Pair) + (1|Date), data = cali.comb.div)
summary(cali.div)
Anova(cali.div)
emmip(cali.div, Castilleja|Year ~ Site)
emmeans(cali.div, pairwise ~ Castilleja|Site)

cali.rich <- lmer(rich ~ Castilleja*Site*Year + (1|Pair)+ (1|Date), data = cali.comb.div)
summary(cali.rich)
Anova(cali.rich)
emmip(cali.rich, Castilleja ~ Site)
emmeans(cali.rich, pairwise ~ Castilleja|Site)

cali.even <- lmer(even ~ Castilleja*Site*Year + (1|Pair)+ (1|Date), data = cali.comb.div)
summary(cali.even)
Anova(cali.even)
emmip(cali.even, Castilleja/Year ~ Site)
emmeans(cali.even, pairwise ~ Castilleja|Site)

#-------plotting------#
cali.diversity <- ggplot(cali.comb.div, aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Shannon Diversity") +
  ylim(0,3)

cali.diversity

cali.richness <- ggplot(cali.comb.div, aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Richness") +
  ylim(0,20)

cali.richness

cali.evenness <- ggplot(cali.comb.div, aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Evenness") +
  ylim(.4,1)

cali.evenness

linariifolia.plot <- ggarrange(cali.diversity, cali.richness, cali.evenness,
                         labels = c("A", "B","C"), 
                         nrow = 1, common.legend = TRUE, legend = "bottom")
linariifolia.plot



