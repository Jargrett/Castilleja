#Combined diversity analysis
#Initiated: 2/7/24


#setwd
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Makena Data")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # linear regression
library(lme4) # for linear mixed effect model
library(ggpubr)
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)

#import files
case <- read.csv("CASE Diversity.csv")
cali <- read.csv("CALI Diversity.csv")
case <- case[ -c(1)]
cali<- cali[ -c(1)]
cali <- cali %>% 
  rename("div" = "div.cov",
         "even" = "even.cov")

#merge dataframes
cd <- rbind(case,cali)
cd <- cd %>% 
  rename(Castilleja = Treatment)
cd$Castilleja[cd$Castilleja == 'Castilleja'] <- 'Present'
cd$Castilleja[cd$Castilleja == 'Control'] <- 'Absent'

cd <- cd %>%
  mutate(species = case_when(
    (Site == "Deer Creek 1") ~ "Castilleja linarifolia",
    (Site == "Deer Creek 2") ~ "Castilleja linarifolia",
    (Site == "Avery") ~ "Castilleja septentrionalis",
    (Site == "Emerald Lake") ~ "Castilleja septentrionalis",
  ))

#Export cd 
write.csv(cd, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Makena Data/Combined Diversity.csv", row.names=FALSE)

#this is the final master file
cd.pair <- read.csv("Combined Diversity Pair.csv")

#--------------------Models--------------------#
#Analysis for Cali
cali.div <- lmer(div ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja linarifolia"))
summary(cali.div)
Anova(cali.div)

cali.rich <- lmer(rich ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja linarifolia"))
summary(cali.rich)
Anova(cali.rich)

cali.even <- lmer(even ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja linarifolia"))
summary(cali.even)
Anova(cali.even)


#Analysis for Case
case.div <- lmer(div ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja septentrionalis"))
summary(case.div)
Anova(case.div)

case.rich <- lmer(rich ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja septentrionalis"))
summary(case.rich)
Anova(case.rich)

case.even <- lmer(even ~ Castilleja*Site + (1|Pair), data=subset(cd.pair, species == "Castilleja septentrionalis"))
summary(case.even)
Anova(case.even)

#--------------------graphs--------------------#
#Diversity
cali.div.plot <- ggplot(subset(cd.pair, species %in% "Castilleja linarifolia"), aes(x = Site, y = div, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Shannon Diversity") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,2.5)
cali.div.plot

case.div.plot <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Site, y = div, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Shannon Diversity") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,2.5)
case.div.plot

div.plot <- ggarrange(cali.div.plot, case.div.plot, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
div.plot

#richness
cali.rich.plot <- ggplot(subset(cd.pair, species %in% "Castilleja linarifolia"), aes(x = Site, y = rich, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species Richness") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,15)
cali.rich.plot
case.rich.plot <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Site, y = rich, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species Richness") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,15)
case.rich.plot

rich.plot <- ggarrange(cali.rich.plot, case.rich.plot, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
rich.plot

#evenness
cali.even.plot <- ggplot(subset(cd.pair, species %in% "Castilleja linarifolia"), aes(x = Site, y = even, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species Evenness") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,1)
cali.even.plot

case.even.plot <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Site, y = even, fill = Castilleja)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species Evenness") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,1)
case.even.plot

even.plot <- ggarrange(cali.even.plot, case.even.plot, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
even.plot

