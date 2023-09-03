#Combined diversity analysis
#Initiated: 7828/23


#setwd
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # linear regression

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
  mutate(species = case_when(
    (Site == "Deer Creek 1") ~ "Castilleja linarifolia",
    (Site == "Deer Creek 2") ~ "Castilleja linarifolia",
    (Site == "Avery") ~ "Castilleja septentrionalis",
    (Site == "Emerald Lake") ~ "Castilleja septentrionalis",
  ))

#Analysis
#F-statistic: 3.544 on 1 and 158 DF,  p-value: 0.06161
even.lm <- lm(even ~ Treatment, data = cd)
summary(even.lm)
Anova(even.lm)

#F-statistic: 24.62 on 1 and 158 DF,  p-value: 1.792e-06
rich.lm <- lm(rich ~ Treatment*species, data = cd)
summary(rich.lm)
Anova(rich.lm)

#F-statistic: 13.84 on 1 and 158 DF,  p-value: 0.000276
div.lm <- lm(div ~ Treatment, data = cd)
summary(div.lm)
Anova(div.lm)


rich.plot <- ggplot(cd, aes(x = Treatment, y = rich, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species Richness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2() +
  ylim(2.5,15) +
  facet_wrap(~species)


rich.plot


div.plot <- ggplot(cd, aes(x = Treatment, y = div, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Shannon") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2() +
  ylim(0,3) +
  facet_wrap(~species)

div.plot


even.plot <- ggplot(cd, aes(x = Treatment, y = even, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Species Eveness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2() +
  ylim(0,1.3) +
  facet_wrap(~species)
even.plot

