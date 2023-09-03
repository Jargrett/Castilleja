#Makena Septentrionalis Analysis

# Species Diversity
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)

case.cover <- read.csv("Combined septentrionalis plant data - Cover.csv")
case.env <- subset(case.cover, select=c(2,4,5,6)) # holding our environmental stuff
case.cover.nobg <- case.cover[-c(1:7)] #remove bare ground and annuals
case.cover.nobg[is.na(case.cover.nobg)] <- 0

rich <- specnumber(case.cover.nobg)
div <- diversity(case.cover.nobg, index = "shannon")
even <- diversity(case.cover.nobg, index = "shannon") / log(specnumber(case.cover.nobg))     
case.div <- cbind(case.env,div,even,rich)
case.div$Pair <- as.factor(case.div$Pair)

#exporting for combined analyisis
case.cov.div <- case.div[-c(2)] # removing pairs
write.csv(case.cov.div, "C:\\Users\\jargr\\Dropbox\\PC\\Desktop\\Data Analysis\\RMBL\\Makena Data\\CASE Diversity.csv")

aggregate(div ~ Treatment,
          data = case.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)
aggregate(rich ~ Treatment,
          data = case.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

aggregate(even ~ Treatment,
          data = case.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

#---------Analysis of the effect of treatment on diversity--------#

#Linear model -> Anova Type II 
div.lm <- lm(div ~ Treatment*Site*Pair, data = case.div)
summary(div.lm)
Anova(div.lm)

div.plot <- ggplot(case.div, aes(x =Site , y = div, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Shannon") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0.5,3)

div.plot

ggplot(case.div, aes(x = Treatment, y = div, color = Treatment)) +
  theme_bw() + 
  geom_point() + 
  stat_summary(fun.y=mean, position = "dodge", size = 1) + 
  # Lines by species using grouping
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="grey57") +
  ylab("Shannon Diversity")
#---------Analysis of the effect of treatment on richness--------#

#Linear model -> Anova Type II 
rich.lm <- lm(rich ~ Treatment*Site*Pair, data = case.div)
summary(rich.lm)
Anova(rich.lm)

rich.plot <- ggplot(case.div, aes(x = Site, y = rich, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species richness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0,20)
rich.plot

rich.plot <- ggplot(case.div, aes(x = Site, y = rich, color = Treatment)) +
  geom_point() +
  labs(x = "Treatment", y = "Species richness") +
  facet_wrap(~Pair)+
  #geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0,20)
rich.plot

ggplot(case.div, aes(x = Treatment, y = rich, color = Treatment)) +
  theme_bw() + 
  geom_point() + 
  stat_summary(fun.y=mean, position = position_dodge(0.3), width = 1,size = 1) +
  # Lines by species using grouping
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="grey57") +
  ylab("Species Richness") +
  facet_wrap(~Site) +
  ylim(0,20)
#---------Analysis of the effect of treatment on Evenness--------#

#Linear model -> Anova Type II 
even.lm <- lm(even ~ Treatment*Site*Pair, data = case.div)
summary(even.lm)
Anova(even.lm)

even.plot <- ggplot(case.div, aes(x = Site, y = even, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species evenness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0.6,1)

even.plot
