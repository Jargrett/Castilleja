#Makena Scratch

# Species Diversity
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)

deer.creek.cover <- read.csv("Combined deer creek data - Cover.csv")
dc.env <- subset(deer.creek.cover, select=c(2,4,5,6)) # holding our environmental stuff
dc.cover <- deer.creek.cover[ -c(1:7,27,28)]
dc.cover.nobg <- deer.creek.cover[ -c(1:7,27,28)] #remove bare ground and annuals
dc.cover.nobg[is.na(dc.cover.nobg)] <- 0

# Calculating Shannon diversity for plots using cover data
rich <- specnumber(dc.cover.nobg)
div <- diversity(dc.cover.nobg, index = "shannon")
even <- diversity(dc.cover.nobg, index = "shannon") / log(specnumber(dc.cover.nobg))     
dc.div <- cbind(dc.env,rich,div,even)

#---------Analysis of the effect of treatment on diversity--------#
#Linear model -> Anova Type II 
div.lm <- lm(div ~ Treatment*Pair, data = dc.div)
summary(div.lm)
Anova(div.lm)

aggregate(div ~ Treatment,
          data = dc.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)
aggregate(rich ~ Treatment,
          data = dc.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

aggregate(even ~ Treatment,
          data = dc.div,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)


#Plot
div.plot <- ggplot(dc.div, aes(x = Site, y = div, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Shannon") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0.5,3)

div.plot

ggplot(dc.div, aes(x = Treatment, y = div, color = Treatment)) +
  theme_bw() + 
  geom_point() + 
  stat_summary(fun.y=mean, position = position_dodge(0.3), width = 1,size = 1) +
  # Lines by species using grouping
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="grey57") +
  ylab("Shannon Diversity") +
  facet_wrap(~Site) +
  ylim(0.5,3)

#---------Analysis of the effect of treatment on evenness--------#

#linear regression and Anova type II
even.lm <- lm(even ~ Treatment*Pair, data = dc.div)
summary(even.lm)
Anova(even.lm)

even.plot <- ggplot(dc.div, aes(x = Site, y = even, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species eveness") +
  geom_jitter(shape=16, alpha = 0.3) +
  theme_classic2()

even.plot

#---------Analysis of the effect of treatment on richness--------#

#linear regression and Anova type II
#Strong significance of the effect of treatment on species richness
rich.lm <- lm(rich ~ Treatment*Pair, data = dc.div)
summary(rich.lm)
Anova(rich.lm)


#This could just be due to the presence of Castilleja so lets remove it
dc.rich.cast <- deer.creek.cover[ -c(1:7,17,27,28)]
rich.nocast <- specnumber(dc.rich.cast)
div.nocast <- diversity(dc.rich.cast, index = "shannon")
even.nocast <- diversity(dc.rich.cast, index = "shannon") / log(specnumber(dc.rich.cast))     
dc.div.nocast <- cbind(dc.env,rich.nocast,div.nocast,even.nocast)

#linear regression and Anova type II
nocast.rich.lm <- lm(rich.nocast ~ Treatment*Pair, data = dc.div.nocast)
summary(nocast.rich.lm)
Anova(nocast.rich.lm)


even.plot <- ggplot(dc.div, aes(x = Site, y = rich, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species richness") +
  geom_jitter(shape=16, alpha = 0.3) +
  theme_classic2()

even.plot

ggplot(dc.div, aes(x = Treatment, y = div, color = Treatment)) +
  theme_bw() + 
  geom_point() + 
  stat_summary(fun.y=mean, position = position_dodge(0.3), width = 1,size = 1) +
  # Lines by species using grouping
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="grey57") +
  ylab("Shannon") +
  facet_wrap(~Site) +
  ylim(0,3)

nocast.even.plot <- ggplot(dc.div.nocast, aes(x = Treatment, y = rich.nocast, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Species richness") +
  geom_jitter(shape=16, alpha = 0.3) +
  theme_classic2()

nocast.even.plot

#---------Analysis of the effect of treatment on bare ground--------#
#Does the presence of Castilleja have an impact on the proportion of bare ground
dc.bare <- subset(deer.creek.cover, select=c(1:7))
#linear regression and Anova type II
bare.lm <- lm(Bare.ground ~ Treatment*Pair, data = dc.bare)
summary(bare.lm)
Anova(bare.lm)

bare.plot <- ggplot(dc.bare, aes(x = Treatment, y = Bare.ground, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Bare ground") +
  geom_jitter(shape=16, alpha = 0.3) +
  theme_classic2()

bare.plot

#----Functional Group Analysis----#




