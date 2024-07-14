
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Makena Data")
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(emmeans)
library(rstatix)
library(ggrepel)
library(devtools)
library(pairwiseAdonis)
library(indicspecies)
library(simboot)
library(lme4)

# importing dataframe

dc24 <- read.csv("2024 Deer Creek 1 Counts.csv")
av24 <- read.csv("2024 Avery Counts.csv")
dc24.c <- dc24[ -c(1:6)] # removing unnecessary columns 
av24.c <- av24[ -c(1:6)]
dc24.env <- subset(dc24, select=c(2,4:6))
av24.env <- subset(av24, select=c(2,4:6))

rich <- specnumber(av24.c)

# Calculating Shannon diversity for plots
div <- diversity(av24.c, index = "shannon")

# Calculating species evenness for plots
even <- diversity(av24.c, index = "shannon") / log(specnumber(av24.c))     
  
#combined data set with environmental and calculated values
av24.div <- cbind(av24.env,rich,div,even)
av24.div <- av24.div[-c(41:81), ]


av24.diver <- lmer(div ~ Treatment + (1|Pair), data = av24.div)
summary(av24.diver)
Anova(av24.diver)

av24.rich <- lmer(rich ~ Treatment + (1|Pair), data = av24.div)
summary(av24.rich)
Anova(av24.rich)

av24.even <- lmer(even ~ Treatment + (1|Pair), data = av24.div)
summary(av24.even)
Anova(av24.even)




av24.diversity <- ggplot(av24.div, aes(x = Treatment, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Treatment), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("coral2","burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Shannon Diversity") +
  ylim(1,2.8)

av24.diversity

av24.rich <- ggplot(av24.div, aes(x = Treatment, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Treatment), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("coral2","burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Richness") +
  ylim(0,17)

av24.rich

dc24.even <- ggplot(av24.div, aes(x = Treatment, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Treatment), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("coral2","burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Evenness") +
  ylim(0.5,1)

dc24.even
