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
library(ggthemes)

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

#----BSA Graphs-----#
bsa.div <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun.y=mean, geom = "crossbar", position = position_dodge(1), size = 1, width = 0.25, col = "grey34") +
  facet_wrap(~Site) +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "khaki3")) +
  labs(x = "Castilleja septentrionalis", y = "Shannon Diversity") +
  ylim(1,2.5)

bsa.div

bsa.rich <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun.y=mean, geom = "crossbar", position = position_dodge(1), size = 1, width = 0.25, col = "grey34") +
  facet_wrap(~Site) +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "khaki3")) +
  labs(x = "Castilleja septentrionalis", y = "Species Richness") +
  ylim(3,16)

bsa.rich

bsa.even <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun.y=mean, geom = "crossbar", position = position_dodge(1), size = 1, width = 0.25, col = "grey34") +
  facet_wrap(~Site) +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "khaki3")) +
  labs(x = "Castilleja septentrionalis", y = "Species Evenness") +
  ylim(0.6,1)

bsa.even

bsa.plot <- ggarrange(bsa.div, bsa.rich, bsa.even, labels = c("A", "B","C"),nrow = 1, common.legend = TRUE, legend = "bottom")
bsa.plot
#--------------------graphs--------------------#
#Diversity
cali.div.plot <- ggplot(subset(cd.pair, species %in% "Castilleja linarifolia"), aes(x = Site, y = div, fill = Castilleja)) +
  geom_point() +
  labs(x = "Population", y = "Shannon Diversity") +
  facet_wrap(~species) +
  theme_classic2()+
  ylim(0,2.5)
cali.div.plot

cd.pair %>% filter(species=="Castilleja linarifolia") %>%
  mutate(population=paste(Site,Castilleja,sep="_")) %>%
  ggplot(aes(x=population,y=div,fill=Castilleja)) +
  geom_point() +
  theme_classic() 
  
cd_sum <- cd.pair %>% group_by(Site,Castilleja,species) %>%
  summarize(mean_div = mean(div),
            mean_rich = mean(rich),
            mean_even = mean(even))

cd.pair %>% filter(species=="Castilleja linarifolia") %>%
  mutate(population=paste(Site,Castilleja,sep="_")) %>%
  ggplot(aes(x=population,y=div,fill=Castilleja)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun="mean",color="lightblue",size=8,geom="point") +
  theme_classic() +
  theme(legend.position = "none")

case.div.plot <- ggplot(subset(cd.pair, species %in% "Castilleja septentrionalis"), aes(x = Site, y = div, fill = Castilleja)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = .3, position = position_jitter(seed = 1, width = .2)) +
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


#---------------error plots--------#

