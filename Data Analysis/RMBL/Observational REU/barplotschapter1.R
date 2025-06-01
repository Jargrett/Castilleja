setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")
library(ggplot2)
library(devtools)
library(remotes)
library(ggpubr)
library(ggpattern)
library(tidyverse)#for data wrangling and restructuring
castilleja.cover <- read.csv("castilleja cover complete.csv")
castilleja.cover$castilleja[castilleja.cover$castilleja == "Control"] <- "Absent"
castilleja.cover$castilleja[castilleja.cover$castilleja == "Castilleja"] <- "Present"
castilleja.cover <- as.data.frame(unclass(castilleja.cover),stringsAsFactors=TRUE)
castilleja.cover$year = as.factor(castilleja.cover$year)
#Standard error calculations
castilleja.div <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(div),
                   se = sd(div)/sqrt(n()))
castilleja.rich <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(rich),
                   se = sd(rich)/sqrt(n()))
castilleja.even <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(even),
                   se = sd(even)/sqrt(n()))

cols<-ifelse(COroots$Invasion=="Invaded", "#E69F00", "#009E73")
alfa<-as.factor(0.5)
COroots$grp = paste(COroots$Invasion, COroots$Treatment)


ggplot(data = castilleja.div, aes(x = year, y = mean, fill = castilleja)) +
  geom_col(position = "dodge", width = 0.9, linewidth = 0.5, alpha = 0.5, color = "black", size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.9), width = 0.07) +
  theme_pubr() +
  scale_fill_manual(values=c("#e4b98c", "#80a784")) +
  labs(x = "Sampling year", y = "Shannon Diversity of co-occuring species") +
  ylim(0,2)

  
  
