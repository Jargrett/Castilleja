#Georgia Parasites
#Initiated: 10/18/24

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie data")

#load-in Packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(gridExtra)
library(webr)

#Load-in Data
georgia.parasites <- read.csv("Georgia Parasites - Master.csv")
georgia.parasites$taxa <- interaction(georgia.parasites$genus, georgia.parasites$species, sep = " ")
georgia.parasites <- georgia.parasites %>% relocate(taxa)
#Removing unnecessary coloumns
georgia.data <- georgia.parasites[ -c(31:35)]
#organize data strucutre, remove link columns, changing to factors
str(georgia.data)
georgia.data <- as.data.frame(unclass(georgia.data),stringsAsFactors=TRUE)
summary(georgia.data)
#Question 1:
# What is the current breakdown of status for Georgia's root hemiprasitic plants
# We only have data for 20 out of 40 Root HP:
# of which S1/S1? = 11, S2/S2? = 5, S3 = 1, SH = 1
filter.data <- georgia.data %>% 
  filter((parasitic_habit == "Root Hemiparasite")) %>% 
  filter((native == "Yes"))

summary(filter.data)

ggplot(filter.data, aes(x = ecosystem)) + 
  geom_bar() +
  geom_bar(aes(x = DNR_status_GA))

PieDonut(georgia.data, aes(x = parasitic_habit, y = genus), r0=0.3, explode=0,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)


PieDonut(georgia.data, aes(x = DNR_status_lump, y = DNR_monitor), r0=0.3, explode=0,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)

PieDonut(georgia.data, aes(x = ecosystem_lump, y = parasitic_habit), r0=0.3, explode=0,showRatioThreshold =.001, labelpositionThreshold=.21,
         ratioByGroup = TRUE)

S1.data <- filter.data %>% 
  filter((DNR_status_GA == "S1"))

agalinis.data <- filter.data %>% 
  filter((genus == "Agalinis"))

PieDonut(S1.data, aes(x = south_inat, y =species ),
         ratioByGroup = F)
#-----------------Example Analysis-------------------#
# lets look into creating a graph showing the...
ggplot(data = filter.data, aes(x = genus, y = inat_GA)) +
  geom_bar(stat ="identity", position=position_dodge(1)) +
  labs(x = "Genera", y = "Proportion of observatoins in GA") +
  scale_fill_manual(values=c('black','lightgray')) +
  theme_classic()


