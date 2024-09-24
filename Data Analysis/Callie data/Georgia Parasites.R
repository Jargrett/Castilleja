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


#Load-in Data
georgia.parasites <- read.csv("Georgia Parasites - Master.csv")
#Removing unnecessary coloumns
georgia.data <- georgia.parasites[ -c(18:22)]
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
  filter((DNR_status_GA != "N/A"))

summary(filter.data)
 
filter.data$DNR_status_GA[filter.data$DNR_status_GA=="S1?"] <- "S1" 
filter.data$DNR_status_GA[filter.data$DNR_status_GA=="S2?"] <- "S2"
filter.data$DNR_status_GA[filter.data$DNR_status_GA=="S3?"] <- "S3"

ggplot(filter.data, aes(x = DNR_status_GA)) + 
  geom_bar()


#-----------------Example Analysis-------------------#
# lets look into creating a graph showing the...
ggplot(data = georgia.data, aes(x = genus, y = inat_GA, fill = native)) +
  geom_bar(stat ="identity", position=position_dodge(1)) +
  labs(x = "Genera", y = "Elemental Occurances") +
  scale_fill_manual(values=c('black','lightgray')) +
  theme_classic()


