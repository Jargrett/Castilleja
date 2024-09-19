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
#Removing unneccesary coloumns
georgia.data <- georgia.parasites[ -c(18:22)]
#organize data strucutre, remove link columns, changing to factors
str(georgia.data)
georgia.data <- as.data.frame(unclass(georgia.data),stringsAsFactors=TRUE)

#-----------------Example Analysis-------------------#
# lets look into creating a graph showing the...
ggplot(data = georgia.data, aes(x = genus, y = inat_GA, fill = native)) +
  geom_bar(stat ="identity", position=position_dodge(1)) +
  labs(x = "Genera", y = "Elemental Occurances") +
  scale_fill_manual(values=c('black','lightgray')) +
  theme_classic()
