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

#organize data strucutre, remove link columns, changing to factors

str(georgia.parasites)
georgia.parasites$genus <- as.factor(georgia.parasites$genus)
georgia.parasites$species <- as.factor(georgia.parasites$species)
georgia.parasites$species <- as.factor(georgia.parasites$species)