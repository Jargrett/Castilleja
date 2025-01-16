#Initiated:1/2/2025
#Completed:
#For the analysis of Fall Biomass harvest

#Guide
#Anything after the # is a comment meant to explain or support your understanding and purpose of the code

#Here we set our working directory
#This establishes where on the computer we will be pulling data from
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Chelsea Data")

#Here we will load in packages for our research
#Packages allow for us to perform functions not in Base R
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis

#To install packages use this code template:
#install.packages("package_name")

#------------------------Data QA-QC-------------------------#
#We will now import our text files
#note how we no longer need to specify the directory due to (setwd())
chelsea <- read.csv("Host Specificity.csv")
str(chelsea) #allows us to check the structure of the data and see variable coding

#We will want to change some of our integers and characters to factors
#this will allow us to group and compare by these factor groupings
chelsea <- as.data.frame(unclass(chelsea),stringsAsFactors=TRUE)
summary(chelsea) #This allows us to view a summary of the dataframe

#For this analysis we are only concerned with Biomass data
#We are subsetting the dataset to only include the columns useful for this analysis
biomass <- subset(chelsea, select=c(1:8))

#------------------------Preliminary Visualizations-------------------------#
host.biomass <- ggerrorplot(biomass, x = "host", y = "Biomass_host",color = "parasite")
host.biomass

hemi.biomass <- ggerrorplot(biomass, x = "host", y = "biomass_parasite")
hemi.biomass
#------------------------Statistical Analysis-------------------------#

biomass.host <- lm(Biomass_host ~ host*parasite, data = biomass)
summary(biomass.host)
Anova(biomass.host) #this is the meat of what we want
emmip(biomass.host, parasite ~ host)
emmeans(biomass.host, pairwise ~ parasite|host)

biomass.hemi <- lm(biomass_parasite ~ host, data = biomass)
summary(biomass.hemi)
Anova(biomass.hemi)

