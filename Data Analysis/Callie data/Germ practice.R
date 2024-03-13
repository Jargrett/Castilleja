# March Seed Germination Data Exploration
# initated: 3/13/24
# Last edited: 3/13/24

# Set working directory (this tells the computer which filezone we will be working from)
# change this every time you work from a new computer or filepath
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie data")

# Load-in packages
# for packages you do not have you can first type install.packages(package)

library(tidyverse) # for data working
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)
library(lme4)

#--------------Data working-------------------#

# First step is to load in data
germ <- read.csv("March Germination.csv")

# We need to clean the data and isolate the relavant columns for our analysis
# Remove columns = Description, notes (4,16) 
germ_clean <- germ [ -c(4,16)]

# Remove all rows with N/A or blanks present
# removing NA
germ_clean_NA <- germ_clean %>% drop_na()

#---------------------Visualization--------------------#

# Visualize Average Germ for LUTE (lupinus texensis) based off of seed source
Lute.plot <- ggplot(subset(germ_clean_NA, Species %in% "Lupinus texensis"), aes(x = Stratification.length , y = percent.germ, fill = Sterilized)) +
  geom_point() + 
  facet_wrap(~Sterilized)

Lute.plot
  
  
  
  
  
  
  