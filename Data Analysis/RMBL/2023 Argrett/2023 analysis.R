# Emerald plant cover 2023
# Initiated: 8/6/23
# Completed: TBD

#####################
#                   # - Cover analysis of paired plots
#      Question     # - Total Castilleja cover
#                   # - NN analysis and richness
##################### - 

setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/2023 Argrett")

# Load-in packages
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(ggpubr)
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)

plants <- read.csv("Plant Data - Emerald Lake.csv")

# visualize
ggplot(data=plants, aes(x= Code, y= Cover)) +
  geom_bar(stat="identity")

