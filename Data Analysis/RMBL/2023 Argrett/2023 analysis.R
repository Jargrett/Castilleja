# Emerald plant cover 2023
# Initiated: 8/6/23
# Completed: TBD

#####################
#                   # - Cover analysis of paired plots
#      Question     # - Total Castilleja cover
#                   # - NN analysis and richness
##################### - 



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

el.plants <- read.csv("Plant Data - Emerald Lake.csv")
sp.plants <- read.csv("Plant Data - Schofield Park.csv")
el.plot <- read.csv("")
sp.plot <- read.csv("")

