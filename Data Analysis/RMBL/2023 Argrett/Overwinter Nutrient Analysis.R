#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2023 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis

soil <- read.csv("Soil Nutrients - Overwinter 2023.csv")


nitrate <- ggplot(soil, aes(x = Litter, y = Ammonium)) +
  geom_point(aes(color = (Removal))) +
  labs(x = "Litter Treatment", y = "Ammonium")
  
nitrate

