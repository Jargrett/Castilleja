setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Soil Moisture")

# Load-in packages
library(tidyverse) # for data working
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # Regression Analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models

sm23 <- read.csv("Emerald Lake Soil Moisture Data - 2023.csv")
sm23 <- as.data.frame(unclass(sm23),stringsAsFactors=TRUE)
sm24 <- read.csv("Emerald Lake Soil Moisture Data - 2024.csv")
sm24 <- as.data.frame(unclass(sm24),stringsAsFactors=TRUE)
sm25 <- read.csv("Emerald Lake Soil Moisture Data - 2025.csv")
sm25 <- as.data.frame(unclass(sm25),stringsAsFactors=TRUE)

#Create a new DF that groups by averages per plot
august.sm23 <- sm23 %>%
  group_by(plot) %>%
  summarise_at(vars(average_soil_moisture), list(asm23 = mean))

august.sm24 <- sm24 %>%
  group_by(plot) %>%
  summarise_at(vars(average_soil_moisture), list(asm24 = mean))

august.sm25 <- sm25 %>%
  group_by(plot) %>%
  summarise_at(vars(average_soil_moisture), list(asm25 = mean))

soil.moisture <- merge(august.sm23, august.sm24, by = "plot")
soil.moisture <- merge(soil.moisture, august.sm25, by = "plot")
write.csv(soil.moisture, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Soil Moisture/Moisture.csv", row.names=FALSE)



