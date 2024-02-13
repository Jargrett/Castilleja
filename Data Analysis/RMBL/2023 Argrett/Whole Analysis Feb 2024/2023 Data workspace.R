# plant/plot data from 2023
# Initiated: 2/7/24
# Completed: TBD

#set workind directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2023 Argrett/Whole Analysis Feb 2024")

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
library(lme4)

el.plants <- read.csv("Plant Data - Emerald Lake.csv")
sp.plants <- read.csv("Plant Data - Schofield Park.csv")
el.plot <- read.csv("Plot Data - Emerald Lake.csv")
sp.plot <- read.csv("Plot Data - Schofield Park.csv")

#Combine the plot and plant data for analysis
el.plants <- el.plants [ -c(6,7,9,11,15)]
sp.plants <- sp.plants [ -c(6,7,8,10,12)]
el.plot <- el.plot[ -c(2,4,6,8:11)]
sp.plot <- sp.plot[ -c(2,5,6,8:10)]
comb.plot <- rbind(el.plot,sp.plot)
comb.plants <- rbind(el.plants,sp.plants)
comb.plant <- subset(comb.plants, comb.plants$Collection.Point=='Post')

#-----------------Checking data structure-----------------#
str(comb.plant)
comb.plant$Code <- as.factor(comb.plant$Code)
comb.plant$Functional.Group <- as.factor(comb.plant$Functional.Group)
comb.plant$Removal <- as.factor(comb.plant$Removal)
comb.plant$Life.History <- as.factor(comb.plant$Life.History)


#------------Diversity Analysis------------#
# here we will be working with the count column
# we will need to conververt data to a matrix format
comb.c <- subset(comb.plant, select = c('Plot','Code','Count'))
comb.count<- matrify(comb.c)

# Calculating the species richness for plots
rich <- specnumber(comb.count)

# Calculating Shannon diversity for plots
div <- diversity(comb.count, index = "shannon")

# Calculating species evenness for plots
even <- diversity(comb.count, index = "shannon") / log(specnumber(comb.count))  
comb.div <- cbind(comb.plot,div,rich,even)


