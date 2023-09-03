#Nearest Neighbor Analysis
#Initiated: 8/29/23

#setwd
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/2023 Argrett")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # linear regression

#import EL and SP Data
EL <- read.csv("Emerald Lake Plant Data.csv")
SP <- read.csv("Schofield Park Plant Data.csv")

#First We need to see if EL sites nearest neighbor differ significantly enough
#this is to check as I have both pre and post collection data
#I want to confirm that NN data does not differ seasonally

#Lets begin by subsetting the necessary columns (NN,collection point,code,functionalgroup)
EL.NN<- subset(EL, select=c(3,5,10,16,17,19,21))
#I will do the same for Schofield however I will include only necessary coloums for now
SP.NN<- subset(SP, select=c(3,16,11,17,19,21))

#now lets look at NN by collection point for EL
str(EL.NN)
EL.NN$Collection.Point <- as.factor(EL.NN$Collection.Point)
EL.NN$Removal <- as.factor(EL.NN$Removal)
EL.NN$Code <- as.factor(EL.NN$Code)
EL.NN.Post <- EL.NN[EL.NN$Collection.Point == 'Post',]
EL.NN.Post <- EL.NN.Post[ -c(2)]

#after cleaning and shaping data frame we have:
#combined NN data for both sites (Post Removal)
NN <- rbind(EL.NN.Post,SP.NN)
NN$Removal[NN$Removal == 'C'] <- 'Control'
NN$Removal[NN$Removal == 'R'] <- 'Removal'
NN$Functional.Group[NN$Functional.Group == 'sedge'] <- 'Sedge'
NN$Functional.Group[NN$Functional.Group == 'environmental'] <- 'Environmental'
NN$Functional.Group[NN$Functional.Group == 'forb'] <- 'Forb'
NN$Functional.Group[NN$Functional.Group == 'grass'] <- 'Grass'
NN$Functional.Group[NN$Functional.Group == 'legum'] <- 'Legume'
NN$Functional.Group[NN$Functional.Group == 'hemiparasite'] <- 'Hemiparasite'
NN$Functional.Group[NN$Functional.Group == 'shrub'] <- 'Shrub'
#it looks like NN data changed between collection points.
#Not surprising as we have 
ggplot(NN, aes(x=Functional.Group, y=Nearest.Neighbor)) +
  geom_bar(stat="identity") +
  facet_wrap(~Site)



