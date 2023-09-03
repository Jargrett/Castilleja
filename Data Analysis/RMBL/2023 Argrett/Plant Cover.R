#Percent Cover Analysis
#Initiated: 9/2/23

#setwd
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/2023 Argrett")

library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # linear regression
library(labdsv) 

#import EL and SP Data
EL <- read.csv("Emerald Lake Plant Data.csv")
SP <- read.csv("Schofield Park Plant Data.csv")

#Here we will attempt to summarize cover data per species for our collected data

#Emerald Lake
EL.post <- EL[EL$Collection.Point == 'Post',]
EL.matrix <- subset(EL.post, select = c('Code','Plot','Cover'))
EL.matrix <- matrify(EL.matrix)
EL.matrix$Average.cover <-apply(EL.matrix,1,mean)
EL.average <- subset(EL.matrix, select = c(41))

#
SP.post <- SP[SP$Collection.Point == 'Post',]
SP.matrix <- subset(SP.post, select = c('Code','Plot','Cover'))
SP.matrix <- matrify(SP.matrix)
SP.matrix$Average.cover <-apply(SP.matrix,1,mean)
SP.average <- subset(SP.matrix, select = c(41))
