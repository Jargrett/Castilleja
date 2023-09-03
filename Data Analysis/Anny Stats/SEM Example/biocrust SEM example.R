#Biocrust SEM example
#AC 210411 adapted from ES

#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/SEM Example")
#Load packages (both new packages)
library(lavaan)
library(semPlot)


#Import data and wrangle--------------------
SEM.dat <- read.csv("biocrust soil properties.csv")

head(SEM.dat)

#Make new "Trt" variable
SEM.dat[SEM.dat$Treatment=="C","Trt"] <- 0 #Control is coded "0"
SEM.dat[SEM.dat$Treatment=="NC","Trt"] <- 1 #Stomp (Not Control) is coded "1"

#Split data into shrubland and grassland only datasets
SEMS <- SEM.dat[which(SEM.dat$Site == "S"),]
SEMG <- SEM.dat[which(SEM.dat$Site == "G"),]

#Specify model to test direct and indirect effects on infiltration depth----------
mod.chl <-" 
Depth ~ c*Trt

Chloro ~ a*Trt
Depth ~  b*Chloro

# indirect effect (a*b)
indirect := a*b
total := c + (a*b)"

#Test SEM model on shrubland and grassland plots separately-------------
# this code only gives parameter estimates, not model test results, 
# because it is saturated (no df remain to test model fit)
 
# shrubland
mod.chlS <- sem(mod.chl, data = SEMS)
summary(mod.chlS,  rsquare = T, standardized = T)
# only "a" is significant path
#graph the sem model
semPaths(mod.chlS,"std")

# grassland
mod.chlG <- sem(mod.chl, data = SEMG)
summary(mod.chlG,  rsquare = T, standardized = T)
# "a" and "b" are significant paths, "c" marginally insignificant
# total indirect effect is significant and positive 
#graph the sem model
semPaths(mod.chlG,"std")


