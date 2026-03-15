setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")
#load in relevant packages
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(ggpubr)
library(rstatix)
library(indicspecies)
castilleja.cover <- readRDS("Processed Data/Total Castilleja Cover.rds")


#-------------------------Indicator Species Analysis---------------------------#
#In this analysis we will be assessing species that are found more often in one treatment group compared to another.
#This package takes into account both the relative abundance in a given plot as well as presence absence so we will use cover data
#we will assess this by species and sites but across years

#We will need a species matrix (castilleja removed) and vector that contains our presence absence plot info

#Avery
case.avery <- filter(castilleja.cover, site == "Avery")#filtering for a specific site
case.avery.matrix <- case.avery %>%#this segment selects just our species matrix and removes castilleja
  select(17:158)

case.avery.cast = case.avery$castilleja
case.avery.pair = case.avery$pair

case.avery.inv = multipatt(case.avery.matrix, case.avery.cast, func = "r.g",
                           control = how(blocks = case.avery.pair, nperm = 9999))

summary(case.avery.inv, alpha = 0.1)

#Emerald Lake
case.emerald <- filter(castilleja.cover, site == "Emerald Lake")
case.emerald.matrix <- case.emerald %>%
  select(17:158)

case.emerald.cast = case.emerald$castilleja
case.emerald.pair = case.emerald$pair

case.emerald.inv = multipatt(case.emerald.matrix, case.emerald.cast, func = "r.g",
                             control = how(blocks = case.emerald.pair, nperm=9999))

summary(case.emerald.inv, alpha = 0.1)

#Copper Creek
case.copper <- filter(castilleja.cover, site == "Copper Creek")
case.copper.matrix <- case.copper %>%
  select(17:158)

case.copper.cast = case.copper$castilleja
case.copper.pair = case.copper$pair

case.copper.inv = multipatt(case.copper.matrix, case.copper.cast, func = "r.g",
                            control = how(blocks = case.copper.pair, nperm=9999))

summary(case.copper.inv, alpha = 0.1)#Nothing significant for Copper Creek
#Absent: Cymopterus.lemmonii p = 0.0893 Present: Thalictrum.fendleri p = 0.0682, Campanula.petiolata p = 0.0849

#Now we do linariifolia
#Deer Creek 1
cali.dc1 <- filter(castilleja.cover, site == "Deer Creek 1")
cali.dc1.matrix <- cali.dc1 %>%
  select(17:158)

cali.dc1.cast = cali.dc1$castilleja
cali.dc1.pair = cali.dc1$pair

cali.dc1.inv = multipatt(cali.dc1.matrix, cali.dc1.cast, func = "r.g",
                         control = how(blocks = cali.dc1.pair, nperm=9999))

summary(cali.dc1.inv, alpha = 0.1)

#Deer Creek 2
cali.dc2 <- filter(castilleja.cover, site == "Deer Creek 2")
cali.dc2.matrix <- cali.dc2 %>%
  select(17:158)

cali.dc2.cast = cali.dc2$castilleja
cali.dc2.pair = cali.dc2$pair

cali.dc2.inv = multipatt(cali.dc2.matrix, cali.dc2.cast, func = "r.g",
                         control = how(blocks = cali.dc2.pair, nperm=9999))

summary(cali.dc2.inv, alpha = 0.1)#Castilleja Group: Delphinum.nuttalliianum p = 0.0050, Koeleria.macrantha = 0.0437
#Absent: Gayophytum.sp. p = 0.0512, Erigeron.speciosus p = 0.0661

#Johnson Hill
cali.johnson <- filter(castilleja.cover, site == "Johnson Hill")
cali.johnson.matrix <- cali.johnson %>%
  select(17:158)

cali.johnson.cast = cali.johnson$castilleja
cali.johnson.pair = cali.johnson$pair

cali.johnson.inv = multipatt(cali.johnson.matrix, cali.johnson.cast, func = "r.g",
                             control = how(blocks =cali.johnson.pair, nperm=9999))

summary(cali.johnson.inv, alpha = 0.1)

#Almont
cacr.Almont <- filter(castilleja.cover, site == "Almont")
cacr.cover.matrix <- cacr.Almont %>%
  select(17:158)

cacr.cover.cast = cacr.Almont$castilleja
cacr.cover.pair = cacr.Almont$pair

cacr.cover.inv = multipatt(cacr.cover.matrix, cacr.cover.cast, func = "r.g", 
                           control = how(blocks =cacr.cover.pair, nperm=9999))

summary(cacr.cover.inv, alpha = 0.1)
