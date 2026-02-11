setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Pollenator Project")

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
library(magrittr)
library(rsthemes)

community <- read.csv("plantsurvey_castilleja - plantsurvey_Raw.csv")
network <- read.csv("network_castilleja - network_castilleja.csv")

#coding for castileja
network %<>% 
  group_by(transect,segment) %>% 
  mutate(cas_yn = str_detect(plant,"Castilleja"),
         cast_presence = ifelse(any(cas_yn==TRUE),1,0))

community %<>% 
  group_by(transect,segment) %>% 
  mutate(cas_yn = str_detect(plant,"Castilleja"),
         cast_presence = ifelse(any(cas_yn==TRUE),1,0))

