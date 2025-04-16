setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")
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

sm.23 <- read.csv("Emerald Lake Plant Data - 2023.csv")
EL.24 <- read.csv("Emerald Lake Plant Data - 2024.csv")