setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
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
library(magrittr)#for data wrangling and restructuring
library(plyr)
library(ggnewscale)

#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicted::conflicts_prefer(plyr::mutate)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::summarise)



#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))