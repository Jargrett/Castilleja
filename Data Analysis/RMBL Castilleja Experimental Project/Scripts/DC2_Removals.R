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

#import data
removal_cover <- read.csv("Raw Data/EL Removal Data - Cover.csv")
removal_meta <- read.csv("Raw Data/EL Removal Data - Meta.csv")

#remove castilleja and environmental rows for analysis
removal_cover %<>% 
  mutate(percent_cover = cover/100) %>% 
  subset (removed !="y") %>% 
  drop_na(count)



removal_clean <- subset(removal_cover, select = c('year','removal_species','code','percent_cover'))


removal_24 <- removal_clean %>% 
  filter(year == "2024") %>%
  select(-c(year))

removal_25 <- removal_clean %>% 
  filter(year == "2025") %>%
  select(-c(year))

#convert to matrix format for diversity calculations
library(labdsv)#enables restructuring for ecological analysis
removal.24.matrix <- matrify(removal_24)
removal.25.matrix <- matrify(removal_25)

#---------------Diversity Calculations---------------#
library(vegan)#for calculating diversity
# Calculating Shannon diversity for plots
div.24 <- diversity(removal.24.matrix, index = "shannon")
div.25 <- diversity(removal.25.matrix, index = "shannon")


# Calculating species richness for plots
rich.24 <- specnumber(removal.24.matrix)
rich.25 <- specnumber(removal.25.matrix)

# Calculating species evenness for plots 
even.24 <- diversity(removal.24.matrix, index = "shannon") / log(specnumber(removal.24.matrix))
even.25 <- diversity(removal.25.matrix, index = "shannon") / log(specnumber(removal.25.matrix))


rem_24 <- cbind(removal_meta,div.24,rich.24,even.24)
rem_24 %<>% 
  plyr::mutate(year = '2024') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.24, rich = rich.24, even = even.24)

rem_25 <- cbind(removal_meta,div.25,rich.25,even.25)
rem_25 %<>% 
  plyr::mutate(year = '2025') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.25, rich = rich.25, even = even.25)

rem_div <- rbind(rem_24, rem_25)


rich.lm <- lm(div ~ year*biomass_removed, data = rem_div)
summary(rich.lm)
Anova(rich.lm)#:removal:year p = 0.0008, Chisq = 16.8171, df = 3


rem_mean <- rem_div %>% 
  group_by(year) %>% 
  dplyr::summarise(div_mean = mean(div), div_se = sd(div)/sqrt(n()),
                   rich_mean = mean(rich), rich_se = sd(rich)/sqrt(n()),
                   even_mean = mean(even), even_se = sd(even)/sqrt(n()))

ggplot(rem_mean, aes(x = year, y = rich_mean)) +
  geom_point(data = rem_div, aes(x = year, y = rich, color = species_removed),
             position = position_jitterdodge(0.2, dodge.width = .3), size = 2, alpha = 0.5) +
  geom_point(shape = 18, size = 4.5, position = position_dodge(width = 0.2))+
  geom_errorbar(aes(ymin = rich_mean - rich_se, ymax = rich_mean + rich_se),
              position = position_dodge(width = 0.2), width = 0.07)


# Mean change in richness for R plots 2024 to 2025
r_change <- rem_div %>%
  pivot_wider(id_cols = plot, names_from = year, values_from = rich) %>%
  mutate(rich_change = `2025` - `2024`) %>%
  summarise(mean_change = mean(rich_change),
            sd_change   = sd(rich_change))

# Compare to mean change in your main removal plots over the same window

