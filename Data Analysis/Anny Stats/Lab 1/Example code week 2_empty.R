#Week 2 Example Code#
#AC 201230

#A few housekeeping things-----------

#Set working directory
setwd("C:/Users/jargr/Desktop/Anny Stats/Lab 1")

#Make sure clean environment


#Load libraries
#If you don't have these libraries, install first
library(ggplot2)
library(tidyverse)
library(ggpubr)

#Import data and check--------------
penguins <- read.csv("penguins.csv")

#Check data structure
str(penguins)

#Notice that 'year' is imported as integer, which technically is correct, but unhelpful for us
#Change 'year' to factor
penguins$year <- as.factor(penguins$year)
penguins1 <- na.omit (penguins)
penguins$body_mass_g <- as.numeric(penguins$body_mass_g)

#Request frequencies for factors--------------
#How many of each penguin?
table(penguins$species)
#Tidyverse alternative
count(penguins.species)
#Tidyverse piped alternative
penguins %>% count(species)

#How many of each island?
table(penguins$island)

#How many of each penguin each year?
table(penguins$species,penguins$year)


#Data manipulation-------------------
#Subsetting: Make 3 dataframes each containing one of the 3 penguin species
#There are multiple ways to do this
adelie_sp <- subset(penguins, species == "Adelie")
chinstrap_sp <- subset(penguins, species == "Chinstrap")
Gentoo_sp <- subset(penguins, species == "Gentoo")
#Tidyverse


adelie_sp_female_dream <- filter(adelie_sp, sex == "female", island == "Dream")
count(adelie_sp_female_dream) #answer to question 1

#What if we only want female adelie penguins that weigh >3500g?
adelie_sp_large <- subset(adelie_sp, body_mass_g > 3500)
#Creat new variables
#Let's say bill area = bill length * bill depth

penguins$bill_area_mm <- penguins$bill_length_mm*penguins$bill_depth_mm
#Summary statistics-----------------
#Can get overall summary fast

#Note NA's, how 'year' is handled

#Simple summary stats functions
mean()
sd()
quantile(penguins$flipper_length_mm)#why error?
quantile(penguins$flipper_length_mm,na.rm = T)

mean(adelie_sp$bill_depth_mm, na.rm = T) # answer to 2
mean(chinstrap_sp$bill_depth_mm, na.rm = T) # answer to 2
mean(Gentoo_sp$bill_depth_mm, na.rm = T)# answer to 2

sd(penguins$body_mass_g, na.rm = T)
mean(penguins$body_mass_g, na.rm = T)


penguins2008 <- filter(penguins, year == "2008")

penguins2008 %>%
  ggplot( aes(x=bill_area_mm, fill=species)) +
  geom_histogram(color="black", alpha=0.5, position = 'identity') +
  scale_fill_grey() +
  theme_classic2()# answer to 3 


sd <- sd(penguins1$body_mass_g)
##pnegu <- summary(penguins1, measurevar="body_mass_g", groupvars=c("species","sex"))

Q4<-aggregate(body_mass_g~sex+species,data=penguins1,FUN=mean_sd)

ag <- aggregate(penguins, function(body_mass_g) c(mean = mean(body_mass_g), sd = sd(body_mass_g)))

Q4 %>%
  ggplot (aes(x = species, y = body_mass_g, fill = species ))

penguins1 %>%
  ggplot (aes(x = species, y = body_mass_g, fill = species )) +
  scale_fill_brewer(palette="Dark2") + 
  geom_bar(position = "dodge", stat = "summary", fun = "mean_sd", na.rm = T ) +
  geom_errorbar(aes(ymin=body_mass_g-sd, ymax=body_mass_g+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~sex) +
  theme_minimal() # answer to 4 

penguins2 <- penguins1 %>%
  group_by(sex,species) %>%
  summarize(mean = mean(body_mass_g), sd = sd(body_mass_g))



  ggplot (penguins2, aes(x= species, y = mean, fill = species)) + # answer to 4
    scale_fill_brewer(palette="Dark2") + 
    geom_bar(position = "dodge", stat = "summary", fun = "mean", na.rm = T ) +
    print(labs(y = "Mean Body Mass(g)", x = "Species")) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9)) +
    facet_grid (~sex) + 
    theme_minimal()
  
#Summarize by group?
#average mass of male vs. female penguins

#break down further by species

#tidyverse alternative (group, then summarize)


#Plot data distributions------------------
#Basic histograms

#Basic box plots


#Multiple histograms
#Together

#In multiple panels


#Nice and easy error plots with ggpubr

#Note the y axis truncation!

