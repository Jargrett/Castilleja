#########################
##      Packages       ##
#########################
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/In Lab Examples/Cristy Lab")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("nlme")

library(dplyr)
library(tidyr)
library(nlme)
library(ggplot2)

#########################
##   Data + Cleaning   ##
#########################

## Read in data as downloaded from the web
## doi:10.6073/pasta/ec929b4fec86eaf676d2a293694175d2
## https://www.cedarcreek.umn.edu/research/data/dataset?ple141
## I deleted meta data from csv files

rawdata <- read.csv("e141_Plant aboveground biomass data.csv")

## I like keeping the raw data df we imported around as I clean the data
## So we will simply assign it to another df we will work on
## This way, I can double check easily if I second guess myself at any point

biocon0 <- rawdata %>%
  tidyr::separate(Date, c("month", "Day", "year"), sep = "/") %>%
  select(year, month, Plot, Ring, 
         CO2_Treatment, Nitrogen_Treatment, 
         CountOfSpecies, CountOfGroup, Species, 
         Aboveground_Biomass)

# Here I will paste the 19 or 20 before the 2 other date digits

b <- biocon0$year
b <- if_else(b > 50, paste0(19, b), paste0(20, b))
biocon0$Year <- as.integer(b)
rm(b) # remove things we don't need!

str(biocon0)


## Issue with CO2 having more than 2 levels
unique(biocon0$CO2_Treatment)
# biocon0$CO2_Treatment <- if_else(biocon0$CO2_Treatment == "Cenriched", paste0("Cenrich"), paste0(biocon0$CO2_Treatment))


biocon0$Nitrogen_Treatment <- as.factor(biocon0$Nitrogen_Treatment)
biocon0$CO2_Treatment <- as.factor(biocon0$CO2_Treatment)


## The absolute value of year can affect your results when doing time series analysis
## We will add a column that has years listed starting from 1

ExpYr <- c(1:21)  # for experimental year
yr <- c(1998:2018) # for calendar year

df <- data.frame(yr, ExpYr) # data frame that we will join with big data frame
str(df)

biocon00 <- left_join(biocon0, df, by = c("Year" = "yr"))
  
## At the start of the experiment, biomass was sampled 2 times per year (June/August)
## We want to have one measurement per year
## Therefore, we will only use the August sampling for all years

biocon_long0 <- biocon00 %>%
  filter(month != 6 ) %>%
  select(-month, -year)   # Also getting rid of variables we will not use

biocon_long <- biocon_long0 %>%
  filter(!Species %in%  c("Miscellaneous litter" , "Miscellaneous Litter", 
           "Miscellaneous woody plants" , "Seed Contaminant" , 
           "Seed contaminant - forb" , "16 species Weeds", "16 Species Weeds",
           "Oak Leaves", "Moss", "Mosses & lichens", 
           "Seed contaminant - grass", "Bare ground"))  # these are species we don't want to consider
str(biocon_long)


biocon_wide0 <- biocon_long %>%
  group_by(Plot, Year) %>%
  summarise(TotalBio = sum(Aboveground_Biomass)) ## Making a column for Total Biomass

biocon_wide00 <- biocon_long %>%
  group_by(Species, Year, Plot, Ring) %>%
  mutate(grouped_id = row_number()) %>%  # Included this because of repeated IDs
  spread(Species, Aboveground_Biomass, fill=0) 

biocon_wide <- left_join(biocon_wide00, biocon_wide0, by = c("Year" = "Year", "Plot" = "Plot"))

str(biocon_wide)
  
## Nice, now we can get rid of clutter!
rm(biocon_long0, biocon_wide0, biocon_wide00, biocon0, biocon00, df, rawdata, ExpYr, yr) # remove things we don't need!

## Now we can get to work!



#########################
##    Data Analysis    ##
#########################

str(biocon_wide)

totbio_mod1<-lme(TotalBio~CountOfSpecies*ExpYr, random=~1|Plot, data=biocon_wide)
summary(totbio_mod1)
anova(totbio_mod1)


totbio_mod2<-lme(TotalBio~CountOfSpecies*ExpYr, random=~1|Plot, correlation=corAR1(form=~1|Plot), data=biocon_wide)
summary(totbio_mod2)
anova(totbio_mod2)

## 1st AR1
residuals.1<-residuals(totbio_mod1)
hist(residuals.1)
qqnorm(residuals.1)
qqline(residuals.1)

## 2nd with Corr structure
residuals.2<-residuals(totbio_mod2)
hist(residuals.2)
qqnorm(residuals.2)
qqline(residuals.2)

## How do the 2 correlation structures compare?
anova(totbio_mod1,totbio_mod2)
## How do these two compare?
## Should we try a different correlation structure?

## Other tests we can run?
ACF(totbio_mod1, maxLag = 4)

## Here you can find other correlation structure classes
## https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/corClasses.html

totbio_mod3<-lme(TotalBio~CountOfSpecies*ExpYr, random=~1|Plot, correlation=corARMA(form=~1|Plot, p=1, q=2), data=biocon_wide)
summary(totbio_mod3)
anova(totbio_mod3)

anova(totbio_mod2,totbio_mod3)

## We still have a not so great QQ plot
## Remember we manipulated more variables, but haven't used them in the model

totbio_modfull<-lme(TotalBio~CountOfSpecies*CO2_Treatment*Nitrogen_Treatment*ExpYr, random=~1|Plot, data=biocon_wide)
summary(totbio_modfull)
anova(totbio_modfull)

## *Thankfully* no 4 way interaction
## Let's get rid of non-sig interactions

totbio_modfull0<-lme(TotalBio~CountOfSpecies*CO2_Treatment*Nitrogen_Treatment*Year
                 -CountOfSpecies:CO2_Treatment:Nitrogen_Treatment:Year
                 -CO2_Treatment:Nitrogen_Treatment:Year 
                 -CO2_Treatment:CountOfSpecies:Nitrogen_Treatment,
                 random=~1|Ring, data=biocon_wide)
summary(totbio_modfull0)
anova(totbio_modfull0)


## Should we plot the data? 
plot <- ggplot()+
  geom_line(data=biocon_wide, 
               aes(x=Year,y=TotalBio, group=factor(CountOfSpecies), color=factor(CountOfSpecies)),
            stat = "summary", fun = mean) +
  facet_grid(~Nitrogen_Treatment) +
  labs(y = expression(paste("Aboveground biomass (g ", m^{-2},")"))) +
  theme_classic()

plot

