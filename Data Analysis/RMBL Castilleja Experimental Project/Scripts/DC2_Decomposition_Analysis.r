setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data import, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)
library(litterfitter)#for k-curve fitting


conflicted::conflicts_prefer(dplyr::filter)

#load-in the data
within1 <- read.csv("Raw Data/Litter Decomposition - 2024 Within Season.csv")
within2 <- read.csv("Raw Data/Litter Decomposition - 2025 Within Season.csv")
overwinter1 <- read.csv("Raw Data/Litter Decomposition - 2024 Overwinter.csv")
overwinter2 <- read.csv("Raw Data/Litter Decomposition - 2025 Overwinter.csv")
year1 <- read.csv("Raw Data/Litter Decomposition - 2024 Full Year.csv")
year2 <- read.csv("Raw Data/Litter Decomposition - 2025 Full Year.csv")
twoyear <- read.csv("Raw Data/Litter Decomposition - 2023-2025 Two Year.csv")
found <- read.csv("Raw Data/Litter Decomposition - Found Bags.csv")

#Data Cleaning and Restructuring
 within1 %<>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing)))

 within2 %<>% 
   dplyr::select(-c((missing)))
 
overwinter1 %<>% 
  filter(missing != "Yes") %>% 
  dplyr::select(-c((missing))) %>%
  dplyr::select(-c((redo_coin_litter_dry_weight))) %>% 
  dplyr::select(-c((diff))) 

overwinter2 %<>% 
  filter(ups_missing != "yes") %>% 
  dplyr::select(-c((ups_missing))) %>% 
  filter(missing != "yes") %>% 
  dplyr::select(-c((missing)))

year1 %<>% 
  filter(missing != "Yes") %>%
  dplyr::select(-c((redo_coin_litter_dry_weight))) %>% 
  dplyr::select(-c((missing)))

year2 %<>% 
  filter(missing != "yes") %>% 
  dplyr::select(-c((missing)))

twoyear %<>% 
  filter(missing != "yes") %>% 
  dplyr::select(-c((missing)))

#Combine datasets
decomp <- bind_rows(within1, within2 , overwinter1, overwinter2, year1, year2, twoyear)
decomp <- as.data.frame(unclass(decomp),stringsAsFactors=TRUE)
#Calculate Mass Remaining and time
o <- 0.91
decomp %<>% mutate(mass_remaining = final_dry_weight/initial_dry_weight) %>% 
  mutate(time = deployment_duration/365) %>% 
  drop_na(time) %>% 
  filter(mass_remaining <= o)

#Filter by litter type
#Mixed
decomp.mixed <- filter(decomp, litter == "Mixed")
decomp.mixed.removal <- filter(decomp, litter == "Mixed") %>% 
  filter(removal == "R")
decomp.mixed.control <- filter(decomp, litter == "Mixed") %>% 
  filter(removal == "C")

#Castilleja
decomp.castilleja <- filter(decomp, litter == "Castilleja")
decomp.castilleja.removal <- filter(decomp, litter == "Castilleja") %>% 
  filter(removal == "R")
decomp.castilleja.control <- filter(decomp, litter == "Castilleja") %>% 
  filter(removal == "C")

#Community
decomp.community <- filter(decomp, litter == "Community")
decomp.community.removal <- filter(decomp, litter == "Community") %>% 
  filter(removal == "R")
decomp.community.control <- filter(decomp, litter == "Community") %>% 
  filter(removal == "C")

#Originally from Melanie
#This is from Cornwell and Weedon, 2013 supplement. 
#Weibull fitting function, half-life and mrt calculation.
fit.weibull.nls = function(time_data, mass_data){
  fit = nls(mass_data ~ exp(- (time_data/beta)ˆalpha), 
            start = list(beta = 1, alpha = 1), 
            algorithm = "port", 
            lower = c(0.0001, 0.0001))
  return(fit)
}

#create dataframe to add predictions of models (for plotting)
s <-expand.grid(time_dat = seq(0, 2.5, 0.1))

#half-life calculation:
half.life.calc = function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(2))ˆ(1/pars[2])
  names(hl) ="half.life"
  return(hl)
}
#mean residence time calculation:
mrt.calc = function(nls.mod){
  pars= coef(nls.mod)
  mrt=pars[1] * gamma(1+(1/pars[2]))
  names(mrt)="mrt"
  return(mrt)
}

#Subset data into plot/litter type combinations () and fit Weibull.
#Take a look at fit for each. Calculate halflife and mrt.

#Mixed Control
mxd.control <- fit.weibull.nls(time_data = decomp.mixed.control$mass_remaining,
                                mass_data = decomp.mixed.control$proportion_remaining)
summary(mxd.control)

#Mixed removal
mxd.removal <- fit.weibull.nls(time_data = decomp.mixed.removal$mass_remaining,
                                  mass_data = decomp.mixed.removal$proportion_remaining)
summary(mxd.removal)

#Castilleja Control
cst.control <- fit.weibull.nls(time_data = decomp.castilleja.control$mass_remaining,
                               mass_data = decomp.castilleja.control$proportion_remaining)
summary(cst.control)

#Mixed removal
cst.removal <- fit.weibull.nls(time_data = decomp.castilleja.removal$mass_remaining,
                               mass_data = decomp.castilleja.removal$proportion_remaining)
summary(cst.removal)

#Mixed Control
com.control <- fit.weibull.nls(time_data = decomp.community.control$mass_remaining,
                               mass_data = decomp.community.control$proportion_remaining)
summary(com.control)

#Mixed removal
com.removal <- fit.weibull.nls(time_data = decomp.community.removal$mass_remaining,
                               mass_data = decomp.community.removal$proportion_remaining)
summary(com.removal)



s$pred_AM1B1<-predict(AM1B1_weibull, s)



ggplot(AM1B1, mapping= aes(x=years_deployed, y= proportion_remaining)) +
  theme_bw() +
  theme(panel.border = element_rect(colour="black", size=.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line =
          element_line(colour = "black", size=.4)) +
  xlab("Years deployed")+
  ylab(expression(proportion~litter~remaining)) +
  geom_point(shape = 1) +
  geom_line(data = s, aes(x= time_data, y = pred_AM1B1), color = "blue") +
  ggtitle("1B1 AM")
