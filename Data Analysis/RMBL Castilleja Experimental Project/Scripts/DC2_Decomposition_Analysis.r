setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data import, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)
library(litterfitter)#for k-curve fitting
library(ggplot2)
library(ggpubr)#extended functions for plottinglibrary(remotes)
library(ggpattern)
library(purrr)

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
  filter(mass_remaining <= o) %>% 
  mutate(across(c("time"), ~ round(.x, 2)))

#function to fit weibull
fit.weibull.nls = function(time_data, mass_data){
  fit = nls(mass_data ~ exp(- (time_data/beta)^alpha), 
            start = list(beta = 1, alpha = 1), 
            algorithm = "port", 
            lower = c(0.0001, 0.0001))
  return(fit)
}

#create dataframe to add predictions of models (for plotting)
t <- expand.grid(time_data = seq(0, 2.5, 0.01))

#half-life calculation:
half.life.calc = function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(2))^(1/pars[2])
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

#funtion to fit all data within the dataframe
fit_extract_predict <- function(df, t) {
  mod <- fit.weibull.nls(
    time_data = df$time,
    mass_data = df$mass_remaining)
  hl  <- half.life.calc(nls.mod = mod)
  mrt <- mrt.calc(nls.mod = mod)
  pred <- predict(mod, t)
  list(half_life = hl, mrt = mrt, pred = pred)
}

#create a results dataframe for each plot as well as a list object for predicitive line
weibull_results <- decomp %>%
  group_by(plot, litter, removal) %>%
  group_modify(~{
    out <- fit_extract_predict(.x, t)
    tibble(
      half_life = out$half_life,
      mrt = out$mrt,
      pred = list(out$pred)
    )
  }) %>%
  ungroup()

#create a new dataframe in long for for the predictive lines
pred_weibull <- weibull_results %>%
  select(plot, litter, removal, pred) %>%
  tidyr::unnest(pred) %>%
  mutate(time_data = rep(t$time_data, times = nrow(weibull_results))) %>%
  select(plot, litter, removal, time_data, pred)
write.csv(pred_weibull, "Processed Data/Weibull Predictions.csv", row.names=FALSE)

#subset for plotting
pred_cas <- pred_weibull %>% 
  filter(litter == "Castilleja") %>% 
  group_by(time_data, removal) %>% 
  summarise(mean = mean(pred), se = sd(pred)/sqrt(n()))
  
pred_com <- pred_weibull %>% 
  filter(litter == "Community") %>% 
  group_by(time_data, removal) %>% 
  summarise(mean = mean(pred), se = sd(pred)/sqrt(n()))

pred_mix <- pred_weibull %>% 
  filter(litter == "Mixed") %>% 
  group_by(time_data, removal) %>% 
  summarise(mean = mean(pred), se = sd(pred)/sqrt(n()))

#remove predictive lines to create final dataset
#cleaned dataframe for mrt and half_life plotting
summary_decomp <- weibull_results %>%
  select(plot, litter, removal, half_life, mrt)
write.csv(summary_decomp, "Processed Data/Litter Decomp Values.csv", row.names=FALSE)

#Analysis


#plotting
decomp_mrt <- summary_decomp %>% 
  group_by(litter, removal) %>% 
  summarise(mean = mean(mrt),
                   se = sd(mrt)/sqrt(n())) %>% 
  mutate(pat = ifelse(removal == "R", "stripe", "none"))


decomp_hl <- summary_decomp %>% 
  group_by(litter, removal) %>% 
  dplyr::summarise(mean = mean(half_life),
                   se = sd(half_life)/sqrt(n())) %>% 
  mutate(pat = ifelse(removal == "R", "stripe", "none"))

#
ggplot() +
  geom_line(data = pred_cas, aes(time_data, mean), color = "#4b3b40") +
  geom_line(data = pred_com, aes(time_data, mean), color = "#9b7e46") +
  geom_line(data = pred_mix, aes(time_data, mean), color = "#808f87") +
  facet_wrap(~ removal) +
  theme_pubr() +
  labs(x = "time since deployment", y = "proportion of litter mass remaining")



#Mean Residence Time
ggplot(decomp_mrt, aes(x = litter , y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  ylim(0, 2.5) +
  labs(x = "Litter Type", y = "Mean Resdience Time (yrs)")


#Half Life
ggplot(decomp_hl, aes(x = litter , y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "top", panel.grid = element_blank()) +
  guides(pattern = "none") +
  ylim(0, 0.5) +
  labs(x = "Litter Type", y = "Half Life (yrs)")




#Analysis
hl.lme <- lm(half_life ~ litter*removal, data = summary_decomp)
summary(hl.lme)
Anova(hl.lme) #No significance

mrt.lme <- lm(mrt ~ litter*removal, data = summary_decomp)
summary(hl.lme)
Anova(hl.lme) #No significance

