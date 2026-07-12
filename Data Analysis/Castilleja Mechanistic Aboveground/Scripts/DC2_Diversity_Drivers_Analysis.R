setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(rstatix)
library(magrittr)#for data wrangling and restructuring
library(conflicted)
library(broom)
library(broom.mixed)
library(forcats)

conflicts_prefer(lme4::lmer)
# --- overall model ---------------------------------------------------------
drivers <- readRDS("Processed Data/Drivers.rds")
mod_all <- lmer(rich ~ removal + litter + elevation + soil_moisture +
                  disturbance_distance + factor(year) + (1|block/pair/plot),
                data = drivers)

overall <- broom.mixed::tidy(mod_all, effects = "fixed", conf.int = TRUE) %>%
  filter(!grepl("factor\\(year\\)", term)) %>%
  mutate(panel = "Overall")

per_year <- drivers %>%
  group_by(year) %>%
  group_modify(~ broom.mixed::tidy(
    lmer(rich ~ removal + litter + elevation + soil_moisture +
           disturbance_distance + (1|block/pair),
         data = .x),
    effects = "fixed", conf.int = TRUE)) %>%
  ungroup() %>%
  mutate(panel = paste("Year", year))

per_year

coefs <- bind_rows(per_year, overall) %>%
  filter(term != "(Intercept)") %>%
  filter(!grepl("^litter", term)) %>%
  mutate(
    type = ifelse(grepl("removal", term), "Experimental", "Environmental"),
    term = recode(term,
                  "removalRemoved"       = "Castilleja removal",
                  "elevation"            = "Elevation",
                  "soil_moisture"        = "Soil moisture",
                  "disturbance_distance" = "Meadow edge"),
    panel = factor(panel, levels = c("Year 1", "Year 2", "Year 3", "Overall")),
    # fixed row order across all panels
    term  = factor(term, levels = c("Elevation", "Meadow edge",
                                    "Soil moisture", "Castilleja removal"))
  )


saveRDS(coefs, "Processed Data/Drivers Coefficients.rds")
