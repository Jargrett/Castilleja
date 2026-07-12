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

conflicted::conflicts_prefer(lme4::lmer)

#----------------------------------------------------------#
#---------------------TEMPORAL DIVERSITY-------------------#
#----------------------------------------------------------#

diversity <- readRDS("Processed Data/Plant Diversity Full.rds")
diversity %<>% filter(year != "0")

#diversity
div.lmm <- lmer(div ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(div.lmm)
qqnorm(resid(div.lmm))#Check passed
qqline(resid(div.lmm))#Check passed
plot(div.lmm)#Check passed
Anova(div.lmm)#removal:year chisq = 10.546, df = 3, p = 0.005
emmeans(div.lmm, pairwise ~ removal|year)
emmip(div.lmm, removal ~ year)

#richness
rich.lmm <- lmer(rich ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(rich.lmm)
qqnorm(resid(rich.lmm))#Check passed
qqline(resid(rich.lmm))#Check passed
plot(rich.lmm)#Check passed
Anova(rich.lmm)#removal:year chisq = 14.927, df = 2, p < 0.001
emmeans(rich.lmm, pairwise ~ removal|year)
emmip(rich.lmm, removal ~ year)

#eveness
even.lmm <- lmer(even ~ removal*litter*year + (1|block) + (1|pair), data = diversity)
summary(even.lmm)
qqnorm(resid(even.lmm))#Check passed
qqline(resid(even.lmm))#Check passed
plot(even.lmm)#Check passed
Anova(even.lmm)#n.s.  
emmeans(even.lmm, pairwise ~ removal|year)
emmip(even.lmm, removal ~ year)

#----------------------------------------------------------#
#---------------------CHANGE IN DIVERSITY------------------#
#----------------------------------------------------------#

#---------delta diversity Analysis---------#
delta_diversity <- readRDS("Processed Data/Delta Diversity.rds")

#Change in diversity from year 1 to year 3
delta.rich <- lmer(delta_rich ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.rich)
Anova(delta.rich)#removal chisq = 16.9246, df = 1, p < 0.001
emmip(delta.rich, litter~removal)
emmeans(delta.rich, pairwise ~ litter|removal)

#Change in richness from year 1 to year 3
delta.div <- lmer(delta_div ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.div)
Anova(delta.div)#removal chisq = 17.0773, df = 1, p < 0.001
emmip(delta.div, litter~removal)
emmeans(delta.div, pairwise ~ litter|removal)

#Change in evenness from year 1 to year 3
delta.even <- lmer(delta_even ~ litter*removal + (1|block) + (1|pair), data = delta_diversity)
summary(delta.even)
Anova(delta.even)#litter:removal chisq = 6.5534, df = 3, p = 0.088
emmip(delta.even, litter~removal)
emmeans(delta.even, pairwise ~ litter|removal)

#----------------------------------------------------------#
  #-------RICHNESS: CONTROLLING FOR SAMPLING EFFORT----------#
  #----------------------------------------------------------#
  # Pre-empts the critique that higher richness in Castilleja-present plots is a
  # sampling artifact of greater total plant cover. If the removal effect holds
  # with total cover in the model, the richness difference is not a cover artifact.
  #
  # NOTE: counts in the raw data represent ramets rather than genets for clonal
  # species, so individual-based rarefaction is not appropriate here; cover is
  # used as the abundance measure throughout.
  
  cover <- readRDS("Processed Data/Plant Cover.rds") %>% filter(year != "0")
  drivers <- readRDS("Processed Data/Drivers.rds")
cover_tot <- cover %>%
  group_by(year, plot) %>%
  dplyr::summarise(total_cover = sum(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year)),
         plot = as.character(plot))

# envi_div_z is built in the figures/analysis script (rich + predictors, z-scaled)
check_df <- drivers %>%
  mutate(plot = as.character(plot)) %>%
  left_join(cover_tot, by = c("year", "plot")) %>%
  mutate(total_cover_z = as.numeric(scale(total_cover)))

# richness ~ removal, WITHOUT cover
mod_nocov <- lmer(rich ~ removal + factor(year) + (1|block/pair/plot),
                  data = check_df)

# richness ~ removal, CONTROLLING for total plot cover
mod_cov <- lmer(rich ~ removal + total_cover_z + factor(year) + (1|block/pair/plot),
                data = check_df)

summary(mod_nocov)
summary(mod_cov)
anova(mod_cov)

# side-by-side comparison of the removal coefficient
rbind(
  broom.mixed::tidy(mod_nocov, effects = "fixed", conf.int = TRUE) %>%
    filter(grepl("removal", term)) %>% mutate(model = "Without cover"),
  broom.mixed::tidy(mod_cov, effects = "fixed", conf.int = TRUE) %>%
    filter(grepl("removal", term)) %>% mutate(model = "With cover")
)

