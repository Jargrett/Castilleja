



#----------------------------------------------------------#
#------------------ENVIRONMENTAL GRADIENTS-----------------#
#----------------------------------------------------------#
site_cover <- readRDS("Processed Data/Castilleja Cover Yearly.rds")

# one row per plot: env + Castilleja cover, soil moisture averaged across years
plot_env <- site_cover %>%
  group_by(plot, block) %>%
  summarise(elevation        = first(elevation),
            disturbance_type = first(disturbance_type),
            soil_moisture    = mean(soil_moisture, na.rm = TRUE),
            cas_cover_mean   = mean(cas_cover, na.rm = TRUE),
            removal          = first(removal),
            .groups = "drop") %>%
  mutate(block = factor(block), disturbance_type = factor(disturbance_type))

#---------do the blocks capture real heterogeneity?---------#
#soil moisture by block
moist.block <- lm(soil_moisture ~ block, data = plot_env)
Anova(moist.block)#block F = 0.85, df = 4,35, p = 0.503 

#elevation by block
elev.block <- lm(elevation ~ block, data = plot_env)
Anova(elev.block)#block F = 26.73, df = 4,35, p < 0.001

#----------------------------------------------------------#
#--------------------CASTILLEJA NICHE----------------------#
#----------------------------------------------------------#
site_cover <- readRDS("Processed Data/Castilleja Cover Yearly.rds")

# one row per plot: env + mean Castilleja cover.  removal is coded R / C here.
plot_env <- site_cover %>%
  group_by(plot, block) %>%
  summarise(elevation        = first(elevation),
            disturbance_type = first(disturbance_type),
            soil_moisture    = mean(soil_moisture, na.rm = TRUE),
            cas_cover_mean   = mean(cas_cover, na.rm = TRUE),
            removal          = first(removal),
            .groups = "drop") %>%
  mutate(block = factor(block), disturbance_type = factor(disturbance_type))

# niche from control plots only (C = Castilleja present); n = 20
niche_dat <- plot_env %>% filter(removal == "C")

cas.niche <- lm(cas_cover_mean ~ soil_moisture + elevation + disturbance_type, data = niche_dat)
summary(cas.niche)
qqnorm(resid(cas.niche))#Check
qqline(resid(cas.niche))#Check
plot(cas.niche)#Check
Anova(cas.niche)#fill in once run

#----------------------------------------------------------#
#----------------SITE-LEVEL COVER DESCRIPTION--------------#
#----------------------------------------------------------#
site_cover <- readRDS("Processed Data/Castilleja Cover Yearly.rds")
site_cover_pre <- readRDS("Processed Data/Castilleja Cover Pre.rds")

#overall and yearly means of the cover components
site_cover %>%
  group_by(year) %>%
  summarise(across(c(plant_cover, envi_cover, cas_cover, cas_count),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x,   na.rm = TRUE))), .groups = "drop")

#----------baseline equivalence (pre-treatment, R vs C; n = 40)----------#
base.cas <- lmer(case_cover ~ removal*litter*year + (1|block) + (1|pair), data = site_cover_pre)
Anova(base.cas)

base.plant <- lm(plant_cover ~ removal, data = site_cover_pre)
Anova(base.plant)

base.envi  <- lm(envi_cover  ~ removal, data = site_cover_pre)
Anova(base.envi)
#----------treatment response of total cover and bare ground----------#
site_cover$year <- as.factor(site_cover$year)

plant.lmm <- lmer(plant_cover ~ removal*litter*year + (1|block) + (1|pair), data = site_cover)
Anova(plant.lmm)
emmeans(plant.lmm, pairwise ~ litter | removal) # unpack the removal:litter interaction
emmeans(plant.lmm, pairwise ~ removal)  

envi.lmm <- lmer(envi_cover ~ removal*litter*year + (1|block) + (1|pair), data = site_cover)
Anova(envi.lmm)
emmeans(envi.lmm,  pairwise ~ removal)          # which way bare ground moved
