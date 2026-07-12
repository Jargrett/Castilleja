setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")

#packages
library(plyr)
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(vegan)#for calculating diversity
library(labdsv)#enables restructuring for ecological analysis
library(conflicted)
library(pwr)
library(MuMIn)
library(car)
library(indicspecies)
library(codyn)
library(broom)
library(broom.mixed)
library(forcats)

#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(purrr::set_names)
conflict_prefer("select","dplyr")
conflict_prefer("count","dplyr")

#----------------------------------------------------------#
#-------------------PLOT DATA PROCESSING-------------------#
#----------------------------------------------------------#
plot_data <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot_data$plot <- as.factor(plot_data$plot)
plot_data$pair <- as.factor(plot_data$pair)
plot_data$block <- as.factor(plot_data$block)
plot_data %<>% mutate(removal = recode(removal,"C" = "Present","R" = "Removed"))
saveRDS(plot_data, "Processed Data/Plot Data.rds")

#----------------------------------------------------------#
#---------------PLANT COVER DATA PROCESSING----------------#
#----------------------------------------------------------#
#importing scripts
cover.pre <- read.csv("Raw Data/Emerald Lake Plant Data - pre.csv")
cover.23 <- read.csv("Raw Data/Emerald Lake Plant Data - 2023.csv")
cover.24 <- read.csv("Raw Data/Emerald Lake Plant Data - 2024.csv")
cover.25 <- read.csv("Raw Data/Emerald Lake Plant Data - 2025.csv")

#combine datasets
cover.comb <- rbind.fill(cover.pre,cover.23,cover.24,cover.25)
cover.comb <- as.data.frame(unclass(cover.comb),stringsAsFactors=TRUE)
cover.comb$plot <- as.factor(cover.comb$plot)
cover.comb$pair <- as.factor(cover.comb$pair)
cover.comb$block <- as.factor(cover.comb$block)
cover.comb %<>% mutate(removal = recode(removal,"C" = "Control","R" = "Removal"))
cover.comb <- cover.comb %>% distinct()
cover.comb %<>% 
  mutate(year = recode(year,
                       "Pre" = "0",
                       "2023" = "1",
                       "2024" = "2",
                       "2025" = "3",)) %>% 
  mutate(removal = recode(removal,
                          Control = "Present",
                          Removal = "Removed"))
saveRDS(cover.comb, "Processed Data/Dirty Cover.rds")

#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
cover.comb.clean$plot <- as.factor(cover.comb.clean$plot)
cover.comb.clean$pair <- as.factor(cover.comb.clean$pair)
cover.comb.clean$block <- as.factor(cover.comb.clean$block)
saveRDS(cover.comb.clean, "Processed Data/Plant Cover.rds")

comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))
saveRDS(comb.cov, "Processed Data/Cleaned Cover.rds")

#----------------------------------------------------------#
#--------------DIVERSITY ANALYSIS PROCESSING---------------#
#----------------------------------------------------------#

#filter for year/pre and calculate
emerald.pre <- comb.cov %>% 
  filter(year == "0") %>%
  select(-c(year))
emerald.23 <- comb.cov %>% 
  filter (year == "1") %>%
  select(-c(year))
emerald.24 <- comb.cov %>% 
  filter (year == "2") %>%
  select(-c(year))
emerald.25 <- comb.cov %>% 
  filter (year == "3") %>%
  select(-c(year))

#convert to matrix format for diversity calculations
emerald.pre.matrix <- matrify(emerald.pre)
emerald.23.matrix <- matrify(emerald.23)
emerald.24.matrix <- matrify(emerald.24)
emerald.25.matrix <- matrify(emerald.25)

#---------------Diversity Calculations---------------#
# Calculating Shannon diversity for plots
div.pre <- diversity(emerald.pre.matrix, index = "shannon")
div.23 <- diversity(emerald.23.matrix, index = "shannon")
div.24 <- diversity(emerald.24.matrix, index = "shannon")
div.25 <- diversity(emerald.25.matrix, index = "shannon")
# Calculating species richness for plots
rich.pre <- specnumber(emerald.pre.matrix)
rich.23 <- specnumber(emerald.23.matrix)
rich.24 <- specnumber(emerald.24.matrix)
rich.25 <- specnumber(emerald.25.matrix)
# Calculating species evenness for plots 
even.pre <- diversity(emerald.pre.matrix, index = "shannon") / log(specnumber(emerald.pre.matrix))
even.23 <- diversity(emerald.23.matrix, index = "shannon") / log(specnumber(emerald.23.matrix))
even.24 <- diversity(emerald.24.matrix, index = "shannon") / log(specnumber(emerald.24.matrix))
even.25 <- diversity(emerald.25.matrix, index = "shannon") / log(specnumber(emerald.25.matrix))

#---------------Combining results and exporting---------------#
el.pre <- cbind(plot_data,div.pre,rich.pre,even.pre)
el.pre %<>% 
  mutate(year = '0') %>% 
  relocate(year) %>% 
  rename(div = div.pre, rich = rich.pre, even = even.pre)
el.23 <- cbind(plot_data,div.23,rich.23,even.23)
el.23 %<>% 
  plyr::mutate(year = '1') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.23, rich = rich.23, even = even.23)
el.24 <- cbind(plot_data,div.24,rich.24,even.24)
el.24 %<>% 
  plyr::mutate(year = '2') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.24, rich = rich.24, even = even.24)
el.25 <- cbind(plot_data,div.25,rich.25,even.25)
el.25 %<>% 
  plyr::mutate(year = '3') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.25, rich = rich.25, even = even.25)

diversity <- rbind.fill(el.pre,el.23,el.24,el.25)

saveRDS(diversity, "Processed Data/Plant Diversity Full.rds")

#---------------Delta Diversity---------------#
diversity <- readRDS("Processed Data/Plant Diversity Full.rds")

delta_diversity <- diversity %>%
  filter(year != "0") %>% 
  pivot_wider(names_from  = year, 
              values_from = c(div, rich, even),
              names_sep   = "_yr") %>% 
  mutate(
    delta_div  = div_yr3  - div_yr1,
    delta_div_32  = div_yr3  - div_yr2,
    div_change_21  = div_yr2  - div_yr1,
    
    delta_rich = rich_yr3 - rich_yr1,
    delta_rich_32 = rich_yr3 - rich_yr2,
    delta_rich_21 = rich_yr2 - rich_yr1,
    
    delta_even = even_yr3 - even_yr1,
    delta_even_32 = even_yr3 - even_yr2,
    delta_even_21 = even_yr2 - even_yr1
  )

saveRDS(delta_diversity, "Processed Data/Delta Diversity.rds")
#----------------------------------------------------------#
#--------------ENVIRONMENTAL DATA PROCESSING---------------#
#----------------------------------------------------------#
cover.comb <- readRDS("Processed Data/Dirty Cover.rds")
envi.cover <- cover.comb %>% 
  filter(functional_group == 'environmental') %>% 
  subset(select = c('year','plot','pair','block','removal','litter','code','percent_cover'))

envi.cover$pair <- as.factor(envi.cover$pair)
envi.cover$plot <- as.factor(envi.cover$plot)
envi.cover$block <- as.factor(envi.cover$block)

split_and_name <- function(df, column) {
  # Get the dataframe name as a string
  df_name <- deparse(substitute(df))
  # Ensure column exists
  if (!column %in% colnames(df)) {
    stop("Column not found in dataframe")
  }
  # Split the dataframe by the column
  df_list <- split(df, df[[column]])
  # Remove dataframes with 0 rows
  df_list <- df_list[sapply(df_list, nrow) > 0]
  # If nothing remains, stop
  if (length(df_list) == 0) {
    warning("All split dataframes have 0 rows; nothing to assign.")
    return(invisible(NULL))
  }
  # Make valid R names for safety
  names(df_list) <- paste0(df_name, "_", make.names(names(df_list)))
  # Assign each dataframe into the global environment
  list2env(df_list, envir = .GlobalEnv)
  # Return the list invisibly
  invisible(df_list)
}
split_and_name(envi.cover, "code")

total.envi <- envi.cover %>%
  group_by(year,plot,pair,block,removal,litter) %>%
  summarise(
    total_cover = sum(percent_cover, na.rm = TRUE)
  )

saveRDS(total.envi, "Processed Data/Environmental Cover.rds")

#----------------------------------------------------------#
#-------------------DRIVERS OF DIVERSITY-------------------#
#----------------------------------------------------------#

site <- read.csv("Raw Data/Site Level Data - EL.csv")

site %<>% rename(year_raw = Year) %>%
  mutate(year = recode(year_raw, `2023` = 1, `2024` = 2, `2025` = 3),
         removal = recode(removal, "C" = "Present", "R" = "Removed"),
         removal = factor(removal, levels = c("Present", "Removed")))

diverse <- readRDS("Processed Data/Plant Diversity Full.rds")

envi_div <- diverse %>%
  filter(year != "0") %>%
  mutate(year = as.numeric(as.character(year)),
         plot = as.character(plot)) %>%
  select(year, plot, rich, div, even) %>%
  left_join(site %>%
              mutate(plot = as.character(plot)) %>%
              select(year, plot, block, pair, removal, litter, elevation,
                     disturbance_distance, disturbance_type, soil_moisture),
            by = c("year", "plot")) %>%
  na.omit()

envi_div_z <- envi_div %>%
  mutate(across(c(elevation, soil_moisture, disturbance_distance),
                ~ as.numeric(scale(.))))

saveRDS(envi_div_z, "Processed Data/Drivers.rds")
#----------------------------------------------------------#
#--------------------PLANT BIOMASS DATA--------------------#
#----------------------------------------------------------#
#import and restructure  biomass data (raw) 
biomass <- read.csv("Raw Data/EL Biomass - Biomass.csv")
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)
biomass %<>% mutate(removal = recode(removal,"C" = "Present","R" = "Removed"))
biomass$pair <- as.factor(biomass$pair)
biomass$plot <- as.factor(biomass$plot)
biomass$block <- as.factor(biomass$block)
str(biomass)
saveRDS(biomass, "Processed Data/Plant Biomass.rds")

#----------------------------------------------------------#
#------------------SPECIES TURNOVER DATA-------------------#
#----------------------------------------------------------#
clean_cover <- readRDS("Processed Data/Plant Cover.rds")
clean_cover %<>% filter (year != "0")
clean_cover$year <- as.integer(clean_cover$year)

#Differences in community composition over the treatment duration
comp_change <- rate_change_interval(clean_cover,
                                    time.var = "year",
                                    species.var = "code",
                                    abundance.var = "percent_cover",
                                    replicate.var = "plot")
comp_change <- as.data.frame(unclass(comp_change),stringsAsFactors=TRUE)

total_comp_change <- comp_change  %>%  
  filter (interval != "1") %>% 
  left_join(plot_data, comp_change, by = "plot")
saveRDS(total_comp_change, "Processed Data/Compositional Change.rds")


#Compare whether the rank abundance structure differs between removal at end of study
between_cover <- clean_cover %>% filter(year %in% c("1", "3"))   # adjust to your year codes

total_turn <- turnover(between_cover,
                       time.var = "year", species.var = "code",
                       abundance.var = "percent_cover", replicate.var = "plot")

app_turn <- turnover(between_cover,
                     time.var = "year", species.var = "code",
                     abundance.var = "percent_cover", replicate.var = "plot",
                     metric = "appearance")

diss_turn <- turnover(between_cover,
                      time.var = "year", species.var = "code",
                      abundance.var = "percent_cover", replicate.var = "plot",
                      metric = "disappearance")

total_turn <- plot_data %>%
  left_join(total_turn %>% select(plot, total), by = "plot") %>%
  left_join(app_turn %>% select(plot, appearance), by = "plot") %>%
  left_join(diss_turn %>% select(plot, disappearance), by = "plot")

saveRDS(total_turn, "Processed Data/Species Turnover.rds")

#----------------------------------------------------------#
#-----------------COMMUNITY COMPOSITION--------------------#
#----------------------------------------------------------#
cover <- readRDS("Processed Data/Plant Cover.rds")
cover %<>% filter (year != "0")
## 1. One community matrix, rows = plot-year ----------------------------
com_cover <- cover %>%
  select(year, plot, code, percent_cover) %>%
  mutate(sample = paste(plot, year, sep = "_")) %>%   # unique plot-year ID
  select(sample, code, percent_cover) %>%
  as.data.frame()

comm.matrix <- matrify(com_cover)

meta_raw <- read.csv("Raw Data/Site Level Data - EL.csv")

plot_meta <- meta_raw %>%
  mutate(year = recode(Year, "2023" = "1",
                         "2024" = "2", "2025" = "3"),
         sample = paste(plot, year, sep = "_")) %>%
  select(sample, year, plot, pair, block, litter, removal)    # match your columns

#align metadata to matrix row order
plot_meta <- plot_meta[match(rownames(comm.matrix), plot_meta$sample), ]

plot_meta %<>% mutate(across(c(year, pair, block, litter, removal), factor))

#Bray-Curtis on the full matrix 
dist.all <- vegdist(comm.matrix, method = "bray")
saveRDS(dist.all, "Processed Data/dist_all_bray.rds")

#One NMDS across all plot-years (ordination figure)
set.seed(20)
nmds.all <- metaMDS(dist.all, distance = "bray", k = 3, try = 1000)
nmds.all #k = 0.17

nmds.scores <- as.data.frame(vegan::scores(nmds.all, display = "sites"))
NMDS <- cbind(plot_meta, nmds.scores)

saveRDS(NMDS, "Processed Data/Community NMDS.rds")

#----------------------------------------------------------#
#--------------------NEAREST NEIGHBOR----------------------#
#----------------------------------------------------------#
plant <- readRDS("Processed Data/Plant Cover.rds")
cover <- subset(plant, select = c('year','plot','code','percent_cover', 'nearest_neighbor'))

#Sum for each species the number of times they appear as a NN
nearest <- cover %>%
  group_by(code, year) %>%
  summarise(
    total_cover = sum(percent_cover, na.rm = TRUE),
    nn_count    = sum(nearest_neighbor, na.rm = TRUE),
    .groups = "drop") %>%
  group_by(year) %>%
  mutate(rel_abund_cover = total_cover / sum(total_cover, na.rm = TRUE),
    nn_freq = nn_count/sum(nn_count, na.rm = TRUE)) %>%
  ungroup()

nn_yearly <- nearest %>%
  group_by(year) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(nn_freq ~ rel_abund_cover, data = .x)),
    results = map2(model, data, ~ {bind_cols(.y, as.data.frame(
          predict(.x, newdata = .y, interval = "prediction", level = 0.95)) 
      ) %>%
        mutate(NN_association = case_when(
            nn_freq > upr ~ "Frequent Neighbor",
            nn_freq < lwr ~ "Infrequent Neighbor",
            TRUE ~ "As expected"))
      })) %>%
  select(year, results) %>%
  unnest(results)

split_and_name(nn_yearly, "year")

saveRDS(nn_yearly_X0, "Processed Data/NN Pre.rds")
saveRDS(nn_yearly_X1, "Processed Data/NN 23.rds")
saveRDS(nn_yearly_X2, "Processed Data/NN 24.rds")
saveRDS(nn_yearly_X3, "Processed Data/NN 25.rds")

#----------------------------------------------------------#
#----------------CASTILLEJA SUMMARY DATA-------------------#
#----------------------------------------------------------#
dirty_cover <- readRDS("Processed Data/Dirty Cover.rds")

# per plot-year summaries
cover_sum <- dirty_cover %>%
  group_by(year, plot) %>%
  summarise(
    cas_cover = sum(cover[code == "CASE"], na.rm = TRUE),
    cas_count = sum(count[code == "CASE"], na.rm = TRUE),
    envi_cover = sum(cover[functional_group == "environmental"], na.rm = TRUE),
    plant_cover = sum(cover[functional_group != "environmental" & code != "CASE"], na.rm = TRUE),
    .groups = "drop")

cover_sum <- cover_sum %>%
  mutate(year = as.numeric(year),
         plot = as.integer(as.character(plot)))


# yearly + mean + delta (yr3 - yr1) for every cover metric, one row per plot
cover_metrics <- cover_sum %>%
  filter(year %in% c(1, 2, 3)) %>%
  select(plot, year, cas_cover, cas_count, plant_cover, envi_cover) %>%
  pivot_wider(names_from  = year,
              values_from = c(cas_cover, cas_count, plant_cover, envi_cover),
              names_sep   = "_yr") %>%
  mutate(
    cas_cover_mean = rowMeans(across(c(cas_cover_yr1, cas_cover_yr2, cas_cover_yr3)), na.rm = TRUE),
    cas_cover_delta = cas_cover_yr3   - cas_cover_yr1,
    cas_count_mean = rowMeans(across(c(cas_count_yr1, cas_count_yr2, cas_count_yr3)), na.rm = TRUE),
    cas_count_delta = cas_count_yr3   - cas_count_yr1,
    plant_cover_mean = rowMeans(across(c(plant_cover_yr1, plant_cover_yr2, plant_cover_yr3)), na.rm = TRUE),
    plant_cover_delta = plant_cover_yr3 - plant_cover_yr1,
    envi_cover_mean = rowMeans(across(c(envi_cover_yr1, envi_cover_yr2, envi_cover_yr3)), na.rm = TRUE),
    envi_cover_delta = envi_cover_yr3  - envi_cover_yr1
  )

site <- read.csv("Raw Data/Site Level Data - EL.csv") %>%
  mutate(year = case_when(Year == 2023 ~ 1,
                          Year == 2024 ~ 2,
                          Year == 2025 ~ 3)) %>% 
  mutate(year = as.numeric(year)) 

# per-plot design info (constant within a plot)
site_plot <- site %>%
  filter(year %in% c(1, 2, 3)) %>%
  distinct(plot, field_plot, pair, block, litter, removal)

site_cover <- site %>%
  left_join(cover_sum, by = c("plot", "year"))

saveRDS(site_cover, "Processed Data/Castilleja Cover Yearly.rds")

# PRE: year 0, attached to the 2023 site metadata (since pre was collected then)
pre_cover <- cover_sum %>% filter(year == 0) %>% select(-year)

site_cover_pre <- site %>%
  filter(year == 1) %>%                                  # 2023 metadata rows = 40 plots
  distinct(plot, pair, block, removal, litter,
           elevation, soil_moisture, disturbance_type) %>%
  left_join(cover_sum %>% filter(year == 0) %>% select(-year), by = "plot") %>%
  mutate(removal = factor(removal))

saveRDS(site_cover_pre, "Processed Data/Castilleja Cover Pre.rds")

#----------------------------------------------------------#
#-------------SPECIES RESPONSE TO TREATMENT----------------#
#----------------------------------------------------------#
#----------Importing Data----------#
plant_data <- readRDS("Processed Data/Plant Cover.rds")
rac_cover <- subset(plant_data, select = c('year','plot','pair', 'removal',
                                           'litter','code','count',
                                           'functional_group', 'percent_cover'))
rac_cover %<>% filter (year != "0")

rac_cover <- as.data.frame(unclass(rac_cover),stringsAsFactors=TRUE)
species_info <- read.csv("Raw Data/EL Species List - EL.csv")
species_info <- as.data.frame(unclass(species_info),stringsAsFactors=TRUE)
cover <- readRDS("Processed Data/Plant Cover.rds")
plot_meta <- readRDS("Processed Data/Plot Data.rds")

#----------Spatial Rarity (SR)----------#
# Proportional rank of mean cover across all control plots and years
# SR near 1 = rare, SR near 0 = dominant
spatial_rarity <- rac_cover %>%
  filter(removal == "Present") %>%
  group_by(code) %>%
  summarize(mean_cover_control = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(SR = 1 - percent_rank(mean_cover_control))

#----------Temporal Rarity (TR)----------#
# Proportion of years absent from control plots
# TR near 1 = temporally rare, TR near 0 = temporally persistent

total_years <- rac_cover %>%
  summarize(n = n_distinct(year)) %>%
  pull(n)

temp_rarity <- rac_cover %>%
  filter(removal == "Present") %>%
  group_by(code) %>%
  summarize(present = n_distinct(year), .groups = "drop") %>%
  mutate(TR = 1 - (present / total_years))

#----------Abundance group----------#
# Spatial:  SR >= 0.25 = Sparse,       SR < 0.25 = Common
# Temporal: TR >= 0.50 = Intermittent, TR < 0.50 = Persistent

rarity <- spatial_rarity %>%
  left_join(temp_rarity, by = "code") %>%
  mutate(
    spatial  = if_else(SR >= 0.25, "Sparse", "Common"),
    temporal = if_else(TR >= 0.50, "Intermittent", "Persistent"),
    class = case_when(
      spatial == "Common" & temporal == "Persistent" ~ "CP",
      spatial == "Common" & temporal == "Intermittent" ~ "CI",
      spatial == "Sparse" & temporal == "Persistent" ~ "SP",
      spatial == "Sparse" & temporal == "Intermittent" ~ "SI"
    )
  )

#----------Response Ratio (RR)----------#
# RR,i = (C_removal - C_control) / (C_removal + C_control)
# Pooled across all years and plots — one value per species
# Species absent from control plots excluded

response_ratio <- rac_cover %>%
  group_by(code, removal) %>%
  summarize(mean_cover = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = removal, values_from = mean_cover) %>%
  filter(!is.na(Removed) & !is.na(Present)) %>%
  mutate(RR = (Removed - Present) / (Removed + Present))

#----------Merge RR with Rarity and Species Info----------#
rarity_response <- response_ratio %>%
  rename(mean_cover_control = Present,
         mean_cover_removal = Removed) %>%
  left_join(rarity, by = "code") %>%
  left_join(species_info %>% select(code, functional_group, life_history, growth_form), by = "code") %>%
  filter(!is.na(SR) & !is.na(TR))

#----------Sensitivity Check — type Cutoffs ± 0.10----------#
make_rarity <- function(sr_cut, tr_cut) {
  spatial_rarity %>%
    left_join(temp_rarity, by = "code") %>%
    mutate(
      spatial  = if_else(SR >= sr_cut, "Sparse",       "Common"),
      temporal = if_else(TR >= tr_cut, "Intermittent", "Persistent"),
      class = case_when(
        spatial == "Common" & temporal == "Persistent"   ~ "CP",
        spatial == "Common" & temporal == "Intermittent" ~ "CI",
        spatial == "Sparse" & temporal == "Persistent"   ~ "SP",
        spatial == "Sparse" & temporal == "Intermittent" ~ "SI"
      ))}
rarity_low  <- make_rarity(0.15, 0.40)
rarity_high <- make_rarity(0.35, 0.60)

sensitivity_check <- function(rarity_df, label) {
  df   <- response_ratio %>%
    left_join(rarity_df %>% select(code, SR, TR, class), by = "code") %>%
    filter(!is.na(SR) & !is.na(TR))
  m.SR <- lm(RR ~ SR, data = df)
  m.TR <- lm(RR ~ TR, data = df)
  cat("\n---", label, "---\n")
  cat("SR p:", coef(summary(m.SR))["SR", "Pr(>|t|)"], "\n")
  cat("TR p:", coef(summary(m.TR))["TR", "Pr(>|t|)"], "\n")
}
sensitivity_check(rarity_low,  "Cutoffs - 0.10")
sensitivity_check(rarity_high, "Cutoffs + 0.10")

#----------Sensitivity Check — Exclude Transients----------#
# Transient = not present in year 1 of control plots
rac_cover$year <- as.integer(rac_cover$year)
year_one <- rac_cover %>% filter(removal == "Present", year == min(year)) %>% distinct(code)

no_transients <- rarity_response %>% filter(code %in% year_one$code)

model.SR.nt <- lm(RR ~ SR, data = no_transients)
model.TR.nt <- lm(RR ~ TR, data = no_transients)
model.nt    <- lm(RR ~ SR * TR, data = no_transients)

summary(model.SR.nt); Anova(model.SR.nt, type = "II"); confint(model.SR.nt)
summary(model.TR.nt); Anova(model.TR.nt, type = "II"); confint(model.TR.nt)
summary(model.nt)

pwr.f2.test(u = 1, v = 32, f2 = 0.096 / (1 - 0.096), sig.level = 0.05)# current power 0.45
pwr.f2.test(u = 1, f2 = 0.096 / (1 - 0.096), sig.level = 0.05, power = 0.80) #calulating number of species needed

#----------Delta Occurrence (ΔOcc)----------#
# Proportion of plot-year combinations present per treatment
# ΔOcc = Removed - Present (positive = more common in removal)

delta_occ <- cover %>%
  filter(year != "0") %>%
  complete(plot, year, code, fill = list(cover = 0, removal = NA)) %>%
  group_by(plot) %>%
  mutate(removal = first(na.omit(removal))) %>%
  ungroup() %>%
  filter(!code %in% c("bare", "litter", "rock", "CASE")) %>%
  mutate(present = ifelse(cover > 0, 1, 0)) %>%
  group_by(code, removal) %>%
  summarise(occ_rate = sum(present) / n(), .groups = "drop") %>%
  filter(!is.na(removal)) %>%
  pivot_wider(names_from = removal, values_from = occ_rate, values_fill = 0) %>%
  mutate(delta_occ = Removed - Present)

delta_rarity <- delta_occ %>%
  left_join(rarity_response %>% select(code, SR, TR, RR, class), by = "code") %>%
  filter(!is.na(SR))

#----------Indicator Species Analysis----------#
# Mean cover per plot averaged across years to avoid pseudoreplication
# Permutations constrained within paired plot blocks

ind_pooled_mean <- cover %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot, code) %>%
  summarise(cover = mean(cover, na.rm = TRUE), .groups = "drop")

pooled_matrix <- ind_pooled_mean %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")

ind_meta <- cover %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot) %>%
  slice(1) %>%
  arrange(plot)

case.removal_pooled <- ind_meta %>% pull(removal)
case.pair_pooled    <- ind_meta %>% pull(pair)

pooled_inv <- multipatt(pooled_matrix, case.removal_pooled, func = "r.g",
                        control = how(blocks = case.pair_pooled, nperm = 9999))

summary(pooled_inv, alpha = 0.1)

pooled_ind_results <- as.data.frame(pooled_inv$sign) %>%
  rownames_to_column("code") %>%
  select(code, stat, p.value) %>%
  rename(indval = stat) %>%
  mutate(significant = p.value < 0.1) %>%
  arrange(p.value)

sig_species <- pooled_ind_results %>% filter(significant) %>% pull(code) %>% as.character()

indval_rarity <- pooled_ind_results %>%
  left_join(delta_rarity %>% select(code, SR, TR, RR, delta_occ, class), by = "code")

#----------Plot-Level Response Classification----------#
# Identifies species that gained or lost plots in one treatment relative to the other
# Combines colonizer and extirpation analyses into a single asymmetry-based classification
# Threshold: response in one treatment must be at least 1 plot AND at least 2x the response in the other

plot_changes <- cover %>%
  filter(year %in% c(1, 3), !code %in% c("bare", "litter", "rock", "CASE")) %>%
  mutate(present = ifelse(cover > 0, 1, 0)) %>%
  group_by(code, removal, year) %>%
  summarise(n_plots = sum(present), .groups = "drop") %>%
  pivot_wider(names_from = c(removal, year), values_from = n_plots, values_fill = 0) %>%
  rename(control_y1 = Present_1, control_y3 = Present_3,
         removal_y1 = Removed_1, removal_y3 = Removed_3) %>%
  mutate(
    change_control = control_y3 - control_y1,
    change_removal = removal_y3 - removal_y1,
    turnover_magnitude = change_control - change_removal,
    turnover_class = case_when(
      change_control >=  1 & change_control >=  2 * abs(change_removal) & change_removal <= 0 ~ "Control colonizer",
      change_removal >=  1 & change_removal >=  2 * abs(change_control) & change_control <= 0 ~ "Removal colonizer",
      change_control <= -1 & abs(change_control) >=  2 * abs(change_removal) & change_removal >= 0 ~ "Control extirpation",
      change_removal <= -1 & abs(change_removal) >=  2 * abs(change_control) & change_control >= 0 ~ "Removal extirpation",
      TRUE                                                                                          ~ "No clear response"
    )
  )
# Tag species in indval_rarity with their response type and plot-level changes
indval_rarity <- indval_rarity %>%
  left_join(
    plot_changes %>% select(code, change_control, change_removal,
                            turnover_magnitude, turnover_class),
    by = "code"
  )

#---------Summary Dataframe---------#
robinhood <- indval_rarity %>%
  left_join(
    species_info %>% select(code, family, genus, species,
                            functional_group, life_history, growth_form),
    by = "code"
  ) %>%
  select(code, family, genus, species, functional_group, life_history, growth_form,
         spatial_rarity = SR, temporal_rarity = TR, response_ratio = RR, dominance_class = class, occupancy_shift = delta_occ,
         indicator = significant,
         turnover_class, turnover_magnitude)
saveRDS(robinhood, "Processed Data/Robinhood Summary.rds")
write.csv(robinhood, "Processed Data/Robinhood Summary.csv", row.names = FALSE)

# 1 — Rarity and Response Ratio
robin_rr <- robinhood %>%
  select(code, spatial_rarity, temporal_rarity, response_ratio, dominance_class) %>%
  filter(!is.na(spatial_rarity) & !is.na(temporal_rarity) & !is.na(response_ratio)) %>%
  arrange(spatial_rarity)

saveRDS(robin_rr, "Processed Data/Species Response Ratio.rds")

# 2 — Delta Occurrence
robin_delta <- delta_occ %>%
  filter(!is.na(delta_occ)) %>%
  select(code, occ_present = Present, occ_removed = Removed, delta_occ) %>%
  arrange(delta_occ)

saveRDS(robin_delta, "Processed Data/Species Change in Occurance.rds")

# 3 — Indicator Species
robin_indicator <- indval_rarity %>%
  select(code, indval, p.value) %>%
  arrange(p.value)

saveRDS(robin_indicator,  "Processed Data/Indicator Species Mean.rds")

# 4 — Turnover (colonizers and extirpated species combined)
robin_turnover <- indval_rarity %>%
  filter(turnover_class != "No clear response") %>%
  select(code, change_control, change_removal, turnover_magnitude, turnover_class) %>%
  arrange(turnover_class, desc(abs(turnover_magnitude)))

saveRDS(robin_turnover, "Processed Data/Robinhood Turnover.rds")




#----------------------------------------------------------#
#-----------------------SPECIES TRAITS---------------------#
#----------------------------------------------------------#
species <- readRDS("Processed Data/Robinhood Summary.rds")

