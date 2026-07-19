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
library(purrr)

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

#Functions
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
  # All rarity and response metrics are calculated from COVER only, following
  # the approach in Watkins et al. (2026). Metrics:
  #   SR   - spatial rarity   (baseline cover rank in parasite-present plots)
  #   TR   - temporal rarity  (years present in any parasite-present plot)
  #   RR   - response ratio   (cover change under removal)
  #   dOcc - occupancy change (proportion of plot-years present)
  #   NN   - nearest-neighbour association (frequent / infrequent neighbour)
  #
  # Rarity metrics use parasite-present plots only, so they describe the
  # baseline community independent of treatment. Trait-based predictors of
  # species response are deferred to a separate analysis using external trait


#----------Importing Data----------#
plant_data   <- readRDS("Processed Data/Plant Cover.rds")
species_info <- read.csv("Raw Data/EL Species List - EL.csv")
species_info <- as.data.frame(unclass(species_info), stringsAsFactors = TRUE)

#---- taxonomic corrections (apply to species_info so they propagate) -------
species_info <- species_info %>%
  mutate(
    genus   = recode(as.character(genus),
                     "Taraxicum"     = "Taraxacum",
                     "Chamaenuerion" = "Chamaenerion"),
    species = recode(as.character(species),
                     "Hoodii" = "hoodii"),
    family  = recode(as.character(family),
                     "Thalictraceae" = "Ranunculaceae")
  )

#---- build cover working frame (records only; absences are missing rows) ---
rac_cover <- plant_data %>%
  select(year, plot, pair, removal, litter, code, count,
         functional_group, percent_cover) %>%
  filter(year != "0") %>%
  mutate(across(c(plot, pair, removal, litter, code, functional_group), as.factor))

# numeric year kept in a SEPARATE column (never overwrite the factor 'year')
rac_cover   <- rac_cover %>% mutate(year_num = as.integer(as.character(year)))
total_years <- n_distinct(rac_cover$year)

#----------Spatial Rarity (SR)----------#
# 1 - proportional cover rank in parasite-present plots. Near 1 = rare.
spatial_rarity <- rac_cover %>%
  filter(removal == "Present", percent_cover > 0) %>%
  group_by(code) %>%
  summarize(mean_cover_control = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(SR = 1 - percent_rank(mean_cover_control))

#----------Temporal Rarity (TR)----------#
# 1 - (years present in ANY parasite-present plot / total years). Near 1 = rare.
temp_rarity <- rac_cover %>%
  filter(removal == "Present", percent_cover > 0) %>%
  group_by(code) %>%
  summarize(present = n_distinct(year), .groups = "drop") %>%
  mutate(TR = 1 - (present / total_years))

#----------Abundance class----------#
rarity <- spatial_rarity %>%
  left_join(temp_rarity, by = "code") %>%
  mutate(
    spatial  = if_else(SR >= 0.25, "Sparse", "Common"),
    temporal = if_else(TR >= 0.50, "Intermittent", "Persistent"),
    class = case_when(
      spatial == "Common" & temporal == "Persistent"   ~ "CP",
      spatial == "Common" & temporal == "Intermittent" ~ "CI",
      spatial == "Sparse" & temporal == "Persistent"   ~ "SP",
      spatial == "Sparse" & temporal == "Intermittent" ~ "SI"))

#----------Response Ratio (RR) — cover----------#
# RR = (C_removal - C_present) / (C_removal + C_present), pooled over plots/years.
# Only species present in BOTH treatments retained.
response_ratio <- rac_cover %>%
  group_by(code, removal) %>%
  summarize(mean_cover = mean(percent_cover, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = removal, values_from = mean_cover) %>%
  filter(!is.na(Removed) & !is.na(Present)) %>%
  mutate(RR = (Removed - Present) / (Removed + Present)) %>%
  rename(mean_cover_control = Present, mean_cover_removal = Removed)

#----------Merge RR with rarity + species info----------#
rarity_response <- response_ratio %>%
  left_join(rarity, by = "code") %>%
  left_join(species_info %>% select(code, functional_group, life_history, growth_form),
            by = "code") %>%
  filter(!is.na(SR) & !is.na(TR))

#----------Predictors of the cover response----------#
Anova(lm(RR ~ SR, data = rarity_response), type = "II")
Anova(lm(RR ~ TR, data = rarity_response), type = "II")
summary(lm(RR ~ SR, data = rarity_response))
summary(lm(RR ~ TR, data = rarity_response))
Anova(lm(RR ~ functional_group, data = rarity_response), type = "II")
Anova(lm(RR ~ life_history,      data = rarity_response), type = "II")

#----------Sensitivity — class cutoffs +/- 0.10----------#
make_rarity <- function(sr_cut, tr_cut) {
  spatial_rarity %>%
    left_join(temp_rarity, by = "code") %>%
    mutate(spatial  = if_else(SR >= sr_cut, "Sparse", "Common"),
           temporal = if_else(TR >= tr_cut, "Intermittent", "Persistent"))
}
sensitivity_check <- function(rarity_df, label) {
  df <- response_ratio %>%
    left_join(rarity_df %>% select(code, SR, TR), by = "code") %>%
    filter(!is.na(SR) & !is.na(TR))
  cat("\n---", label, "---\n")
  cat("SR p:", coef(summary(lm(RR ~ SR, data = df)))["SR", "Pr(>|t|)"], "\n")
  cat("TR p:", coef(summary(lm(RR ~ TR, data = df)))["TR", "Pr(>|t|)"], "\n")
}
sensitivity_check(make_rarity(0.15, 0.40), "Cutoffs - 0.10")
sensitivity_check(make_rarity(0.35, 0.60), "Cutoffs + 0.10")

#----------Sensitivity — exclude transients (present in year 1)----------#
year_one <- rac_cover %>%
  filter(removal == "Present", year_num == min(year_num), percent_cover > 0) %>%
  distinct(code)
no_transients <- rarity_response %>% filter(code %in% year_one$code)
summary(lm(RR ~ SR, data = no_transients)); Anova(lm(RR ~ SR, data = no_transients), type = "II")
summary(lm(RR ~ TR, data = no_transients)); Anova(lm(RR ~ TR, data = no_transients), type = "II")

#----------Delta Occurrence (dOcc)----------#
# Proportion of plot-years present per treatment; absences are true zeros.
delta_occ <- plant_data %>%
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

Anova(lm(delta_occ ~ SR, data = delta_rarity), type = "II")
Anova(lm(delta_occ ~ TR, data = delta_rarity), type = "II")

#----------Indicator Species Analysis----------#
# Mean cover per plot across years (avoids pseudoreplication); permutations
# constrained within paired-plot blocks.
ind_pooled_mean <- plant_data %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot, code) %>%
  summarise(cover = mean(cover, na.rm = TRUE), .groups = "drop")

pooled_matrix <- ind_pooled_mean %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")

ind_meta <- plant_data %>%
  filter(year != "0", !code %in% c("bare", "litter", "rock", "CASE")) %>%
  group_by(plot) %>% slice(1) %>% arrange(plot)

pooled_inv <- multipatt(pooled_matrix, ind_meta$removal, func = "r.g",
                        control = how(blocks = ind_meta$pair, nperm = 9999))
summary(pooled_inv, alpha = 0.1)

pooled_ind_results <- as.data.frame(pooled_inv$sign) %>%
  rownames_to_column("code") %>%
  select(code, stat, p.value) %>%
  rename(indval = stat) %>%
  mutate(significant = p.value < 0.1) %>%
  arrange(p.value)

indval_rarity <- pooled_ind_results %>%
  left_join(delta_rarity %>% select(code, SR, TR, RR, delta_occ, class), by = "code")

#----------Nearest Neighbour----------#
# Computed on ALL years (year 0 retained for the figure objects); the summary
# column for robinhood uses post-treatment years only.
nn_cover <- plant_data %>%
  select(year, plot, code, percent_cover, nearest_neighbor)   # no year filter

nearest <- nn_cover %>%
  group_by(code, year) %>%
  summarise(total_cover = sum(percent_cover, na.rm = TRUE),
            nn_count    = sum(nearest_neighbor, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(year) %>%
  mutate(rel_abund_cover = total_cover / sum(total_cover, na.rm = TRUE),
         nn_freq         = nn_count / sum(nn_count, na.rm = TRUE)) %>%
  ungroup()

nn_yearly <- nearest %>%
  group_by(year) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(nn_freq ~ rel_abund_cover, data = .x)),
         results = map2(model, data, ~ bind_cols(.y, as.data.frame(
           predict(.x, newdata = .y, interval = "prediction", level = 0.95))) %>%
             mutate(NN_association = case_when(
               nn_freq > upr ~ "Frequent Neighbor",
               nn_freq < lwr ~ "Infrequent Neighbor",
               TRUE          ~ "As expected")))) %>%
  select(year, results) %>%
  unnest(results)

# per-year figure objects (year 0 included)
split_and_name(nn_yearly, "year")
saveRDS(nn_yearly_X0, "Processed Data/NN Pre.rds")
saveRDS(nn_yearly_X1, "Processed Data/NN 23.rds")
saveRDS(nn_yearly_X2, "Processed Data/NN 24.rds")
saveRDS(nn_yearly_X3, "Processed Data/NN 25.rds")

# summary column for robinhood — POST-TREATMENT years only
nn_summary <- nn_yearly %>%
  filter(year != "0") %>%
  group_by(code) %>%
  summarise(nn_freq_years   = sum(NN_association == "Frequent Neighbor"),
            nn_infreq_years = sum(NN_association == "Infrequent Neighbor"),
            .groups = "drop") %>%
  mutate(nn_association = case_when(
    nn_freq_years > 0 & nn_infreq_years > 0 ~ "Mixed",
    nn_freq_years > 0                        ~ "Frequent Neighbor",
    nn_infreq_years > 0                      ~ "Infrequent Neighbor",
    TRUE                                     ~ "As expected"))
#----------Plot-level colonizer / extirpation classification----------#
plot_changes <- plant_data %>%
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
      change_control <= -1 & abs(change_control) >= 2 * abs(change_removal) & change_removal >= 0 ~ "Control extirpation",
      change_removal <= -1 & abs(change_removal) >= 2 * abs(change_control) & change_control >= 0 ~ "Removal extirpation",
      TRUE ~ "No clear response"))

indval_rarity <- indval_rarity %>%
  left_join(plot_changes %>% select(code, change_control, change_removal,
                                    turnover_magnitude, turnover_class),
            by = "code")

#----------Summary dataframe (robinhood)----------#
# growth_form retained as raw reference data for later merging with database
# traits; no derived groupings computed here.
robinhood <- indval_rarity %>%
  left_join(species_info %>% select(code, family, genus, species,
                                    functional_group, life_history, growth_form),
            by = "code") %>%
  left_join(nn_summary, by = "code") %>%
  select(code, family, genus, species,
         functional_group, life_history, growth_form,
         spatial_rarity = SR, temporal_rarity = TR,
         response_ratio = RR, occupancy_shift = delta_occ,
         dominance_class = class, indicator = significant,
         nn_association, nn_freq_years, nn_infreq_years,
         turnover_class, turnover_magnitude)

saveRDS(robinhood, "Processed Data/Robinhood Summary.rds")
write.csv(robinhood, "Processed Data/Robinhood Summary.csv", row.names = FALSE)

robinhood <- readRDS("Processed Data/Robinhood Summary.rds")

#----------------------------------------------------------#
#---------------------SIZE OR NUMBER-----------------------#
#----------------------------------------------------------#
# Cover conflates the number of plant units with their size. Decomposing it
# tests whether Castilleja alters community structure through establishment/
# persistence (density) or growth (size).
#
# NOTE: counts are ramets/stems, not genets (many species are clonal), so
# "density" = stem density and "size" = cover per counted unit.

library(dplyr); library(lme4); library(lmerTest); library(car); library(emmeans)

cover_all <- readRDS("Processed Data/Plant Cover.rds")

plot_meta <- site %>%
  mutate(plot = as.character(plot)) %>%
  select(plot, block, pair, litter, removal) %>%
  distinct(plot, .keep_all = TRUE)

struct_all <- cover_all %>%
  group_by(year, plot) %>%
  dplyr::summarise(total_cover = sum(percent_cover, na.rm = TRUE),
                   total_count = sum(count, na.rm = TRUE),
                   .groups = "drop") %>%
  mutate(mean_size = total_cover / total_count,
         plot = as.character(plot)) %>%
  left_join(plot_meta, by = "plot")

#---- baseline equivalence (pre-treatment, year 0) -------------------------
pre <- struct_all %>% filter(year == "0")

Anova(lmer(total_cover    ~ removal + (1|block) + (1|pair), data = pre))
Anova(glmer(total_count   ~ removal + (1|block) + (1|pair), data = pre, family = poisson))
Anova(lmer(log(mean_size) ~ removal + (1|block) + (1|pair), data = pre))

#---- treatment years ------------------------------------------------------
struct_df <- struct_all %>% filter(year != "0")

cover.lmm <- lmer(total_cover ~ removal*litter*year + (1|block) + (1|pair),
                  data = struct_df)
count.glmm <- glmer(total_count ~ removal*litter*year + (1|block) + (1|pair),
                    data = struct_df, family = poisson)
size.lmm <- lmer(log(mean_size) ~ removal*litter*year + (1|block) + (1|pair),
                 data = struct_df)

# diagnostics
qqnorm(resid(cover.lmm)); qqline(resid(cover.lmm)); plot(cover.lmm)
qqnorm(resid(size.lmm));  qqline(resid(size.lmm));  plot(size.lmm)

Anova(cover.lmm)
Anova(count.glmm)
Anova(size.lmm)

# removal contrast within each year
emmeans(cover.lmm,  pairwise ~ removal|year)
emmeans(count.glmm, pairwise ~ removal|year, type = "response")
emmeans(size.lmm,   pairwise ~ removal|year)

emmip(cover.lmm,  removal ~ year)
emmip(count.glmm, removal ~ year)
