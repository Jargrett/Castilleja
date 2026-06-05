setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")

#packages
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(vegan)#for calculating diversity
library(labdsv)#enables restructuring for ecological analysis
library(conflicted)
#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(purrr::set_names)

#----------------------------------------------------------#
#-------------------PLOT DATA PROCESSING-------------------#
#----------------------------------------------------------#
plot_data <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot_data$plot <- as.factor(plot_data$plot)
plot_data$pair <- as.factor(plot_data$pair)
plot_data$block <- as.factor(plot_data$block)
plot_data %<>% mutate(removal = recode(removal,"C" = "Present","R" = "Removed"))

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
#--------------------PLANT BIOMASS DATA--------------------#
#----------------------------------------------------------#
#import and restructure  biomass data (raw) 
biomass <- read.csv("Raw Data/EL Biomass - Biomass.csv")
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)
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
    cas_cover_mean = rowMeans(across(c(cas_cover_yr1, cas_cover_yr2, cas_cover_yr3)),       na.rm = TRUE),
    cas_cover_delta = cas_cover_yr3   - cas_cover_yr1,
    cas_count_mean = rowMeans(across(c(cas_count_yr1, cas_count_yr2, cas_count_yr3)),       na.rm = TRUE),
    cas_count_delta = cas_count_yr3   - cas_count_yr1,
    plant_cover_mean = rowMeans(across(c(plant_cover_yr1, plant_cover_yr2, plant_cover_yr3)), na.rm = TRUE),
    plant_cover_delta = plant_cover_yr3 - plant_cover_yr1,
    envi_cover_mean = rowMeans(across(c(envi_cover_yr1, envi_cover_yr2, envi_cover_yr3)),    na.rm = TRUE),
    envi_cover_delta = envi_cover_yr3  - envi_cover_yr1
  )

cast_cov <- cover_metrics %>% filter(cas_cover_mean != 0)
ggplot(cast_cov, aes(cas_cover_delta, plant_cover_delta)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(x = "Castilleja cover Delta (%)", y = "Co-occurring plant cover delta (%)") +
  theme_bw()

m <- lm(plant_cover_delta ~ cas_cover_delta, data = cast_cov)
Anova(m)

# per-plot design info (constant within a plot)
site_plot <- site %>%
  filter(year %in% c(1, 2, 3)) %>%
  distinct(plot, field_plot, pair, block, litter, removal)

# merge: one row per plot, design + all cover metrics
plot_summary <- site_plot %>%
  left_join(cover_metrics, by = "plot")

site <- read.csv("Raw Data/Site Level Data - EL.csv") %>%
  mutate(year = case_when(Year == 2023 ~ 1,
                          Year == 2024 ~ 2,
                          Year == 2025 ~ 3)) %>% 
  mutate(year = as.numeric(year)) 


site_cover <- site %>%
  left_join(cover_sum, by = c("plot", "year"))

saveRDS(site_cover, "Processed Data/Castilleja Cover Yearly.rds")
# PRE: year 0, attached to the 2023 site metadata (since pre was collected then)
pre_cover <- cover_sum %>% filter(year == 0) %>% select(-year)

site_cover_pre <- site %>%
  filter(year == 1) %>%                  # the 2023 site rows
  left_join(pre_cover, by = "plot")


saveRDS(site_cover, "Processed Data/Castilleja Cover Yearly.rds")
#----------------------------------------------------------#
#-------------SPECIES RESPONSE TO TREATMENT----------------#
#----------------------------------------------------------#
