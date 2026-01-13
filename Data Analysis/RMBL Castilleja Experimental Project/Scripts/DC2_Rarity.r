setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)
library(forcats) 

#-----Assessing species Dominance---------#
#Importing Species Data
rac.cover <- read.csv("Processed Data/Subset.csv")
species_info <- read.csv("Raw Data/EL Species List - EL.csv")

species_code <- species_info %>%
  distinct(code, functional_group)

rac.cover.full <- rac.cover %>%
  group_by(year, plot, pair, block, removal, litter) %>%
  complete(
    code = unique(species_info$code),
    fill = list(percent_cover = 0, count = 0)
  ) %>%
  ungroup() %>%
  left_join(species_code, by = "code") %>% 
  mutate(functional_group = coalesce(functional_group.y, functional_group.x)) %>%
  select(-functional_group.x, -functional_group.y)


#Calculate relative cover and dominance per plot-year
plot_dom <- rac.cover.full %>%
  group_by(plot, year) %>%
  dplyr::mutate(
    total_cover = sum(percent_cover, na.rm = TRUE),
    rel_cover = if_else(
      percent_cover > 0 & total_cover > 0,
      percent_cover / total_cover,
      0),
    q66 = quantile(rel_cover[rel_cover > 0], 0.66, na.rm = TRUE),
    q33 = quantile(rel_cover[rel_cover > 0], 0.33, na.rm = TRUE),
    dominance_class = case_when(
      percent_cover == 0 ~ NA_character_,
      rel_cover >= q66 | rel_cover >= 0.10 ~ "dominant",
      rel_cover <= q33 | rel_cover <= 0.01  ~ "rare",
      percent_cover == 0 ~ NA_character_,
      TRUE ~ "intermediate"
    )
  ) %>%
  ungroup() %>%
  group_by(plot, code) %>%
  summarise(
    mean_percent_cover_plot = mean(percent_cover, na.rm = TRUE),
    se_percent_cover_plot = sd(percent_cover, na.rm = TRUE) /
      sqrt(sum(!is.na(percent_cover))),
    mean_rel_cover_plot = mean(rel_cover, na.rm = TRUE),
    prop_dom_plot = mean(dominance_class == "dominant", na.rm = TRUE),
    prop_rare_plot = mean(dominance_class == "rare", na.rm = TRUE),
    prop_int_plot = mean(dominance_class == "intermediate", na.rm = TRUE),
    occur = mean(percent_cover > 0),
    .groups = "drop"
  )

site_dom <- plot_dom %>%
  group_by(code) %>%
  summarise(
    freq = mean(occur),
    mean_rel_cover_site = mean(mean_rel_cover_plot[mean_rel_cover_plot > 0], na.rm = TRUE),
    se_percent_cover_site = sd(mean_rel_cover_plot[mean_rel_cover_plot > 0], na.rm = TRUE) /
      sqrt(sum(!is.na(mean_rel_cover_plot[mean_rel_cover_plot > 0]))),
    prop_dom_site = mean(prop_dom_plot, na.rm = TRUE),
    prop_rare_site = mean(prop_rare_plot, na.rm = TRUE),
    prop_int_site = mean(prop_int_plot, na.rm = TRUE),
    rarity = case_when( #The proportion of times that a given species occurs as one of the most abundant species in a give plot-year
      prop_dom_site > 0.5 ~ "Dominant",
      prop_rare_site > 0.5 ~ "Rare",
      prop_int_site > 0.10 ~ "Intermediate",
      TRUE ~ "Intermediate" ),
    occurance = case_when(
      freq > 0.66 ~ "Frequent",
      freq <=0.66 & freq >= 0.33 ~ "Intermediate",
      freq < 0.33 ~ "Infrequent",
      TRUE ~ "Intermediate" ),
    .groups = "drop"
  ) %>%
  filter(freq > 0)


species.rarity <- full_join(species_info, site_dom, by = "code") %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>% 
  drop_na()

species.rarity <- as.data.frame(unclass(species.rarity),stringsAsFactors=TRUE)
write.csv(species.rarity, "Processed Data/Species Rarity.csv", row.names=FALSE)

library(webr)
species.pie <- species.rarity %>%
  group_by(rarity, occurance) %>%
  mutate(count = n()) %>% 
  summarise(count = n(), .groups = 'drop')

PieDonut(species.pie, aes(occurance, rarity, count=count, fill = occurance )) +
  scale_fill_manual(values=c("#582f0e","#7f4f24", "#a68a64","#b6ad90","#c2c5aa", "#656d4a","#414833", "#333d29" ))

