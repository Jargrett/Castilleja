#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)

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
  mutate(
    total_cover = sum(percent_cover[percent_cover > 0], na.rm = TRUE),
    rel_cover = if_else(
      percent_cover > 0 & total_cover > 0,
      percent_cover / total_cover * 100,
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
    mean_rel_cover_site = mean(mean_rel_cover_plot, na.rm = TRUE),
    se_percent_cover_site = sd(mean_rel_cover_plot, na.rm = TRUE) /
      sqrt(sum(!is.na(mean_rel_cover_plot))),
    prop_dom_site = mean(prop_dom_plot, na.rm = TRUE),
    prop_rare_site = mean(prop_rare_plot, na.rm = TRUE),
    prop_int_site = mean(prop_int_plot, na.rm = TRUE),
    Dominance = case_when(
      #Commonly occuring and most often of a specific class
      freq >= 0.50 & prop_dom_site > 0.5 ~ "DC",
      freq >= 0.50 & prop_rare_site > 0.5 ~ "RC",
      #less commonly occurring but most often of a specific class
      freq < 0.50 & prop_dom_site > 0.5 ~ "DU",
      freq < 0.50 & prop_rare_site > 0.5 ~ "RU",
      #Intermediate speices that are more often dominant 
      freq >= 0.50 & prop_int_site > 0.10 ~ "IC",
      freq < 0.50  & prop_int_site > 0.10 ~ "IR",
      TRUE ~ "I" ), 
    .groups = "drop"
  ) %>%
  filter(freq > 0)

ggplot(site_dom, aes(x = rel_cover, y = reorder(code, rel_cover), fill = Dominance)) +
  geom_col() +                       # horizontal bar
  scale_fill_brewer(palette = "Set2") +  # nice colors for dominance
  labs(
    x = "Relative Percent Cover",
    y = "Species",
    fill = "Dominance Class"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )