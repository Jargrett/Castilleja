setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
library(plyr)
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
library(car)#for regression analysis
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(codyn)#compositional analysis
library(labdsv)#enables restructuring for ecological analysis
library(ggeffects)
library(ggcharts)
library(ggthemes)
library(plyr)
library(tidytext)

cover <- readRDS("Processed Data/Cleaned Cover.rds")
plot <- readRDS("Processed Data/Plot Data.rds")

#------Treatment Affinity Calculations-------#
#gain — the difference in gain rates between treatments (Removed - Present). 
#A positive value means that species colonized a greater proportion of plots under parasite removal than in the control.
#loss — the difference in loss rates between treatments (Removed - Present). 
#A positive value means that species disappeared from a greater proportion of plots under parasite removal.
#affinity — calculated as gain - loss. This combines both dynamics into a single value:

#Positive affinity means the species colonized more and was lost less when CASE was removed
#A negative affinity would mean the species is associated with parasite presence
# Near zero means treatment had little net effect on that species
#gain_diff (0.05) = colonized 5% more plots under parasite removal
#loss_diff = -0.10 → lost from 10% fewer plots under parasite removal

cover %<>% filter (year != "0")
cover$year <- as.integer(cover$year)

#filter extra info and isolate
pa_cover <- cover %>%
  complete(plot, year, code, fill = list(cover = 0)) %>%
  mutate(present = ifelse(cover > 0, 1, 0)) %>% 
  select(year, plot, code, present) %>% 
  filter(code != "bare") %>% 
  filter(code != "litter") %>% 
  filter(code != "rock") %>% 
  filter(code != "CASE") %>% 
  filter(year != "2") %>% 
  pivot_wider(names_from = year, values_from = present) %>% 
  rename("y1" = "1", "y3" = "3")


persist <- pa_cover %>%
  mutate(change = case_when(
    y1 == 0 & y3 == 1 ~ "gain",
    y1 == 1 & y3 == 0 ~ "loss",
    y1 == 1 & y3 == 1 ~ "persist",
    TRUE ~ "absent")) %>% 
  filter(change != "absent")

meta <- plot %>% 
  select(plot, pair, block, removal)

persist %<>%
  left_join(meta, by = "plot") 

# Step 1: plots per treatment as denominator
n_plots <- meta %>%
  group_by(removal) %>%
  summarise(total_plots = n_distinct(plot))

# Step 2: rate of gain/loss per species per treatment
rate <- persist %>%
  group_by(code, removal, change) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(n_plots, by = "removal") %>%
  mutate(rate = n / 20) %>%
  select(-n, -total_plots) %>%
  pivot_wider(names_from = change, values_from = rate, values_fill = 0)

# Step 3: calculate difference between treatments
diff <- rate %>%
  pivot_longer(cols = c(gain, loss, persist), names_to = "change", values_to = "rate") %>%
  pivot_wider(names_from = removal, values_from = rate, values_fill = 0) %>%
  mutate(diff = Removed - Present) %>%  # positive = more common in Removed
  arrange(desc(abs(diff)))

# How many plots does each species appear in at all?
freq <- persist %>%
  group_by(code) %>%
  summarise(
    n_plots = n_distinct(plot),
    n_gains = sum(change == "gain"),
    n_losses = sum(change == "loss"),
    n_persist = sum(change == "persist")
  ) %>%
  arrange(n_plots)


# Only keep species present in at least 5 plots
common_species <- freq %>%
  filter(n_plots >= 5) %>%
  pull(code)

# Total gain/loss events per treatment and change type
total_changes <- persist %>%
  group_by(removal, change) %>%
  summarise(total_change = n(), .groups = "drop")

summary_df <- persist %>%
  group_by(code, removal, change) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(total_changes, by = c("removal", "change")) %>%
  mutate(proportion = n / total_change) %>%
  select(-n, -total_change) %>%
  pivot_wider(names_from = change, values_from = proportion, values_fill = 0)

species_data <- read.csv("Raw Data/EL Species List - EL.csv")

affinity <- diff %>%
  filter(change != "persist") %>%
  select(code, change, diff) %>%
  pivot_wider(names_from = change, values_from = diff, values_fill = 0) %>%
  rename(gain_diff = gain, loss_diff = loss) %>%  # rename to make clear these are already differences
  mutate(
    affinity = gain_diff - loss_diff,
    abs_affinity = abs(affinity)
  ) %>%
  left_join(species_data %>% select(code, functional_group), by = "code") %>%
  arrange(desc(abs_affinity))%>% 
  mutate_if(is.numeric, round, digits = 3) 

saveRDS(affinity, "Processed Data/Removal Affinity.rds")
common_affinity <- affinity %>% filter(code %in% common_species)

ral_abund <- cover %>%  # or whatever your main df is named
  group_by(code) %>%
  summarise(total_cover = sum(cover, na.rm = TRUE)) %>%
  mutate(relative_abundance = total_cover / sum(total_cover)) %>%
  select(code, relative_abundance)

# Merge into affinity dataframe
affinity_merged <- affinity %>%
  left_join(ral_abund, by = "code")

ggplot(affinity_merged, aes(x = relative_abundance, y = affinity)) +
  geom_point(aes(color = functional_group), size = 3, alpha = 0.7) +
  geom_text(aes(label = code), vjust = -0.8, size = 3, check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("#4B3B40", "#9B7E46", "#CBBBAA", "#A8C7BB", "#808F87")) +
  labs(x = "Species Relative Abundance",y = "Affinity",
    color = "Functional Group"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

#----------Indicator Species Analysis---------#
#Avery

ind.y3 <- filter(cover, year == 3)#filtering for a specific site
ind.y2 <- filter(cover, year == 2)#filtering for a specific site
ind.y1 <- filter(cover, year == 1)#filtering for a specific site
y3.matrix <- ind.y3 %>% select(plot, code, cover) %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")
y2.matrix <- ind.y2 %>% select(plot, code, cover) %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")
y1.matrix <- ind.y1 %>% select(plot, code, cover) %>%
  pivot_wider(names_from = code, values_from = cover, values_fill = 0) %>%
  column_to_rownames("plot")

case.removal <- ind.y2 %>%
  group_by(plot) %>%
  slice(1) %>%
  arrange(plot) %>%
  pull(removal)

case.pair <- ind.y2 %>% group_by(plot) %>%
  slice(1) %>%
  arrange(plot) %>%
  pull(pair)

y3.inv = multipatt(y2.matrix, case.removal, func = "r.g",
                           control = how(blocks = case.pair, nperm = 9999))

summary(y3.inv, alpha = 0.1)
