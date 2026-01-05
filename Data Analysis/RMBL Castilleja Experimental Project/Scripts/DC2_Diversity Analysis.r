setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)

#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::filter)

#important cover data (raw)
cover.pre <- read.csv("Raw Data/Emerald Lake Plant Data - pre.csv")
cover.23 <- read.csv("Raw Data/Emerald Lake Plant Data - 2023.csv")
cover.24 <- read.csv("Raw Data/Emerald Lake Plant Data - 2024.csv")
cover.25 <- read.csv("Raw Data/Emerald Lake Plant Data - 2025.csv")

#combine datasets
cover.comb <- rbind.fill(cover.pre,cover.23,cover.24,cover.25)
cover.comb <- as.data.frame(unclass(cover.comb),stringsAsFactors=TRUE)
cover.comb %<>%
  mutate(removal = recode(removal,
                          "C" = "Control",
                          "R" = "Removal"))
#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))
#filter for year/pre and calculate
emerald.pre <- comb.cov %>% 
  filter(year == "Pre") %>%
  select(-c(year))
emerald.23 <- comb.cov %>% 
  filter (year == "2023") %>%
  select(-c(year))
emerald.24 <- comb.cov %>% 
  filter (year == "2024") %>%
  select(-c(year))
emerald.25 <- comb.cov %>% 
  filter (year == "2025") %>%
  select(-c(year))

#convert to matrix format for diversity calculations
library(labdsv)#enables restructuring for ecological analysis
emerald.pre.matrix <- matrify(emerald.pre)
emerald.23.matrix <- matrify(emerald.23)
emerald.24.matrix <- matrify(emerald.24)
emerald.25.matrix <- matrify(emerald.25)

#---------------Diversity Calculations---------------#
library(vegan)#for calculating diversity
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
plot <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot %<>%
  mutate(removal = recode(removal,
                          "C" = "Control",
                          "R" = "Removal"))

el.pre <- cbind(plot,div.pre,rich.pre,even.pre)
el.pre %<>% 
  plyr::mutate(year = 'pre') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.pre, rich = rich.pre, even = even.pre)
el.23 <- cbind(plot,div.23,rich.23,even.23)
el.23 %<>% 
  plyr::mutate(year = '2023') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.23, rich = rich.23, even = even.23)
el.24 <- cbind(plot,div.24,rich.24,even.24)
el.24 %<>% 
  plyr::mutate(year = '2024') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.24, rich = rich.24, even = even.24)
el.25 <- cbind(plot,div.25,rich.25,even.25)
el.25 %<>% 
  plyr::mutate(year = '2025') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.25, rich = rich.25, even = even.25)

diversity.pre <- rbind.fill(el.pre,el.23,el.24,el.25)
diversity <- rbind.fill(el.23,el.24,el.25)
write.csv(diversity, "Processed Data/Plant Diversity.csv", row.names=FALSE)

#----------Rank Abundance Curve Analysis----------#
library(codyn)
rac.cover <- cover.comb.clean %>% 
  filter (year != "Pre") %>% 
  subset(select = c('year','plot','pair','block','removal','litter','functional_group', 'code','count', 'percent_cover'))
rac.cover$year <- as.integer(rac.cover$year)
write.csv(rac.cover, "Processed Data/Subset.csv", row.names=FALSE)

rac.diff <- RAC_difference(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  treatment.var = "removal",
  block.var = "pair")

rac.change <- RAC_change(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  reference.time = "1")

rac.change %<>% 
  filter (year2 != "2")

merged.rac.change <- left_join(plot, rac.change, by = "plot")

rs <- rank_shift(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot")


cs <- community_stability(rac.cover,
                    time.var = "year",
                    abundance.var = "percent_cover",
                    replicate.var = "plot")
turn <- turnover(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  metric = "total")


#-----Assessing species Dominance---------#
#Importing Species Data
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
#----------RAC Analysis----------#
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis

rank.lm<- lmer(rank_change ~ litter*removal + (1|block), data = merged.rac.change)
summary(rank.lm)
Anova(rank.lm)
emmip(rank.lm, litter ~ removal)
emmeans(rank.lm, pairwise ~ removal|litter)

hist(rac.cover$percent_cover, breaks = 50,
     main = "Histogram of Percent Cover", xlab = "Percent Cover")
summary(rac.cover$percent_cover)

quantile(rac.cover$percent_cover, probs = c(0.25, 0.5, 0.75, 0.90))

#----------Nearest Neighbor Analysis----------#
comb.cov.NN <- subset(cover.comb, select = c('year','plot','code','functional_group','percent_cover','nearest_neighbor'))

