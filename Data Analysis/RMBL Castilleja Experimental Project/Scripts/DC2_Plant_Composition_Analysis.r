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



#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::filter)
conflicts_prefer(lme4::lmer)
conflicts_prefer(dplyr::summarise)

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
                          "C" = "Present",
                          "R" = "Removed")) %>% 
  mutate(year = recode(year,
                       "Pre" = "0",
                       "2023" = "1",
                       "2024" = "2",
                       "2025" = "3",))

#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
cover.comb.clean$plot <- as.factor(cover.comb.clean$plot)
cover.comb.clean$pair <- as.factor(cover.comb.clean$pair)
cover.comb.clean$block <- as.factor(cover.comb.clean$block)
saveRDS(cover.comb.clean, "Processed Data/Cleaned Cover.rds")


#import plot data
plot_data <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot_data$plot <- as.factor(plot_data$plot)
plot_data$pair <- as.factor(plot_data$pair)
plot_data$block <- as.factor(plot_data$block)
plot_data %<>%
  mutate(removal = recode(removal,
                          "C" = "Present",
                          "R" = "Removed"))

#----Codyn Work-----#
clean_cover <- readRDS("Processed Data/Cleaned Cover.rds")
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

comp.lmm <- lmer(distance ~ litter*removal + (1|block) + (1|pair), data = total_comp_change)
summary(comp.lmm)
Anova(comp.lmm)#removal: chisq = 4.0428, Df= 1, p = 0.0445*
emmeans(comp.lmm, pairwise ~ removal|litter)
emmip(comp.lmm, removal ~ litter)

ggplot(total_comp_change, aes(x = removal, y = distance, fill = removal)) +
  geom_boxplot(alpha = 0.8, color = "black", width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
    dodge.width = 0.10),
    alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c("#333d29", "#4A3D21")) +
  scale_fill_manual(values = c("#c5c6af", "#D3BC8D")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Castilleja",
       y = "Rate of compositional change") +
  theme(legend.position = "none")


#Community Stability
cs <- community_stability(clean_cover,
                          time.var = "year",
                          abundance.var = "percent_cover",
                          replicate.var = "plot")
com_stability <- cs  %>%  
  left_join(plot_data, comp_change, by = "plot")

stab.lmm <- lmer(stability ~ litter*removal + (1|block) + (1|pair), data = com_stability)
summary(stab.lmm)
Anova(stab.lmm)#litter: chisq = 7.3853, Df= 3, p = 0.06058.
emmeans(stab.lmm, pairwise ~ removal|litter)
emmip(stab.lmm, removal ~ litter)

#Compare whether the rank abundance structure differs between removal at end of study
between_cover <- clean_cover %>% 
  filter (year != "2")
  
  
total_turn <- turnover(between_cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot")
  
app_turn <- turnover(between_cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  metric = "appearance")

diss_turn <- turnover(between_cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  metric = "disappearance")

total_turn <- left_join(plot_data, total_turn, by = "plot")
total_turn <- left_join(total_turn, app_turn, by = "plot")
total_turn <- left_join(total_turn, diss_turn, by = "plot")

total_turn %<>% select(-c("year.x", "year.y", "year"))


turn.lmm <- lmer(total ~ litter*removal + (1|block) + (1|pair), data = total_turn)
summary(turn.lmm)
Anova(turn.lmm)#No significance 
emmeans(turn.lmm, pairwise ~ removal|litter)
emmip(turn.lmm, litter ~ removal)

gain.lmm <- lmer(appearance ~ litter*removal + (1|block) + (1|pair), data = total_turn)
summary(gain.lmm)
Anova(gain.lmm)#removal: Chisq = 11.9239, Df= 1, p = 0.0005542***
emmeans(gain.lmm, pairwise ~ removal|litter)
emmip(gain.lmm, removal ~ litter)

loss.lmm <- lmer(disappearance ~ litter*removal + (1|block) + (1|pair), data = total_turn)
summary(loss.lmm)
Anova(loss.lmm)#removal: Chisq = 12.6935, Df= 1, p = 0.0003669***
emmeans(loss.lmm, pairwise ~ removal|litter)
emmip(loss.lmm, removal ~ litter)

turn_total <- total_turn %>% 
  group_by(removal) %>% 
  dplyr::summarise(mean = mean(total),
                   se = sd(total)/sqrt(n()))

turn_gain <- total_turn %>% 
  group_by(removal) %>% 
  dplyr::summarise(mean = mean(appearance),
                   se = sd(appearance)/sqrt(n()))

turn_loss <- total_turn %>% 
  group_by(removal) %>% 
  dplyr::summarise(mean = mean(disappearance),
                   se = sd(disappearance)/sqrt(n()))

turn_total %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "Removed", "stripe", "none"))
turn_gain %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "Removed", "stripe", "none"))
turn_loss %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "Removed", "stripe", "none"))


turn.plot <- ggplot(turn_total, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.75, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Total species turnover") +
  ylim(0,0.75)

gain.plot <- ggplot(turn_gain, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.75, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Proportion of species gained") +
  ylim(0,0.75)

loss.plot <- ggplot(turn_loss, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.75, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Proportion of species lost") +
  ylim(0,0.75)


turnover.plots <- ggarrange(turn.plot, gain.plot, loss.plot,
                          labels = c("A", "B", "C"), 
                          nrow = 1, ncol = 3)
turnover.plots

#---species identity changes---#
pa_cover <- clean_cover %>%
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

meta <- plot_data %>% 
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
  mutate(rate = n / total_plots) %>%
  select(-n, -total_plots) %>%
  pivot_wider(names_from = change, values_from = rate, values_fill = 0)

# Step 3: calculate difference between treatments
diff <- rate %>%
  pivot_longer(cols = c(gain, loss, persist), names_to = "change", values_to = "rate") %>%
  pivot_wider(names_from = removal, values_from = rate, values_fill = 0) %>%
  mutate(diff = Removed - Present) %>%  # positive = more common in Removed
  arrange(desc(abs(diff)))

# How many plots does each species appear in at all?
species_frequency <- persist %>%
  group_by(code) %>%
  summarise(
    n_plots = n_distinct(plot),
    n_gains = sum(change == "gain"),
    n_losses = sum(change == "loss"),
    n_persist = sum(change == "persist")
  ) %>%
  arrange(n_plots)


# Only keep species present in at least 5 plots
common_species <- species_frequency %>%
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

common_affinity <- affinity %>% filter(code %in% common_species)

# join functional group info
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


#-----species abundance changes----#

#import plot data
plot_data <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot_data$plot <- as.factor(plot_data$plot)
plot_data$pair <- as.factor(plot_data$pair)
plot_data$block <- as.factor(plot_data$block)
plot_data %<>%
  mutate(removal = recode(removal,
                          "C" = "Present",
                          "R" = "Removed"))

plot_d <- plot_data %>% 
  select(plot,removal)

abund.cover <- abundance_change(clean_cover,
                     time.var = "year",
                     species.var = "code",
                     abundance.var = "percent_cover",
                     replicate.var = "plot",
                     reference.time = "1")

abund.cover %<>%
  filter(year2 != "2") %>% 
  select(plot, code, change) %>% 
  rename(percent_change = change)


abund_merge_cover <- merge(abund.cover, plot_d, by = "plot")

abund_mean_cover <- abund_merge_cover %>% 
  group_by(code, removal) %>%
  summarise(mean_change_cover = mean(percent_change),
            se_change_cover = sd(percent_change)/sqrt(n()))


abund.count <- abundance_change(clean_cover,
                                time.var = "year",
                                species.var = "code",
                                abundance.var = "count",
                                replicate.var = "plot",
                                reference.time = "1")

abund.count %<>%
  filter(year2 != "2") %>% 
  select(plot, code, change) %>% 
  rename(count_change = change)
  

abund_merge_count <- merge(abund.count, meta, by = "plot")
  
abund_mean_cover <- abund_merge_cover %>% 
  group_by(code, removal) %>%
  summarise(mean_change_cover = mean(percent_change),
            se_change_cover = sd(percent_change)/sqrt(n()))

abund_mean_count <- abund_merge_count %>% 
  group_by(code, removal) %>%
  summarise(mean_change_count = mean(count_change),
            se_change_count = sd(count_change)/sqrt(n()))

abund_full <- inner_join(abund_mean_cover, abund_mean_count, by = c("code", 'removal'))


abund.species <- abundance_change(clean_cover,
                                time.var = "year",
                                species.var = "code",
                                abundance.var = "count",
                                replicate.var = "plot",
                                reference.time = "1")


#importing dominance data
species_info <- read.csv("Raw Data/EL Species List - EL.csv")
species_code <- species_info %>%
  distinct(code, functional_group, life_history, growth_form)

dom_year <- readRDS("Processed Data/Dominance y1.rds")
dom_year %<>%
  mutate(removal = recode(removal,
                                 "Control" = "Present",
                                 "Removal" = "Removed")) 

abund_merge_year <- merge(dom_year, abund_full,
                          by = c("code","removal"),
                          all.x = T, all.y = T) %>% 
  drop_na(freq)

yearly.rarity <- full_join(species_code, abund_merge_year, by = "code") %>% 
  drop_na(removal)


dom_total <- readRDS("Processed Data/Dominance total.rds")

dom_total %<>%
  mutate(removal = recode(removal,
                          "Control" = "Present",
                          "Removal" = "Removed")) 

abund_merge_total <- merge(dom_total, abund_full,
                           by = c("code","removal"),
                           all.x = T, all.y = T)

total.rarity <- full_join(species_code, abund_merge_total, by = "code") %>% 
  drop_na(removal) %>% 
  drop_na(mean_change_cover)

ggplot(total.rarity, aes(x = growth_form, y = mean_change_cover, color = dominance)) +
  geom_point() +
  ylim(-0.3, 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~removal)

#-------Abundance change Analysis-----#
#year 1 abund only data
ac.lm <- lm(mean_change_count ~ removal*dominance, data = yearly.rarity)
summary(ac.lm)
Anova(ac.lm)#No significance 

#Functional Group Change
ac.lm <- lm(mean_change_cover ~ growth_form, data = total.rarity)
summary(ac.lm)
Anova(ac.lm)

#Full data
ac.lm <- lm(mean_change_cover ~ functional_group, data = total.rarity)
summary(ac.lm)
Anova(ac.lm)


pa_data <- plant_data %>%
  mutate(present = ifelse(cover > 0, 1, 0))

species_change <- pa_data %>%
  group_by(species, year) %>%
  summarise(present = sum(present)>0) %>%
  arrange(species, year) %>%
  group_by(species) %>%
  mutate(change = present - lag(present))

#convert to matrix format for diversity calculations
library(labdsv)#enables restructuring for ecological analysis
library(vegan)
comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))
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

emerald.pre.matrix <- matrify(emerald.pre)
emerald.23.matrix <- matrify(emerald.23)
emerald.24.matrix <- matrify(emerald.24)
emerald.25.matrix <- matrify(emerald.25)
plot <- read.csv("Raw Data/Site Level Data - EL.csv") #importing metadata
plot.25 <- plot %>% 
  filter(Year == "2025") %>%
  select(-c(Year))
plot.23 <- plot %>% 
  filter(Year == "2023") %>%
  select(-c(Year))
plot.24 <- plot %>% 
  filter(Year == "2024") %>%
  select(-c(Year))
plot.pre <- plot %>% 
  filter(Year == "pre") %>%
  select(-c(Year))

#calculate distance matrix
dist.pre <-vegdist(emerald.pre.matrix, method="bray")
beta.pre <- betadisper(dist.pre, plot$removal)
dist.23 <-vegdist(emerald.23.matrix, method="bray")
beta.23 <- betadisper(dist.23, plot$removal)
dist.24 <-vegdist(emerald.24.matrix, method="bray")
beta.24 <- betadisper(dist.24, plot$removal)
dist.25 <-vegdist(emerald.25.matrix, method="bray")
beta.25 <- betadisper(dist.25, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
nmds.25 <- metaMDS(dist.25, distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=1000) #for publication I recommend 500)
nmds.25#stress value 0.14 which is below .2 so we need to investigate


ordiplot(nmds.25, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds.25))

NMDS <- cbind(plot.25, nmds.scores) #final dataset

cap.mod <- capscale(dist.25 ~ removal*litter + Condition(block), 
                     data = NMDS)

perms <- how(blocks = NMDS$block, nperm = 9999)

anova_res <- anova(cap.mod, permutations = perms, by = "terms")

perm <- adonis2(dist ~ removal*litter, data = NMDS, permutations=9999)
perm

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=litter , shape = removal), size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c("#4B3B40", "#9B7E46", "#A8C7BB", "#808F87")) +
  coord_equal() +
  theme_bw()
