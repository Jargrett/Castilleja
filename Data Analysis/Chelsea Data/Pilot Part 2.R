#This sets our working directory so that we have access to our files
setwd("~/Desktop/Castilleja/Data Analysis/Chelsea Data")

#The goal of this is to assess how light influences biomas growth for host and parasite

#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
library(car)
library(lme4)
library(emmeans)# for comparison of means
library(ggpubr)
library(ggplot2)
library(ggpattern)

#conflicts
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

#import data
full_data <- read.csv("Pilot - Total.csv")
full_data <- as.data.frame(unclass(full_data),stringsAsFactors=TRUE)
full_data %<>% drop_na(biomass)

#Filter datasets by species so that we can conduct analysis and graph
cain_data <- filter(full_data, plant == "CAIN")#parasite data
cain_data %<>% mutate(pot_type = recode(pot_type,
                          "HP" = "with host",
                          "P" = "alone"))
cain_above <- filter(cain_data, bio_type == "AGB")#parasite data
cain_below <- filter(cain_data, bio_type == "BGB")#parasite data


scsc_data <- filter(full_data, plant == "SCSC")#host data
scsc_data %<>% mutate(pot_type = recode(pot_type,
                                        "HP" = "with parasite",
                                        "H" = "alone"))
scsc_above <- filter(scsc_data, bio_type == "AGB")#parasite data
scsc_below <- filter(scsc_data, bio_type == "BGB")#parasite data
  

#------Analysis-------#
#SCSC Final Height
scsc.height.lm <- lmer(final_height ~ light*pot_type + (1|replicate), data = scsc_data)
summary(scsc.height.lm)
Anova(scsc.height.lm) #pot_type p < 0.001
emmeans(scsc.height.lm, pairwise ~ light*pot.type)
emmip(scsc.height.lm, ~ light ~ pot.type)

#CAIN Final Height
cain.height.lm <- lmer(final_height ~ light*pot_type + (1|replicate), data = cain_data)
summary(cain.height.lm)
Anova(cain.height.lm) #pot_type p < 0.001
emmeans(cain.height.lm, pairwise ~ light*pot.type)
emmip(cain.height.lm, ~ light ~ pot.type)


#SCSC Above
scsc.above.lm <- lmer(biomass ~ light*pot_type + (1|replicate), data = scsc_above)
summary(scsc.above.lm)
Anova(scsc.above.lm) #light p = 0.03, pot_type p < 0.001
emmeans(scsc.above.lm, pairwise ~ light*pot_type)
emmip(scsc.above.lm, ~ light ~ pot_type)

#SCSC Below
scsc.below.lm <- lmer(biomass ~ light*pot_type + (1|replicate), data = scsc_below)
summary(scsc.below.lm)
Anova(scsc.below.lm) #light:pot_type p < 0.0001
emmeans(scsc.below.lm, pairwise ~ light*pot_type)
emmip(scsc.below.lm, ~ light ~ pot_type)

#CAIN Above
cain.above.lm <- lmer(biomass ~ light*pot_type + (1|replicate), data = cain_above)
summary(cain.above.lm)
Anova(cain.above.lm) #light p < 0.001, pot_type p < 0.001
emmeans(cain.above.lm, pairwise ~ light*pot_type)
emmip(cain.above.lm, ~ light ~ pot_type)

#CAIN Below
cain.below.lm <- lmer(biomass ~ light*pot_type + (1|replicate), data = cain_below)
summary(cain.below.lm)
Anova(cain.below.lm) #light:pot_type p < 0.035
emmeans(cain.below.lm, pairwise ~ light*pot_type)
emmip(cain.below.lm, ~ light ~ pot_type)


#remove NA values from rows and columns

#----------Graphing----------#
cain_above <- filter(cain_data, bio_type == "AGB")#parasite data
cain_below <- filter(cain_data, bio_type == "BGB")#parasite data

scsc_above <- filter(scsc_data, bio_type == "AGB")#parasite data
scsc_below <- filter(scsc_data, bio_type == "BGB")#parasite data

cain_above_mean <- cain_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

cain_height_mean <- cain_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(final_height), se = sd(final_height)/sqrt(n()))

cain_below_mean <- cain_below %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

scsc_above_mean <- scsc_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))

scsc_height_mean <- scsc_above %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(final_height), se = sd(final_height)/sqrt(n()))

scsc_below_mean <- scsc_below %>% 
  group_by(light, pot_type) %>% 
  summarise(mean = mean(biomass), se = sd(biomass)/sqrt(n()))


ggplot(data = cain_above, aes(x = pot_type, y = biomass, color = light)) + 
  geom_point() +
  scale_color_manual(values=c("#4b3b40","#b6ad90")) +
  labs(x = "Pot combination", y = "Parasite Aboveground Biomass") +
  ylim(0,0.6)
  

#We will graph Paraite biomass by our two treatments (pot_type & light)
cain.above <- ggplot(cain_above_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Aboveground Biomass") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr() +
  ylim(0,0.4)
cain.above

cain.below <- ggplot(cain_below_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Belowground Biomass") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr() +
  ylim(0,0.4)
cain.below

scsc.above <- ggplot(scsc_above_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Aboveground Biomass") +
  scale_fill_manual(values=c("#45463e","#b3b1be")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr() +
  ylim(0,0.75)
scsc.above

scsc.below <- ggplot(scsc_below_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Belowground Biomass") +
  scale_fill_manual(values=c("#45463e","#b3b1be")) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_pubr() +
  ylim(0,0.75)
scsc.below

biomass.plots <- ggarrange(cain.above, scsc.above, cain.below, scsc.below,
                           labels = c("A", "B","C", "D"), 
                           nrow = 2, ncol = 2)
biomass.plots

#---------Height----------#
cain.height <- ggplot(cain_height_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Parasite Height") +
  scale_fill_manual(values=c("#4b3b40","#b6ad90")) +
  theme_pubr() +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 15))+
  ylim(0,10)
cain.height

scsc.height <- ggplot(scsc_height_mean, aes(x=pot_type, y=mean, fill=light)) + 
  geom_col(position = position_dodge(0.7), width = .6, linewidth = 0.75, alpha = 0.9, size = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.7), width = 0.15) +
  labs(x = "Pot combination", y = "Host Height") +
  scale_fill_manual(values=c("#45463e","#b3b1be")) +
  theme_pubr() +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 15))+
  ylim(0,50)
scsc.height


height.plots <- ggarrange(cain.height, scsc.height,
                           labels = c("A", "B"), 
                           nrow = 1, ncol = 2)
height.plots

ggsave(plot = height.plots, filename = 'height pilot.png',
       width = 10 , height = 8, units = "in", dpi = 600)

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
fill_colors <- c(
  "Host.Alone"     = "#45463e",
  "Host.HP"        = "#b3b1be",
  "Parasite.Alone" = "#4b3b40",
  "Parasite.HP"    = "#b6ad90"
)

# ---------------------------------------------------------------------------
# Shared y limits so zero lines align
# ---------------------------------------------------------------------------
both <- full_data %>%
  filter((plant == "SCSC" & pot_type %in% c("H",  "HP")) | 
           (plant == "CAIN" & pot_type %in% c("P",  "HP"))) %>%
  group_by(plant, light, pot_type, bio_type) %>%
  summarise(mean_mass = mean(biomass, na.rm = TRUE),
            se = sd(biomass,   na.rm = TRUE) / sqrt(n()),
            .groups   = "drop") %>%
  mutate(mean_mass = ifelse(bio_type == "BGB", -mean_mass, mean_mass))

y_max <- max(both$mean_mass + both$se, na.rm = TRUE) * 1.1
y_min <- min(both$mean_mass - both$se, na.rm = TRUE) * 1.1


# ---------------------------------------------------------------------------
# CAIN (Parasite) — summarize
# ---------------------------------------------------------------------------
cain_pd <- full_data %>%
  filter(plant == "CAIN", pot_type %in% c("P", "HP")) %>%
  group_by(light, pot_type, bio_type) %>%
  summarise(mean_mass = mean(biomass, na.rm = TRUE),
    se        = sd(biomass,   na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  mutate(mean_mass = ifelse(bio_type == "BGB", -mean_mass, mean_mass),
    tissue = factor(ifelse(bio_type == "AGB", "Shoot", "Root"),
                         levels = c("Shoot", "Root")),
    light = factor(light, levels = c("ambient", "low"),
                         labels = c("Ambient", "Low")),
    pot_label   = factor(ifelse(pot_type == "P", "Alone", "HP"),
                         levels = c("Alone", "HP")),
    color_group = interaction("Parasite", pot_label),
    bar_fill    = fill_colors[as.character(color_group)],
    x_num = case_when(
      light == "Ambient" & pot_label == "Alone" ~ 1.0,
      light == "Ambient" & pot_label == "HP"    ~ 1.5,
      light == "Low"     & pot_label == "Alone" ~ 2.2,
      light == "Low"     & pot_label == "HP"    ~ 2.7
    )
  )

cain_shoot <- cain_pd %>% filter(tissue == "Shoot")
cain_root  <- cain_pd %>% filter(tissue == "Root")

# ---------------------------------------------------------------------------
# CAIN (Parasite) — plot
# ---------------------------------------------------------------------------
cain.biomass <- ggplot(cain_pd, aes(x = x_num, y = mean_mass,
                                    fill = bar_fill, pattern = tissue)) +
  geom_col_pattern(position = "stack", width = 0.45, color = "black",
                   linewidth = 0.3, pattern_fill = "transparent",
                   pattern_colour = "white", pattern_angle = 45,
                   pattern_density  = 0.06, pattern_spacing  = 0.018,
                   pattern_key_scale_factor = 0.4) +
  geom_errorbar(data = cain_shoot, 
                aes(x = x_num, ymin = mean_mass - se, ymax = mean_mass + se), 
                width = 0.12, linewidth = 0.35) +
  geom_errorbar(data = cain_root, 
                aes(x = x_num, ymin = mean_mass - se, ymax = mean_mass + se), 
                width = 0.12, linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_identity() +
  scale_pattern_manual(name = "Tissue", 
                       values = c("Shoot" = "none", "Root" = "stripe"),
                       guide  = guide_legend(override.aes = list(fill = "grey60", 
                                            pattern_colour = "white"))) +
  scale_fill_identity(name = "Pot type", 
                      guide  = guide_legend(override.aes = list(fill = c("#4b3b40", "#b6ad90"),
                      pattern = c("none", "none"))), labels = c("P", "HP"), breaks = c("#4b3b40", "#b6ad90")) +
  scale_x_continuous(breaks = c(1.25, 2.45), labels = c("Ambient", "Low"),limits = c(0.6, 3.1)) +
  scale_y_continuous(labels = abs, limits = c(y_min, y_max)) +
  labs(title = "Parasite", x = "Light", y = "Plant dry mass (g)") +
  theme_pubr() +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5), legend.position = "top")

cain.biomass
# ---------------------------------------------------------------------------
# SCSC (Host) — summarize
# ---------------------------------------------------------------------------
scsc_pd <- full_data %>%
  filter(plant == "SCSC", pot_type %in% c("H", "HP")) %>%
  group_by(light, pot_type, bio_type) %>%
  summarise(
    mean_mass = mean(biomass, na.rm = TRUE),
    se        = sd(biomass,   na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  mutate(
    mean_mass   = ifelse(bio_type == "BGB", -mean_mass, mean_mass),
    tissue      = factor(ifelse(bio_type == "AGB", "Shoot", "Root"),
                         levels = c("Shoot", "Root")),
    light       = factor(light, levels = c("ambient", "low"),
                         labels = c("Ambient", "Low")),
    pot_label   = factor(ifelse(pot_type == "H", "Alone", "HP"),
                         levels = c("Alone", "HP")),
    color_group = interaction("Host", pot_label),
    bar_fill    = fill_colors[as.character(color_group)],
    x_num = case_when(
      light == "Ambient" & pot_label == "Alone" ~ 1.0,
      light == "Ambient" & pot_label == "HP"    ~ 1.5,
      light == "Low"     & pot_label == "Alone" ~ 2.2,
      light == "Low"     & pot_label == "HP"    ~ 2.7))

scsc_shoot <- scsc_pd %>% filter(tissue == "Shoot")
scsc_root  <- scsc_pd %>% filter(tissue == "Root")

# ---------------------------------------------------------------------------
# SCSC (Host) — plot
# ---------------------------------------------------------------------------
scsc.biomass <- ggplot(scsc_pd, aes(x = x_num, y = mean_mass,
                                    fill = bar_fill, pattern = tissue)) +
  geom_col_pattern(position = "stack", width = 0.45, color = "black", 
                   linewidth = 0.3, pattern_fill = "transparent", pattern_colour = "white",
                   pattern_angle = 45, pattern_density = 0.06, pattern_spacing = 0.018,
                   pattern_key_scale_factor = 0.4) +
  geom_errorbar(data = scsc_shoot, aes(x = x_num, ymin = mean_mass - se, ymax = mean_mass + se),
                width = 0.12, linewidth = 0.35) +
  geom_errorbar(data = scsc_root, aes(x = x_num, ymin = mean_mass - se, ymax = mean_mass + se),
                width     = 0.12, linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_identity(name = "Pot type", 
                      guide  = guide_legend(override.aes = list(fill = c("#45463e", "#b3b1be"),
        pattern = c("none", "none"))), labels = c("H", "HP"), breaks = c("#45463e", "#b3b1be")) +
  scale_pattern_manual(
    name   = "Tissue",
    values = c("Shoot" = "none", "Root" = "stripe"),
    guide  = guide_legend(override.aes = list(fill = "grey60", pattern_colour = "white"))) +
  scale_x_continuous(breaks = c(1.25, 2.45), labels = c("Ambient", "Low"), limits = c(0.6, 3.1)) +
  scale_y_continuous(labels = abs, limits = c(y_min, y_max)) +
  labs(title = "Host", x = "Light", y = "Plant dry mass (g)") +
  theme_pubr() +
  theme(axis.text=element_text(size = 12), axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5), legend.position = "top")

scsc.biomass
# ---------------------------------------------------------------------------
# Combine — legend pulled from scsc.biomass
# ---------------------------------------------------------------------------
shared_legend <- get_legend(scsc.biomass)

biomass.plots <- ggarrange(cain.biomass, scsc.biomass,
  labels  = c("A", "B"),
  nrow    = 1,
  ncol    = 2,
  widths  = c(1, 1, 0.3), common.legend = T)

biomass.plots



ggsave(plot = biomass.plots, filename = 'AGBBGB final3.png',
       width = 10, height = 10, units = "in", dpi = 600)

