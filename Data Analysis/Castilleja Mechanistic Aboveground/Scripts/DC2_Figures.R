setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")

#packages
library(ggeffects)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(ggnewscale)
library(ggthemes)
library(ggcharts)
library(ggpattern)

#----------------------------------------------------------#
#---------------------DIVERSITY PLOTS----------------------#
#----------------------------------------------------------#
#load in data:
diversity <- readRDS("Processed Data/Plant Diversity Full.rds")
diversity %<>% filter(year != "0")

#Calculating Mean and SE
div.mean <- diversity %>% 
  group_by(year,removal) %>% 
  summarise(mean = mean(div), se = sd(div)/sqrt(n()))

rich.mean <- diversity %>% 
  group_by(year,removal) %>% 
  summarise(mean = mean(rich), se = sd(rich)/sqrt(n()))

even.mean <- diversity %>% 
  group_by(year,removal) %>% 
  summarise(mean = mean(even), se = sd(even)/sqrt(n()))


rich.plot <- ggplot(data = rich.mean, aes(x = year, y = mean)) +
  geom_point(data = diversity, aes(x = year, y = rich, color = removal),
             position = position_jitterdodge(0.2, dodge.width = .3), size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#909256", "#C1A575")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(shape = removal, color = removal), shape = 18, size = 4.5, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = removal),
                position = position_dodge(width = 0.2), width = 0.07) +
  geom_line(aes(color = removal, group = removal), 
            position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#333d29", "#A17D5D")) +  # ← your new colors here
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent",
                                    color = "gray", linewidth = 0.12)) +
  ylim(5, 25) +
  labs(x = "Growing Season", y = "Species Richness")

rich.plot

even.plot <- ggplot(data = even.mean, aes(x = year, y = mean)) +
  geom_point(data = diversity, aes(x = year, y = even, color = removal),
             position = position_jitterdodge(0.2, dodge.width = .3), size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#909256", "#C1A575")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(shape = removal, color = removal), shape = 18, size = 4.5, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = removal),
                position = position_dodge(width = 0.2), width = 0.07) +
  geom_line(aes(color = removal, group = removal), 
            position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#333d29", "#A17D5D")) +  # ← your new colors here
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent",
                                    color = "gray", linewidth = 0.12)) +
  ylim(0.5, 1) +
  labs(x = "Growing Season", y = "Species Eveness")

even.plot

div.plot <- ggplot(data = div.mean, aes(x = year, y = mean)) +
  geom_point(data = diversity, aes(x = year, y = div, color = removal),
             position = position_jitterdodge(0.2, dodge.width = .3), size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#909256", "#C1A575")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(shape = removal, color = removal), shape = 18, size = 4.5, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = removal),
                position = position_dodge(width = 0.2), width = 0.07) +
  geom_line(aes(color = removal, group = removal), 
            position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#333d29", "#A17D5D")) +  # ← your new colors here
  theme_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent",
                                    color = "gray", linewidth = 0.12)) +
  ylim(1, 3) +
  labs(x = "Growing Season", y = "Shannon Diversity")

div.plot

div.time.plots <- ggarrange(div.plot, rich.plot, even.plot,
                            labels = c("A", "B","C"), 
                            nrow = 3, ncol = 1,
                            common.legend = TRUE)

div.time.plots

ggsave(plot = div.time.plots, filename = 'Figures and Tables/Diversity_Season.png',
       width = 5, height = 12, units = "in", dpi = 800)


#---------------Delta Diversity---------------#
delta_diversity <- readRDS("Processed Data/Delta Diversity.rds")

nudge_value <- 0.6 

delta_long <- delta_diversity %>%
  mutate(
    field_plot2 = as.numeric(gsub("A|B", "", field_plot)),
    rich_before = rich_yr1,                               # year-1 value, kept for the arrow start
    arrow_draw  = ifelse(delta_rich == 0, "no arrow", "arrow"),
    rich_max    = pmax(rich_yr1, rich_yr3)
  ) %>%
  select(field_plot, field_plot2, removal, rich_before, delta_rich,
         arrow_draw, rich_max, rich_yr1, rich_yr3) %>%
  pivot_longer(c(rich_yr1, rich_yr3), names_to = "year", values_to = "rich") %>%
  mutate(
    year         = recode(year, rich_yr1 = "Year 1", rich_yr3 = "Year 3"),
    text_nudge_x = ifelse(rich == rich_max, 1.2, -1.2)
  ) %>%
  arrange(field_plot, year)

delta.graph <- delta_long %>%
  ggplot(aes(x = rich, y = field_plot2)) +
  geom_segment(data = . %>% filter(year == "Year 3" & arrow_draw == "arrow"),
    aes(x = ifelse(delta_rich > 0, rich_before + 0.35, rich_before - 0.35),
        xend = ifelse(delta_rich > 0, rich - 0.35,         rich + 0.35)),
    arrow = arrow(angle = 25, type = "closed", length = unit(0.24, "cm")),
    color = "gray55") +
  geom_point(aes(color = year), size = 3.5) +
  geom_text(aes(label = rich, x = rich + text_nudge_x), size = 3.25) +
  scale_color_manual(values = c("#cbbbaa", "#45463e")) +
  theme_classic() +
  labs_pubr() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent",
                                    color = "gray23", linewidth = 0.05)) +
  xlim(6, 28) +
  scale_y_continuous(breaks = seq(1, 20, by = 1)) +
  facet_wrap(~ removal) +
  labs(x = "Species Richness", y = "Paired Plot")

ggsave(plot = delta.graph, filename = 'Figures and Tables/Delta_Diversity.png',
       width = 10, height = 12, units = "in", dpi = 600)

#----------------------------------------------------------#
#---------------PHYSICAL STRUCTURE PLOTS-------------------#
#----------------------------------------------------------#

envi <- readRDS("Processed Data/Environmental Cover.rds")
envi %<>% filter(year != "0")

envi.total <- envi %>% 
  group_by(year,removal) %>% 
  summarise(mean = mean(total_cover),
                   se = sd(total_cover)/sqrt(n()))

envi.total %<>% mutate(pat = ifelse(removal == "Removal", "stripe", "none"))

ggplot(envi.total, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  facet_wrap(~year) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Environmental Cover (Bareground + Litter + Rock)")


#----------------------------------------------------------#
#------------------PRODUCTIVITY PLOTS----------------------#
#----------------------------------------------------------#
biomass <- readRDS("Processed Data/Plant Biomass.rds")

ggplot(biomass, aes(x = removal, y = total_no_cas,fill = removal, color = removal)) +
  geom_boxplot( lwd = 0.7, outlier.shape = NA,
                position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.10),
             alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c("#333d29", "#4A3D21")) +
  scale_fill_manual(values = c("#c5c6af", "#D3BC8D")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  ylim(0,425) +
  labs(x = "Parasite", y = "Total Biomass (No Castilleja)")

#----------------------------------------------------------#
#--------------------TURNOVER PLOTS------------------------#
#----------------------------------------------------------#
total_turn <- readRDS("Processed Data/Species Turnover.rds")

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

#----------------------------------------------------------#
#-------------COMMUNITY COMPOSITION PLOTS------------------#
#----------------------------------------------------------#
NMDS <- readRDS("Processed Data/Community NMDS.rds")
NMDS %<>% mutate(removal = recode(removal,
                        C = "Present",
                        R = "Removed"))

sig_labels <- data.frame(
  litter = c("Castilleja", "Mixed", "Community", "Control"),
  label  = c("p = 0.048", "p = 0.081", "n.s.", "n.s."))

ggplot(NMDS, aes(NMDS1, NMDS2, color = removal, fill = removal)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, color = NA) +
  stat_ellipse(linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ litter) +
  geom_text(data = sig_labels, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3.2) +
  scale_color_manual(values = c(Present = "#333d29", Removed = "#b6ad90")) +
  scale_fill_manual(values  = c(Present = "#333d29", Removed = "#b6ad90")) +
  coord_equal() + 
   stheme_bw() +
  labs(color = "Parasite", fill = "Parasite")

#----------------------------------------------------------#
#------------------NEAREST NEIGHBOR PLOTS------------------#
#----------------------------------------------------------#

nn_yearly_Pre <- readRDS("Processed Data/NN Pre.rds")

pre.nearest <- ggplot(nn_yearly_Pre, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance (Cover)", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.101, 0.188811189), label = "ELGL", color = "grey22", nudge_y = - 0.007, size = 3) +
  geom_text(aes(0.107, 0.020979021), label = "LIPO", color = "grey22", nudge_y = - 0.007, size = 3) +
  theme_minimal() +
  ggtitle("2023 Pre") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
pre.nearest

nn_yearly_X2023 <- readRDS("Processed Data/NN 23.rds")

y1.nearest <- ggplot(nn_yearly_X2023, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance (Cover)", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  geom_text(aes(0.0985, 0.24590164), label = "ELGL", color = "grey22", nudge_y = - 0.0075, size = 3) +
  geom_text(aes(0.0475, 0.16393443), label = "LUAR", color = "grey22", nudge_y = - 0.0075, size = 3) +
  ggtitle("2023 Post") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
y1.nearest

nn_yearly_X2024 <- readRDS("Processed Data/NN 24.rds")

y2.nearest <- ggplot(nn_yearly_X2024, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance (Cover)", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("2024") +
  theme(legend.position = "none") +
  geom_text(aes(0.0435, 0.13253012), label = "ELGL", color = "grey22", nudge_y = - 0.005, size = 3) +
  geom_text(aes(0.058, 0.09638554), label = "MESP", color = "grey22", nudge_y = 0.005, size = 3) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))

y2.nearest

nn_yearly_X2025 <- readRDS("Processed Data/NN 25.rds")

y3.nearest <- ggplot(nn_yearly_X2025, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance (Cover)", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  geom_text(aes(0.1069609991, 0.04210526), label = "HEQU", color = "grey22", nudge_y = - 0.004, size = 3) +
  geom_text(aes(0.0668506244, 0.11578947), label = "ELGL", color = "grey22", nudge_y = - 0.004, size = 3) +
  geom_text(aes(0.0391147271, 0.08421053), label = "POGR", color = "grey22", nudge_y = - 0.004, size = 3) +
  ggtitle("2025") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
y3.nearest


nearestplots <- ggarrange(pre.nearest, y1.nearest, y2.nearest, y3.nearest,
                          labels = c("A", "B","C", "D"), 
                          nrow = 2, ncol = 2)

nearestplots 

ggsave(plot = nearestplots, filename = 'Figures and Tables/Nearest Neighbor.png',
       width = 12, height = 8, units = "in", dpi = 800)

#----------------------------------------------------------#
#----------------------CASTILLEJA DATA---------------------#
#----------------------------------------------------------#
cast <- readRDS("Processed Data/Castilleja Summary.rds")

cast_rem <- cast %>% filter(removal == "C")

ggplot(cast_rem, aes(cas_cover, plant_cover, color = year)) +
  geom_point(aes(color = factor(year)), alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_manual( values=c("#335c67", "#e09f3e","#540b0e")) +
  labs(x = "Castilleja cover (%)", y = "Co-occurring plant cover (%)", color = "Year") +
  theme_bw()



