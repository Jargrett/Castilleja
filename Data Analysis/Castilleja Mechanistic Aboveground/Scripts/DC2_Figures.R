setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")

#packages
library(tidyverse)#for data wrangling and restructuring
library(ggeffects)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(ggnewscale)
library(ggthemes)
library(ggcharts)
library(ggpattern)
library(ggpmisc)
library(magrittr)#for data wrangling and restructuring
library(showtext)
library(ggrepel)
library(ggstar)
library(sysfonts)
library(sf)
library(ggspatial)


conflicts_prefer(ggplot2::annotate)

#graph themes
font_add(family = "Times", regular = "Times New Roman.ttf",
         bold = "Times New Roman Bold.ttf",
         italic = "Times New Roman Italic.ttf")
showtext_auto()

theme_pub <- theme_bw(base_size = 12, base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = unit(-1.4, "mm"),                    # ticks point INWARD
    axis.text.x  = element_text(colour = "black", margin = margin(t = 4)),
    axis.text.y = element_text(colour = "black", margin = margin(r = 4)),
    axis.title = element_text(colour = "black"),
    legend.key = element_blank(),
    legend.background = element_blank()
  )

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

div.plot <- ggplot(data = div.mean, aes(x = year, y = mean)) +
  geom_line(aes(group = removal), position = position_dodge(width = 0.15),
            colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = removal),
                position = position_dodge(width = 0.15), width = 0.1,
                colour = "black", linewidth = 0.5) +
  geom_point(aes(fill = removal), shape = 21, size = 2, stroke = 0.6,
             colour = "black", position = position_dodge(width = 0.15)) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  ylim(1.5, 2.5) +
  labs(x = "Growing season", y = "Shannon diversity", fill = "Castilleja") +
  theme_pub +
  theme(legend.position = c(0.05, 0.94),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        legend.margin = margin(3, 4, 3, 4),
        legend.title = element_text(size = 10, family = "Times", hjust = 0.5),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key = element_blank(),
        legend.key.size = unit(4, "mm"),
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 11, family = "Times"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10, family = "Times"),
        strip.background = element_blank())
div.plot


rich.plot <- ggplot(data = rich.mean, aes(x = year, y = mean)) +
  geom_line(aes(group = removal), position = position_dodge(width = 0.15),
            colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = removal),
                position = position_dodge(width = 0.15), width = 0.1,
                colour = "black", linewidth = 0.5) +
  geom_point(aes(fill = removal), shape = 21, size = 2, stroke = 0.6,
             colour = "black", position = position_dodge(width = 0.15)) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  ylim(10,18) +
  labs(x = "Growing season", y = "Species richness", fill = "Castilleja") +
  theme_pub +
  theme(legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 11, family = "Times"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10, family = "Times"),
        strip.background = element_blank())
rich.plot


even.plot <- ggplot(data = even.mean, aes(x = year, y = mean)) +
  geom_line(aes(group = removal), position = position_dodge(width = 0.15),
            colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = removal),
                position = position_dodge(width = 0.15), width = 0.1,
                colour = "black", linewidth = 0.5) +
  geom_point(aes(fill = removal), shape = 21, size = 2, stroke = 0.6,
             colour = "black", position = position_dodge(width = 0.15)) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  ylim(0.7, 0.9) +
  labs(x = "Growing season", y = "Species evenness", fill = "Castilleja") +
  theme_pub +
  theme(legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 11, family = "Times"),
        strip.text = element_text(size = 10, family = "Times"),
        strip.background = element_blank())
even.plot


div.time.plots <- ggarrange(div.plot, rich.plot, even.plot,
                            labels = c("A", "B", "C"),
                            font.label = list(family = "Times", size = 12),
                            nrow = 3, ncol = 1)
div.time.plots

showtext_opts(dpi = 800)
ggsave(plot = div.time.plots, filename = 'Figures and Tables/Diversity_Season.png',
       width = 3, height = 6.5, units = "in", dpi = 800)

ggsave(plot = div.time.plots, filename = 'Figures and Tables/Diversity_Season.pdf',
       width = 3, height = 6.5, units = "in", device = cairo_pdf)


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
               aes(x = ifelse(delta_rich > 0, rich_before + 0.4, rich_before - 0.4),
                   xend = ifelse(delta_rich > 0, rich - 0.4, rich + 0.4)),
               arrow = arrow(angle = 28, type = "open", length = unit(0.15, "cm")),
               color = "black",
               linewidth = 0.5,
               lineend = "round") +
  geom_point(aes(fill = year), shape = 21, size = 2, stroke = 0.5, colour = "black") +
  geom_text(aes(label = rich, x = rich + text_nudge_x -0.05), size = 2.5, family = "Times") +
  scale_fill_manual(values = c("Year 1" = "white", "Year 3" = "black")) +
  scale_x_continuous(breaks = seq(0, 25, by = 5)) +
  scale_y_continuous(breaks = seq(1, 20, by = 1)) +
  facet_wrap(~ removal) +
  labs(x = "Species richness", y = "Paired plot", fill = "Year") +
  theme_pub +
  theme(legend.position = c(0.985, 0.97),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        legend.margin = margin(4, 5, 4, 5),
        legend.title = element_text(size = 11,hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        text = element_text(size = 10, family = "Times"),
        axis.text = element_text(size = 11, family = "Times"),
        axis.title = element_text(size = 15, family = "Times"),
        axis.ticks.length = unit(1.4, "mm"),
        strip.text = element_text(size = 11, family = "Times"),
        strip.background = element_blank())

delta.graph

showtext_opts(dpi = 800)
ggsave(plot = delta.graph, filename = 'Figures and Tables/Delta_Diversity.png',
       width = 4.5, height = 5, units = "in", dpi = 800)
ggsave(plot = delta.graph, filename = 'Figures and Tables/Delta_Diversity.pdf',
       width = 4.5, height = 5, units = "in", device = cairo_pdf)

#----------------------------------------------------------#
#-----------------------COMPOSITION------------------------#
#----------------------------------------------------------#
NMDS <- readRDS("Processed Data/Community NMDS.rds")
NMDS %<>% mutate(removal = recode(removal,
                                  C = "Present",
                                  R = "Removed"))
#--------per litter--------#
sig_labels <- data.frame(
  litter = c("Castilleja", "Mixed", "Community", "Control"))

community.comp <- ggplot(NMDS, aes(NMDS1, NMDS2, fill = removal)) +
  geom_point(shape = 21, size = 2, stroke = 0.6, colour = "black") +
  xlim(-1, 1) +
  scale_fill_manual(values = c(Present = "black", Removed = "white")) +
  theme_bw(base_size = 10, base_family = "Times") +
  facet_wrap(~litter) +
  labs(fill = "Parasite") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.text.x = element_text(colour = "black", family = "Times", margin = margin(t = 4)),
        axis.text.y = element_text(colour = "black", family = "Times", margin = margin(r = 4)),
        axis.title = element_text(colour = "black", family = "Times"),
        strip.text = element_text(size = 10, family = "Times"),
        strip.background = element_blank(),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        legend.margin = margin(4, 5, 4, 5),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key = element_blank(),
        legend.key.size = unit(4, "mm"))
community.comp

showtext_opts(dpi = 800)
ggsave(plot = community.comp, filename = 'Figures and Tables/Composition.png',
       width = 5.5, height = 5.5, units = "in", dpi = 800)

ggsave(plot = community.comp, filename = 'Figures and Tables/Composition.pdf',
       width = 5.5, height = 5.5, units = "in", device = cairo_pdf)

#--------total------#
NMDS <- readRDS("Processed Data/Community NMDS.rds")
NMDS <- readRDS("Processed Data/Community NMDS.rds") %>%
  mutate(removal = recode(removal, "C" = "Present", "R" = "Removed"),
         removal = factor(removal, levels = c("Present", "Removed")))

community.comp <- ggplot(NMDS, aes(NMDS1, NMDS2, fill = removal, linetype = removal)) +
  stat_ellipse(aes(colour = removal), type = "t", level = 0.95,
               linewidth = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, size = 2, stroke = 0.6, colour = "black") +
  scale_fill_manual(values = c(Present = "black", Removed = "white")) +
  scale_colour_manual(values = c(Present = "black", Removed = "grey45")) +
  scale_linetype_manual(values = c(Present = "solid", Removed = "dashed")) +
  labs(x = "NMDS1", y = "NMDS2") +
  guides(fill = guide_legend(title = NULL, override.aes = list(linetype = 0)),
         linetype = "none") +
  theme_pub +
  theme(legend.position = c(0.97, 0.97),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        legend.margin = margin(4, 5, 4, 5),
        legend.key = element_blank(),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key.size = unit(4, "mm"),
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"))
community.comp

showtext_opts(dpi = 800)
ggsave(plot = community.comp, filename = 'Figures and Tables/Full Composition.png',
       width = 4.5, height = 4.5, units = "in", dpi = 800)
ggsave(plot = community.comp, filename = 'Figures and Tables/Full Composition.pdf',
       width = 4.5, height = 4.5, units = "in", device = cairo_pdf)

#----------------------------------------------------------#
#------------------PRODUCTIVITY PLOTS----------------------#
#----------------------------------------------------------#
biomass <- readRDS("Processed Data/Plant Biomass.rds")
cast    <- readRDS("Processed Data/Castilleja Summary.rds")

# Panel A: total biomass
bio_sum <- biomass %>%
  group_by(removal) %>%
  dplyr::summarise(mean = mean(total_no_cas, na.rm = TRUE),
                   se   = sd(total_no_cas, na.rm = TRUE) / sqrt(sum(!is.na(total_no_cas))),
                   .groups = "drop")

# Panel B: total plant cover, 2025 only
cov_sum <- cast %>%
  filter(Year == 2025) %>%
  mutate(removal = recode(removal, "C" = "Present", "R" = "Removed"),
         removal = factor(removal, levels = c("Present", "Removed"))) %>%
  group_by(removal) %>%
  dplyr::summarise(mean = mean(plant_cover, na.rm = TRUE),
                   se   = sd(plant_cover, na.rm = TRUE) / sqrt(sum(!is.na(plant_cover))),
                   .groups = "drop")

bar_theme <- theme_pub +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(margin = margin(r = 4)),
        legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"),
        strip.background = element_blank(),
        plot.margin = margin(2, 4, 2, 2))

prod <- ggplot(bio_sum, aes(x = removal, y = mean, fill = removal)) +
  geom_col(colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("Present" = "white", "Removed" = "black")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Castilleja", y = "Total biomass (g)") +
  bar_theme

cover.plot <- ggplot(cov_sum, aes(x = removal, y = mean, fill = removal)) +
  geom_col(colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("Present" = "white", "Removed" = "black")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 118), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Castilleja", y = "Plant cover (%)") +
  bar_theme +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.margin = margin(4, 5, 4, 5),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key.spacing.y = unit(4, "pt"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(5, "mm"))
cover.plot

biomass.cover <- ggarrange(prod, cover.plot,
                           labels = c("A", "B"),
                           font.label = list(family = "Times", size = 12),
                           align = "h",
                           ncol = 2, nrow = 1)
biomass.cover

showtext_opts(dpi = 800)
ggsave(plot = biomass.cover, filename = 'Figures and Tables/Biomass_Cover.png',
       width = 5, height = 3, units = "in", dpi = 800)

ggsave(plot = biomass.cover, filename = 'Figures and Tables/Biomass_Cover.pdf',
       width = 5, height = 3, units = "in", device = cairo_pdf)
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

turn.plot <- ggplot(turn_total, aes(x = removal, y = mean, fill = removal)) +
  geom_col(colour = "black", width = 0.85) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2),
                     expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  labs(x = "Castilleja", y = "Total species turnover") +
  theme_pub +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(margin = margin(r = 4)),
        legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(2, 1, 2, 1))

gain.plot <- ggplot(turn_gain, aes(x = removal, y = mean, fill = removal)) +
  geom_col(colour = "black", width = 0.85) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2),
                     expand = expansion(mult = c(0, 0.05))) +
  guides(fill = "none") +
  labs(x = "Castilleja", y = "Species gained") +
  theme_pub +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(margin = margin(r = 4)),
        legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(2, 1, 2, 1))

loss.plot <- ggplot(turn_loss, aes(x = removal, y = mean, fill = removal)) +
  geom_col(colour = "black", width = 0.85) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("Present" = "black", "Removed" = "white")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2),
                     expand = expansion(mult = c(0, 0.05))) +
  guides(fill = guide_legend(title = NULL)) +
  labs(x = "Castilleja", y = "Species lost") +
  theme_pub +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-1.4, "mm"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(margin = margin(r = 4)),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.margin = margin(4, 5, 4, 5),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key.spacing.y = unit(4, "pt"),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(2, 0, 2, 0))

turnover.plots <- ggarrange(turn.plot, gain.plot, loss.plot,
                            labels = c("A", "B", "C"),
                            font.label = list(family = "Times", size = 12),
                            nrow = 1, ncol = 3)

turnover.plots <- annotate_figure(turnover.plots,
                                  bottom = text_grob("Castilleja",
                                                     family = "Times", size = 12))
turnover.plots

showtext_opts(dpi = 800)
ggsave(plot = turnover.plots, filename = 'Figures and Tables/Turnover.png',
       width = 5.5, height = 3, units = "in", dpi = 800)

ggsave(plot = turnover.plots, filename = 'Figures and Tables/Turnover.pdf',
       width = 5.5, height = 3, units = "in", device = cairo_pdf)
#----------------------------------------------------------#
#------------------NEAREST NEIGHBOR PLOTS------------------#
#----------------------------------------------------------#

nn_yearly_Pre <- readRDS("Processed Data/NN Pre.rds")
nn_yearly_X2023 <- readRDS("Processed Data/NN 23.rds")
nn_yearly_X2024 <- readRDS("Processed Data/NN 24.rds")
nn_yearly_X2025 <- readRDS("Processed Data/NN 25.rds")

# --- one editable label table per panel (adjust nx/ny to place each) ------
lab_pre <- data.frame(code = c("ELGL","LIPO"),
                      nx = c(-0.01, -0.006), 
                      ny = c(-0.000, -0.015))
lab_y1  <- data.frame(code = c("ELGL","LUAR"),
                      nx = c(-0.011, -0.011), 
                      ny = c(-0.000, 0.000))
lab_y2  <- data.frame(code = c("ELGL","MESP"),
                      nx = c(0.01, 0.00), 
                      ny = c(-0.000, 0.012))
lab_y3  <- data.frame(code = c("HEQU","ELGL","POGR"),
                      nx = c(-0.006, -0.01, -0.01), 
                      ny = c(-0.015, 0.000, -0.000))

pre_lab <- merge(nn_yearly_Pre,   lab_pre, by = "code")
y1_lab  <- merge(nn_yearly_X2023, lab_y1,  by = "code")
y2_lab  <- merge(nn_yearly_X2024, lab_y2,  by = "code")
y3_lab  <- merge(nn_yearly_X2025, lab_y3,  by = "code")

pre.nearest <- ggplot(nn_yearly_Pre, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method = lm, color = "black", fill = "grey85", se = FALSE, linewidth = 0.6) +
  geom_point(shape = 16, size = 1, colour = "black") +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black", linewidth = 0.35) +
  geom_line(aes(y = upr), linetype = "dashed", col = "black", linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = 0.05),
                     breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(limits = c(-0.05, 0.25),
                     breaks = seq(0, 0.2, length.out = 3),
                     expand = expansion(mult = 0.05)) +
  labs(x = "Relative abundance", y = "Neighbor frequency") +
  geom_text_repel(data = pre_lab, aes(label = code),
                  family = "Times", fontface = "italic", size = 2,
                  nudge_x = pre_lab$nx, nudge_y = pre_lab$ny,
                  box.padding = 0, point.padding = 0.2,
                  force = 0, force_pull = 0,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.colour = "black", segment.size = 0.3, seed = 1) +
  annotate("text", x = -Inf, y = Inf, label = "2023 Pre",
           family = "Times", fontface = "bold", size = 3,
           hjust = -0.15, vjust = 1.6) +
  theme_pub +
  theme(legend.position = "none",
        aspect.ratio = 0.7,
        text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 9, family = "Times"),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
        strip.background = element_blank())
pre.nearest

y1.nearest <- ggplot(nn_yearly_X2023, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method = lm, color = "black", fill = "grey85", se = FALSE, linewidth = 0.6) +
  geom_point(shape = 16, size = 1, colour = "black") +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black", linewidth = 0.35) +
  geom_line(aes(y = upr), linetype = "dashed", col = "black", linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = 0.05),
                     breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(limits = c(-0.06, 0.3),
                     breaks = seq(0, 0.25, length.out = 3),
                     expand = expansion(mult = 0.05)) +
  labs(x = "Relative abundance", y = "Neighbor frequency") +
  geom_text_repel(data = y1_lab, aes(label = code),
                  family = "Times", fontface = "italic", size = 2,
                  nudge_x = y1_lab$nx, nudge_y = y1_lab$ny,
                  box.padding = 0, point.padding = 0.2,
                  force = 0, force_pull = 0,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.colour = "black", segment.size = 0.3, seed = 1) +
  annotate("text", x = -Inf, y = Inf, label = "2023",
           family = "Times", fontface = "bold", size = 3,
           hjust = -0.15, vjust = 1.6) +
  theme_pub +
  theme(legend.position = "none",
        aspect.ratio = 0.7,
        text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 9, family = "Times"),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
        axis.title.y = element_blank(),
        strip.background = element_blank())
y1.nearest

y2.nearest <- ggplot(nn_yearly_X2024, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method = lm, color = "black", fill = "grey85", se = FALSE, linewidth = 0.6) +
  geom_point(shape = 16, size = 1, colour = "black") +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black", linewidth = 0.35) +
  geom_line(aes(y = upr), linetype = "dashed", col = "black", linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = 0.05),
                     breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(limits = c(-0.05, 0.2),
                     breaks = seq(0, 0.15, length.out = 3),
                     expand = expansion(mult = 0.05)) +
  labs(x = "Relative abundance", y = "Neighbor frequency") +
  geom_text_repel(data = y2_lab, aes(label = code),
                  family = "Times", fontface = "italic", size = 2,
                  nudge_x = y2_lab$nx, nudge_y = y2_lab$ny,
                  box.padding = 0, point.padding = 0.2,
                  force = 0, force_pull = 0,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.colour = "black", segment.size = 0.3, seed = 1) +
  annotate("text", x = -Inf, y = Inf, label = "2024",
           family = "Times", fontface = "bold", size = 3,
           hjust = -0.15, vjust = 1.6) +
  theme_pub +
  theme(legend.position = "none",
        aspect.ratio = 0.7,
        text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 10, family = "Times"),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
        strip.background = element_blank())
y2.nearest

y3.nearest <- ggplot(nn_yearly_X2025, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method = lm, color = "black", fill = "grey85", se = FALSE, linewidth = 0.6) +
  geom_point(shape = 16, size = 1, colour = "black") +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black", linewidth = 0.35) +
  geom_line(aes(y = upr), linetype = "dashed", col = "black", linewidth = 0.35) +
  scale_x_continuous(expand = expansion(mult = 0.05),
                     breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(limits = c(-0.06, 0.2),
                     breaks = seq(0, 0.15, length.out = 3),
                     expand = expansion(mult = 0.05)) +
  labs(x = "Relative abundance", y = "Neighbor frequency") +
  geom_text_repel(data = y3_lab, aes(label = code),
                  family = "Times", fontface = "italic", size = 2,
                  nudge_x = y3_lab$nx, nudge_y = y3_lab$ny,
                  box.padding = 0, point.padding = 0.2,
                  force = 0, force_pull = 0,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.colour = "black", segment.size = 0.3, seed = 1) +
  annotate("text", x = -Inf, y = Inf, label = "2025",
           family = "Times", fontface = "bold", size = 3,
           hjust = -0.15, vjust = 1.6) +
  theme_pub +
  theme(legend.position = "none",
        aspect.ratio = 0.7,
        text = element_text(size = 8, family = "Times"),
        axis.title = element_text(size = 9, family = "Times"),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
        axis.title.y = element_blank(),
        strip.background = element_blank())
y3.nearest

nearestplots <- (pre.nearest | y1.nearest) / (y2.nearest | y3.nearest) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(family = "Times", size = 12, face = "bold"))
nearestplots

showtext_opts(dpi = 800)
ggsave(plot = nearestplots, filename = 'Figures and Tables/Nearest Neighbor.png',
       width = 6, height = 4, units = "in", dpi = 800)

ggsave(plot = nearestplots, filename = 'Figures and Tables/Nearest Neighbor.pdf',
       width = 6, height = 4, units = "in", device = cairo_pdf)

#----------------------------------------------------------#
#----------------------CASTILLEJA DATA---------------------#
#----------------------------------------------------------#
cast <- readRDS("Processed Data/Castilleja Summary.rds")
cast_rem <- cast %>% filter(removal == "C")

cas.cover <- ggplot(cast_rem, aes(cas_cover, plant_cover)) +
  geom_smooth(method = "lm", color = "black", fill = "grey85", se = TRUE, linewidth = 0.6) +
  geom_point(aes(shape = factor(year), fill = factor(year)),
             size = 2, stroke = 0.5, colour = "black") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("white", "grey55", "black")) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  labs(x = "Castilleja cover (%)", y = "Co-occurring plant cover (%)",
       shape = "Year", fill = "Year") +
  theme_pub +
  theme(legend.position = c(0.94, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        legend.margin = margin(4, 6, 4, 6),
        legend.text = element_text(size = 9, family = "Times"),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(4, "pt"),
        legend.key.size = unit(3, "mm"),
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 13, family = "Times"))

cas.count <- ggplot(cast_rem, aes(cas_count, plant_cover)) +
  geom_smooth(method = "lm", color = "black", fill = "grey85", se = TRUE, linewidth = 0.6) +
  geom_point(aes(shape = factor(year), fill = factor(year)),
             size = 2, stroke = 0.5, colour = "black") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("white", "grey55", "black")) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  labs(x = "Castilleja Count", y = "Co-occurring plant cover (%)",
       shape = "Year", fill = "Year") +
  theme_pub +
theme(legend.position = "none",
      axis.title.y = element_blank())

cas.plots <- ggarrange(cas.cover, cas.count,
                          labels = c("A", "B"),
                          font.label = list(family = "Times", size = 12),
                          nrow = 1, ncol = 2)

cas.plots
showtext_opts(dpi = 800)
ggsave(plot = cas.plots, filename = 'Figures and Tables/Castilleja metrics.pdf',
       width = 8, height = 4, units = "in", device = cairo_pdf)
ggsave(plot = cas.plots, filename = 'Figures and Tables/Castilleja metrics.pdf',
       width = 8, height = 4, units = "in", device = cairo_pdf)

#----------------------------------------------------------#
#---------------------SPECIES RESPONSE---------------------#
#----------------------------------------------------------#
robinhood <- readRDS("Processed Data/Robinhood Summary.rds")
robinhood %<>% filter(!is.na(occupancy_shift), !is.na(response_ratio))

# --- fill: indicator status (open = indicator, filled = non-indicator) -----
robinhood %<>%
  mutate(ind_fill = factor(ifelse(indicator, "Indicator", "Non-indicator"),
                           levels = c("Indicator", "Non-indicator")))

# --- labels: indicators only -----------------------------------------------
lab_codes <- c("RIMO", "AGGL", "BAVU", "LIPO", "HYCA")
ind_labels <- data.frame(
  code = c("RIMO", "AGGL", "BAVU", "LIPO", "HYCA"),
  nx = c(-0.015, -0.012, -0.015,  0.000, -0.005),
  ny = c( 0.09,  0.000, -0.080, -0.120,  0.100))

robinhood %<>%
  left_join(ind_labels, by = "code") %>%
  mutate(lab = ifelse(code %in% lab_codes, code, NA),
         nx = ifelse(is.na(nx), 0, nx),
         ny = ifelse(is.na(ny), 0, ny))

species <- ggplot(robinhood, aes(x = occupancy_shift, y = response_ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_star(aes(fill = ind_fill), starshape = 15,
            size = 2.6, starstroke = 0.7, colour = "black") +
  scale_fill_manual(values = c("Indicator" = "white", "Non-indicator" = "black")) +
  geom_text_repel(aes(label = lab), family = "Times", size = 3,
                  nudge_x = robinhood$nx, nudge_y = robinhood$ny,
                  box.padding = 0.5, point.padding = 0.4,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.colour = "grey40", segment.size = 0.3, seed = 1) +
  guides(fill = guide_legend(position = "top", nrow = 1,
                             override.aes = list(starshape = 15))) +
  annotate("text", x =  0.055, y =  1, hjust = 0.5, vjust = 1, size = 3.5,
           fontface = "bold", colour = "grey15", family = "Times", lineheight = 0.9,
           label = "Suppressed,\nlost") +
  annotate("text", x = -0.055, y =  1, hjust = 0.5, vjust = 1, size = 3.5,
           fontface = "bold", colour = "grey15", family = "Times", lineheight = 0.9,
           label = "Suppressed,\ngained") +
  annotate("text", x = -0.055, y = -1, hjust = 0.5, vjust = 0, size = 3.5,
           fontface = "bold", colour = "grey15", family = "Times", lineheight = 0.9,
           label = "Facilitated,\ngained") +
  annotate("text", x =  0.055, y = -1, hjust = 0.5, vjust = 0, size = 3.5,
           fontface = "bold", colour = "grey15", family = "Times", lineheight = 0.9,
           label = "Facilitated,\nlost") +
  labs(x = "Change in occupancy", y = "Response ratio (Removal)", fill = NULL) +
  coord_fixed(ratio = 0.1, xlim = c(-0.1, 0.1), ylim = c(-1, 1)) +
  theme_pub +
  theme(legend.position = "top", legend.key = element_blank(),
        legend.key.width = unit(10, "pt"), legend.background = element_blank(),
        legend.key.spacing.x = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        text = element_text(family = "Times"),
        axis.title = element_text(size = 12, family = "Times"))
species

showtext_opts(dpi = 800)
ggsave(plot = species, filename = 'Figures and Tables/Species Response.png',
       width = 5, height = 5, units = "in", dpi = 800)
ggsave(plot = species, filename = 'Figures and Tables/Species Response.pdf',
       width = 5, height = 5, units = "in", device = cairo_pdf)

#----------------------------------------------------------#
#-------------------DRIVERS OF DIVERSITY-------------------#
#----------------------------------------------------------#

coefs<- readRDS( "Processed Data/Drivers Coefficients.rds")

envi_div_z$litter <- relevel(factor(envi_div_z$litter), ref = "Control")

effect.plot <- ggplot(coefs, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.15, colour = "black", linewidth = 0.45) +
  geom_point(aes(fill = type), shape = 21, size = 2.2, stroke = 0.6, colour = "black") +
  scale_fill_manual(values = c("Experimental" = "black", "Environmental" = "white")) +
  guides(fill = "none") +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Standardized effect on species richness", y = NULL) +
  theme_pub +
  theme(legend.position = "none",
        text = element_text(size = 10, family = "Times"),
        axis.title = element_text(size = 12, family = "Times"),
        axis.text.y = element_text(size = 9, family = "Times"),
        strip.text = element_text(size = 10, family = "Times", face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(6, "pt"))
effect.plot
showtext_opts(dpi = 800)
ggsave(plot = effect.plot, filename = 'Figures and Tables/Effect_Sizes.png',
       width = 9, height = 3.5, units = "in", dpi = 800)

ggsave(plot = effect.plot, filename = 'Figures and Tables/Effect_Sizes.pdf',
       width = 9, height = 3.5, units = "in", device = cairo_pdf)
