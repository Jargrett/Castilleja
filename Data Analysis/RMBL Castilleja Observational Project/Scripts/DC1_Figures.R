setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")

#Packages
library(ggplot2)
library(ggpubr)
library(hrbrthemes)

castilleja.cover <- readRDS("Processed Data/Total Castilleja Cover.csv")

#Diversity Graphs
#Standard error calculations
castilleja.div <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(div),
                   se = sd(div)/sqrt(n()))

castilleja.rich <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(rich),
                   se = sd(rich)/sqrt(n()))
castilleja.even <- castilleja.cover %>% 
  group_by(castilleja, year) %>% 
  dplyr::summarise(mean= mean(even),
                   se = sd(even)/sqrt(n()))

castilleja.div$colcast = castilleja.div$castilleja
castilleja.rich$colcast = castilleja.rich$castilleja
castilleja.even$colcast = castilleja.even$castilleja
#Diveristy
div.plot <- ggplot(data = castilleja.div, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.6, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.6, colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Shannon Diversity of co-occuring species") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(1.4,2)

div.plot

#Richness
rich.plot <- ggplot(data = castilleja.rich, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.7, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.6, colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species Richness of co-occuring species") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(6,12)

rich.plot

#Evenness
even.plot <- ggplot(data = castilleja.even, aes(x = reorder(castilleja, -mean), y = mean, color = castilleja)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.8, width = 0.09) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), size = 0.7, width = 0.09, color = "grey38") +
  geom_point(shape = 18 ,size = 5.2,colour = "grey38") +
  geom_point(aes(colour=castilleja),shape = 18, size = 4) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species evenness of co-occuring species") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  ylim(0.6,1)

even.plot

diversity.plots <- ggarrange(div.plot, rich.plot, even.plot,
                             labels = c("A", "B","C"), 
                             nrow = 1, common.legend = TRUE)
diversity.plots  

ggsave(plot = diversity.plots, filename = 'Figures and Tables/Diversity.png',
       width = 12 ,height = 5, units = "in", dpi = 600, 
       bg = "transparent")

#Composition
NMDS <- read_rds("Processed Data/NMDS.rds") 

ordination.plot <- ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=castilleja , shape = species), size = 2.2, alpha = 0.8) +
  scale_color_manual( values=c("#D6A839", "#71A4A0")) +
  stat_ellipse(geom = "polygon", segments = 20, linetype = 2, alpha = 0.1, aes(group = site)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = site)) +
  geom_text(aes(-1.53,0.72), label = "EL", color = "grey22", size = 3.5) +
  geom_text(aes(-0.9,0.1), label = "AP", color = "grey22", size = 3.5) +
  geom_text(aes(-0.3,0.58), label = "CC", color = "grey22", size = 3.5) +
  geom_text(aes(0.85,1.85), label = "AL", color = "grey22", size = 3.5) +
  geom_text(aes(0.85,0.9), label = "DC1", color = "grey22", size = 3.5) +
  geom_text(aes(0.23,0), label = "DC2", color = "grey22", size = 3.5) +
  geom_text(aes(-0.14,-0.5), label = "JH", color = "grey22", size = 3.5) +
  coord_equal() +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme_bw()
ordination.plot

ggsave(plot = ordination.plot, filename = 'Figures and Tables/Ordination.png',
       width = 10 ,height = 10, units = "in", dpi = 600, 
       bg = "transparent")

#Nearest Neighbor
#Case
EL.nn <- read_rds("Processed Data/EL Nearest Neighbor.rds") 
AP.nn <- read_rds("Processed Data/AP Nearest Neighbor.rds") 
CC.nn <- read_rds("Processed Data/CC Nearest Neighbor.rds") 

emerald.nearest <- ggplot(EL.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.1071574642, 0.24050633), label = "Deschampsia cespitosa", color = "grey22", nudge_y = 0.015, size = 3.5) +
  geom_text(aes(0.1251533742, 0.18987342), label = "Fragaria virginiana", color = "grey22", nudge_y = -0.01, size = 3.5) +
  geom_text(aes(0.1480333652, 0.20430108), label = "Fragaria virginiana", color = "grey22", nudge_y = 0.01, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
emerald.nearest

avery.nearest <- ggplot(AP.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.1531226486, 0.29411765), label = "Poa pratensis", color = "grey22", nudge_y = 0.015, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
avery.nearest

copper.nearest <- ggplot(CC.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0778150748, 0.13333333), label = "Thalictrum fendleri", color = "grey22", nudge_y = 0.01, size = 3.5) +
  geom_text(aes(0.0234589564, 0.06666667), label = "Mertensia brevistyla", color = "grey22", nudge_y = 0.01, size = 3.5) +
  theme_minimal() +
  ggtitle("Castilleja septentrionalis") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
copper.nearest

#Cali
DC.nn <- read_rds("Processed Data/DC Nearest Neighbor.rds") 
JH.nn <- read_rds("Processed Data/JH Nearest Neighbor.rds") 

deercreek.nearest <- ggplot(DC.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.0746228926, 0.155038760), label = "Eremogone congesta", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0143414410, 0.121212121), label = "Bromus inermis", color = "grey22",nudge_x = 0.007, nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0402707664, 0.111111111), label = "Carex sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0383203304, 0.090909091), label = "Achnatherum sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0463176575, 0.085271318), label = "Carex sp.", color = "grey22", nudge_y = -0.005, size = 3.5) +
  theme_minimal() +
  ggtitle("Castilleja linariifolia") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
deercreek.nearest

johnson.nearest <- ggplot(JH.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.0738952297, 0.14492754), label = "Lathyrus lanszwertii", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0496049166, 0.10144928), label = "Viola praemorsa", color = "grey22", nudge_y = 0.007, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja linariifolia") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
johnson.nearest

#Cach
AL.nn <- read_rds("Processed Data/AL Nearest Neighbor.rds") 

almont.nearest <- ggplot(AL.nn, aes(x = rel_abund_cover, y = NN_freq)) +
  geom_smooth(method=lm , color="darkred", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#e07a5f", "#3d405b")) +
  scale_shape_manual(values = c(20, 18)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  geom_text(aes(0.0574450029, 0.16417910), label = "Stipeae", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0168418304, 0.08955224), label = "Crepis sp.", color = "grey22", nudge_y = 0.007, size = 3.5) +
  geom_text(aes(0.0577061166, 0.01492537), label = "Artemisia arbuscula", color = "grey22", nudge_y = -0.006, size = 3.5) +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("Castilleja chromosa") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
almont.nearest

#complete
nearestplots <- ggarrange(emerald.nearest, avery.nearest, copper.nearest, deercreek.nearest, johnson.nearest, almont.nearest,
                          labels = c("A", "B", "C", "D", "E", "F"), 
                          nrow = 2, ncol = 3, common.legend = TRUE, legend = "top")
nearestplots

ggsave(plot = nearestplots, filename = 'Figures and Tables/Nearest Neighbor.png',
       width = 15 ,height = 10, units = "in", dpi = 600, 
       bg = "transparent")
