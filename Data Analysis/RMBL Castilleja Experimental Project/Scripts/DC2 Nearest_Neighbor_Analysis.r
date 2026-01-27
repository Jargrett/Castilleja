setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
library(car)
library(broom)


conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)
#-----Assessing Nearest Neighbor data---------#
#Importing Species Data
cover <- read.csv("Processed Data/Neighbor Cover.csv")

#Sum for each species the number of times they apeard as a NN
nearest <- cover %>%
  group_by(code, year) %>%
  summarise(
    total_cover = sum(percent_cover, na.rm = TRUE),
    nn_count    = sum(nearest_neighbor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  mutate(
    rel_abund_cover = total_cover / sum(total_cover, na.rm = TRUE),
    nn_freq    = nn_count / sum(nn_count, na.rm = TRUE)
  ) %>%
  ungroup()


nn_yearly <- nearest %>%
  group_by(year) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(nn_freq ~ rel_abund_cover, data = .x)),
    results = map2(model, data, ~ {
      bind_cols(
        .y,
        as.data.frame(
          predict(.x, newdata = .y, interval = "prediction", level = 0.95)
        ) 
      ) %>%
        mutate(
          NN_association = case_when(
            nn_freq > upr ~ "Frequent Neighbor",
            nn_freq < lwr ~ "Infrequent Neighbor",
            TRUE           ~ "As expected"
          )
        )
    })
  ) %>%
  select(year, results) %>%
  unnest(results)

split_and_name <- function(df, column) {
  # Get the dataframe name as a string
  df_name <- deparse(substitute(df))
  # Ensure column exists
  if (!column %in% colnames(df)) stop("Column not found in dataframe")
  # Split the dataframe by the column
  df_list <- split(df, df[[column]])
  # Make valid R names for safety
  names(df_list) <- paste0(df_name, "_", make.names(names(df_list)))
  # Assign each dataframe into the global environment
  list2env(df_list, envir = .GlobalEnv)
  # Return the list invisibly (optional)
  invisible(df_list)
}

split_and_name(nn_yearly, "year")

pre.nearest <- ggplot(nn_yearly_Pre, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  geom_text(aes(0.101, 0.188811189), label = "ELGL", color = "grey22", nudge_y = - 0.007, size = 3) +
  geom_text(aes(0.107, 0.020979021), label = "LIPO", color = "grey22", nudge_y = - 0.007, size = 3) +
  theme_minimal() +
  ggtitle("2023 Pre") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
pre.nearest

y1.nearest <- ggplot(nn_yearly_X2023, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  geom_text(aes(0.0985, 0.24590164), label = "ELGL", color = "grey22", nudge_y = - 0.0075, size = 3) +
  geom_text(aes(0.0475, 0.16393443), label = "LUAR", color = "grey22", nudge_y = - 0.0075, size = 3) +
  ggtitle("2023 Post") +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
y1.nearest

y2.nearest <- ggplot(nn_yearly_X2024, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
  theme_minimal() +
  ggtitle("2024") +
  theme(legend.position = "none") +
  geom_text(aes(0.0435, 0.13253012), label = "ELGL", color = "grey22", nudge_y = - 0.005, size = 3) +
  geom_text(aes(0.058, 0.09638554), label = "MESP", color = "grey22", nudge_y = 0.005, size = 3) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
y2.nearest

y3.nearest <- ggplot(nn_yearly_X2025, aes(x = rel_abund_cover, y = nn_freq)) +
  geom_smooth(method=lm , color="#582f0e", fill = "cornsilk3", se=TRUE) + 
  geom_point(aes(color = year, shape = year), size = 2.3) +
  scale_color_manual( values=c("#936639")) +
  scale_shape_manual(values = c(20)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = lwr), linetype = "dashed", col = "black") +
  geom_line(aes(y = upr), linetype = "dashed", col = "black") +
  labs(x = "Relative Abundance", y = "Nearest Neighbor Frequency") +
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
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) + #transparent legend pane
  theme(plot.title = element_text(hjust = 0.92, vjust= -0.12))
y3.nearest


nearestplots <- ggarrange(pre.nearest, y1.nearest, y2.nearest, y3.nearest,
                               labels = c("A", "B","C", "D"), 
                               nrow = 2, ncol = 2)

nearestplots 
