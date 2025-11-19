#Walder script

setwd("/Users/jargrett/Downloads")
#Packages
library(tidyverse)#for data wrangling and restructuring
library(flextable)#table creation and formatting
library(conflicted)#settle package conflicts
library(gtsummary)#to make summary tables for models

#Create theme for All Tables
my_theme <- function(x, ...) {
  x <- colformat_double(x, big.mark = "'", decimal.mark = ".", digits = 3)
  x <- set_table_properties(x, layout = "fixed")
  x <- border_remove(x)
  thick <- fp_border_default(width = 1.5)
  x <- font(x, part = "all", fontname = "Times New Roman")
  x <- bold(x, part = "header")
  x <- align(x, part = "all", align = "center")
  #x <- border_inner_h(x, part="all", border = thin)
  x <- hline_top(x, part = "body", border = thick)
  x <- hline_bottom(x, part = "body", border = thick)
  x <- vline(x, j = 1, border = thick)
  autofit(x)
}
thin <- fp_border_default(width = 0.5)
thick <- fp_border_default(width = 1.5)


partnership <- read.csv("Collaborators and Partnerships - Final.csv",check.names = FALSE)
partners <- flextable(partnership) %>% 
  merge_at(i = 1:2, j = 1) %>% 
  merge_at(i = 3:5, j = 1) %>% 
  merge_at(i = 6:9, j = 1) %>% 
  my_theme(partners) 
partners <- width(partners, width = 6, unit = "cm")
partners <- hline(partners, i = 2, border = thin)
partners <- hline(partners, i = 5, border = thin)
partners <- hline(partners, i = 9, border = thin)
partners <- hline(partners, i = 10, border = thin)
partners

save_as_image(x = partners, path = "Biota Partnerships.png", res = 600)
