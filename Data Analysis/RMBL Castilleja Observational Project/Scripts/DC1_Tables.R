setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")
#Packages
library(tidyverse)#for data wrangling and restructuring
library(flextable)#table creation and formatting
library(conflicted)#settle package conflicts
library(gtsummary)#to make summary tables for models

conflicts_prefer(flextable::font)
conflicts_prefer(plyr::mutate)

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

#Create Site Data Table
site.main <- read.csv("Processed Data/Castilleja Observational Project Tables - Site.csv",check.names = FALSE)
site <- flextable(site.main) %>% 
  italic(j = "Castilleja species", part = "body") %>% 
  merge_at(i = 1:3, j = 1) %>% 
  merge_at(i = 4:6, j = 1) %>% 
  my_theme(site)
site <- hline(site, i = 3, border = thin)
site <- hline(site, i = 6, border = thin)
# site.title <- as_paragraph(
#   as_chunk("Table I. ", props = fp_text_default(bold = TRUE, font.family = "Times New Roman")),
#   as_chunk("Site and ", props = fp_text_default(font.family = "Times New Roman")),
#   as_chunk("Castilleja", props = fp_text_default(italic = TRUE, font.family = "Times New Roman")),
#   as_chunk(" species for observational paired-plots.", props = fp_text_default(font.family = "Times New Roman")))
# site <- set_caption(site, site.title)
site

#Indicator Data Table
ind.main <-  read.csv("Processed Data/Castilleja Observational Project Tables - Indicator real.csv",check.names = FALSE)
ind.main <- ind.main %>%
  mutate(
    pval_numeric = `P-value`,  # keep numeric copy for conditions
    `P-value` = case_when(
      pval_numeric < 0.001 ~ paste0(round(pval_numeric, 3), "***"),
      pval_numeric < 0.01  ~ paste0(round(pval_numeric, 3), "**"),
      pval_numeric < 0.05  ~ paste0(round(pval_numeric, 3), "*"),
      pval_numeric < 0.1   ~ paste0(round(pval_numeric, 3), "."),
      TRUE                 ~ as.character(round(pval_numeric, 3))
    )
  )

ind <- flextable(ind.main) %>%  
  italic(j = "Castilleja species", part = "body") %>% 
  italic(j = "Indicator species", part = "body") %>% 
  my_theme(site)

ind <- bold(ind, i = ~ pval_numeric < 0.05, j = "Indicator species")  
ind <- bold(ind, i = ~ pval_numeric < 0.05, j = "P-value")

ind <- delete_columns(ind, j = "pval_numeric")

ind <- merge_at(ind, i = 1:6, j = 1) %>% 
  merge_at(i = 7:14, j = 1) %>%
  merge_at(i = 15:16, j = 1) %>% 
  merge_at(i = 2:4, j = 2) %>% 
  merge_at(i = 5:6, j = 2) %>% 
  merge_at(i = 7:8, j = 2) %>%
  merge_at(i = 9:13, j = 2) %>% 
  merge_at(i = 15:16, j = 2) %>% 
  hline(i = 6, border = thick) %>%
  hline(i = 14, border = thick) %>%
  hline(i = 1, border = thin) %>%
  hline(i = 4, border = thin) %>%
  hline(i = 8, border = thin) %>%
  hline(i = 13, border = thin)
ind

#Host Data Table
host.sup <-  read.csv("Processed Data/Castilleja Observational Project Tables - Host.csv",check.names = FALSE)
host <- flextable(host.sup) %>% 
  italic(j = "Castilleja species", part = "body") %>% 
  italic(j = "Excavated species", part = "body") %>% 
  merge_at(i = 1:7, j = 1) %>% 
  merge_at(i = 8:11, j = 1) %>% 
  merge_at(i = 12:24, j = 1) %>% 
  merge_at(i = 25:29, j = 1) %>% 
  merge_at(i = 1:7, j = 2) %>% 
  merge_at(i = 8:9, j = 2) %>% 
  merge_at(i = 12:24, j = 2) %>% 
  merge_at(i = 28:29, j = 2) %>% 
  my_theme(host)
host <- hline(host, i = 7, border = thin)
host <- hline(host, i = 11, border = thin)
host <- hline(host, i = 24, border = thin)
host <- width(host, width = 2.5)
host
#Environmental Data Table
envi.sup <-  read.csv("Processed Data/Castilleja Observational Project Tables - Environmental.csv",check.names = FALSE)
envi.sup$`Sampling year` <- as.factor(envi.sup$`Sampling year`)

my_theme2 <- function(x, ...) {
  x <- colformat_double(x, big.mark = "'", decimal.mark = ".", digits = 2)
  x <- set_table_properties(x, layout = "fixed")
  x <- border_remove(x)
  thick <- fp_border_default(width = 1.5)
  x <- font(x, part = "all", fontname = "Times New Roman")
  x <- bold(x, part = "header")
  x <- align(x, part = "all", align = "center")
  #x <- border_inner_h(x, part="all", border = thin)
  x <- hline_top(x, part = "body", border = thick)
  x <- hline_bottom(x, part = "body", border = thick)
  x <- vline(x, j = 3, border = thick)
  x <- vline(x, j = 11, border = thick)
  autofit(x)
}
envi <- flextable(envi.sup) %>%
  separate_header() %>% 
  my_theme2(envi)
envi <- vline(envi, j = 5, part = "header", border = thick)
envi <- vline(envi, j = 7, part = "header", border = thick)
envi <- vline(envi, j = 9, part = "header", border = thick)
envi

#Diversity stats output Table

#Exporting Tables

#Main tables
save_as_image(x = site, path = "Figures and Tables/Site Table.png", res = 600)
save_as_image(x = ind, path = "Figures and Tables/Indicator Species Table.png", res = 600)
#Supplement
save_as_image(x = envi, path = "Figures and Tables/Environmental Table.png", res = 600)
save_as_image(x = host, path = "Figures and Tables/Excavation Specimen Table.png", width = 8, height = 11, res = 600)
save_as_image(x = diversity.model, path = "Figures and Tables/Diversity LMM Table.png", res = 600)
