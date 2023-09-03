library(metagear)
library(tidyverse)

#set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Parasitic Plant Review/Host")

#Web of science search with just hemiparsaite+conservation
reading_list <- read.csv("hemiparasitehostpreference.csv")
reading_list <- reading_list %>% rename("ABSTRACT" = "Abstract",
                                        "TITLE" = "Article.Title")

effort_distribute(reading_list, initialize = TRUE, reviewers = "jordan", save_split = TRUE)
# initialize screener GUI
abstract_screener("effort_jordan.csv",
                  aReviewer = "jordan",
                  theButtons = c("YES","Maybe","NO","Background"),
                  keyBindingToButtons = c("y","m","n","b"))
