library(metagear)
library(tidyverse)

#set working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Parasitic Plant Review/Castilleja")

#Web of science search with just Castilleja
reading_list <- read.csv("Castillejareview.csv")
reading_list <- reading_list %>% rename("ABSTRACT" = "Abstract",
                                        "TITLE" = "Article.Title")

effort_distribute(reading_list, initialize = TRUE, reviewers = "jordan", save_split = TRUE)
# initialize screener GUI
abstract_screener("effort_jordan.csv",
                  aReviewer = "jordan",
                  theButtons = c("YES","Maybe","NO","Background"),
                  keyBindingToButtons = c("y","m","n","b"))
