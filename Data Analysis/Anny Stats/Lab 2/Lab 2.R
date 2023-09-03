
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 2")
surveys <- read_csv("portal_data_joined.csv.csv")

#inspects data
str(surveys)

# preview the data
view(surveys)

# selects specific parts of dataframe to keep
select(surveys, plot_id, species_id, weight)

# - next to record_id makes it so that select all but specific columns that are marked
select(surveys, -record_id, -species_id)

# filters by specific criteria
filter(surveys, year == 1995)

#nested functions are possible and done in this syntax (elavuated from inside out)
surveys_sml <- select(filter(surveys, weight < 5), species_id, sex, weight)

#pipped form... same command order as aboce though (read it as the word "then")
surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)

#mutating creates new colomns based on values in existing
surveys %>% 
  mutate (weight_kg = weight / 1000)

# removes first few NA rows. 
surveys %>%
  filter(!is.na(weight)) %>%
  mutate(weight_kg = weight / 1000) %>%
  head()

#Challenge
surveysC <- surveys %>% 
  filter(!is.na(hindfoot_length)) %>% 
  mutate(hindfoot_cm = hindfoot_length / 10) %>% 
  filter (hindfoot_cm < 3 ) %>% 
  select (species_id, hindfoot_cm)

surveyG <- surveys %>% 
  filter(!is.na(sex)) %>%
  group_by(sex) %>% 
  summarize(mean_weight = mean(weight, na.rm = T))


surveyD <- surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight)) %>%
  print(n = 15)


# To export:

survey_complete <- surveys %>% 
  filter(!is.na(weight), !is.na(hindfoot_length), !is.na(sex))

species_count <- survey_complete %>% 
  count(species_id) %>% 
  filter (n >= 50)

survey_complete <- survey_complete %>% 
  filter (species_id %in% species_count$species_id)

dim(survey_complete)

write.csv(survey_complete, file = "survey_complete.csv")
