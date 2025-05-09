
## Load in libraries for processing and visualizing data
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

#### Read in mature tree data
adults = read_xlsx("./data/L0/MNF_matureTree_allSites_2024.xlsx")

# List of basal area columns
ba_columns <- c("ba_fall99", "ba_aug05", "ba_aug07", "ba_oct09", 
                "ba_aug11", "ba_sep13", "ba_fall16", "ba_fall19", "ba_22")

# Summarize total BA per species per census year
species_summary <- adults %>%
  group_by(Sp) %>%
  summarise(
    across(all_of(ba_columns), ~ sum(as.numeric(.x), na.rm = TRUE), 
           .names = "total_{.col}")  # Sum total BA
  )

# Calculate average BA per species across all 9 census years
species_summary$ba_avg <- rowMeans(species_summary[, 2:10], na.rm = TRUE)

# Add column that is the sum of the values
species_summary$total_ba_avg = rep(sum(species_summary$ba_avg), 
                                   nrow(species_summary))

# Calculate the proportion of total BA
species_summary$sp_perc_BA = 
      (species_summary$ba_avg/species_summary$total_ba_avg)*100

# Round values
species_summary$sp_perc_BA = round(species_summary$sp_perc_BA, 2)

View(species_summary)

sp_summary_percBA = species_summary[,c(1,13)]


