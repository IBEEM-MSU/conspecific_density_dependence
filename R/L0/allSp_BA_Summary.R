
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

# Summarize total BA per species/site per census year
speciesSite_summary <- adults %>%
  group_by(Site,Sp) %>%
  summarise(
    across(all_of(ba_columns), ~ sum(as.numeric(.x), na.rm = TRUE), 
           .names = "total_{.col}")  # Sum total BA
  )
# Calculate average BA per species/site across all 9 census years
speciesSite_summary$spBA_avg <- rowMeans(speciesSite_summary[, 3:11], na.rm = TRUE)

# Subset out to just include relevant columns for exploratory plots
speciesSite_sum <- speciesSite_summary[,c(1,2,12)]
speciesSite_sum = as.data.frame(speciesSite_sum)

# Step 1: Create a complete grid of all Site × Sp combinations
site_sp_grid <- expand.grid(
  Site = unique(speciesSite_sum$Site),
  Sp = unique(speciesSite_sum$Sp),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# Step 2: Left join the original data onto this grid
speciesSite_full <- site_sp_grid %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) %>%
  mutate(spBA_avg = ifelse(is.na(spBA_avg), 0, spBA_avg))  # fill missing with 0

# Add 1 to each BA sum to account for missing data (/0 errors)
speciesSite_full$spBA_corr = speciesSite_full$spBA_avg + 1

# Write csv
write.csv(speciesSite_full, "./data/L1/mean_totalSpSiteBA.csv", row.names = FALSE)

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


