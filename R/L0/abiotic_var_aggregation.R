
#################################################
####### Derivation of environmental covariates needed for abiotic models of survivorship
####### Author: Bailey McNichol 
####### Date: 20 May 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(dplyr)
library(ggplot2)
library(lubridate)

##### Read in csv files with currently derived environmental covariates
## Light, VWC, baseCEC
envVars <- read.csv("./data/L0/quad_envVars_allSites.csv", row.names = 1)

# Subset environmental variables to only include the relevant columns
envVars_mod = envVars[,c(1:3,5:6,13)]
# Multiple percOpen by 100 so it's actually a percentage
envVars_mod$percOpen = envVars_mod$percOpen*100
envVars_mod$Site = as.character(envVars_mod$Site)

## Daily climate data - for temperature
clim_site1 <- read.csv("./data/L0/clim_49.01_1997-2023.csv")
clim_site2 <- read.csv("./data/L0/clim_59.12_1997-2023.csv")
clim_site3 <- read.csv("./data/L0/clim_8.21_1997-2023.csv")
clim_site4 <- read.csv("./data/L0/clim_39.21_1997-2023.csv")
clim_site5 <- read.csv("./data/L0/clim_16.35_1997-2023.csv")
clim_site6 <- read.csv("./data/L0/clim_60.35_1997-2023.csv")
clim_site7 <- read.csv("./data/L0/clim_7.40_1997-2023.csv")
clim_site8 <- read.csv("./data/L0/clim_40.40_1997-2023.csv")
clim_site9 <- read.csv("./data/L0/clim_21.43_1997-2023.csv")
clim_site10 <- read.csv("./data/L0/clim_25.43_1997-2023.csv")
clim_site11 <- read.csv("./data/L0/clim_6.45_1997-2023.csv")
clim_site12 <- read.csv("./data/L0/clim_24.45_1997-2023.csv")
## Monthly growing season climate data - for PDSI
pdsi_site1 <- read.csv("./data/L0/monthlyPDSI_49.01_1997-2023.csv")
pdsi_site2 <- read.csv("./data/L0/monthlyPDSI_59.12_1997-2023.csv")
pdsi_site3 <- read.csv("./data/L0/monthlyPDSI_8.21_1997-2023.csv")
pdsi_site4 <- read.csv("./data/L0/monthlyPDSI_39.21_1997-2023.csv")
pdsi_site5 <- read.csv("./data/L0/monthlyPDSI_16.35_1997-2023.csv")
pdsi_site6 <- read.csv("./data/L0/monthlyPDSI_60.35_1997-2023.csv")
pdsi_site7 <- read.csv("./data/L0/monthlyPDSI_7.40_1997-2023.csv")
pdsi_site8 <- read.csv("./data/L0/monthlyPDSI_40.40_1997-2023.csv")
pdsi_site9 <- read.csv("./data/L0/monthlyPDSI_21.43_1997-2023.csv")
pdsi_site10 <- read.csv("./data/L0/monthlyPDSI_25.43_1997-2023.csv")
pdsi_site11 <- read.csv("./data/L0/monthlyPDSI_6.45_1997-2023.csv")
pdsi_site12 <- read.csv("./data/L0/monthlyPDSI_24.45_1997-2023.csv")


################ Climate variable summarizing for each site:
#### Site 1
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site1 <- clim_site1 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site1 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
    max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("49.01"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site1 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
    survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site1 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site1) = c("TransectID", "Site", "clim_year",
            "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 2
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site2 <- clim_site2 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site2 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("59.12"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site2 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site2 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site2) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 3
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site3 <- clim_site3 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site3 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("8.21"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site3 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site3 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site3) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 4
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site4 <- clim_site4 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site4 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("39.21"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site4 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site4 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site4) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 5
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site5 <- clim_site5 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site5 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("16.35"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site5 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site5 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site5) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 6
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site6 <- clim_site6 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site6 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("60.35"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site6 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site6 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site6) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 7
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site7 <- clim_site7 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site7 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("7.40"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site7 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site7 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site7) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 8
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site8 <- clim_site8 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site8 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("40.40"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site8 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site8 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site8) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 9
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site9 <- clim_site9 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site9 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("21.43"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site9 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site9 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site9) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 10
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site10 <- clim_site10 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site10 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("25.43"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site10 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site10 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site10) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 11
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site11 <- clim_site11 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site11 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("6.45"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site11 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site11 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site11) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

#### Site 12
# 1) find min/max temp between June 1st of recruit year and survival year, and average 
# Create a 'clim_year' column that starts in June and ends in May next year
clim_site12 <- clim_site12 %>%
  mutate(clim_year = if_else(month(date) >= 6,
                             paste0(year, "-", year + 1),
                             paste0(year - 1, "-", year)))
# Summarize by climate year
clim_summ <- clim_site12 %>%
  group_by(clim_year) %>%
  summarise(min_tmin_c = min(tmin_c, na.rm = TRUE),
            max_tmax_c = max(tmax_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(survYear = as.integer(sub(".*-", "", clim_year)))  %>%
  filter(., survYear != 1997) %>%
  mutate(Site = rep("24.45"))
clim_summ$Site = as.character(clim_summ$Site)
# 2) find mean PDSI between May-September of recruit year and survival year, and average
pdsi_summ <- pdsi_site12 %>%
  mutate(survYear1 = year,       # Year affects survival in same year
         survYear2 = year + 1 ) %>%   # And also affects survival in next year
  # Reshape to long format so each row contributes to two survYears
  tidyr::pivot_longer(cols = c(survYear1, survYear2), 
                      names_to = "source", values_to = "survYear") %>%
  group_by(Site, Transect_id, survYear) %>%
  summarise(mean_pdsi = mean(pdsi, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(clim_year = paste0(survYear - 1, "-", survYear)) %>%
  filter(., survYear != 1997) 
pdsi_summ$Site = as.character(pdsi_summ$Site)
# 3) Merge via a left join
clim_merge <- pdsi_summ %>%
  left_join(clim_summ, by = c("Site", "survYear", "clim_year"))
# 4) Reorder and rename columns
clim_site12 = clim_merge[,c(2,1,5,3,4,6,7)]
colnames(clim_site12) = c("TransectID", "Site", "clim_year",
                         "survYear", "mean_pdsi", "tmin_c","tmax_c")

####### Join together climate variables across all sites
clim_allSites = rbind(clim_site1,clim_site2,clim_site3,
                      clim_site4,clim_site5,clim_site6,
                      clim_site7,clim_site8,clim_site9,
                      clim_site10,clim_site11,clim_site12)

# Merge with environmental data via a left join
clim_env_allSites <- clim_allSites %>%
  left_join(envVars_mod, by = c("TransectID","Site"))
clim_env_allSites = as.data.frame(clim_env_allSites)

#### Write csv to data/L1 folder
write.csv(clim_env_allSites, "./data/L1/clim_envVars_allSites.csv", row.names = FALSE)
