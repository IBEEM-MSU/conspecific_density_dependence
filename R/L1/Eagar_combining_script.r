# Joining species slope coefficients and other data for one master dataframe of biotic explanatory variables and responses.

# A. Eagar 5-2025

library(tidyverse)
					
AB_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/AB_slope_modeling_data.csv") %>%
					mutate(Sp = "ab") %>%
					select(-c(4, 6:14))

# AE_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/AE_slope_modeling_data.csv") %>%
					# mutate(Sp = "ae") %>%
					# select(-c(4, 6:14))

BC_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BC_slope_modeling_data.csv")  %>%
					mutate(Sp = "bc") %>%
					select(-c(4, 6:14))

BONRO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BONRO_slope_modeling_data.csv")  %>%
					mutate(Sp = "bonro") %>%
					select(-c(4, 6:14))

BS_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BS_slope_modeling_data.csv")  %>%
					mutate(Sp = "bs") %>%
					select(-c(4, 6:14))

IW_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/IW_slope_modeling_data.csv")  %>%
					mutate(Sp = "iw") %>%
					select(-c(4, 6:14))

NRO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/NRO_slope_modeling_data.csv")  %>%
					mutate(Sp = "nro") %>%
					select(-c(4, 6:14))

RM_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/RM_slope_modeling_data.csv")  %>%
					mutate(Sp = "rm") %>%
					select(-c(4, 6:14))

SM_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/SM_slope_modeling_data.csv")  %>%
					mutate(Sp = "sm") %>%
					select(-c(4, 6:14))

WA_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WA_slope_modeling_data.csv")  %>%
					mutate(Sp = "wa") %>%
					select(-c(4, 6:14))

WO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WO_slope_modeling_data.csv")  %>%
					mutate(Sp = "wo") %>%
					select(-c(4, 6:14))

# WP_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WP_slope_modeling_data.csv")  %>%
					# mutate(Sp = "wp") %>%
					# select(-c(4, 6:14))

Combined_slope_data <- rbind(AB_slope_modeling_data, BC_slope_modeling_data, BONRO_slope_modeling_data, BS_slope_modeling_data, IW_slope_modeling_data, NRO_slope_modeling_data, RM_slope_modeling_data, SM_slope_modeling_data, WA_slope_modeling_data, WO_slope_modeling_data) %>%
					mutate(Joiner = paste0(SiteID, quadrat, Cohort.Year, Sp)) %>%
					select(c(4,5,7)) %>%
					unique()

Demog_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/allSp_surv_bioticVars.csv") %>%
					mutate(Joiner = paste0(SiteID, quadrat, First.Year, Sp)) %>%
					unique() %>%
					left_join(Combined_slope_data, by = "Joiner")