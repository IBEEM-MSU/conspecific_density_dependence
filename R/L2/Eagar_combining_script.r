# Joining species slope coefficients and other data for one master dataframe of biotic explanatory variables and responses.

# A. Eagar 5-2025

library(tidyverse)

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")

Site_BA_quads_species_BA <- Site_BA %>%
					mutate(quadrat = as.numeric(X) + 0.5, ba_joiner = paste(TransectID,".",Year)) %>%
					select(-c(X.1, X, Y, Midpoint_X, Midpoint_Y)) %>%
					replace(is.na(.), 0) %>%
					select(c(2, 40:57, 59)) %>%
					gather(c(2:19), key = "Species_BA", value = "BA_15m") %>%
					group_by(ba_joiner, Species_BA) %>%
					mutate(Total_BA = sum(BA_15m), Myc_type = if("ae_15m" %in% Species_BA){"AM"} else{if("bc_15m" %in% Species_BA){"AM"} else{if("rm_15m" %in% Species_BA){"AM"} else{if("sm_15m" %in% Species_BA){"AM"} else{if("svb_15m" %in% Species_BA){"AM"} else{if("wa_15m" %in% Species_BA){"AM"} else{"ECM"}}}}}}) %>%
					ungroup() %>%
					group_by(ba_joiner, Myc_type) %>%
					mutate(Myc_type_BA = sum(Total_BA)) %>%
					ungroup() %>%
					select(ba_joiner, Species_BA, Total_BA) %>%
					unique() %>%
					spread(key = Species_BA, value = Total_BA) %>%
					unique()
	
	head(Site_BA_quads_species_BA)
	
Site_BA_quads_myc_BA <- Site_BA %>%
					mutate(quadrat = as.numeric(X) + 0.5, ba_joiner = paste(TransectID,".",Year)) %>%
					select(-c(X.1, X, Y, Midpoint_X, Midpoint_Y)) %>%
					replace(is.na(.), 0) %>%
					select(c(2, 40:57, 59)) %>%
					gather(c(2:19), key = "Species_BA", value = "BA_15m") %>%
					group_by(ba_joiner, Species_BA) %>%
					mutate(Total_BA = sum(BA_15m), Myc_type = if("ae_15m" %in% Species_BA){"AM"} else{if("bc_15m" %in% Species_BA){"AM"} else{if("rm_15m" %in% Species_BA){"AM"} else{if("sm_15m" %in% Species_BA){"AM"} else{if("svb_15m" %in% Species_BA){"AM"} else{if("wa_15m" %in% Species_BA){"AM"} else{"ECM"}}}}}}) %>%
					ungroup() %>%
					group_by(ba_joiner, Myc_type) %>%
					mutate(Myc_type_BA = sum(Total_BA)) %>%
					select(ba_joiner, Myc_type, Myc_type_BA) %>%
					unique() %>%
					spread(key = Myc_type, value = Myc_type_BA) %>%
					unique()

	head(Site_BA_quads_myc_BA)
	
AB_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/ab_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(AB_prior_densities)
	
AB_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/AB_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "ab") %>%
					unique()

	head(AB_slope_modeling_data)
	
AB_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/ab_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(AB_counts)

Joined_AB_modeling_data <- AB_slope_modeling_data %>%
					left_join(AB_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(AB_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9))

	head(Joined_AB_modeling_data)
	
# AE_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/AE_slope_modeling_data.csv") %>%
					# mutate(Sp = "ae") %>%
					# select(-c(4, 6:14))

BC_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bc_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(BC_prior_densities)
	
BC_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BC_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "bc") %>%
					unique()

	head(BC_slope_modeling_data)
	
BC_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bc_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(BC_counts)

Joined_BC_modeling_data <- BC_slope_modeling_data %>%
					left_join(BC_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(BC_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9))
					
	head(Joined_BC_modeling_data)		

BONRO_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bonro_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(BONRO_prior_densities)
	
BONRO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BONRO_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "bonro") %>%
					unique()

	head(BONRO_slope_modeling_data)
	
BONRO_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bonro_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(BONRO_counts)

Joined_BONRO_modeling_data <- BONRO_slope_modeling_data %>%
					left_join(BONRO_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(BONRO_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9))

	head(Joined_BONRO_modeling_data)

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

Combined_demography_modeling_data <- rbind(Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data, Joined_AB_modeling_data)

Demog_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/allSp_surv_bioticVars.csv") %>%
					mutate(Joiner = paste0(SiteID, quadrat, First.Year, Sp)) %>%
					unique() %>%
					left_join(Combined_slope_data, by = "Joiner")

Trb <- Combined_slope_data[!(Combined_slope_data$Joiner %in% Demog_data$Joiner), ]
					
	write.csv(Combined_slope_data, "C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/Combined.slope.data.csv", row.names = FALSE)