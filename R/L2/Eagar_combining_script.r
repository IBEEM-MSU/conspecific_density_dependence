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
					mutate(Total_BA = sum(BA_15m), Myc_type = if("ae_15m" %in% Species_BA){"AM_15m"} else{if("bc_15m" %in% Species_BA){"AM_15m"} else{if("rm_15m" %in% Species_BA){"AM_15m"} else{if("sm_15m" %in% Species_BA){"AM_15m"} else{if("svb_15m" %in% Species_BA){"AM_15m"} else{if("wa_15m" %in% Species_BA){"AM_15m"} else{"ECM_15m"}}}}}}) %>%
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
					mutate(Total_BA = sum(BA_15m), Myc_type = if("ae_15m" %in% Species_BA){"AM_15m"} else{if("bc_15m" %in% Species_BA){"AM_15m"} else{if("rm_15m" %in% Species_BA){"AM_15m"} else{if("sm_15m" %in% Species_BA){"AM_15m"} else{if("svb_15m" %in% Species_BA){"AM_15m"} else{if("wa_15m" %in% Species_BA){"AM_15m"} else{"ECM_15m"}}}}}}) %>%
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
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

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
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "AM")
					
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
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

	head(Joined_BONRO_modeling_data)

BS_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bs_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(BS_prior_densities)
	
BS_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/BS_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "bs") %>%
					unique()

	head(BS_slope_modeling_data)
	
BS_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/bs_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(BS_counts)

Joined_BS_modeling_data <- BS_slope_modeling_data %>%
					left_join(BS_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(BS_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

	head(Joined_BS_modeling_data)

IW_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/iw_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(IW_prior_densities)
	
IW_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/IW_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "iw") %>%
					unique()

	head(IW_slope_modeling_data)
	
IW_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/iw_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(IW_counts)

Joined_IW_modeling_data <- IW_slope_modeling_data %>%
					left_join(IW_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(IW_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

	head(Joined_IW_modeling_data)

NRO_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/nro_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(NRO_prior_densities)
	
NRO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/NRO_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "nro") %>%
					unique()

	head(NRO_slope_modeling_data)
	
NRO_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/nro_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(NRO_counts)

Joined_NRO_modeling_data <- NRO_slope_modeling_data %>%
					left_join(NRO_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(NRO_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

	head(Joined_NRO_modeling_data)

RM_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(RM_prior_densities)
	
RM_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/RM_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "rm") %>%
					unique()

	head(RM_slope_modeling_data)
	
RM_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(RM_counts)

Joined_RM_modeling_data <- RM_slope_modeling_data %>%
					left_join(RM_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(RM_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "AM")

	head(Joined_RM_modeling_data)

SM_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/sm_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(SM_prior_densities)
	
SM_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/SM_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "sm") %>%
					unique()

	head(SM_slope_modeling_data)
	
SM_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/sm_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(SM_counts)

Joined_SM_modeling_data <- SM_slope_modeling_data %>%
					left_join(SM_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(SM_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "AM")

	head(Joined_SM_modeling_data)

Debug <- Joined_SM_modeling_data %>%
					filter(SiteID == 7 & quadrat == 25)
					
WA_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/wa_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(WA_prior_densities)
	
WA_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WA_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "wa") %>%
					unique()

	head(WA_slope_modeling_data)
	
WA_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/wa_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(WA_counts)

Joined_WA_modeling_data <- WA_slope_modeling_data %>%
					left_join(WA_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(WA_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "AM")

	head(Joined_WA_modeling_data)

WO_prior_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/wo_prior_cohort_densities.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					gather(c(4:6), key = "Prior_count_categories", value = "Prior_Counts") %>%
					group_by(Count_joiner) %>%
					mutate(Avg_3_year_prior_counts = mean(Prior_Counts)) %>%
					spread(key = Prior_count_categories, value = Prior_Counts) %>%
					select(4:8) %>%
					unique()
	
	head(WO_prior_densities)
	
WO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WO_slope_modeling_data.csv") %>%
					select(SiteID, quadrat, Cohort.Year, ba_joiner, GLM.Slope.Estimate) %>%
					unique() %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  Cohort.Year)) %>%
					mutate(Sp = "wo") %>%
					unique()

	head(WO_slope_modeling_data)
	
WO_counts <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/wo_count_demography.csv") %>%
					mutate(Count_joiner = paste0(SiteID, ".", quadrat, ".",  First.Year)) %>%
					select(c(5,6,7)) %>%
					unique()
	head(WO_counts)

Joined_WO_modeling_data <- WO_slope_modeling_data %>%
					left_join(WO_counts, by = "Count_joiner") %>%
					spread(key = Time.Since.First.Year, value = Count) %>%
					left_join(WO_prior_densities, by = "Count_joiner") %>%
					left_join(Site_BA_quads_myc_BA, by = "ba_joiner") %>%
					left_join(Site_BA_quads_species_BA, by = "ba_joiner") %>%
					select(-c(4,6)) %>%
					rename("Year_0_counts" = c(6), "Year_1_counts" = c(7), "Year_2_counts" = c(8), "Year_3_counts" = c(9)) %>%
					mutate(Seedling_myc = "ECM")

	head(Joined_WO_modeling_data)

# WP_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/WP_slope_modeling_data.csv")  %>%
					# mutate(Sp = "wp") %>%
					# select(-c(4, 6:14))

Combined_demography_modeling_data <- rbind(Joined_AB_modeling_data, Joined_BC_modeling_data, Joined_BONRO_modeling_data, Joined_BS_modeling_data, Joined_IW_modeling_data, Joined_NRO_modeling_data, Joined_RM_modeling_data, Joined_SM_modeling_data, Joined_WA_modeling_data, Joined_WO_modeling_data) %>%
replace(is.na(.), 0)

Debug <- Combined_demography_modeling_data %>%
					filter(SiteID == 7 & quadrat == 25 & Sp == "sm")
								
	write.csv(Combined_demography_modeling_data, "C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L3/Combined_demography_modeling_data.csv", row.names = FALSE)