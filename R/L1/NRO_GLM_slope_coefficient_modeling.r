# Authors: A Eagar & B McNichol 5-2025

# GLM slope coefficient modeling script. This conducts a non-parametric test on the quadrat-level slope values representing a seedling cohort's mortality rate over three years. This method follows Anderson, M. and Braak, C.T., 2003. Permutation tests for multi-factorial analysis of variance. Journal of statistical computation and simulation, 73(2), pp.85-113.

library(tidyverse) # For data manipulation & graphing. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.

# Import data for LM slope coefficient modeling

NRO_GLM_slopes <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/nro_GLM_slopes.csv") %>%
					mutate(Joiner = paste(SiteID,".", quadrat,".",Cohort.Year))

NRO_prior_cohort_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/nro_prior_cohort_densities.csv") %>%
					mutate(Joiner = paste(SiteID,".", quadrat,".",Cohort.Year))

NRO_demography <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/nro_count_demography.csv") %>%
					mutate(Joiner = paste(SiteID,".", quadrat,".",First.Year))

All_cohort_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/allSp_prior_cohort_densities.csv") %>% 
					mutate(Transect_joiner = paste(TransectID,".",Obs.Year))

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")

Site_BA_quads <- Site_BA %>%
					mutate(quadrat = as.numeric(X) + 0.5, ba_joiner = paste(TransectID,".",Year)) %>%
					select(-c(X.1, X, Y, Midpoint_X, Midpoint_Y)) %>%
					replace(is.na(.), 0)

## Calculate conmycorrhizal, heterospecific, and heteromycorrhizal BA for RM
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (nro), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH rm
{ Site_BA_quads$conMyc_5m = rowSums(Site_BA_quads[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","nro_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
Site_BA_quads$conMyc_10m = rowSums(Site_BA_quads[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                        "nro_10m","nro_10m", "bo_10m","bonro_10m",
                                        "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
Site_BA_quads$conMyc_15m = rowSums(Site_BA_quads[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                        "nro_15m","nro_15m", "bo_15m","bonro_15m",
                                        "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - NO rm
Site_BA_quads$conMycNS_5m = rowSums(Site_BA_quads[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_10m = rowSums(Site_BA_quads[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                        "iw_10m", "bo_10m","bonro_10m",
                                        "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_15m = rowSums(Site_BA_quads[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                        "iw_15m", "bo_15m","bonro_15m",
                                        "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
Site_BA_quads$hetMyc_5m = rowSums(Site_BA_quads[, c("sm_5m", "rm_5m", "ae_5m", "wa_5m", 
                                         "bc_5m", "svb_5m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_10m = rowSums(Site_BA_quads[, c("sm_10m", "rm_10m", "ae_10m", "wa_10m", 
                                         "bc_10m", "svb_10m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_15m = rowSums(Site_BA_quads[, c("sm_15m", "rm_15m", "ae_15m", "wa_15m", 
                                         "bc_15m", "svb_15m")],  na.rm = TRUE)
}
										
## Use an ifelse statement to add a column with the closest BA census year -a B. McNichol					
{ NRO_GLM_slopes$baYear = ifelse(NRO_GLM_slopes$Cohort.Year == "1999", "1999", 
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2000", "1999",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2001", "1999",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2002", "1999",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2003", "2005",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2004", "2005",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2005", "2005",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2006", "2005",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2007", "2007",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2008", "2007",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2009", "2009",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2010", "2009",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2011", "2011",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2012", "2011",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2013", "2013",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2014", "2013",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2015", "2016",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2016", "2016",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2017", "2016",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2018", "2019",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2019", "2019",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2020", "2019",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2021", "2022",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2022", "2022",
                        ifelse(NRO_GLM_slopes$Cohort.Year == "2023", "2022",
                               "2022")))))))))))))))))))))))))
} 

# Join all data together. Will export .csv file for final upload.

NRO_slope_modeling_data <- NRO_GLM_slopes %>%
					mutate(ba_joiner = paste(TransectID,".",baYear), Transect_joiner = paste(TransectID,".",Cohort.Year)) %>%
					left_join(NRO_prior_cohort_densities, by = "Joiner") %>%
					left_join(NRO_demography, by = "Joiner") %>%
					left_join(Site_BA_quads, by = "ba_joiner") %>%
					left_join(All_cohort_densities, by = "Transect_joiner") %>%
					select(SiteID.x, quadrat.x, Cohort.Year.x, ba_joiner, GLM.Slope.Estimate, Count, Prior_3_year_cohorts_counts, Prior_2_year_cohorts_and_new_1st_year_counts,  Prior_1_year_cohorts_and_1st_year_and_new_2nd_year_counts, nro_15m, conMycNS_15m, hetMyc_15m) %>%
					rename("SiteID" = c(1), "quadrat" = c(2), "Cohort.Year" = c(3)) %>%
					unique() %>% 
					mutate(ECM_prop_15 = (nro_15m + conMycNS_15m)/(nro_15m + conMycNS_15m + hetMyc_15m), nro_prop_15 = nro_15m/(nro_15m + conMycNS_15m + hetMyc_15m), Avg_3_yr_prior_cohort_counts = replace_na((Prior_3_year_cohorts_counts + Prior_2_year_cohorts_and_new_1st_year_counts + Prior_1_year_cohorts_and_1st_year_and_new_2nd_year_counts)/3, 0)) %>%
					replace(is.na(.), 0)

	# write.csv(NRO_slope_modeling_data, "C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/NRO_slope_modeling_data.csv", row.names = FALSE)

# Non-parametric penroutation testing

NRO_slope_modeling_data <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/NRO_slope_modeling_data.csv")

nperm <- 999 # To test this code, drop this number to 9 or 99 - at 999 penroutations, this code takes about 20-30 minutes to compile on a PC with an i7 core and 64 gb RAM.

Rand_results_con_den <- matrix(nrow = nperm, ncol = 1, 0) # we have 7 terms!!! to test in this model, so we'll need 4 dataframes to store the randomized F-value test results.
	Rand_results_con_myc <- matrix(nrow = nperm, ncol = 1, 0)
	Rand_results_het_myc <- matrix(nrow = nperm, ncol = 1, 0)
	Rand_results_con_den_het_myc <- matrix(nrow = nperm, ncol = 1, 0)
	Rand_results_con_myc_het_myc <- matrix(nrow = nperm, ncol = 1, 0)
	Rand_results_con_seedlings <- matrix(nrow = nperm, ncol = 1, 0)
	Rand_results_other_seedlings <- matrix(nrow = nperm, ncol = 1, 0)

for(i in 1:nperm) # For loop for randomized F-values
{

	NRO_slope_modeling_data <- NRO_slope_modeling_data %>%
					group_by(Cohort.Year, SiteID) %>%
					mutate(Rand.GLM.Slope.Estimate = sample(GLM.Slope.Estimate, n())) # This randomizes the slope coefficients for each quadrat within a site and cohort year without replacing values.

		# Test_data <- NRO_slope_modeling_data %>%
							# filter(Cohort.Year == 2002, SiteID == 3) %>%
							# select(GLM.Slope.Estimate, Rand.GLM.Slope.Estimate)
			
			# as.data.frame(Test_data)
	
	nro_mod_BA_rand <- lmer(Rand.GLM.Slope.Estimate ~ nro_15m + conMycNS_15m + 
                    hetMyc_15m + nro_15m:hetMyc_15m + conMycNS_15m:hetMyc_15m + log10(Count + 1) + log(Avg_3_yr_prior_cohort_counts + 1) +
                    (1 | SiteID) + (1 | Cohort.Year),
					data = NRO_slope_modeling_data)
					
	Rand_aov <- anova(nro_mod_BA_rand) # Save results to store F-values below
					
		Rand_results_con_den[i] <- Rand_aov$"F value"[1]				
		Rand_results_con_myc[i] <- Rand_aov$"F value"[2]
		Rand_results_het_myc[i] <- Rand_aov$"F value"[3]
		Rand_results_con_den_het_myc[i] <- Rand_aov$"F value"[4]
		Rand_results_con_myc_het_myc[i] <- Rand_aov$"F value"[5]
		Rand_results_con_seedlings[i] <- Rand_aov$"F value"[6]
		Rand_results_other_seedlings[i] <- Rand_aov$"F value"[7]
					
}

	# LM test with real data. We will compare these F-values to those from the penroutations above by ranking their magnitude. If our observations are not statistically significant, the true F-values will rank somewhere within the "normal" simulated distribution of F-values. If our observations are significant, our F-values will fall towards one of the tails of the distribution and a critical threshold (p-value at 0.05) can be estimated.
			
	nro_mod_BA_true <- lmer(GLM.Slope.Estimate ~ nro_15m + conMycNS_15m + 
                    hetMyc_15m + nro_15m:hetMyc_15m + conMycNS_15m:hetMyc_15m + log10(Count + 1) + log(Avg_3_yr_prior_cohort_counts + 1) +
                    (1 | SiteID) + (1 | Cohort.Year),
					data = NRO_slope_modeling_data)
					
	Obs_aov <- anova(nro_mod_BA_true)

		Obs_results_con_den <- as.vector(Obs_aov[1,5])
			All_results_con_den <- rbind(Rand_results_con_den, Obs_results_con_den)
			Ranks_con_den <- rank(All_results_con_den)
			P_values_con_den <- 1-Ranks_con_den[nperm+1]/(nperm+1)

		Obs_results_con_myc <- as.vector(Obs_aov[2,5])
			All_results_con_myc <- rbind(Rand_results_con_myc, Obs_results_con_myc)
			Ranks_con_myc <- rank(All_results_con_myc)
			P_values_con_myc <- 1-Ranks_con_myc[nperm+1]/(nperm+1)
		
		Obs_results_het_myc <- as.vector(Obs_aov[3,5])
			All_results_het_myc <- rbind(Rand_results_het_myc, Obs_results_het_myc)
			Ranks_het_myc <- rank(All_results_het_myc)
			P_values_het_myc <- 1-Ranks_het_myc[nperm+1]/(nperm+1)
		
		Obs_results_con_den_het_myc <- as.vector(Obs_aov[4,5])
			All_results_con_den_het_myc <- rbind(Rand_results_con_den_het_myc, Obs_results_con_den_het_myc)
			Ranks_con_den_het_myc <- rank(All_results_con_den_het_myc)
			P_values_con_den_het_myc <- 1-Ranks_con_den_het_myc[nperm+1]/(nperm+1)		
		Obs_results_con_myc_het_myc <- as.vector(Obs_aov[5,5])
			All_results_con_myc_het_myc <- rbind(Rand_results_con_myc_het_myc, Obs_results_con_myc_het_myc)
			Ranks_con_myc_het_myc <- rank(All_results_con_myc_het_myc)
			P_values_con_myc_het_myc <- 1-Ranks_con_myc_het_myc[nperm+1]/(nperm+1)
		
		Obs_results_con_seedlings <- as.vector(Obs_aov[6,5])
			All_results_con_seedlings <- rbind(Rand_results_con_seedlings, Obs_results_con_seedlings)
			Ranks_con_seedlings <- rank(All_results_con_seedlings)
			P_values_con_seedlings <- 1-Ranks_con_seedlings[nperm+1]/(nperm+1)
			
		Obs_results_other_seedlings <- as.vector(Obs_aov[7,5])
			All_results_other_seedlings <- rbind(Rand_results_other_seedlings, Obs_results_other_seedlings)
			Ranks_other_seedlings <- rank(All_results_other_seedlings)
			P_values_other_seedlings <- 1-Ranks_other_seedlings[nperm+1]/(nperm+1)

	P_values_con_den
		hist(All_results_con_den) # vs
		Obs_results_con_den
	
	P_values_con_myc
		hist(All_results_con_myc) # vs
		Obs_results_con_myc
		
	P_values_het_myc
		hist(All_results_het_myc) # vs
		Obs_results_het_myc

	P_values_con_den_het_myc
		hist(All_results_con_den_het_myc)
		Obs_results_con_den_het_myc

	P_values_con_myc_het_myc
		hist(All_results_con_myc_het_myc)
		Obs_results_con_myc_het_myc

	P_values_con_seedlings
		hist(All_results_con_seedlings)
		Obs_results_con_seedlings
		
	P_values_other_seedlings
		hist(All_results_other_seedlings)
		Obs_results_other_seedlings
	
	r.squaredGLMM(nro_mod_BA_true)
	
	# Store modeling results & export for plotting.

	summary(nro_mod_BA_true)
	
	NRO_GLM_slope_coefficient_modeling_results <- data.frame(species = c("nro", "nro", "nro", "nro", "nro", "nro", "nro"), 
															term = c("Conspecific_BA", "Conmycorrhizal_BA", "Heteromycorrhizal_BA", "Log_Recruit_Density", "Log_Prior_Density", "Consp_x_Heteromycorrhizal_BA", "Conmyc_x_Heteromycorrhizal_BA"), 
															estimate = c("0.7775", "1.9199", "-2.2478", "21.9498", "5.3454", "-0.1363", "0.9898"), 
															std.error = c("1.0830", "2.0415", "2.1879", "1.1318", "1.0040", "0.9890", "2.1562"), 
															statistic = c("0.718", "0.940", "-1.027", "19.394", "5.324", "-0.138", "0.459"), 
															p.value = c(P_values_con_den, P_values_con_myc, P_values_het_myc, P_values_con_seedlings, P_values_other_seedlings, P_values_con_den_het_myc, P_values_con_myc_het_myc),
															r.squared = c("0.1239852", "0.1239852", "0.1239852", "0.1239852", "0.1239852", "0.1239852", "0.1239852"))

	# write.csv(NRO_GLM_slope_coefficient_modeling_results, "C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L2/NRO_GLM_slope_coefficient_modeling_results.csv", row.names = FALSE)