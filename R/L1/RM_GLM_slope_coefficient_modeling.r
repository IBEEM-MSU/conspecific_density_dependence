# Andrew C. Eagar 2025

# GLM slope coefficient modeling

library(tidyverse) # For data manipulation & graphing. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.

# Import data for GLM slope coefficient modeling

RM_GLM_slopes <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_GLM_slopes.csv")

RM_prior_cohort_densities <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_prior_cohort_densities.csv")

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")

Site_BA_quads <- Site_BA %>%
					mutate(quadrat = as.numeric(X) + 0.5) %>%
					mutate(Transect_quad_year_unifier = paste(TransectID, quadrat, Year)) %>%
					select(-c(1, 2, 4, 5, 6, 7, 63))

## Calculate conmycorrhizal, heterospecific, and heteromycorrhizal BA for RM
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH rm
{ Site_BA_quads$conMyc_5m = rowSums(Site_BA_quads[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
Site_BA_quads$conMyc_10m = rowSums(Site_BA_quads[, c("sm_10m", "ae_10m", "wa_10m", 
                                        "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
Site_BA_quads$conMyc_15m = rowSums(Site_BA_quads[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - NO rm
Site_BA_quads$conMycNS_5m = rowSums(Site_BA_quads[, c("sm_5m", "ae_5m", "wa_5m", 
                                         "bc_5m", "svb_5m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_10m = rowSums(Site_BA_quads[, c("sm_10m", "ae_10m", "wa_10m", 
                                          "bc_10m", "svb_10m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_15m = rowSums(Site_BA_quads[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
Site_BA_quads$hetMyc_5m = rowSums(Site_BA_quads[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_10m = rowSums(Site_BA_quads[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                          "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                          "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_15m = rowSums(Site_BA_quads[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
}
										
## Use an ifelse statement to add a column with the closest BA census year -a B. McNichol					
{ RM_GLM_slopes$baYear = ifelse(RM_GLM_slopes$Cohort.Year == "1999", "1999", 
                        ifelse(RM_GLM_slopes$Cohort.Year == "2000", "1999",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2001", "1999",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2002", "1999",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2003", "2005",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2004", "2005",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2005", "2005",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2006", "2005",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2007", "2007",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2008", "2007",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2009", "2009",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2010", "2009",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2011", "2011",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2012", "2011",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2013", "2013",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2014", "2013",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2015", "2016",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2016", "2016",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2017", "2016",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2018", "2019",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2019", "2019",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2020", "2019",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2021", "2022",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2022", "2022",
                        ifelse(RM_GLM_slopes$Cohort.Year == "2023", "2022",
                               "2022")))))))))))))))))))))))))
} 
							   
RM_GLM_slopes_spread <- RM_GLM_slopes %>%
				spread(key = "Time.Since.First.Year", value = "Count")	
					
RM_GLM_slopes_spread$Starting_counts <- RM_GLM_slopes_spread$'0' 

RM_GLM_slopes_BA <- RM_GLM_slopes_spread %>%
				gather(c(6:9), key = "Time.Since.First.Year", value = "Count") %>%
				mutate(Transect_quad_year_unifier = paste(TransectID, quadrat, baYear)) %>%
				left_join(Site_BA_quads, by = "Transect_quad_year_unifier")

(#1 yr / #0 yr) ~ Conspecific BA * % EMF hosts + Con-myc BA + Hetero-myc BA + (total # RM present at 0 yr - #0 yr seedlings) + #0 yr seedlings + (1| Site + survYear)
