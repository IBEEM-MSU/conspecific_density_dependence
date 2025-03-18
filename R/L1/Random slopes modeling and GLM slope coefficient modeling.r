# Andrew C. Eagar 2025

# Random slopes modeling and GLM slope coefficient modeling

library(tidyverse) # For data manipulation & graphing. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.

# Import data for random slopes modeling

RM_demography <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_count_demography.csv")

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")

	# Reformat data for GLM syntax
	{
	RM_demography_spread <- RM_demography %>%
					spread(key = "Time.Since.First.Year", value = "Count")
					
	RM_demography_spread$Starting_counts <- RM_demography_spread$'0'
	
	RM_demography_spread <- RM_demography_spread %>%
					gather(c(5:8), key = "Time.Since.First.Year", value = "Count")
					
	RM_demography_spread$quad_sites <- paste(RM_demography_spread$quadrat, RM_demography_spread$SiteID)
	}

RM_rand_slopes_glm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + (as.numeric(Time.Since.First.Year) | SiteID), family = binomial(link = "logit"), data = RM_demography_spread) # Uncorrelated relationship between random slope and random intercept. Why not just test for an interaction between Time x Myc BA, with a significant interaction telling us that the slopes are variable? Will we get a better fit model? OR would I include the BA terms, then have the BA be the random slopes component instead of time? Probably the latter, I'd think. But again, why not use it as a direct interactive test?

# RM_rand_slopes_glm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + ( 1 | SiteID), family = binomial(link = "logit"), data = RM_demography_spread)

	anova(RM_rand_slopes_glm)
	summary(RM_rand_slopes_glm)
	r.squaredGLMM(RM_rand_slopes_glm)
	
	# liklihood ratio test anova(Model 1, model 2).
	
# Import data for GLM slope coefficient modeling

RM_GLM_slopes <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_GLM_slopes.csv")

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")
