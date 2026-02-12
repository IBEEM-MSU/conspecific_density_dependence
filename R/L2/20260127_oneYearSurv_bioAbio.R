
#################################################
####### Integrating biotic and abiotic covariates - mixed models
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 27 January 2026
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting lme models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)
library(ggplot2)
library(DHARMa)
library(broom.mixed)
library(tibble)
library(performance)
library(janitor)
library(forcats)
library(scales)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")
allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

# Add X column (quadrat by the 0.5)
allSpDat$survYear = as.factor(allSpDat$survYear)
allSpDat$X = allSpDat$quadrat - 0.5

# Add true site number
allSpDat$SiteNo = ifelse(allSpDat$SiteID == 1, "49.01",
                  ifelse(allSpDat$SiteID == 2, "59.12",
                  ifelse(allSpDat$SiteID == 3, "8.21",
                  ifelse(allSpDat$SiteID == 4, "39.21",
                  ifelse(allSpDat$SiteID == 5, "16.35",
                  ifelse(allSpDat$SiteID == 6, "60.35",
                  ifelse(allSpDat$SiteID == 7, "7.40",
                  ifelse(allSpDat$SiteID == 8, "40.40",
                  ifelse(allSpDat$SiteID == 9, "21.43",
                  ifelse(allSpDat$SiteID == 10, "25.43",
                  ifelse(allSpDat$SiteID == 11, "6.45", "24.45")))))))))))

# Create TransectID
allSpDat$TransectID = paste0(allSpDat$SiteNo, allSpDat$X)

### Left join to add climate variables to data
allSpDat_surv_clim <- allSpDat %>%
  left_join(clim_env_allSites, by = c("TransectID", "survYear", "X"))

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

##### Subset out only the relevant columns
dat <- allSpDat_surv_clim[,c(2,6,8,35:38,41:43,49,52,54,55,56)]
# Rename columns
colnames(dat) <- c("SiteID","Sp","Year_1_counts","Seedling_myc", 
                   "ConspBA","Log_Recruits","Log_Prior","survYear",
                   "AM_BA","ECM_BA", "Mean_PDSI","perc_Openness",
                   "Base_CEC", "Min_Temp", "Max_Temp") 

## Subset for AM
AM_dat = filter(dat, Seedling_myc == "AM")
## Subset for ECM
ECM_dat = filter(dat, Seedling_myc == "ECM")

#####  Helper function: VIF for glmmTMB fixed effects
compute_vif <- function(formula, data) {
  mm <- model.matrix(formula, data)
  # Remove intercept
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) < 2)
    stop("Need at least 2 predictors to compute VIF.")
  dummy <- rnorm(nrow(mm))
  lm_obj <- lm(dummy ~ ., data = as.data.frame(mm))
  car::vif(lm_obj)
}

############################################
############ AM Biotic/Abiotic Model

m_AM_seedlings <- glmmTMB(Year_1_counts ~
    # Biotic main effects
      ConspBA + AM_BA + ECM_BA +
    # Biotic interactions (from best model)
      ConspBA:AM_BA +
      ConspBA:ECM_BA +
    # Abiotic direct effects
      Mean_PDSI + Min_Temp + Max_Temp +
      perc_Openness + Base_CEC +
    # Controls
      Log_Prior + Log_Recruits +
    # Random effects
    (1|Sp) + (1|survYear) + (1|SiteID),
  family = nbinom2, data = AM_dat)
summary(m_AM_seedlings)

# VIF 
vif_AM <- compute_vif(~ ConspBA + AM_BA + ECM_BA +
                        ConspBA:AM_BA + ConspBA:ECM_BA +
                        Mean_PDSI + Min_Temp + Max_Temp +
                        perc_Openness + Base_CEC +
                        Log_Prior + Log_Recruits, AM_dat)
vif_AM # All <6

# R2
# Marginal prediction R² (log-scale)
pred_marg <- predict(m_AM_seedlings, type = "response", re.form = NA)
obs <- AM_dat$Year_1_counts
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(m_AM_seedlings, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond

### Diagnostics
sim_AM <- DHARMa::simulateResiduals(m_AM_seedlings, n = 1000)
DHARMa::testDispersion(sim_AM) # No overdispersion
DHARMa::testZeroInflation(sim_AM) # No zero-inflation
DHARMa::testUniformity(sim_AM) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_AM)




############################################
############ ECM Biotic/Abiotic Model

m_ECM_seedlings <- glmer(Year_1_counts ~
    # Biotic main effects
      ConspBA + AM_BA + ECM_BA +
    # Biotic interactions (supported)
      ConspBA:AM_BA +
      ConspBA:ECM_BA +
    # Abiotic direct effects (best model)
      perc_Openness + Base_CEC +
      # Controls
      Log_Prior + Log_Recruits +
      # Random effects
      (1|Sp) + (1|survYear) + (1|SiteID),
  family = poisson(link = "log"), data = ECM_dat)
summary(m_ECM_seedlings)

# R2
# Marginal prediction R² (log-scale)
pred_marg <- predict(m_ECM_seedlings, type = "response", re.form = NA)
obs <- ECM_dat$Year_1_counts
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(m_ECM_seedlings, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond
  

### Diagnostics
sim_ECM <- DHARMa::simulateResiduals(m_ECM_seedlings, n = 1000)
DHARMa::testDispersion(sim_ECM) # No overdispersion
DHARMa::testZeroInflation(sim_ECM) # No zero-inflation
DHARMa::testUniformity(sim_ECM) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_ECM)



