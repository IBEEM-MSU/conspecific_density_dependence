
#################################################
####### Integrating biotic and abiotic covariates - mixed models
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 27 October 2026
#################################################

## Read in packages
library(tidyverse) 
library(lme4) 
library(MuMIn)
library(lmerTest)
library(emmeans)
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
#allSpDat <- read.csv("final_seedling_demogr_modeling_data.csv")

allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

### Create combined variable that is conspecific 0-yr and prior 3 yr counts
allSpDat$Year_0_3_consp_counts <- allSpDat$Year_0_counts + allSpDat$Prior_3_year_cohorts_counts

## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
#clim_env_allSites <- read.csv("clim_envVars_allSites.csv")
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
# dat <- allSpDat_surv_clim[,c(2,6:8,35:37,45,41:43,50,53,55:57)]
# # Rename columns
# colnames(dat) <- c("SiteID","Sp","Year_0_counts","Year_1_counts",
#                    "Seedling_myc","ConspBA","Log_Recruits",
#                    "Prior_Count","survYear","AM_BA","ECM_BA", 
#                    "Mean_PDSI","perc_Openness",
#                    "Base_CEC", "Min_Temp", "Max_Temp") 
### Proportional model!!
dat <- allSpDat_surv_clim[,c(2,6:8,14,35,37,39,44,45,41,50,53,55:57)]
# Rename columns
colnames(dat) <- c("SiteID","Sp",
                   "Year_0_counts","Year_1_counts",
                   "Prior_Count","Seedling_myc",
                   "Log_Recruits","propConsp",
                   "propHetAM", "Prior_0_3_Count",
                   "survYear", "Mean_PDSI",
                   "perc_Openness", "Base_CEC", 
                   "Min_Temp", "Max_Temp") 

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

m_AM_seedlings <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
    # Biotic main effects
      propConsp + propHetAM +
    # Biotic interactions (from best model)
      propConsp:propHetAM + 
      # Abiotic direct effects (climate + environment)
      Mean_PDSI + Min_Temp + Max_Temp +
      perc_Openness + Base_CEC +
    # Controls
    Log_Recruits + log(Prior_Count + 1) +
    # Random effects
    diag(1 + propConsp | Sp) + (1|survYear) + (1|SiteID),
  family = betabinomial, data = AM_dat)
summary(m_AM_seedlings)

# VIF 
vif_AM <- compute_vif(~ propConsp + propHetAM + propConsp:propHetAM +
                        Mean_PDSI + Min_Temp + Max_Temp +
                        perc_Openness + Base_CEC +
                        Log_Recruits + log(Prior_Count + 1), AM_dat)
vif_AM # All < 3

# R2 — marginal and conditional (survival-proportion scale)
mf   <- model.frame(m_AM_seedlings)     # rows actually used in the fit
resp <- mf[[1]]                         # the cbind() matrix
n0   <- rowSums(resp)                   # cohort size
obs_prop <- resp[, 1] / n0              # observed survival proportion

pred_marg <- predict(m_AM_seedlings, type = "response", re.form = NA)
pred_cond <- predict(m_AM_seedlings, type = "response")

R2_pred_marg <- cov.wt(cbind(obs_prop, pred_marg), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_cond <- cov.wt(cbind(obs_prop, pred_cond), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_marg
R2_pred_cond

### Diagnostics
sim_AM <- DHARMa::simulateResiduals(m_AM_seedlings, n = 1000)
DHARMa::plotQQunif(sim_AM) # QQ plot with tests
DHARMa::plotResiduals(sim_AM) # residuals vs. predicted
DHARMa::testDispersion(sim_AM) # No overdispersion
DHARMa::testZeroInflation(sim_AM) # No zero-inflation
DHARMa::testUniformity(sim_AM) # Non-uniform, normal with large datasets 
DHARMa::testOutliers(sim_AM)

############################################
############ ECM Biotic/Abiotic Model

m_ECM_seedlings <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
    # Biotic main effects
      propConsp + propHetAM +
    # Biotic interactions (from best model)
      propConsp:propHetAM + 
    # Abiotic direct effects (environment only)
      perc_Openness + Base_CEC +
      # Controls
      Log_Recruits + log(Prior_Count + 1) +
      # Random effects
      (diag(1 + propConsp | Sp)) + (1|survYear) + (1|SiteID),
    family = betabinomial, data = ECM_dat)
summary(m_ECM_seedlings)

# R2 — marginal and conditional (survival-proportion scale)
mf   <- model.frame(m_ECM_seedlings)     # rows actually used in the fit
resp <- mf[[1]]                         # the cbind() matrix
n0   <- rowSums(resp)                   # cohort size
obs_prop <- resp[, 1] / n0              # observed survival proportion

pred_marg <- predict(m_ECM_seedlings, type = "response", re.form = NA)
pred_cond <- predict(m_ECM_seedlings, type = "response")

R2_pred_marg <- cov.wt(cbind(obs_prop, pred_marg), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_cond <- cov.wt(cbind(obs_prop, pred_cond), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_marg
R2_pred_cond
  

### Diagnostics
sim_ECM <- DHARMa::simulateResiduals(m_ECM_seedlings, n = 1000)
DHARMa::testDispersion(sim_ECM) # No overdispersion
DHARMa::testZeroInflation(sim_ECM) # No zero-inflation
DHARMa::testUniformity(sim_ECM) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_ECM)


