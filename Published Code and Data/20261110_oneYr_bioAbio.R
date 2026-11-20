
#################################################
####### Integrating biotic and abiotic covariates - mixed models
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 10 November 2026
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
allSpData <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")

allSpData$Seedling_myc <- factor(allSpData$Seedling_myc, 
                                 levels = c("AM", "ECM"))

allSpDat <- subset(allSpData, !is.na(propHetAM)) # 35,791 rows

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

### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

### Guild subsets and standardization
std_vars <- c("consp_15m", "hetAM_15m", "hetECM_15m",
              "mean_pdsi", "percOpen", "Base_CEC", "tmin_c", "tmax_c")

z <- function(x) as.numeric(scale(x))

AM_dat  <- subset(allSpDat_surv_clim, Seedling_myc == "AM")
ECM_dat <- subset(allSpDat_surv_clim, Seedling_myc == "ECM")

for (v in std_vars) {
  AM_dat[[paste0(v, "_z")]]  <- z(AM_dat[[v]])
  ECM_dat[[paste0(v, "_z")]] <- z(ECM_dat[[v]])
}

for (d in c("AM_dat", "ECM_dat")) {
  tmp <- get(d)
  tmp$consp_15m_z <- tmp$consp_15m / sd(tmp$consp_15m)   # scaled, NOT centered
  assign(d, tmp)
}

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
##### AM Biotic/Abiotic Model
m_AM_seedlings <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                            consp_15m_z + hetAM_15m_z + hetECM_15m_z +
                            consp_15m_z:hetAM_15m_z + consp_15m_z:hetECM_15m_z +
                            mean_pdsi_z + tmin_c_z + tmax_c_z +
                            percOpen_z + Base_CEC_z +
                            log0yr + log(Prior_3_year_cohorts_counts + 1) +
                            diag(1 + consp_15m_z | Sp) + (1|survYear) + (1|SiteID),
                            family = betabinomial, data = AM_dat)
summary(m_AM_seedlings)
### VIF 
vif_AM <- compute_vif(~ consp_15m_z + hetAM_15m_z + hetECM_15m_z +
                        consp_15m_z:hetAM_15m_z + consp_15m_z:hetECM_15m_z +
                        mean_pdsi_z + tmin_c_z + tmax_c_z +
                        percOpen_z + Base_CEC_z +
                        log0yr + log(Prior_3_year_cohorts_counts + 1), AM_dat)
vif_AM # All < 5

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
### Plots
DHARMa::plotQQunif(sim_AM) # QQ plot + uniformity/dispersion/outlier tests
DHARMa::testUniformity(sim_AM)
DHARMa::testDispersion(sim_AM)
DHARMa::testOutliers(sim_AM)
DHARMa::testZeroInflation(sim_AM)

############################################
############ ECM Biotic/Abiotic Model
m_ECM_seedlings <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                             consp_15m_z + hetAM_15m_z + hetECM_15m_z +
                             consp_15m_z:hetAM_15m_z + consp_15m_z:hetECM_15m_z +
                             percOpen_z + Base_CEC_z +
                             log0yr + log(Prior_3_year_cohorts_counts + 1) +
                             diag(1 + consp_15m_z | Sp) + (1|survYear) + (1|SiteID),
                           family = betabinomial, data = ECM_dat)
summary(m_ECM_seedlings)
### VIF 
vif_ECM <- compute_vif(~ consp_15m_z + hetAM_15m_z + hetECM_15m_z +
                         consp_15m_z:hetAM_15m_z + consp_15m_z:hetECM_15m_z +
                         percOpen_z + Base_CEC_z +
                         log0yr + log(Prior_3_year_cohorts_counts + 1), ECM_dat)
vif_ECM # All < 5

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
### Plots
DHARMa::plotQQunif(sim_ECM) # QQ plot + uniformity/dispersion/outlier tests
DHARMa::testUniformity(sim_ECM)
DHARMa::testDispersion(sim_ECM)
DHARMa::testOutliers(sim_ECM)
DHARMa::testZeroInflation(sim_ECM)


