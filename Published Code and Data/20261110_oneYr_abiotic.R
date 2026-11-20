
#################################################
####### One-year seedling survival WRT abiotic environment
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
library(glmmTMB)
library(performance)
library(car)

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

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

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

###################################################
################################################### 
## Subset for AM
AM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "AM")
## Subset for ECM
ECM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "ECM")

##### ##### ##### ##### ##### #####
#### Null model (Seedling density only) - AM
AM_abio_null <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                          log0yr + (1|SiteID) + (1|survYear) + (1|Sp),
              family = betabinomial, data = AM_Dat_surv_clim)
summary(AM_abio_null)

##### ##### ##### ##### ##### #####
####  Null model (Seedling density only) - ECM
ECM_abio_null <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                           log0yr + (1|SiteID) + (1|survYear) + (1|Sp),
                         family = betabinomial, data = ECM_Dat_surv_clim)
summary(ECM_abio_null)

##### ##### ##### ##### ##### #####
#### Effects of environment and climate - AM
AM_abio <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                   mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = betabinomial, data = AM_Dat_surv_clim)
summary(AM_abio)
# VIF 
compute_vif(~ mean_pdsi + mintemp + maxtemp +
           percOpen + Base_CEC + log0yr, AM_Dat_surv_clim)

##### ##### ##### ##### ##### #####
#### Effects of environment and climate - ECM
ECM_abio <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                     mean_pdsi + mintemp + maxtemp +
                     percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = betabinomial, data = ECM_Dat_surv_clim)
summary(ECM_abio)
# VIF
compute_vif(~ mean_pdsi + mintemp + maxtemp +
              percOpen + Base_CEC + log0yr, ECM_Dat_surv_clim)

##### ##### ##### ##### ##### #####
#### Effects of climate only - AM
AM_clim <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                     mean_pdsi + mintemp + maxtemp +
                     log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = betabinomial, data = AM_Dat_surv_clim)
summary(AM_clim)
# VIF
compute_vif(~ mean_pdsi + mintemp + maxtemp +
              log0yr, AM_Dat_surv_clim)

##### ##### ##### ##### ##### #####
#### Effects of climate only - ECM
ECM_clim <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                     mean_pdsi + mintemp + maxtemp +
                     log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = betabinomial, data = ECM_Dat_surv_clim)
summary(ECM_clim)
# VIF
compute_vif(~ mean_pdsi + mintemp + maxtemp +
              log0yr, ECM_Dat_surv_clim)

##### ##### ##### ##### ##### #####
#### Effects of environment only - AM
AM_env <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                    percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
            family = betabinomial, data = AM_Dat_surv_clim)
summary(AM_env)
# VIF
compute_vif(~ percOpen + Base_CEC + log0yr, AM_Dat_surv_clim)

##### ##### ##### ##### ##### #####
#### Effects of environment only - ECM
ECM_env <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                    percOpen + Base_CEC + log0yr +
                    (1|SiteID) + (1|survYear) + (1|Sp),
                  family = betabinomial, data = ECM_Dat_surv_clim)
summary(ECM_env)
# VIF
compute_vif(~ percOpen + Base_CEC + log0yr, ECM_Dat_surv_clim)

##### Model comparison
### Compare AICs - AM - Model with climate and env best
aic_table <- AIC(AM_abio_null, AM_abio, AM_clim, AM_env)
aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]

### Compare AICs - ECM - Model with environment only is best
aic_table <- AIC(ECM_abio_null, ECM_abio, ECM_clim, ECM_env)
aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]
 
