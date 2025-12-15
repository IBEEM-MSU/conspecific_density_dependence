
#################################################
####### Analyses of one-year survivorship WRT abiotic environment
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 4 December 2025
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
library(glmmTMB)
library(performance)

### Most complex model to fit (15-m radius)
###  Abiotic model of one-year survivorship 
# (#1 yr) ~ meanGS_PDSI + minDailyTemp + maxDailyTemp + percOpenness + baseCEC + log(#0 yr) + (1 | Site) + (1 | Year)

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
allSpDat$survYear = as.factor(allSpDat$survYear)
allSpDat_surv_clim <- allSpDat %>%
  left_join(clim_env_allSites, by = c("TransectID", "survYear"))

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

#########################################################
#########################################################
################# All Species Model

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

#### Null model
allSp_abio_null <- glmmTMB(Year_1_counts ~ 1 +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                family = nbinom2, data = allSpDat_surv_clim)
summary(allSp_abio_null)
performance::r2(allSp_abio_null)

##### Effects of environment and climate
allSp_abio <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = nbinom2, data = allSpDat_surv_clim)
summary(allSp_abio)
compute_vif(allSp_abio)
performance::r2(allSp_abio)

##### Effects of climate alone
allSp_clim <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr + (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat_surv_clim)
summary(allSp_clim)
compute_vif(allSp_clim)
performance::r2(allSp_clim)

##### Effects of underlying environment alone (site characteristics)
allSp_env <- glmmTMB(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat_surv_clim)
summary(allSp_env)
compute_vif(allSp_env)
performance::r2(allSp_env)

#### Model comparison
### Compare AICs - full model integrating climate and environment is better
AIC(allSp_abio_null, allSp_abio, allSp_clim, allSp_env)

##################################
################# Mycorrhizal type-level Models

## Subset for AM
AM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "AM")
## Subset for ECM
ECM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "ECM")

#### Null model - AM
AM_abio_null <- glmmTMB(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
              family = nbinom2, data = AM_Dat_surv_clim)
summary(AM_abio_null)
performance::r2(AM_abio_null)
#### Null model - ECM
ECM_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_abio_null)
performance::r2(ECM_abio_null)

##### Effects of environment and climate - AM
AM_abio <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = nbinom2, data = AM_Dat_surv_clim)
summary(AM_abio)
compute_vif(AM_abio)
performance::r2(AM_abio)
##### Effects of environment and climate - ECM
ECM_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_abio)
car::vif(ECM_abio)
performance::r2(ECM_abio)

##### Effects of climate alone - AM
AM_clim <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
              family = poisson(link = "log"), data = AM_Dat_surv_clim)
summary(AM_clim)
compute_vif(AM_clim)
performance::r2(AM_clim)
##### Effects of climate alone - ECM
ECM_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_clim)
car::vif(ECM_clim)
performance::r2(ECM_clim)

##### Effects of environment alone (site characteristics) - AM
AM_env <- glmmTMB(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
            family = nbinom2, data = AM_Dat_surv_clim)
summary(AM_env)
compute_vif(AM_env)
performance::r2(AM_env)
##### Effects of environment alone (site characteristics) - ECM
ECM_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear) + (1|Sp),
                family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_env)
car::vif(ECM_env)
performance::r2(ECM_env)

#### Model comparison
### Compare AICs - AM
AIC(AM_abio_null, AM_abio, AM_clim, AM_env) # Model with climate and env best
### Compare AICs - ECM
AIC(ECM_abio_null, ECM_abio, ECM_clim, ECM_env)
anova(ECM_abio, ECM_env) # Model with environment only is best


########## Plotting effect sizes for abiotic models (All, AM, ECM)
# Function to extract and tidy fixed effects
extract_fixef_abio <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(component = case_when(
        term == "mean_pdsi" ~ "Mean PDSI",
        term == "mintemp"  ~ "Min Temperature",
        term == "maxtemp"  ~ "Max Temperature",
        term == "percOpen" ~ "% Canopy Openness",
        term == "Base_CEC" ~ "Base CEC",
        term == "log0yr"   ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type)
}

# ===============================================================
# Extract fixed effects
fixef_all <- extract_fixef_abio(allSp_abio, "All")
fixef_AM  <- extract_fixef_abio(AM_abio, "AM")
fixef_ECM_raw <- extract_fixef_abio(ECM_env, "ECM")

# ===============================================================
# Ensure ECM has explicit NA rows for unsupported climate variables
abiotic_terms <- c(
  "Mean PDSI",
  "Min Temperature",
  "Max Temperature",
  "% Canopy Openness",
  "Base CEC",
  "Log(Recruit Density)")

# Remove term column so rows_update works
fixef_ECM_raw <- fixef_ECM_raw %>%
  select(-term)

fixef_ECM_complete <- tibble(
  component   = abiotic_terms,
  estimate    = NA_real_,
  std.error   = NA_real_,
  p.value     = NA_real_,
  Significant = NA_character_,   # ← IMPORTANT
  Model       = "ECM"
) %>%
  rows_update(fixef_ECM_raw, by = c("component", "Model"))



# ===============================================================
# Combine all models
fixef_df <- bind_rows(fixef_all, fixef_AM, fixef_ECM_complete)

# ===============================================================
# Ordering and plotting positions
desired_order <- abiotic_terms

fixef_df <- fixef_df %>%
  mutate(component   = factor(component, levels = rev(desired_order)),
    Model       = factor(Model, levels = c("All", "AM", "ECM")),
    Significant = factor(Significant, levels = c("yes", "no")),
    y_offset = case_when(
      Model == "All" ~  0.0,
      Model == "AM"  ~ -0.2,
      Model == "ECM" ~  0.2),
    y_pos = as.numeric(component) + y_offset)

# ===============================================================
# Plot
ggplot(fixef_df, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars (only where estimates exist)
  geom_errorbar( data = subset(fixef_df, !is.na(estimate)),
    aes(xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error),
    width = 0.1,color = "black") +
  # Points (only where estimates exist)
  geom_point(data = subset(fixef_df, !is.na(estimate)),
    aes(shape = Model, fill = Significant),
    size = 2.7, color = "black", stroke = 1) + 
  scale_y_continuous(breaks = seq_along(desired_order),
    labels = rev(desired_order)) +
  scale_shape_manual(values = c("All" = 21,
    "AM"  = 24,"ECM" = 22 )) +
  scale_fill_manual(values = c(
    "yes" = "black", "no"  = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(
      override.aes = list(shape = 21, color = "black", size = 3.5)))





