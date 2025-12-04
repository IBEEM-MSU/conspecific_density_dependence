
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
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")
## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1
allSpDat$survYear = as.factor(allSpDat$survYear)

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

##  Helper function: VIF for glmmTMB fixed effects
compute_vif <- function(formula, data) {
  mm <- model.matrix(formula, data = data)
  mm_df <- as.data.frame(mm[, -1])  # Remove intercept
  fake_y <- rnorm(nrow(mm_df))
  lm_fake <- lm(fake_y ~ ., data = mm_df)
  car::vif(lm_fake)
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
                      log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
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
anova(ECM_abio, ECM_clim) # Model with climate and env best

########## Plotting effect sizes for models with climate and environment
# --- Function to extract fixed effects and tidy them
extract_fixef_abio <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      # Convert model terms to nicer labels
      component = case_when(
        term == "mean_pdsi" ~ "Mean PDSI",
        term == "mintemp" ~ "Min Temperature",
        term == "maxtemp" ~ "Max Temperature",
        term == "percOpen" ~ "% Canopy Openness",
        term == "Base_CEC" ~ "Base CEC",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & 
              p.value <= 0.05, "yes", "no"),
      Model = model_type)
}

# --- Extract for all three models
fixef_all  <- extract_fixef_abio(allSp_abio, "All")
fixef_AM   <- extract_fixef_abio(AM_abio, "AM")
fixef_ECM  <- extract_fixef_abio(ECM_abio, "ECM")

# --- Combine
fixef_df <- bind_rows(fixef_all, fixef_AM, fixef_ECM)

# Define vertical order for plotting (top → bottom)
desired_order <- c("Mean PDSI", "Min Temperature", "Max Temperature",
                   "% Canopy Openness", "Base CEC",
                   "Log(Recruit Density)")

# --- Prepare data for plotting
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = rev(desired_order)),
    Model = factor(Model, levels = c("All", "AM", "ECM")),
    Significant = factor(Significant, levels = c("yes", "no")),
    y_offset = case_when(
      Model == "All" ~ 0,
      Model == "AM"  ~ -0.2,
      Model == "ECM" ~  0.2), 
    y_pos = as.numeric(component) + y_offset)
  
# --- Plot
ggplot(fixef_df, aes(x = estimate, y = y_pos)) +   # use numeric y_pos
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error),
                width = 0.1, orientation = "y", color = "black") +
  # Points
  geom_point(aes(shape = Model, fill = Significant),
             size = 2.5, color = "black", stroke = 1) +
  # Y-axis: show factor labels at base positions
  scale_y_continuous(
    breaks = 1:length(desired_order),
    labels = rev(desired_order)) +
  # Shapes for models
  scale_shape_manual(values = c("All" = 21, "AM" = 24, "ECM" = 22)) +
  # Fill for significance
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  # Fix fill legend shape
  guides(fill = guide_legend(override.aes = 
                               list(shape = 21, color = "black", size = 3.5)))

#########################################################
#########################################################
################# Species-Level - RM

##### Subset RM data
rm_surv <- subset(allSpDat_surv_clim, Sp == "rm")

#### Null model
rm_abio_null <- glmmTMB(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear),
                         family = nbinom2, data = rm_surv)
summary(rm_abio_null)
performance::r2(rm_abio_null)

##### Effects of environment and climate
rm_abio <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = rm_surv)
summary(rm_abio)
compute_vif(rm_abio)
performance::r2(rm_abio)

##### Effects of climate alone
rm_clim <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr +
                      (1|SiteID) + (1|survYear),
                    family = nbinom2, data = rm_surv)
summary(rm_clim)
compute_vif(rm_clim)
performance::r2(rm_clim)

##### Effects of underlying environment alone (site characteristics)
rm_env <- glmmTMB(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = nbinom2, data = rm_surv)
summary(rm_env)
compute_vif(rm_env)
performance::r2(rm_env)

#### Model comparison
### Compare AICs
AIC(rm_abio_null, rm_abio, rm_clim, rm_env)
anova(rm_abio, rm_env) #full model integrating climate and environment is better

#########################################################
################# Species-Level - SM

##### Subset sm data
sm_surv <- subset(allSpDat_surv_clim, Sp == "sm")

#### Null model
sm_abio_null <- glmmTMB(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = nbinom2(), data = sm_surv)
summary(sm_abio_null)
performance::r2(sm_abio_null)

##### Effects of environment and climate
sm_abio <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = sm_surv)
summary(sm_abio)
compute_vif(sm_abio)
performance::r2(sm_abio)

##### Effects of climate alone
sm_clim <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = sm_surv)
summary(sm_clim)
compute_vif(sm_clim)
performance::r2(sm_clim)

##### Effects of underlying environment alone (site characteristics)
sm_env <- glmmTMB(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = nbinom2, data = sm_surv)
summary(sm_env)
compute_vif(sm_env)
performance::r2(sm_env)

#### Model comparison
### Compare AICs - model with environment only is better!
AIC(sm_abio_null, sm_abio, sm_clim, sm_env)
anova(sm_abio, sm_env)

#########################################################
################# Species-Level - bc

##### Subset bc data
bc_surv <- subset(allSpDat_surv_clim, Sp == "bc")

#### Null model
bc_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = bc_surv)
summary(bc_abio_null)
performance::r2(bc_abio_null)

##### Effects of environment and climate
bc_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + #maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_abio)
car::vif(bc_abio)
performance::r2(bc_abio)

##### Effects of climate alone
bc_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_clim)
car::vif(bc_clim)
performance::r2(bc_clim)

##### Effects of underlying environment alone (site characteristics)
bc_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = bc_surv)
summary(bc_env)
car::vif(bc_env)
performance::r2(bc_env)

#### Model comparison
### Compare AICs - model with environment only is better!
AIC(bc_abio_null, bc_abio, bc_clim, bc_env)
anova(bc_abio, bc_env)

#########################################################
################# Species-Level - wa

##### Subset wa data
wa_surv <- subset(allSpDat_surv_clim, Sp == "wa")

#### Null model
wa_abio_null <- glmmTMB(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = nbinom2, data = wa_surv)
summary(wa_abio_null)
performance::r2(wa_abio_null)

##### Effects of environment and climate
wa_abio <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = wa_surv)
summary(wa_abio)
compute_vif(wa_abio)
performance::r2(wa_abio)

##### Effects of climate alone
wa_clim <- glmmTMB(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = wa_surv)
summary(wa_clim)
compute_vif(wa_clim)
performance::r2(wa_clim)

##### Effects of underlying environment alone (site characteristics)
wa_env <- glmmTMB(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = nbinom2, data = wa_surv)
summary(wa_env)
compute_vif(wa_env)
performance::r2(wa_env)

#### Model comparison
### Compare AICs - full model integrating climate and environment is better
AIC(wa_abio_null, wa_abio, wa_clim, wa_env)
anova(wa_abio, wa_clim)

#########################################################
################# Species-Level - ab (No RE of site)

##### Subset ab data
ab_surv <- subset(allSpDat_surv_clim, Sp == "ab")

#### Null model
ab_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = ab_surv)
summary(ab_abio_null)
performance::r2(ab_abio_null)

##### Effects of environment and climate
ab_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_abio)
car::vif(ab_abio)
performance::r2(ab_abio)

##### Effects of climate alone
ab_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_clim)
car::vif(ab_clim)
performance::r2(ab_clim)

##### Effects of underlying environment alone (site characteristics)
ab_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|survYear),
                family = poisson(link = "log"), data = ab_surv)
summary(ab_env)
car::vif(ab_env)
performance::r2(ab_env)

#### Model comparison
### Compare AICs - full model integrating climate and environment is better
AIC(ab_abio_null, ab_abio, ab_clim, ab_env)
anova(ab_abio, ab_env)

#########################################################
################# Species-Level - nro (No RE of site)

##### Subset nro data
nro_surv <- subset(allSpDat_surv_clim, Sp == "nro")

#### Null model
nro_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = nro_surv)
summary(nro_abio_null)
performance::r2(nro_abio_null)

##### Effects of environment and climate
nro_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_abio)
car::vif(nro_abio)
performance::r2(nro_abio)

##### Effects of climate alone
nro_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_clim)
car::vif(nro_clim)
performance::r2(nro_clim)

##### Effects of underlying environment alone (site characteristics)
nro_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                 (1|survYear),
                family = poisson(link = "log"), data = nro_surv)
summary(nro_env)
car::vif(nro_env)
performance::r2(nro_env)

#### Model comparison
### Compare AICs - model with environment only is better
AIC(nro_abio_null, nro_abio, nro_clim, nro_env)
anova(nro_abio, nro_env)

#########################################################
################# Species-Level - bonro

##### Subset bonro data
bonro_surv <- subset(allSpDat_surv_clim, Sp == "bonro")

#### Null model
bonro_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = bonro_surv)
summary(bonro_abio_null)
performance::r2(bonro_abio_null)

##### Effects of environment and climate
bonro_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bonro_surv)
summary(bonro_abio)
car::vif(bonro_abio)
performance::r2(bonro_abio)

##### Effects of climate alone
bonro_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr + (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bonro_surv)
summary(bonro_clim)
car::vif(bonro_clim)
performance::r2(bonro_clim)

##### Effects of underlying environment alone (site characteristics)
bonro_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = bonro_surv)
summary(bonro_env)
car::vif(bonro_env)
performance::r2(bonro_env)

#### Model comparison
### Compare AICs - model with environment only is better (but environmental effects NS)
AIC(bonro_abio_null, bonro_abio, bonro_clim, bonro_env)
anova(bonro_clim, bonro_env)

#########################################################
################# Species-Level - wo - No RE of site

##### Subset wo data
wo_surv <- subset(allSpDat_surv_clim, Sp == "wo")

#### Null model
wo_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = wo_surv)
summary(wo_abio_null)
performance::r2(wo_abio_null)

##### Effects of environment and climate
wo_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_abio)
car::vif(wo_abio)
performance::r2(wo_abio)

##### Effects of climate alone
wo_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr + (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_clim)
car::vif(wo_clim)
performance::r2(wo_clim)

##### Effects of underlying environment alone (site characteristics)
wo_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                 (1|survYear),
                family = poisson(link = "log"), data = wo_surv)
summary(wo_env)
car::vif(wo_env)
performance::r2(wo_env)

#### Model comparison
### Compare AICs - model with environment only is better (but environmental effects NS)
AIC(wo_abio_null, wo_abio, wo_clim, wo_env)
anova(wo_abio, wo_env)

#########################################################
################# Species-Level - iw (No RE of Site)

##### Subset iw data
iw_surv <- subset(allSpDat_surv_clim, Sp == "iw")

#### Null model
iw_abio_null <- glmer(Year_1_counts ~ 1 + (1|survYear),
                      family = poisson(link = "log"), data = iw_surv)
summary(iw_abio_null)
performance::r2(iw_abio_null)

##### Effects of environment and climate
iw_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_abio)
car::vif(iw_abio)
performance::r2(iw_abio)

##### Effects of climate alone
iw_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_clim)
car::vif(iw_clim)
performance::r2(iw_clim)

##### Effects of underlying environment alone (site characteristics)
iw_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                 (1|survYear),
                family = poisson(link = "log"), data = iw_surv)
summary(iw_env)
car::vif(iw_env)
performance::r2(iw_env)

#### Model comparison
### Compare AICs - model with climate only is better
AIC(iw_abio_null, iw_abio, iw_clim, iw_env)
anova(iw_abio, iw_clim)


########## Plotting effect sizes for models with climate and environment
# --- Function to extract fixed effects and tidy them
extract_fixef_abio <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      # Convert model terms to nicer labels
      component = case_when(
        term == "mean_pdsi" ~ "Mean PDSI",
        term == "mintemp" ~ "Min Temperature",
        term == "maxtemp" ~ "Max Temperature",
        term == "percOpen" ~ "% Canopy Openness",
        term == "Base_CEC" ~ "Base CEC",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & 
                             p.value <= 0.05, "yes", "no"),
      Model = model_type)
}

## Extract for each species
fixef_RM  <- extract_fixef_abio(rm_abio, "rm")
fixef_SM  <- extract_fixef_abio(sm_abio, "sm")
fixef_BC  <- extract_fixef_abio(bc_abio, "bc")
fixef_WA  <- extract_fixef_abio(wa_abio, "wa")
fixef_AB  <- extract_fixef_abio(ab_abio, "ab")
fixef_NRO  <- extract_fixef_abio(nro_abio, "nro")
fixef_BONRO  <- extract_fixef_abio(bonro_abio, "bonro")
fixef_WO  <- extract_fixef_abio(wo_abio, "wo")
fixef_IW  <- extract_fixef_abio(iw_abio, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

#write.csv(fixef_df, "./R/L1/AbioticSpeciesMods.csv", row.names = FALSE)

# Define vertical order for plotting (top → bottom)
desired_order <- c("Mean PDSI", "Min Temperature", "Max Temperature",
                   "% Canopy Openness", "Base CEC",
                   "Log(Recruit Density)")

## Create a color vector with red shades for AM and blue shades for ECM 
species_cols <- c(
  # AM species first, alphabetically
  rm    = "#FB8A74",  # Acer rubrum
  sm    = "#E63E2D",  # Acer saccharum
  wa    = "#991F1F",  # Fraxinus americana
  bc    = "#FDCBC0",  # Prunus serotina
  # ECM species next, alphabetically
  ab    = "#C6E0EB",  # Fagus grandifolia
  iw    = "#4A90C2",  # Ostrya virginiana
  wo    = "#143D7A",  # Quercus alba
  nro   = "#1F5BA1",  # Quercus rubra
  bonro = "#7FB8D0"   # Quercus velutina
)

## Species labels (italicized scientific names, same order as colors)
species_labels <- c(
  rm    = expression(italic("Acer rubrum")),
  sm    = expression(italic("Acer saccharum")),
  wa    = expression(italic("Fraxinus americana")),
  bc    = expression(italic("Prunus serotina")),
  ab    = expression(italic("Fagus grandifolia")),
  iw    = expression(italic("Ostrya virginiana")),
  wo    = expression(italic("Quercus alba")),
  nro   = expression(italic("Quercus rubra")),
  bonro = expression(italic("Quercus velutina"))
)

## Create a size vector for species to scale points by mean percent ba 
species_sizes <- data.frame(species = c("ab", "bonro", "iw", "nro", 
                                        "wo", "bc", "rm", "sm", "wa"), 
                            percBA = c(2.43, 15.20, 0.11, 29.7, 10.30, 
                                       1.03, 9.40, 15.80, 1.82))

# Representative % BA values for legend
legend_percBA <- c(1, 10, 20)

# --- Prepare plotting data
fixef_df_plot <- fixef_df %>%
  mutate(component = factor(component, levels = desired_order),
         Model = factor(Model, levels = names(species_cols)),
         Significant = factor(Significant, levels = c("yes", "no")),
         Myc_type = case_when(Model %in% c("rm", "sm", "bc", "wa") ~ "AM",
                              TRUE ~ "ECM"),
         percBA = species_sizes$percBA[match(Model, species_sizes$species)],
         color_col = species_cols[Model],
         y_offset = case_when(
           Model == "rm" ~ 0.20,
           Model == "bc" ~ 0.15,
           Model == "sm" ~ 0.10,
           Model == "wa" ~ 0.05,
           Model == "ab" ~ 0,
           Model == "nro" ~ -0.05,
           Model == "bonro" ~ -0.10,
           Model == "wo" ~ -0.15,
           Model == "iw" ~ -0.20,
           TRUE ~ 0),
         y_pos = as.numeric(component) + y_offset,
         fill_col = ifelse(Significant == "yes", species_cols[Model], "white"))

################# Plot - all species
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model),
                width = 0.1, orientation = "y") +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species",
                     values = species_cols,
                     labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA") +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))


################# Plot - AM species only
##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

# Define AM species in alphabetical order for legend
am_species <- c("rm", "sm", "wa", "bc")
am_labels <- species_labels[am_species]
am_colors <- species_cols[am_species]

# Plot
ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model),
                width = 0.1, orientation = "y") +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species",
                     values = am_colors,
                     labels = am_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA") +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

################# Plot - ECM species only
##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

# Define ECM species in alphabetical order for legend
ecm_species <- c("ab", "iw", "wo", "nro", "bonro")
ecm_labels <- species_labels[ecm_species]
ecm_colors <- species_cols[ecm_species]

# Plot
ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model),
                width = 0.1, orientation = "y") +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species",
                     values = ecm_colors,
                     labels = ecm_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA") +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))


