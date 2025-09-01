
#################################################
####### Analyses of one-year survivorship WRT abiotic environment
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 18 June 2025
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

#### Null model
allSp_abio_null <- glmer(Year_1_counts ~ 1 +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_abio_null)
r.squaredGLMM(allSp_abio_null)

##### Effects of environment and climate
allSp_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
               family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_abio)
car::vif(allSp_abio)
r.squaredGLMM(allSp_abio)

##### Effects of climate alone
allSp_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_clim)
car::vif(allSp_clim)
r.squaredGLMM(allSp_clim)

##### Effects of underlying environment alone (site characteristics)
allSp_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_env)
car::vif(allSp_env)
r.squaredGLMM(allSp_env)

#### Model comparison
### Compare AICs - full model integrating climate and environment is better
AIC(allSp_abio_null, allSp_abio, allSp_clim, allSp_env)
anova(allSp_abio, allSp_clim)

##################################
################# Mycorrhizal type-level Models

## Subset for AM
AM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "AM")
## Subset for ECM
ECM_Dat_surv_clim = filter(allSpDat_surv_clim, Seedling_myc == "AM")

#### Null model - AM
AM_abio_null <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
              family = poisson(link = "log"), data = AM_Dat_surv_clim)
summary(AM_Dat_surv_clim)
r.squaredGLMM(AM_Dat_surv_clim)
#### Null model - ECM
ECM_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_Dat_surv_clim)
r.squaredGLMM(ECM_Dat_surv_clim)

##### Effects of environment and climate - AM
AM_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                family = poisson(link = "log"), data = AM_Dat_surv_clim)
summary(AM_abio)
car::vif(AM_abio)
r.squaredGLMM(AM_abio)
##### Effects of environment and climate - ECM
ECM_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_abio)
car::vif(ECM_abio)
r.squaredGLMM(ECM_abio)

##### Effects of climate alone - AM
AM_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
              family = poisson(link = "log"), data = AM_Dat_surv_clim)
summary(AM_clim)
car::vif(AM_clim)
r.squaredGLMM(AM_clim)
##### Effects of climate alone - ECM
ECM_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear) + (1|Sp),
                 family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_clim)
car::vif(ECM_clim)
r.squaredGLMM(ECM_clim)

##### Effects of environment alone (site characteristics) - AM
AM_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
            family = poisson(link = "log"), data = AM_Dat_surv_clim)
summary(AM_env)
car::vif(AM_env)
r.squaredGLMM(AM_env)
##### Effects of environment alone (site characteristics) - ECM
ECM_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear) + (1|Sp),
                family = poisson(link = "log"), data = ECM_Dat_surv_clim)
summary(ECM_env)
car::vif(ECM_env)
r.squaredGLMM(ECM_env)

#### Model comparison
### Compare AICs - AM
AIC(AM_abio_null, AM_abio, AM_clim, AM_env)
anova(AM_abio, AM_clim) # Model with climate and env best
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
  mutate(
    component = factor(component, levels = rev(desired_order)),
    Model = factor(Model, levels = c("All", "AM", "ECM")),
    Significant = factor(Significant, levels = c("yes", "no")),
    y_offset = case_when(
      Model == "All" ~ 0,
      Model == "AM"  ~ -0.2,
      Model == "ECM" ~  0.2), 
    y_pos = as.numeric(component) + y_offset,
    plot_est = ifelse(estimate > 1, estimate - 1, estimate),
    plot_xmin = ifelse(estimate > 1, estimate - 1, estimate) - std.error,
    plot_xmax = ifelse(estimate > 1, estimate - 1, estimate) + std.error)
  
# --- Plot
ggplot(fixef_df, aes(x = plot_est, y = y_pos)) +   # use numeric y_pos
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars
  geom_errorbarh(aes(xmin = plot_xmin, xmax = plot_xmax),
                 height = 0.05, color = "black", linewidth = 0.5)+
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
  # X-axis breaks
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
                     labels = c("-0.5", "0", "0.5", "2.0", "", "2.5")) +
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  # Fix fill legend shape
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 3.5)))

#########################################################
#########################################################
################# Species-Level - RM

##### Subset RM data
rm_surv <- subset(allSpDat_surv_clim, Sp == "rm")

#### Null model
rm_abio_null <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear),
                         family = poisson(link = "log"), data = rm_surv)
summary(rm_abio_null)
r.squaredGLMM(rm_abio_null)

##### Effects of environment and climate
rm_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_abio)
car::vif(rm_abio)
r.squaredGLMM(rm_abio)

##### Effects of climate alone
rm_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      log0yr +
                      (1|SiteID) + (1|survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_clim)
car::vif(rm_clim)
r.squaredGLMM(rm_clim)

##### Effects of underlying environment alone (site characteristics)
rm_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = rm_surv)
summary(rm_env)
car::vif(rm_env)
r.squaredGLMM(rm_env)

#### Model comparison
### Compare AICs - full model integrating climate and environment is better
AIC(rm_abio_null, rm_abio, rm_clim, rm_env)
anova(rm_abio, rm_env)

#########################################################
################# Species-Level - SM

##### Subset sm data
sm_surv <- subset(allSpDat_surv_clim, Sp == "sm")

#### Null model
sm_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = sm_surv)
summary(sm_abio_null)
r.squaredGLMM(sm_abio_null)

##### Effects of environment and climate
sm_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_abio)
car::vif(sm_abio)
r.squaredGLMM(sm_abio)

##### Effects of climate alone
sm_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_clim)
car::vif(sm_clim)
r.squaredGLMM(sm_clim)

##### Effects of underlying environment alone (site characteristics)
sm_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = sm_surv)
summary(sm_env)
car::vif(sm_env)
r.squaredGLMM(sm_env)

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
r.squaredGLMM(bc_abio_null)

##### Effects of environment and climate
bc_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + #maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_abio)
car::vif(bc_abio)
r.squaredGLMM(bc_abio)

##### Effects of climate alone
bc_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_clim)
car::vif(bc_clim)
r.squaredGLMM(bc_clim)

##### Effects of underlying environment alone (site characteristics)
bc_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = bc_surv)
summary(bc_env)
car::vif(bc_env)
r.squaredGLMM(bc_env)

#### Model comparison
### Compare AICs - model with environment only is better!
AIC(bc_abio_null, bc_abio, bc_clim, bc_env)
anova(bc_abio, bc_env)

#########################################################
################# Species-Level - wa

##### Subset wa data
wa_surv <- subset(allSpDat_surv_clim, Sp == "wa")

#### Null model
wa_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = wa_surv)
summary(wa_abio_null)
r.squaredGLMM(wa_abio_null)

##### Effects of environment and climate
wa_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_abio)
car::vif(wa_abio)
r.squaredGLMM(wa_abio)

##### Effects of climate alone
wa_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_clim)
car::vif(wa_clim)
r.squaredGLMM(wa_clim)

##### Effects of underlying environment alone (site characteristics)
wa_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = wa_surv)
summary(wa_env)
car::vif(wa_env)
r.squaredGLMM(wa_env)

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
r.squaredGLMM(ab_abio_null)

##### Effects of environment and climate
ab_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_abio)
car::vif(ab_abio)
r.squaredGLMM(ab_abio)

##### Effects of climate alone
ab_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_clim)
car::vif(ab_clim)
r.squaredGLMM(ab_clim)

##### Effects of underlying environment alone (site characteristics)
ab_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|survYear),
                family = poisson(link = "log"), data = ab_surv)
summary(ab_env)
car::vif(ab_env)
r.squaredGLMM(ab_env)

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
r.squaredGLMM(nro_abio_null)

##### Effects of environment and climate
nro_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_abio)
car::vif(nro_abio)
r.squaredGLMM(nro_abio)

##### Effects of climate alone
nro_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_clim)
car::vif(nro_clim)
r.squaredGLMM(nro_clim)

##### Effects of underlying environment alone (site characteristics)
nro_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                 (1|survYear),
                family = poisson(link = "log"), data = nro_surv)
summary(nro_env)
car::vif(nro_env)
r.squaredGLMM(nro_env)

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
r.squaredGLMM(bonro_abio_null)

##### Effects of environment and climate
bonro_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr + (1|survYear),
                 family = poisson(link = "log"), data = bonro_surv)
summary(bonro_abio)
car::vif(bonro_abio)
r.squaredGLMM(bonro_abio)

##### Effects of climate alone
bonro_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bonro_surv)
summary(bonro_clim)
car::vif(bonro_clim)
r.squaredGLMM(bonro_clim)

##### Effects of underlying environment alone (site characteristics)
bonro_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = bonro_surv)
summary(bonro_env)
car::vif(bonro_env)
r.squaredGLMM(bonro_env)

#### Model comparison
### Compare AICs - model with environment only is better (but environmental effects NS)
AIC(bonro_abio_null, bonro_abio, bonro_clim, bonro_env)
anova(bonro_clim, bonro_env)

#########################################################
################# Species-Level - wo

##### Subset wo data
wo_surv <- subset(allSpDat_surv_clim, Sp == "wo")

#### Null model
wo_abio_null <- glmer(Year_1_counts ~ 1 +
                        (1|SiteID) + (1|survYear),
                      family = poisson(link = "log"), data = wo_surv)
summary(wo_abio_null)
r.squaredGLMM(wo_abio_null)

##### Effects of environment and climate
wo_abio <- glmer(Year_1_counts ~ mean_pdsi + #mintemp + 
                   maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_abio)
car::vif(wo_abio)
r.squaredGLMM(wo_abio)

##### Effects of climate alone
wo_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_clim)
car::vif(wo_clim)
r.squaredGLMM(wo_clim)

##### Effects of underlying environment alone (site characteristics)
wo_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = wo_surv)
summary(wo_env)
car::vif(wo_env)
r.squaredGLMM(wo_env)

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
r.squaredGLMM(iw_abio_null)

##### Effects of environment and climate
iw_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_abio)
car::vif(iw_abio)
r.squaredGLMM(iw_abio)

##### Effects of climate alone
iw_clim <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_clim)
car::vif(iw_clim)
r.squaredGLMM(iw_clim)

##### Effects of underlying environment alone (site characteristics)
iw_env <- glmer(Year_1_counts ~ percOpen + Base_CEC + log0yr +
                 (1|survYear),
                family = poisson(link = "log"), data = iw_surv)
summary(iw_env)
car::vif(iw_env)
r.squaredGLMM(iw_env)

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

# Define vertical order for plotting (top → bottom)
desired_order <- c("Mean PDSI", "Min Temperature", "Max Temperature",
                   "% Canopy Openness", "Base CEC",
                   "Log(Recruit Density)")

## Create a color vector with red shades for AM and blue shades for ECM 
species_cols <- c(ab    = "#C6E0EB", 
                  bonro = "#7FB8D0",iw    = "#4A90C2",
                  nro   = "#1F5BA1",wo    = "#143D7A", 
                  bc    = "#FDCBC0",  rm    = "#FB8A74",
                  sm    = "#E63E2D", wa    = "#991F1F")

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
         # tidy up Model and Significant
         Model = factor(Model, levels = c("ab", "bonro", "iw",
                                          "nro", "wo", 
                                          "bc", "rm", "sm", "wa")),
         Significant = factor(Significant, levels = c("yes", "no")),
         # assign Myc_type based on Model
         Myc_type = case_when(Model %in% c("rm", "sm", "bc", "wa") ~ "AM",
                              TRUE ~ "ECM"),
         # percent BA per species
         percBA = species_sizes$percBA[match(Model, species_sizes$species)],
         # color for species
         color_col = species_cols[Model],
         # species-specific vertical offset
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
         # final y-position for plotting
         y_pos = as.numeric(component) + y_offset,
         # fill color: species color if significant, white otherwise
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"),
         plot_est = ifelse(estimate > 1, estimate - 1, estimate),
         plot_xmin = ifelse(estimate > 1, estimate - 1, estimate) - std.error,
         plot_xmax = ifelse(estimate > 1, estimate - 1, estimate) + std.error)

################# Plot - all species
ggplot(fixef_df_plot, aes(x = plot_est, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = plot_xmin, xmax = plot_xmax, color = Model),
                 height = 0.05, linewidth = 0.5)+
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  # Shape legend for mycorrhizal type
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22) ) +
  # Fill is used only for significance; no separate legend
  scale_fill_identity(guide = "none") +
  # Species colors and legend
  scale_color_manual(
    name = "Species", values = species_cols) +
  # Size scale with interpretable legend
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),          # positions in plot scale
             labels = paste0(legend_percBA, "%"),   # human-readable labels
             name = "% BA" ) +
  # X-axis breaks
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
                     labels = c("-0.5", "0", "0.5", "2.0", "", "2.5")) +
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))


################# Plot - AM species only
##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")
# Plot
ggplot(fixef_df_AM, aes(x = plot_est, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = plot_xmin, xmax = plot_xmax, color = Model),
                 height = 0.05, linewidth = 0.5)+
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  # Shape legend for mycorrhizal type
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22) ) +
  # Fill is used only for significance; no separate legend
  scale_fill_identity(guide = "none") +
  # Species colors and legend
  scale_color_manual(
    name = "Species", values = species_cols) +
  # Size scale with interpretable legend
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),          # positions in plot scale
             labels = paste0(legend_percBA, "%"),   # human-readable labels
             name = "% BA" ) +
  # X-axis breaks
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
                     labels = c("-0.5", "0", "0.5", "2.0", "", "2.5")) +
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

################# Plot - ECM species only
##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")
# Plot
ggplot(fixef_df_ECM, aes(x = plot_est, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = plot_xmin, xmax = plot_xmax, color = Model),
                 height = 0.05, linewidth = 0.5)+
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  # Shape legend for mycorrhizal type
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22) ) +
  # Fill is used only for significance; no separate legend
  scale_fill_identity(guide = "none") +
  # Species colors and legend
  scale_color_manual(
    name = "Species", values = species_cols) +
  # Size scale with interpretable legend
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),          # positions in plot scale
             labels = paste0(legend_percBA, "%"),   # human-readable labels
             name = "% BA" ) +
  # X-axis breaks
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
                     labels = c("-0.5", "0", "0.5", "2.0", "", "2.5")) +
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

