 
#################################################
####### One-year seedling survival WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 10 November 2026
#################################################

## Read in packages
library(tidyverse) 
library(lme4) 
library(emmeans)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggeffects)
library(DHARMa)
library(broom)
library(broom.mixed)
library(tibble)
library(forcats)
library(MuMIn)
library(glmmTMB)
library(performance)
library(car)
library(marginaleffects)
library(patchwork)

#### Read in csv with all survivorship data for all species
allSpData <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")

allSpData$Seedling_myc <- factor(allSpData$Seedling_myc, levels = c("AM", "ECM"))

allSpDat <- subset(allSpData, !is.na(propHetAM)) # 35,791 rows

allSpDatSub <- subset(allSpData, is.na(propHetAM)) # 602 rows

###  Helper function: VIF for glmmTMB fixed effects
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

#########################################################
#########################################################
################# All Species Models - Proportional BA

######### ######### ######### ######### ######### 
######### Null model (seedling density only)
all_sp_prop_null <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                              log0yr + log(Prior_3_year_cohorts_counts + 1) +
                              (1|SiteID) + (1|survYear) + (1|Sp),
                            family = betabinomial, data = allSpDat)
summary(all_sp_15m_null)

######### ######### ######### ######### ######### 
######### Model A — Conspecific Density Dependence (CDD)
all_sp_propA <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                          propConsp*Seedling_myc +
                          log0yr + log(Prior_3_year_cohorts_counts + 1) +
                          (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                        family = betabinomial, data = allSpDat) 
summary(all_sp_propA)
## Significant negative effect of propConsp
# VIF
vif_propA <- compute_vif( ~ propConsp*Seedling_myc +
                            log0yr + 
                            log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_propA # all <4

######### ######### ######### ######### ######### 
######### Model B — Mycorrhizal-Dependent Pathogen Spillover 
all_sp_propB <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                          propConsp*Seedling_myc + 
                          propHetAM + propConsp:propHetAM + 
                          log0yr + log(Prior_3_year_cohorts_counts + 1) +
                          (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                        family = betabinomial, data = allSpDat)
summary(all_sp_propB)
##### Significant negative effects of propConsp and propAM
# VIF
vif_propB <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            log0yr + 
                            log(Prior_3_year_cohorts_counts + 1) , allSpDat)
vif_propB # All <=4

######### ######### ######### ######### ######### 
######### Model C — Conmycorrhizal Density Dependence
all_sp_propC <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                          propConsp*Seedling_myc + 
                          propHetAM + 
                          Seedling_myc:propHetAM + 
                          log0yr + log(Prior_3_year_cohorts_counts + 1) +
                          (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                        family = betabinomial, data = allSpDat)
summary(all_sp_propC)
## Significant negative effects of propConsp and propAM
# VIF
vif_propC <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + Seedling_myc:propHetAM + 
                            log0yr + 
                            log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_propC # All <5

######### ######### ######### ######### ######### 
######### Model D — Conmycorrhizal plus pathogen spillover
all_sp_propD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                          propConsp*Seedling_myc + 
                          propHetAM + 
                          propConsp:propHetAM + 
                          Seedling_myc:propHetAM +
                          log0yr + log(Prior_3_year_cohorts_counts + 1) +
                          (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                        family = betabinomial, data = allSpDat)
summary(all_sp_propD)
## Significant negative effects of propConsp and propAM
# VIF
vif_propD <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            Seedling_myc:propHetAM + 
                            log0yr + 
                            log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_propD # All <=6

######### ######### ######### ######### ######### 
######### Model E - Interactive conmycorrhizal and pathogen spillover
all_sp_propBCD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                            propConsp*Seedling_myc + 
                            propHetAM + 
                            propConsp:propHetAM + 
                            Seedling_myc:propHetAM + 
                            propConsp:Seedling_myc:propHetAM + 
                            log0yr + log(Prior_3_year_cohorts_counts + 1) +
                            (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                          family = betabinomial, data = allSpDat)
summary(all_sp_propBCD)
## Significant negative effects of propConsp and propAM 
## propConsp:Seedling_mycECM interaction sign. and positive
## All other interactions not significant
##### VIF 
vif_propBCD <- compute_vif( ~ propConsp*Seedling_myc + 
                              propHetAM + propConsp:propHetAM + 
                              Seedling_myc:propHetAM + 
                              propConsp:Seedling_myc:propHetAM + 
                              log0yr + 
                              log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_propBCD # four values >10 due to interactions

## Model comparison
### Compute Delta AICs
aic_table <- AIC(all_sp_prop_null, all_sp_propA, 
                 all_sp_propB, all_sp_propC,
                 all_sp_propD, all_sp_propBCD)
aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]
anova(all_sp_propC, all_sp_propB) # But, NOT significantly better than B


#########################################################
#########################################################
#########################################################
##### All Species Models - Absolute BA

######### ######### ######### ######### ######### 
######### Null model (seedling density only)
all_sp_15m_null <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                             log0yr + log(Prior_3_year_cohorts_counts + 1) +
                             (1|SiteID) + (1|survYear) + (1|Sp),
                           family = betabinomial, data = allSpDat)
summary(all_sp_15m_null)

######### ######### ######### ######### ######### 
######### Model A — Conspecific Density Dependence (CDD)
all_sp_15mA <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                         consp_15m * Seedling_myc + 
                         log0yr + log(Prior_3_year_cohorts_counts + 1) +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mA)
## Significant negative effect of conspBA
# VIF 
vif_A <- compute_vif(~ consp_15m * Seedling_myc +
                       log0yr + log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_A # All <4

######### ######### ######### ######### ######### 
######### Model B — Mycorrhizal-Dependent Pathogen Spillover 
all_sp_15mB <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                         consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m +
                         consp_15m:hetECM_15m +
                         log0yr + log(Prior_3_year_cohorts_counts + 1) +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mB)

### Sensitivity test for best supported model
all_sp_15mB2 <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                         consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m +
                         consp_15m:hetECM_15m +
                         log0yr + log(Prior_3_year_cohorts_counts + 1) +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpData)
summary(all_sp_15mB2)
## Significant positive effect of hetECM_15m
## Significant negative consp_15m:hetAM_15m and consp_15m:hetECM_15m interactions 
# VIF 
vif_B <- compute_vif(~ consp_15m * Seedling_myc +
                       hetAM_15m + hetECM_15m +
                       consp_15m:hetAM_15m + consp_15m:hetECM_15m +
                       log0yr + log(Prior_3_year_cohorts_counts + 1),  allSpDat)
vif_B # All < 8 (or < 10 in reduced)

### Species-specific CDD slopes (logit scale) — Model B
fe      <- fixef(all_sp_15mB)$cond
sp_coef <- coef(all_sp_15mB)$cond$Sp
sp_info <- unique(allSpDat[, c("Sp", "Seedling_myc")])
sp_info$Sp <- as.character(sp_info$Sp)

median_hetAM  <- median(allSpDat$hetAM_15m,  na.rm = TRUE)
median_hetECM <- median(allSpDat$hetECM_15m, na.rm = TRUE)

## Slope for species s at given heterospecific AM and ECM BA
sp_slope <- function(s, ecm, hAM, hECM) {
  unname(sp_coef[s, "consp_15m"] +
           ecm  * fe["consp_15m:Seedling_mycECM"] +
           hAM  * fe["consp_15m:hetAM_15m"] +
           hECM * fe["consp_15m:hetECM_15m"])
}

sp_slopes <- do.call(rbind, lapply(seq_len(nrow(sp_info)), function(i) {
  s   <- sp_info$Sp[i]
  ecm <- sp_info$Seedling_myc[i] == "ECM"
  data.frame(Sp = s, 
             Seedling_myc = sp_info$Seedling_myc[i],
             slope_zeroHet = sp_slope(s, ecm, 0, 0),
             slope_medHet  = sp_slope(s, ecm, 
                             median_hetAM, median_hetECM))
}))
sp_slopes[order(sp_slopes$slope_medHet), ]

### Marginal and conditional predictive R² (survival-proportion scale)
mf   <- model.frame(all_sp_15mB)        # rows actually used in the fit
resp <- mf[[1]]                         # the cbind() matrix
n0   <- rowSums(resp)                   # cohort size
obs_prop <- resp[, 1] / n0              # observed survival proportion

pred_marg <- predict(all_sp_15mB, type = "response", re.form = NA)
pred_cond <- predict(all_sp_15mB, type = "response")

R2_pred_marg <- cov.wt(cbind(obs_prop, pred_marg), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_cond <- cov.wt(cbind(obs_prop, pred_cond), 
                       wt = n0, cor = TRUE)$cor[1, 2]^2
R2_pred_marg
R2_pred_cond

######### ######### ######### ######### ######### 
######### Model C — Conmycorrhizal Density Dependence
all_sp_15mC <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                         consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         Seedling_myc:hetAM_15m +
                         Seedling_myc:hetECM_15m +
                         log0yr + log(Prior_3_year_cohorts_counts + 1) +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mC)
## Significant negative effects of conspBA
## Significant positive effect of hetECM_15m
# VIF 
vif_C <- compute_vif( ~ consp_15m * Seedling_myc +
                        hetAM_15m + hetECM_15m +
                        Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
                        log0yr + log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_C # One term with value >10

######### ######### ######### ######### ######### 
######### Model D — Conmycorrhizal plus pathogen spillover
all_sp_15mD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
                         consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m +
                         consp_15m:hetECM_15m +
                         Seedling_myc:hetAM_15m +
                         Seedling_myc:hetECM_15m +
                         log0yr + log(Prior_3_year_cohorts_counts + 1) +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mD)
## Significant positive effect of hetECM_15m
## Significant negative consp_15m:hetAM_15m and consp_15m:hetECM_15m interactions
# VIF 
vif_D <- compute_vif( ~ consp_15m * Seedling_myc +
                        hetAM_15m + hetECM_15m +
                        consp_15m:hetAM_15m + consp_15m:hetECM_15m +
                        Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
                        log0yr + log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_D # one >10

######### ######### ######### ######### ######### 
######### Model E - Interactive conmycorrhizal and pathogen spillover
all_sp_15mBCD <- glmmTMB(
  cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 
    consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    consp_15m:Seedling_myc:hetAM_15m +
    consp_15m:Seedling_myc:hetECM_15m +
    log0yr + log(Prior_3_year_cohorts_counts + 1) +
    (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
  family = betabinomial, data = allSpDat)
summary(all_sp_15mBCD)
## Significant positive main effect of hetECM BA
## Significant negative interaction between consp:hetECM

# VIF and R2
vif_BCD <- compute_vif(~ consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m + consp_15m:hetECM_15m +
                         Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
                         consp_15m:Seedling_myc:hetAM_15m +
                         consp_15m:Seedling_myc:hetECM_15m +
                         log0yr + log(Prior_3_year_cohorts_counts + 1), allSpDat)
vif_BCD # 4 VIFs >10; interaction-induced multicollinearity.

##############################
### Compute Delta AICs
aic_table <- AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
                 all_sp_15mC, all_sp_15mD, all_sp_15mBCD)

aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]
anova(all_sp_15mB, all_sp_15mD)
anova(all_sp_15mB, all_sp_15mBCD)

### Compare Delta AICs - raw BA AND proportional BA
aic_table <- AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
                 all_sp_15mC, all_sp_15mD, all_sp_15mBCD,
                 all_sp_prop_null, all_sp_propA, all_sp_propB,
                 all_sp_propC,all_sp_propD, all_sp_propBCD)
aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]

##### Diagnostics and Post-hoc tests - Model B
simB <- DHARMa::simulateResiduals(all_sp_15mB, n = 1000)
### Plots
DHARMa::plotQQunif(simB) # QQ plot + uniformity/dispersion/outlier tests
DHARMa::testUniformity(simB)
DHARMa::testDispersion(simB)
DHARMa::testOutliers(simB)
DHARMa::testZeroInflation(simB)

################################################## 
################################################## 
################################################## 
##### Plot predictions for Model B

# =========================================================
# Central values and quartiles
# =========================================================
median_consp  <- median(allSpDat$consp_15m,  na.rm = TRUE)
median_hetAM  <- median(allSpDat$hetAM_15m,  na.rm = TRUE)
median_hetECM <- median(allSpDat$hetECM_15m, na.rm = TRUE)
median_log0yr <- median(allSpDat$log0yr,   na.rm = TRUE)
median_Prior <- median(allSpDat$Prior_3_year_cohorts_counts, na.rm = TRUE)

hetAM_quart  <- quantile(allSpDat$hetAM_15m,
                         probs = c(.1,.5,.9), na.rm = TRUE)
hetECM_quart <- quantile(allSpDat$hetECM_15m,
                         probs = c(.1,.5,.9), na.rm = TRUE)

cond_vals <- c(Prior_3_year_cohorts_counts = median_Prior,
               log0yr = median_log0yr)

# =========================================================
## -------- PANEL A — Conspecific Density Dependence
pred_A <- ggpredict(all_sp_15mB,
                    terms = c("consp_15m", "Seedling_myc",
                              paste0("hetAM_15m[",  median_hetAM,  "]"),
                              paste0("hetECM_15m[", median_hetECM, "]")),
                    condition = cond_vals,
                    bias_correction = FALSE)

## Build dataframe
df_A <- data.frame(x = pred_A$x,
                   predicted = pred_A$predicted,
                   conf.low  = pred_A$conf.low,
                   conf.high = pred_A$conf.high,
                   seed_type = factor(pred_A$group, 
                               levels = c("AM","ECM")))

### Create combo variable 
df_A$combo <- factor(df_A$seed_type,levels = c("AM", "ECM"),
                     labels = c("AM_Medium","ECM_Medium")) 

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, 
                     color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0,0.8)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme

plot_A

# =========================================================
## ----PANEL B1 ------ Mycorrhizal-Dependent Pathogen Spillover 
### Consp × hetAM 
## Predict 
pred_B_AM <- ggpredict(all_sp_15mB,
              terms = c("consp_15m",
              paste0("hetAM_15m[", paste(hetAM_quart,
                    collapse = ","), "]"), "Seedling_myc",
                    paste0("hetECM_15m[", median_hetECM, "]")),
              condition = cond_vals,
                       bias_correction = FALSE)

## Build dataframe cleanly
df_B_AM <- data.frame(x = pred_B_AM$x,
                      predicted = pred_B_AM$predicted,
                      conf.low  = pred_B_AM$conf.low,
                      conf.high = pred_B_AM$conf.high,
                      seed_type = factor(pred_B_AM$facet, 
                                         levels = c("AM","ECM")),
                      hetAM_raw = pred_B_AM$group)

## Assign intuitive labels
df_B_AM$hetAM <- factor(df_B_AM$hetAM_raw,
                        levels = unique(pred_B_AM$group),
                        labels = c("Low","Medium","High"))

## Interaction label for color: AM_Low, AM_Medium, etc.
df_B_AM$combo <- interaction(df_B_AM$seed_type, df_B_AM$hetAM, sep="_")
df_B_AM$combo <- factor(df_B_AM$combo, levels = desired_order)

# =========================================================
## ----PANEL B2 ------ Mycorrhizal-Dependent Pathogen Spillover 
### Consp × hetECM 
## Predict 
pred_B_ECM <- ggpredict(all_sp_15mB,
              terms = c("consp_15m", paste0("hetECM_15m[", 
              paste(hetECM_quart, collapse = ","), "]"),
              "Seedling_myc", paste0("hetAM_15m[", median_hetAM, "]")),
              condition = cond_vals,
                        bias_correction = FALSE)

## Build dataframe cleanly
df_B_ECM <- data.frame(x = pred_B_ECM$x,
                       predicted = pred_B_ECM$predicted,
                       conf.low  = pred_B_ECM$conf.low,
                       conf.high = pred_B_ECM$conf.high,
                       seed_type = factor(pred_B_ECM$facet, 
                                          levels = c("AM","ECM")),
                       hetECM_raw = pred_B_ECM$group)

## Assign intuitive labels
df_B_ECM$hetECM <- factor(df_B_ECM$hetECM_raw,
                          levels = unique(pred_B_ECM$group),
                          labels = c("Low","Medium","High"))

## Interaction label for color: AM_Low, AM_Medium, etc.
df_B_ECM$combo <- interaction(df_B_ECM$seed_type, 
                              df_B_ECM$hetECM, sep="_")
df_B_ECM$combo <- factor(df_B_ECM$combo, levels = desired_order)

######### Combined plot
df_all <- bind_rows(df_B_AM %>% mutate(neighbor_type = "AM neighbors"),
                    df_B_ECM %>% mutate(neighbor_type = "ECM neighbors"))

df_all$seed_type <- factor(df_all$seed_type, levels = c("AM", "ECM"))
df_all$neighbor_type <- factor(df_all$neighbor_type,
                               levels = c("AM neighbors", "ECM neighbors"))

### Plot
plot_B <- ggplot(df_all, aes(x = x, y = predicted,
                   color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(seed_type ~ neighbor_type,
             labeller = labeller(seed_type = c(AM = "AM seedlings", 
                                               ECM = "ECM seedlings"))) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0,0.85)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme
plot_B

#####################################
##### Species-specific CDD - Model B
fe      <- fixef(all_sp_15mB)$cond
sp_coef <- coef(all_sp_15mB)$cond$Sp
sp_info <- unique(allSpDat[, c("Sp", "Seedling_myc")])
sp_info$Sp <- as.character(sp_info$Sp)

## Species-specific slope and intercept (logit scale) at given het BA
sp_terms <- function(s, ecm, hAM, hECM) {
  slope <- sp_coef[s, "consp_15m"] +
    ecm  * fe["consp_15m:Seedling_mycECM"] +
    hAM  * fe["consp_15m:hetAM_15m"] +
    hECM * fe["consp_15m:hetECM_15m"]
  int <- sp_coef[s, "(Intercept)"] +
    ecm  * fe["Seedling_mycECM"] +
    hAM  * fe["hetAM_15m"] +
    hECM * fe["hetECM_15m"] +
    fe["log(Prior_3_year_cohorts_counts + 1)"] * log(median_Prior + 1) +
    fe["log0yr"] * median_log0yr
  c(slope = unname(slope), int = unname(int))
}

## Slope table: at zero and at median heterospecific BA
sp_slopes <- do.call(rbind, lapply(seq_len(nrow(sp_info)), function(i) {
  s   <- sp_info$Sp[i]
  ecm <- sp_info$Seedling_myc[i] == "ECM"
  data.frame(Sp = s, Seedling_myc = sp_info$Seedling_myc[i],
             slope_zeroHet = sp_terms(s, ecm, 0, 0)[["slope"]],
             slope_medHet  = sp_terms(s, ecm, median_hetAM, median_hetECM)[["slope"]])
}))
sp_slopes[order(sp_slopes$slope_medHet), ]

## Curves at median heterospecific BA, over each species' observed range
sp_curves <- do.call(rbind, lapply(seq_len(nrow(sp_info)), function(i) {
  s   <- sp_info$Sp[i]
  ecm <- sp_info$Seedling_myc[i] == "ECM"
  b   <- sp_terms(s, ecm, median_hetAM, median_hetECM)
  x   <- seq(min(allSpDat$consp_15m[allSpDat$Sp == s], na.rm = TRUE),
             max(allSpDat$consp_15m[allSpDat$Sp == s], na.rm = TRUE),
             length.out = 100)
  data.frame(Sp = s, myc = ifelse(ecm, "ECM", "AM"), x = x,
             predicted = plogis(b[["int"]] + b[["slope"]] * x))
}))

sp_ends <- do.call(rbind, lapply(split(sp_curves, sp_curves$Sp),
                                 function(d) d[which.max(d$x), ]))

## 4-letter species labels (first two letters of genus + species)
sp_labels <- c(rm = "Acru", sm = "Acsa", wa = "Fram", bc = "Prse",
               ab = "Fagr", iw = "Osvi", wo = "Qual", nro = "Quru",
               bonro = "Quve", bs = "Tiam")

stopifnot(all(sp_curves$Sp %in% names(sp_labels)))   # every species has a label

sp_ends$label <- sp_labels[sp_ends$Sp]

plot_sp <- ggplot(sp_curves, aes(x, predicted, group = Sp, color = myc)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = sp_ends, aes(label = label), hjust = -0.15,
            size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(AM = col_AM, ECM = col_ECM),
                     name = "Seedling type") +
  scale_y_continuous(limits = c(0,0.8)) +
  scale_x_continuous(limits = c(0,4.1)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme +
  theme(legend.position = "none") 
plot_sp
