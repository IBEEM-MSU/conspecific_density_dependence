
#################################################
####### One-year seedling survival WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 23 July 2026
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
#allSpDat <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")
allSpDat <- read.csv("final_seedling_demogr_modeling_data.csv")

allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

### Back-transform log prior seedling densities 
allSpDat$Prior_counts <- 10^allSpDat$logPrior

#########################################################
#########################################################
################# All Species Models - Proportional Basal Area

######### ######### ######### ######### ######### 
######### Null model (Intercept only)
all_sp_prop_null <- 
  glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = betabinomial, data = allSpDat) 
summary(all_sp_prop_null)

######### ######### ######### ######### ######### 
######### Model A — Conspecific Density Dependence (CDD)
all_sp_propA <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                      propConsp*Seedling_myc +
                      Prior_counts + log0yr +
                      (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                    family = betabinomial, data = allSpDat) 
summary(all_sp_propA)
## Significant negative effect of propConsp
# VIF
vif_propA <- compute_vif( ~ propConsp*Seedling_myc +
                            Prior_counts + log0yr, allSpDat)
vif_propA # all <4

######### ######### ######### ######### ######### 
######### Model B — Mycorrhizal-Dependent Pathogen Spillover 
all_sp_propB <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                       propConsp*Seedling_myc + 
                       propHetAM + propConsp:propHetAM + 
                         Prior_counts + log0yr +
                       (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                     family = betabinomial, data = allSpDat)
summary(all_sp_propB)
## Significant negative effects of propConsp and propAM
# VIF and R2
vif_propB <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            Prior_counts + log0yr, allSpDat)
vif_propB # All <4

### Adjust so species-specific slopes are also correct for the ECM spp
fe <- fixef(all_sp_propB)$cond

sp_slopes <- coef(all_sp_propB)$cond$Sp[, c("(Intercept)", "propConsp")]
sp_slopes$Sp <- rownames(sp_slopes)
sp_slopes <- merge(sp_slopes, 
                   unique(allSpDat[, c("Sp", "Seedling_myc")]), by = "Sp")

sp_slopes$slope <- sp_slopes$propConsp +
  (sp_slopes$Seedling_myc == "ECM") * fe["propConsp:Seedling_mycECM"]

sp_slopes[order(sp_slopes$slope), c("Sp", "Seedling_myc", "slope")]

# Marginal prediction R² (log-scale)
obs <- allSpDat$Year_1_counts
pred_marg <- predict(all_sp_propB, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg # 0.081
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_propB, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond # 0.082

######### ######### ######### ######### ######### 
######### Model C — Conmycorrhizal Density Dependence
all_sp_propC <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                        propConsp*Seedling_myc + 
                        propHetAM + 
                        Seedling_myc:propHetAM + 
                        Prior_counts + log0yr +
                        (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                      family = betabinomial, data = allSpDat)
summary(all_sp_propC)
## Significant negative effects of propConsp and propAM
# VIF
vif_propC <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + Seedling_myc:propHetAM + 
                            Prior_counts + log0yr, allSpDat)
vif_propC # All <5

######### ######### ######### ######### ######### 
######### Model D — Conmycorrhizal plus pathogen spillover
all_sp_propD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                        propConsp*Seedling_myc + 
                        propHetAM + 
                        propConsp:propHetAM + 
                        Seedling_myc:propHetAM +
                        Prior_counts + log0yr +
                        (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                      family = betabinomial, data = allSpDat)
summary(all_sp_propD)
## Significant negative effects of propConsp and propAM
# VIF
vif_propD <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            Seedling_myc:propHetAM + 
                            Prior_counts + log0yr, allSpDat)
vif_propD # All <6

######### ######### ######### ######### ######### 
######### Model E - Interactive conmycorrhizal and pathogen spillover
all_sp_propBCD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~
                       propConsp*Seedling_myc + 
                       propHetAM + 
                       propConsp:propHetAM + 
                       Seedling_myc:propHetAM + 
                       propConsp:Seedling_myc:propHetAM + 
                         Prior_counts + log0yr +
                       (1|SiteID) + (1|survYear) + (1 + propConsp|Sp),
                     family = betabinomial, data = allSpDat)
summary(all_sp_propBCD)
## Significant negative effects of propConsp and propAM 
## propConsp:propAMm and propConsp:Seedling_mycECM interactions sign. and positive
## propConsp:Seedling_mycECM:propHetAM significant and negative
## All other interactions not significant
# VIF and R2
vif_propBCD <- compute_vif( ~ propConsp*Seedling_myc + 
                              propHetAM + propConsp:propHetAM + 
                              Seedling_myc:propHetAM + 
                              propConsp:Seedling_myc:propHetAM + 
                              Prior_counts + log0yr, allSpDat)
vif_propBCD # four values >10 due to interactions

### Adjust so species-specific slopes are also correct for the ECM spp
fe <- fixef(all_sp_propBCD)$cond

sp_slopes <- coef(all_sp_propBCD)$cond$Sp[, c("(Intercept)", "propConsp")]
sp_slopes$Sp <- rownames(sp_slopes)
sp_slopes <- merge(sp_slopes, 
                   unique(allSpDat[, c("Sp", "Seedling_myc")]), by = "Sp")

sp_slopes$slope <- sp_slopes$propConsp +
  (sp_slopes$Seedling_myc == "ECM") * fe["propConsp:Seedling_mycECM"]

sp_slopes[order(sp_slopes$slope), c("Sp", "Seedling_myc", "slope")]

# Marginal prediction R² (log-scale)
obs <- allSpDat$Year_1_counts
pred_marg <- predict(all_sp_propBCD, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg # 0.072
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_propBCD, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond # 0.082

## Model comparison
### Compare AICs 
AIC(all_sp_prop_null, all_sp_propA, all_sp_propB, all_sp_propC,
    all_sp_propD, all_sp_propBCD) # lowest is all_sp_propBCD 
anova(all_sp_propB, all_sp_propBCD) # But, NOT significantly better than B
### Delta AICs - raw BA AND proportional BA
aic_table <- AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
                 all_sp_15mC, all_sp_15mD, all_sp_15mBCD,
                 all_sp_prop_null, all_sp_propA, all_sp_propB,
                 all_sp_propC,all_sp_propD, all_sp_propBCD)
aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Diagnostics and Post-hoc tests for Panel E (Best Supported)
simBCDprop <- DHARMa::simulateResiduals(all_sp_propBCD, n = 1000)
DHARMa::plotQQunif(simBCDprop)      # left panel: QQ plot with tests
DHARMa::plotResiduals(simBCDprop)   # right panel: residuals vs. predicted
DHARMa::testDispersion(simBCDprop) # No overdispersion
DHARMa::testZeroInflation(simBCDprop) # No zero-inflation
DHARMa::testUniformity(simBCDprop) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCDprop) # ~0.18% of observations are outliers
# Model diagnostics indicated excellent fit, with no evidence of overdispersion or zero inflation. Although DHARMa’s uniformity and outlier tests showed significant deviations, the magnitude was trivial (<0.18% outliers).

################################################## 
################################################## 
##### Plot predictions for Proportional Model E (Most Supported)

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A and B
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
median_propConsp  <- median(allSpDat$propConsp,  na.rm = TRUE)
median_propHetAM  <- median(allSpDat$propHetAM,  na.rm = TRUE)
median_log0yr   <- median(allSpDat$log0yr,   na.rm = TRUE)
median_Prior <- median(allSpDat$Prior_counts, na.rm = TRUE)

propHetAM_quart  <- quantile(allSpDat$propHetAM,  
                             probs = c(.1,.5,.9), na.rm = TRUE)

### Colors
col_AM  <- "#c0352b"   # AM red
col_ECM <- "#2b4ac0"   # ECM blue

reds  <- colorRampPalette(c("#f7b4ad", col_AM,  "#7a1b14"))
blues <- colorRampPalette(c("#b3c4f7", col_ECM, "#16245c"))

# 3-level gradients for context levels (low, medium, high)
AM_grad  <- reds(3)
ECM_grad <- blues(3)

# Combined color lookup – used in all Panel B/D interaction plots
color_map_combo <- c(
  # AM seedling panel
  "AM_Low"    = AM_grad[1],
  "AM_Medium" = AM_grad[2],
  "AM_High"   = AM_grad[3],
  # ECM seedling panel
  "ECM_Low"   = ECM_grad[1],
  "ECM_Medium"= ECM_grad[2],
  "ECM_High"  = ECM_grad[3])

desired_order <- c("AM_Low","AM_Medium","AM_High",
                   "ECM_Low","ECM_Medium","ECM_High")

### Theme
base_theme <- theme_bw() +
  theme(legend.position = "right",
        axis.text  = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size=15, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# =========================================================
## -------- PANEL A — Conspecific Density Dependence
## Predict 
pred_A <- ggpredict(all_sp_propBCD,
                    terms = c("propConsp", "Seedling_myc",
                              paste0("propHetAM[",  median_propHetAM,  "]")),
                    condition = c(Prior_counts = median_Prior,
                                  log0yr   = median_log0yr),
                    bias_correction = FALSE)

## Build dataframe cleanly
df_A <- data.frame(
  x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high,
  seed_type = factor(pred_A$group, levels = c("AM","ECM")))

### Create combo variable (even though no Low/Med/High here)
df_A$combo <- df_A$seed_type
df_A$combo <- factor(df_A$combo,
                     levels = c("AM", "ECM"),
                     labels = c("AM_Medium","ECM_Medium"))  # middle level color

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme

plot_A

# =========================================================
## ---- PANEL B ------ Mycorrhizal-Dependent Pathogen Spillover 
## Predict 
pred_B_AM <- ggpredict(all_sp_propBCD,
            terms = c("propConsp",
            paste0("propHetAM[", paste(propHetAM_quart, 
            collapse = ","), "]"),"Seedling_myc"),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr   = median_log0yr),
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

## Plot
plot_B_AM <- ggplot(df_B_AM,
                    aes(x = x, y = predicted,
                        color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ seed_type, ncol = 2, 
             labeller = labeller(seed_type = c(
               AM  = "AM seedlings", ECM = "ECM seedlings"))) +
  scale_color_manual(values = color_map_combo, name = "Prop. AM BA",
                     guide = "none") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme
plot_B_AM

# =========================================================
## ---- PANEL C ------ Conmycorrhizal Density Dependence (prop)
pred_C_AM <- ggpredict(all_sp_propBCD,
                       terms = c("propHetAM", "Seedling_myc",
                                 paste0("propConsp[", median_propConsp, "]")),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr       = median_log0yr),
                       bias_correction = FALSE)

df_C_AM <- data.frame(x         = pred_C_AM$x,
                      predicted = pred_C_AM$predicted,
                      conf.low  = pred_C_AM$conf.low,
                      conf.high = pred_C_AM$conf.high,
                      seed_type = factor(pred_C_AM$group, levels = c("AM","ECM")))

## Drop impossible combinations
df_C_AM <- subset(df_C_AM, x <= 1 - median_propConsp)

df_C_AM$combo <- factor(df_C_AM$seed_type,
                        levels = c("AM","ECM"),
                        labels = c("AM_Medium","ECM_Medium"))

plot_C_AM <- ggplot(df_C_AM,
                    aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  labs(x = "Proportion of heterospecific AM basal area",
       y = "Predicted one-year survival probability") +
  base_theme
plot_C_AM

#####################################
fe      <- fixef(all_sp_propBCD)$cond
sp_coef <- coef(all_sp_propBCD)$cond$Sp
sp_info <- unique(allSpDat[, c("Sp", "Seedling_myc")])

sp_curves <- do.call(rbind, lapply(as.character(sp_info$Sp), function(s) {
  ecm <- sp_info$Seedling_myc[sp_info$Sp == s] == "ECM"
  x   <- seq(min(allSpDat$propConsp[allSpDat$Sp == s], na.rm = TRUE),
             max(allSpDat$propConsp[allSpDat$Sp == s], na.rm = TRUE),
             length.out = 100)
  slope <- sp_coef[s, "propConsp"] +
    ecm * fe["propConsp:Seedling_mycECM"] +
    fe["propConsp:propHetAM"] * median_propHetAM
  eta <- sp_coef[s, "(Intercept)"] + ecm * fe["Seedling_mycECM"] +
    slope * x +
    fe["propHetAM"]    * median_propHetAM +
    fe["Prior_counts"] * median_Prior +
    fe["log0yr"]       * median_log0yr
  data.frame(Sp = s, myc = ifelse(ecm, "ECM", "AM"),
             x = x, predicted = plogis(eta), slope = unname(slope))
}))

sp_ends <- do.call(rbind, lapply(split(sp_curves, sp_curves$Sp),
                                 function(d) d[which.max(d$x), ]))

plot_sp <- ggplot(sp_curves, aes(x, predicted, group = Sp, color = myc)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = sp_ends, aes(label = Sp), hjust = -0.15,
            size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(AM = col_AM, ECM = col_ECM),
                     name = "Seedling type") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  scale_x_continuous(limits = c(0,1.05)) +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme +
  theme(legend.position = "none") 
plot_sp

sp_slopes$slope_medHetAM <- sp_slopes$slope +
  fe["propConsp:propHetAM"] * median_propHetAM

sp_slopes[order(sp_slopes$slope_medHetAM),
          c("Sp", "Seedling_myc", "slope", "slope_medHetAM")]


################################################## 
##### Plot predictions for Proportional Model B (2nd Most Supported)
# =========================================================
## -------- PANEL A — Conspecific Density Dependence

## Predict 
pred_A <- ggpredict(all_sp_propB,
                    terms = c("propConsp", "Seedling_myc",
                              paste0("propHetAM[",  median_propHetAM,  "]")),
                    condition = c(Prior_counts = median_Prior,
                                  log0yr   = median_log0yr),
                    bias_correction = FALSE)

## Build dataframe cleanly
df_A <- data.frame(
  x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high,
  seed_type = factor(pred_A$group, levels = c("AM","ECM")))

### Create combo variable (even though no Low/Med/High here)
df_A$combo <- df_A$seed_type
df_A$combo <- factor(df_A$combo,
                     levels = c("AM", "ECM"),
                     labels = c("AM_Medium","ECM_Medium"))  # middle level color

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(0.1,0.8)) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme

plot_A

# =========================================================
## ---- PANEL B ------ Mycorrhizal-Dependent Pathogen Spillover 
## Predict 
pred_B_AM <- ggpredict(all_sp_propB,
                       terms = c("propConsp",
                                 paste0("propHetAM[", paste(propHetAM_quart, 
                                             collapse = ","), "]"),
                                 "Seedling_myc"),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr   = median_log0yr),
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

## Plot
plot_B_AM <- ggplot(df_B_AM,
                    aes(x = x, y = predicted,
                        color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(0.1,0.8)) +
  facet_wrap(~ seed_type, ncol = 2, 
             labeller = labeller(seed_type = c(
               AM  = "AM seedlings", ECM = "ECM seedlings"))) +
  scale_color_manual(values = color_map_combo, name = "Prop. AM BA",
                     guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme
plot_B_AM


# =========================================================
## ---- PANEL C ------ Conmycorrhizal Density Dependence (prop)
pred_C_AM <- ggpredict(all_sp_propB,
                       terms = c("propHetAM", "Seedling_myc",
                                 paste0("propConsp[", median_propConsp, "]")),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr       = median_log0yr),
                       bias_correction = FALSE)

df_C_AM <- data.frame(x         = pred_C_AM$x,
                      predicted = pred_C_AM$predicted,
                      conf.low  = pred_C_AM$conf.low,
                      conf.high = pred_C_AM$conf.high,
                      seed_type = factor(pred_C_AM$group, levels = c("AM","ECM")))

## Drop impossible combinations
df_C_AM <- subset(df_C_AM, x <= 1 - median_propConsp)

df_C_AM$combo <- factor(df_C_AM$seed_type,
                        levels = c("AM","ECM"),
                        labels = c("AM_Medium","ECM_Medium"))

plot_C_AM <- ggplot(df_C_AM,
                    aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  labs(x = "Proportion of heterospecific AM basal area",
       y = "Predicted one-year survival probability") +
  base_theme
plot_C_AM

#####################################
fe      <- fixef(all_sp_propB)$cond
sp_coef <- coef(all_sp_propB)$cond$Sp
sp_info <- unique(allSpDat[, c("Sp", "Seedling_myc")])

sp_curves <- do.call(rbind, lapply(as.character(sp_info$Sp), function(s) {
  ecm <- sp_info$Seedling_myc[sp_info$Sp == s] == "ECM"
  x   <- seq(min(allSpDat$propConsp[allSpDat$Sp == s], na.rm = TRUE),
             max(allSpDat$propConsp[allSpDat$Sp == s], na.rm = TRUE),
             length.out = 100)
  slope <- sp_coef[s, "propConsp"] +
    ecm * fe["propConsp:Seedling_mycECM"] +
    fe["propConsp:propHetAM"] * median_propHetAM
  eta <- sp_coef[s, "(Intercept)"] + ecm * fe["Seedling_mycECM"] +
    slope * x +
    fe["propHetAM"]    * median_propHetAM +
    fe["Prior_counts"] * median_Prior +
    fe["log0yr"]       * median_log0yr
  data.frame(Sp = s, myc = ifelse(ecm, "ECM", "AM"),
             x = x, predicted = plogis(eta), slope = unname(slope))
}))

sp_ends <- do.call(rbind, lapply(split(sp_curves, sp_curves$Sp),
                                 function(d) d[which.max(d$x), ]))

plot_sp <- ggplot(sp_curves, aes(x, predicted, group = Sp, color = myc)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = sp_ends, aes(label = Sp), hjust = -0.15,
            size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(AM = col_AM, ECM = col_ECM),
                     name = "Seedling type") +
  scale_y_continuous(limits = c(0.1,0.8)) +
  scale_x_continuous(limits = c(0,1.05)) +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year survival probability") +
  base_theme +
  theme(legend.position = "none") 
plot_sp

sp_slopes$slope_medHetAM <- sp_slopes$slope +
  fe["propConsp:propHetAM"] * median_propHetAM

sp_slopes[order(sp_slopes$slope_medHetAM),
          c("Sp", "Seedling_myc", "slope", "slope_medHetAM")]

#########################################################
#########################################################
##### All Species Models - Absolute Basal Area

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

######### ######### ######### ######### ######### 
######### Null model (Intercept only)
all_sp_15m_null <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) 
                           ~ 1 + (1|SiteID) + (1|survYear) + (1|Sp),
                           family = betabinomial, data = allSpDat)
summary(all_sp_15m_null)

######### ######### ######### ######### ######### 
######### Model A — Conspecific Density Dependence (CDD)
all_sp_15mA <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ consp_15m * Seedling_myc + Prior_counts + log0yr +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mA)
## Significant negative effect of conspBA
# VIF 
vif_A <- compute_vif(~ consp_15m * Seedling_myc +
                       Prior_counts + log0yr, allSpDat)
vif_A # All <4

######### ######### ######### ######### ######### 
######### Model B — Mycorrhizal-Dependent Pathogen Spillover 
all_sp_15mB <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m +
                         consp_15m:hetECM_15m +
                         Prior_counts + log0yr +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mB)
## Significant positive effect of hetECM_15m
## Significant negative consp_15m:hetAM_15m and consp_15m:hetECM_15m interactions 
# VIF 
vif_B <- compute_vif(~ consp_15m * Seedling_myc +
                       hetAM_15m + hetECM_15m +
                       consp_15m:hetAM_15m + consp_15m:hetECM_15m +
                       Prior_counts + log0yr,  allSpDat)
vif_B # All < 7

######### ######### ######### ######### ######### 
######### Model C — Conmycorrhizal Density Dependence
all_sp_15mC <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         Seedling_myc:hetAM_15m +
                         Seedling_myc:hetECM_15m +
                         Prior_counts + log0yr +
                         (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
                       family = betabinomial, data = allSpDat)
summary(all_sp_15mC)
## Significant negative effects of conspBA and hetAM_15m
## Significant positive effect of hetECM_15m
# VIF 
vif_C <- compute_vif( ~ consp_15m * Seedling_myc +
                        hetAM_15m + hetECM_15m +
                        Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
                        Prior_counts + log0yr, allSpDat)
vif_C # One term with value >10

######### ######### ######### ######### ######### 
######### Model D — Conmycorrhizal plus pathogen spillover
all_sp_15mD <- glmmTMB(cbind(Year_1_counts, Year_0_counts - Year_1_counts) ~ consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m +
                         consp_15m:hetECM_15m +
                         Seedling_myc:hetAM_15m +
                         Seedling_myc:hetECM_15m +
                         Prior_counts + log0yr +
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
                        Prior_counts + log0yr, allSpDat)
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
    Prior_counts + log0yr +
    (1|SiteID) + (1|survYear) + (1 + consp_15m|Sp),
  family = betabinomial, data = allSpDat)
summary(all_sp_15mBCD)
## Significant positive main effect of hetECM BA
## Significant negative interaction between consp:hetECM 
## Non-significant Seedling_mycECM:hetAM_15m and Seedling_mycECM:hetECM_15m
## consp:Seedling_mycECM:hetECM_15m Significant positive; ECM seedlings in ECM-heterogeneous microsites experience less negative or even reversed competitive effects.

# VIF and R2
vif_BCD <- compute_vif(~ consp_15m * Seedling_myc +
                         hetAM_15m + hetECM_15m +
                         consp_15m:hetAM_15m + consp_15m:hetECM_15m +
                         Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
                         consp_15m:Seedling_myc:hetAM_15m +
                         consp_15m:Seedling_myc:hetECM_15m +
                         Prior_counts + log0yr, allSpDat)
vif_BCD # 4 VIFs >10; but high VIFs on interactions do not threaten the validity of the model; they simply reflect interaction-induced multicollinearity.

### Adjust so species-specific slopes are also correct for the ECM spp
fe <- fixef(all_sp_15mBCD)$cond

sp_slopes <- coef(all_sp_15mBCD)$cond$Sp[, c("(Intercept)", "consp_15m")]
sp_slopes$Sp <- rownames(sp_slopes)
sp_slopes <- merge(sp_slopes, 
                   unique(allSpDat[, c("Sp", "Seedling_myc")]), by = "Sp")

sp_slopes$slope <- sp_slopes$consp_15m +
  (sp_slopes$Seedling_myc == "ECM") * fe["consp_15m:Seedling_mycECM"]

sp_slopes[order(sp_slopes$slope), c("Sp", "Seedling_myc", "slope")]

# Marginal prediction R² (log-scale)
obs <- allSpDat$Year_1_counts
pred_marg <- predict(all_sp_15mBCD, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg # 0.058
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mBCD, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond # 0.082

### Compare AICs
### Compute Delta AICs
aic_table <- AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
                 all_sp_15mC, all_sp_15mD, all_sp_15mBCD)

aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
aic_table[order(aic_table$deltaAIC), ]
anova(all_sp_15mB, all_sp_15mBCD) # And just barely not sign. better than B

#### Diagnostics and Post-hoc tests for super-interactive model (best supported)
### Diagnostics - Model E
simBCD <- DHARMa::simulateResiduals(all_sp_15mBCD, n = 1000)
DHARMa::testDispersion(simBCD) # No overdispersion
DHARMa::testZeroInflation(simBCD) # No zero-inflation
DHARMa::testUniformity(simBCD) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCD) # ~0.17% of observations are outliers 
# Model diagnostics indicated no issues with dispersion or zero-inflation. The DHARMa uniformity test was significant, but this is common in large datasets and not indicative of poor fit. Outlier rates were extremely low (~0.21%), confirming that the model adequately captured the distribution of the data

################################################## 
################################################## 
################################################## 
##### Plot predictions for Model E

#### Determine the values for propagule pressure:
10^median(allSpDat$log0yr, na.rm = TRUE) - 1
10^median(allSpDat$logPrior, na.rm = TRUE)

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A–D
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
median_consp  <- median(allSpDat$consp_15m,  na.rm = TRUE)
median_hetAM  <- median(allSpDat$hetAM_15m,  na.rm = TRUE)
median_hetECM <- median(allSpDat$hetECM_15m, na.rm = TRUE)
median_log0yr <- median(allSpDat$log0yr,   na.rm = TRUE)
median_Prior <- median(allSpDat$Prior_counts, na.rm = TRUE)

hetAM_quart  <- quantile(allSpDat$hetAM_15m,  
                         probs = c(.1,.5,.9), na.rm = TRUE)
hetECM_quart <- quantile(allSpDat$hetECM_15m, 
                         probs = c(.1,.5,.9), na.rm = TRUE)

### Colors
col_AM  <- "#c0352b"   # AM red
col_ECM <- "#2b4ac0"   # ECM blue

reds  <- colorRampPalette(c("#f7b4ad", col_AM,  "#7a1b14"))
blues <- colorRampPalette(c("#b3c4f7", col_ECM, "#16245c"))

# 3-level gradients for context levels (low, medium, high)
AM_grad  <- reds(3)
ECM_grad <- blues(3)

# Combined color lookup – used in all Panel B/D interaction plots
color_map_combo <- c(
  # AM seedling panel
  "AM_Low"    = AM_grad[1],
  "AM_Medium" = AM_grad[2],
  "AM_High"   = AM_grad[3],
  # ECM seedling panel
  "ECM_Low"   = ECM_grad[1],
  "ECM_Medium"= ECM_grad[2],
  "ECM_High"  = ECM_grad[3])

desired_order <- c("AM_Low","AM_Medium","AM_High",
                   "ECM_Low","ECM_Medium","ECM_High")

### Theme
base_theme <- theme_bw() +
  theme(legend.position = "right",
        axis.text  = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size=15, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# =========================================================
## -------- PANEL A — Conspecific Density Dependence
pred_A <- ggpredict(all_sp_15mBCD,
                    terms = c("consp_15m", "Seedling_myc",
                              paste0("hetAM_15m[",  median_hetAM,  "]"),
                              paste0("hetECM_15m[", median_hetECM, "]")),
                    condition = c(Prior_counts = median_Prior,
                                  log0yr   = median_log0yr),
                    bias_correction = FALSE)

## Build dataframe cleanly
df_A <- data.frame(x = pred_A$x,
                   predicted = pred_A$predicted,
                   conf.low  = pred_A$conf.low,
                   conf.high = pred_A$conf.high,
                   seed_type = factor(pred_A$group, levels = c("AM","ECM")))

### Create combo variable (even though no Low/Med/High here)
df_A$combo <- df_A$seed_type
df_A$combo <- factor(df_A$combo,levels = c("AM", "ECM"),
                     labels = c("AM_Medium","ECM_Medium"))  # middle level color

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  scale_y_continuous(limits = c(0,0.9)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme

plot_A

# =========================================================
## ----PANEL B1 ------ Mycorrhizal-Dependent Pathogen Spillover 
### Consp × hetAM 
## Predict 
pred_B_AM <- ggpredict(all_sp_15mBCD,
              terms = c("consp_15m",
              paste0("hetAM_15m[", paste(hetAM_quart,
                        collapse = ","), "]"), "Seedling_myc",
                                 paste0("hetECM_15m[", median_hetECM, "]")),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr   = median_log0yr),
                       bias_correction = FALSE)

## Build dataframe cleanly
df_B_AM <- data.frame(x = pred_B_AM$x,
                      predicted = pred_B_AM$predicted,
                      conf.low  = pred_B_AM$conf.low,
                      conf.high = pred_B_AM$conf.high,
                      seed_type = factor(pred_B_AM$facet, levels = c("AM","ECM")),
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
pred_B_ECM <- ggpredict(all_sp_15mBCD,
              terms = c("consp_15m", paste0("hetECM_15m[", 
              paste(hetECM_quart, collapse = ","), "]"),
              "Seedling_myc", paste0("hetAM_15m[", median_hetAM, "]")),
                        condition = c(Prior_counts = median_Prior,
                                      log0yr   = median_log0yr),
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
df_B_ECM$combo <- interaction(df_B_ECM$seed_type, df_B_ECM$hetECM, sep="_")
df_B_ECM$combo <- factor(df_B_ECM$combo, levels = desired_order)

######### Combined plot
df_all <- bind_rows(df_B_AM  %>% 
                      mutate(neighbor_type = "AM neighbors"),
                    df_B_ECM %>% mutate(neighbor_type = "ECM neighbors"))

df_all$seed_type <- factor(df_all$seed_type, levels = c("AM", "ECM"))
df_all$neighbor_type <- factor(df_all$neighbor_type,
                               levels = c("AM neighbors", "ECM neighbors"))

### Plot
ggplot(df_all, aes(x = x, y = predicted,
                   color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(seed_type ~ neighbor_type,
             labeller = labeller(seed_type = c(AM = "AM seedlings", 
                                               ECM = "ECM seedlings"))) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  #scale_y_continuous(limits = c(0,0.9)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme

# =========================================================
## ---- PANEL C1 ------ Conmycorrhizal Density Dependence
## Predict - AM effect
pred_C_AM <- ggpredict(all_sp_15mBCD,
                       terms = c("hetAM_15m","Seedling_myc",
                                 paste0("consp_15m[", median_consp, "]"),
                                 paste0("hetECM_15m[", median_hetECM, "]")),
                       condition = c(Prior_counts = median_Prior,
                                     log0yr   = median_log0yr),
                       bias_correction = FALSE)

## Build dataframe cleanly
df_C_AM <- data.frame(x = pred_C_AM$x,
                      predicted = pred_C_AM$predicted,
                      conf.low  = pred_C_AM$conf.low,
                      conf.high = pred_C_AM$conf.high,
                      seed_type = factor(pred_C_AM$group, levels = c("AM","ECM")))

### Middle-level AM spillover for both seedling types
df_C_AM$combo <- factor(df_C_AM$seed_type,
                        levels = c("AM","ECM"),
                        labels = c("AM_Medium","ECM_Medium"))

## Plot
plot_C_AM <- ggplot(df_C_AM,
                    aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide="none") +
  scale_y_continuous(limits = c(0,0.9)) +
  labs(x = "Heterospecific AM basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme
plot_C_AM

# =========================================================
## ---- PANEL C2 ------ Conmycorrhizal Density Dependence
## Predict - ECM effect
pred_C_ECM <- ggpredict(all_sp_15mBCD,
                        terms = c("hetECM_15m","Seedling_myc",
                                  paste0("consp_15m[", median_consp, "]"),
                                  paste0("hetAM_15m[", median_hetAM, "]")),
                        condition = c(Prior_counts = median_Prior,
                                      log0yr   = median_log0yr),
                        bias_correction = FALSE)

## Build dataframe cleanly
df_C_ECM <- data.frame(x = pred_C_ECM$x,
                       predicted = pred_C_ECM$predicted,
                       conf.low  = pred_C_ECM$conf.low,
                       conf.high = pred_C_ECM$conf.high,
                       seed_type = factor(pred_C_ECM$group, levels = c("AM","ECM")))

### Medium-level ECM spillover for both seedling types
df_C_ECM$combo <- factor(df_C_ECM$seed_type,
                         levels = c("AM","ECM"),
                         labels = c("AM_Medium","ECM_Medium"))

## Plot
plot_C_ECM <- ggplot(df_C_ECM,
                     aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide="none") +
  scale_y_continuous(limits = c(0,0.9)) +
  labs(x = "Heterospecific ECM basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme
plot_C_ECM

#####################################
fe      <- fixef(all_sp_15mBCD)$cond
sp_coef <- coef(all_sp_15mBCD)$cond$Sp
sp_info <- unique(allSpDat[, c("Sp", "Seedling_myc")])
sp_info$Sp <- as.character(sp_info$Sp)

median_hetAM  <- median(allSpDat$hetAM_15m,  na.rm = TRUE)
median_hetECM <- median(allSpDat$hetECM_15m, na.rm = TRUE)

## Species-specific slope and intercept (logit scale) at given het BA
sp_terms <- function(s, ecm, hAM, hECM) {
  slope <- sp_coef[s, "consp_15m"] +
    ecm  * fe["consp_15m:Seedling_mycECM"] +
    hAM  * (fe["consp_15m:hetAM_15m"]  + ecm * fe["consp_15m:Seedling_mycECM:hetAM_15m"]) +
    hECM * (fe["consp_15m:hetECM_15m"] + ecm * fe["consp_15m:Seedling_mycECM:hetECM_15m"])
  int <- sp_coef[s, "(Intercept)"] +
    ecm  * fe["Seedling_mycECM"] +
    hAM  * (fe["hetAM_15m"]  + ecm * fe["Seedling_mycECM:hetAM_15m"]) +
    hECM * (fe["hetECM_15m"] + ecm * fe["Seedling_mycECM:hetECM_15m"]) +
    fe["Prior_counts"] * median_Prior +
    fe["log0yr"]       * median_log0yr
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

plot_sp <- ggplot(sp_curves, aes(x, predicted, group = Sp, color = myc)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = sp_ends, aes(label = Sp), hjust = -0.15,
            size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(AM = col_AM, ECM = col_ECM),
                     name = "Seedling type") +
  scale_y_continuous(limits = c(0.3,0.75)) +
  scale_x_continuous(limits = c(0,4.2)) +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year survival probability") +
  base_theme +
  theme(legend.position = "none") 
plot_sp

sp_slopes$slope_medHetAM <- sp_slopes$slope +
  fe["consp_15m:hetAM_15m"] * median_propHetAM

sp_slopes[order(sp_slopes$slope_medHetAM),
          c("Sp", "Seedling_myc", "slope", "slope_medHetAM")]
