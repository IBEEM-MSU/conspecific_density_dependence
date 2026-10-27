
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### *REVISED VERSIONS BASED ON CHRIS' FEEDBACK*
####### Date: 8 December 2025
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

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")

allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

#########################################################
#########################################################
#########################################################
################# All Species Models - Raw Basal Area

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


######### ######### ######### ######### ######### 
######### Null model
all_sp_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                             (1|SiteID) + (1|survYear) + (1|Sp),
                           family = nbinom2, data = allSpDat)
summary(all_sp_15m_null)
# R2
r2_null <- performance::r2(all_sp_15m_null)
r2_null

######### ######### ######### ######### ######### 
######### PANEL A — Classic Feedback
all_sp_15mA <- glmmTMB(Year_1_counts ~ consp_15m * Seedling_myc +
    logPrior + log0yr +
    (1|SiteID) + (1|survYear) + (1|Sp),
  family = nbinom2, data = allSpDat)
summary(all_sp_15mA)
## Significant negative effect of conspBA
# VIF and R2
vif_A <- compute_vif(~ consp_15m * Seedling_myc +
                        logPrior + log0yr, allSpDat)
vif_A # All <4
r2_A <- performance::r2(all_sp_15mA)
r2_A

######### ######### ######### ######### ######### 
######### PANEL B — Interactive Spillover
all_sp_15mB <- glmmTMB(Year_1_counts ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m +
    consp_15m:hetECM_15m +
    logPrior + log0yr +
    (1|SiteID) + (1|survYear) + (1|Sp),
  family = nbinom2, data = allSpDat)
summary(all_sp_15mB)
## Significant positive effects of conspBA, hetAM_15m, and hetECM_15m
## Significant negative consp_15m:hetAM_15m and consp_15m:hetECM_15m interactions 
# VIF and R2
vif_B <- compute_vif(~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    logPrior + log0yr,  allSpDat)
vif_B # All <8
r2_B <- performance::r2(all_sp_15mB)
r2_B

######### ######### ######### ######### ######### 
######### PANEL C — Additive Spillover
all_sp_15mC <- glmmTMB(
  Year_1_counts ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    Seedling_myc:hetAM_15m +
    Seedling_myc:hetECM_15m +
    logPrior + log0yr +
    (1|SiteID) + (1|survYear) + (1|Sp),
  family = nbinom2,data = allSpDat)
summary(all_sp_15mC)
## Significant negative effects of conspBA and hetAM_15m
# VIF and R2
vif_C <- compute_vif( ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    logPrior + log0yr, allSpDat)
vif_C # One >16
r2_C <- performance::r2(all_sp_15mC)
r2_C

######### ######### ######### ######### ######### 
######### PANEL D — Interactive + Additive (B + C combined)
all_sp_15mD <- glmmTMB(
  Year_1_counts ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m +
    consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m +
    Seedling_myc:hetECM_15m +
    logPrior + log0yr +
    (1|SiteID) + (1|survYear) + (1|Sp),
  family = nbinom2,data = allSpDat)
summary(all_sp_15mD)
## Significant positive effects of conspBA, hetAM_15m, and hetECM_15m
## Significant negative consp_15m:hetAM_15m and consp_15m:hetECM_15m interactions 
# VIF and R2
vif_D <- compute_vif( ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    logPrior + log0yr, allSpDat)
vif_D # one >18
r2_D <- performance::r2(all_sp_15mD)
r2_D

######### ######### ######### ######### ######### 
######### PANEL B+C+D — SUPERINTERACTIVE MODEL
all_sp_15mBCD <- glmmTMB(Year_1_counts ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m +
    consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m +
    Seedling_myc:hetECM_15m +
    consp_15m:Seedling_myc:hetAM_15m +
    consp_15m:Seedling_myc:hetECM_15m +
    logPrior + log0yr +
    (1|SiteID) + (1|survYear) + (1|Sp),
  family = nbinom2,data = allSpDat)
summary(all_sp_15mBCD)
## Significant positive main effects of consp, hetAM, and hetECM BA
## Significant negative consp:Seedling_mycECM interaction
## Significant negative interactions between consp:hetAM and consp:hetECM 
## Non-significant Seedling_mycECM:hetAM_15m and Seedling_mycECM:hetECM_15m 
## consp:Seedling_mycECM:hetECM_15m Significant positive; ECM seedlings in ECM-heterogeneous microsites experience less negative or even reversed competitive effects.
# VIF and R2
vif_BCD <- compute_vif(~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    consp_15m:Seedling_myc:hetAM_15m +
    consp_15m:Seedling_myc:hetECM_15m +
    logPrior + log0yr, allSpDat)
vif_BCD # 3 VIFs >10; but high VIFs on interactions do not threaten the validity of the model; they simply reflect interaction-induced multicollinearity.
r2_BCD <- performance::r2(all_sp_15mBCD)
r2_BCD # 52% of variation explained by fixed effects, additional 14% explained by REs

### Compare AICs
### All sites combined
AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
    all_sp_15mC, all_sp_15mD, all_sp_15mBCD) # BCD is lowest
anova(all_sp_15mD, all_sp_15mBCD) # And is significantly better than B

#### Diagnostics and Post-hoc tests for super-interactive model (best supported)
### Diagnostics
simBCD <- DHARMa::simulateResiduals(all_sp_15mBCD, n = 1000)
DHARMa::testDispersion(simBCD) # No overdispersion
DHARMa::testZeroInflation(simBCD) # No zero-inflation
DHARMa::testUniformity(simBCD) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCD) # ~0.28% of observations are outliers - not concerning
# Model diagnostics indicated no issues with dispersion or zero-inflation, with tests showing residual variation and the frequency of zeros were fully consistent with model expectations. The DHARMa uniformity test was significant, but this is common in large datasets and not indicative of poor fit. Outlier rates were extremely low (~0.28%) and within acceptable bounds, confirming that the model adequately captured the distribution of the data

########## ########## ########## ########## ########## 
########## Plotting effect sizes - for Panel BCD

### Extract fixed effects 
fixef_df <- broom.mixed::tidy(
  all_sp_15mBCD,
  effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)

### Relabel terms
fixef_df <- fixef_df %>%
  mutate(component = case_when(
    term == "consp_15m" ~ "Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "hetAM_15m" ~ "AM BA",
    term == "hetECM_15m" ~ "ECM BA",
    term == "consp_15m:Seedling_mycECM" ~ "Consp. BA:Mycorrhizal type",
    term == "consp_15m:hetAM_15m" ~ "Consp. BA:AM BA",
    term == "consp_15m:hetECM_15m" ~ "Consp. BA:ECM BA",
    term == "Seedling_mycECM:hetAM_15m" ~ "Mycorrhizal type:AM BA",
    term == "Seedling_mycECM:hetECM_15m" ~ "Mycorrhizal type:ECM BA",
    term == "consp_15m:Seedling_mycECM:hetAM_15m" ~ "Consp. BA:Mycorrhizal type:AM BA",
    term == "consp_15m:Seedling_mycECM:hetECM_15m" ~ "Consp. BA:Mycorrhizal type:ECM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))

### Reorder
desired_order <- c(
  "Conspecific BA","AM BA","ECM BA","Mycorrhizal type",
  "Consp. BA:AM BA","Consp. BA:ECM BA","Consp. BA:Mycorrhizal type",
  "Mycorrhizal type:AM BA","Mycorrhizal type:ECM BA",
  "Consp. BA:Mycorrhizal type:AM BA","Consp. BA:Mycorrhizal type:ECM BA",
  "Log(Recruit Density)","Log(Prior Density)")
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = desired_order))

### Add in significance
fixef_df <- fixef_df %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))

### Reorder for plotting
fixef_df$component <- factor(fixef_df$component, levels = rev(desired_order))

### Plot
ggplot(fixef_df, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error), width = 0.1,
    orientation = "y", color = "black")+
  geom_point(aes(shape = Significant, fill = Significant),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"))

############ Estimated marginal means at specific predictor combinations
# Gives predicted survival for AM vs ECM seedlings at low/medium/high levels of all three BA predictors (i.e., predictions under specific values of confounders)
emm_BCD <- emmeans(all_sp_15mBCD,
  ~ Seedling_myc | consp_15m * hetAM_15m * hetECM_15m,
  type = "response",
  at = list(consp_15m   = quantile(allSpDat$consp_15m, c(.1, .5, .9)),
            hetAM_15m   = quantile(allSpDat$hetAM_15m, c(.1, .5, .9)),
            hetECM_15m  = quantile(allSpDat$hetECM_15m, c(.1, .5, .9))))
pairs(emm_BCD, by = c("consp_15m", "hetAM_15m", "hetECM_15m"))

############ marginal effects / average partial effects
# Computes average marginal effect (AME) of each continuous predictor
# Gives the closest analog to “the effect of X on survival,” even in the presence of interactions
ame_continuous <- function(model, var){
  emtrends(model, specs = ~1, var = var)
}
# Marginal effect of conspecific BA
ame_consp <- ame_continuous(all_sp_15mBCD, "consp_15m")
ame_consp %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific AM BA
ame_hetAM <- ame_continuous(all_sp_15mBCD, "hetAM_15m")
ame_hetAM %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific ECM BA
ame_hetECM <- ame_continuous(all_sp_15mBCD, "hetECM_15m")
ame_hetECM %>%
  summary(infer = TRUE)
# Marginal effect of seedling mycorrhizal guild
ame_myc <- emmeans(all_sp_15mBCD, pairwise ~ Seedling_myc)
ame_myc 

################################################## 
##### Plot predictions for superinteractive model

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A–D
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
mid_consp  <- median(allSpDat$consp_15m,  na.rm = TRUE)
mid_hetAM  <- median(allSpDat$hetAM_15m,  na.rm = TRUE)
mid_hetECM <- median(allSpDat$hetECM_15m, na.rm = TRUE)

hetAM_quart  <- quantile(allSpDat$hetAM_15m,  probs = c(.1,.5,.9), na.rm = TRUE)
hetECM_quart <- quantile(allSpDat$hetECM_15m, probs = c(.1,.5,.9), na.rm = TRUE)

### Colors
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
    plot.title = element_text(size=15, face="bold"))

# =========================================================
## -------- PANEL A — Classic feedback (consp × seedling type)
## Predict 
pred_A <- ggpredict(
  all_sp_15mBCD,
  terms = c("consp_15m", "Seedling_myc",
            paste0("hetAM_15m[",  mid_hetAM,  "]"),
            paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = TRUE)

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
                     labels = c("AM","ECM")) +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme

plot_A

# =========================================================
## ---------------------PANEL B1 ------ Consp × hetAM 
## Predict 
pred_B_AM <- ggpredict(all_sp_15mBCD,
  terms = c("consp_15m",
            paste0("hetAM_15m[", paste(hetAM_quart, collapse = ","), "]"),
            "Seedling_myc",
            paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = TRUE)

## Build dataframe cleanly
df_B_AM <- data.frame(x         = pred_B_AM$x,
  predicted = pred_B_AM$predicted,
  conf.low  = pred_B_AM$conf.low,
  conf.high = pred_B_AM$conf.high,
  seed_type = factor(pred_B_AM$facet, levels = c("AM","ECM")),
  hetAM_raw = pred_B_AM$group)

## Assign human-readable labels
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
  scale_color_manual(values = color_map_combo, name = "AM BA Level") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
    y = "Predicted one-year seedling counts") +
  base_theme
plot_B_AM

# =========================================================
## -------- PANEL B2 ----------- Consp × hetECM
## Predict 
pred_B_ECM <- ggpredict(all_sp_15mBCD,
  terms = c("consp_15m",
            paste0("hetECM_15m[", paste(hetECM_quart, collapse = ","), "]"),
            "Seedling_myc",
            paste0("hetAM_15m[", mid_hetAM, "]")),
  bias_correction = TRUE)

## Build dataframe cleanly
df_B_ECM <- data.frame(
  x         = pred_B_ECM$x,
  predicted = pred_B_ECM$predicted,
  conf.low  = pred_B_ECM$conf.low,
  conf.high = pred_B_ECM$conf.high,
  seed_type = factor(pred_B_ECM$facet, levels = c("AM","ECM")),
  hetECM_raw = pred_B_ECM$group)

## Assign human-readable labels
df_B_ECM$hetECM <- factor(
  df_B_ECM$hetECM_raw,
  levels = unique(pred_B_ECM$group),
  labels = c("Low","Medium","High"))

## Interaction label for color: AM_Low, AM_Medium, etc.
df_B_ECM$combo <- interaction(df_B_ECM$seed_type, df_B_ECM$hetECM, sep="_")
df_B_ECM$combo <- factor(df_B_ECM$combo, levels = desired_order)

## Plot
plot_B_ECM <- ggplot(df_B_ECM,
                     aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ seed_type, ncol=2,
             labeller = labeller(seed_type = c(
               AM="AM seedlings", ECM="ECM seedlings"))) +
  scale_color_manual(values = color_map_combo,
                     name = "ECM BA Level", guide="none") +
  scale_fill_manual(values = color_map_combo, guide="none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_B_ECM

# =========================================================
## -------- PANEL C1 -------- Additive AM effect
## Predict
pred_C_AM <- ggpredict(all_sp_15mBCD,
  terms = c("hetAM_15m","Seedling_myc",
            paste0("consp_15m[", mid_consp, "]"),
            paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = TRUE)

## Build dataframe cleanly
df_C_AM <- data.frame(
  x         = pred_C_AM$x,
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
                     labels = c("AM","ECM")) +
  scale_fill_manual(values = color_map_combo, guide="none") +
  labs(x = "Heterospecific AM basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_C_AM

# =========================================================
## -------- PANEL C2 --------- Additive ECM effect
## Predict
pred_C_ECM <- ggpredict(all_sp_15mBCD,
  terms = c("hetECM_15m","Seedling_myc",
            paste0("consp_15m[", mid_consp, "]"),
            paste0("hetAM_15m[", mid_hetAM, "]")),
  bias_correction = TRUE)

## Build dataframe cleanly
df_C_ECM <- data.frame(
  x         = pred_C_ECM$x,
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
                     labels = c("AM","ECM")) +
  scale_fill_manual(values = color_map_combo, guide="none") +
  labs(x = "Heterospecific ECM basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_C_ECM

# # =========================================================
# ## -------- PANEL D1 ------ Consp × hetAM (Panel D logic)
# ## Predict 
# pred_D_AM <- ggpredict(
#   all_sp_15mBCD,
#   terms = c("consp_15m",
#             paste0("hetAM_15m[", paste(hetAM_quart, collapse = ","), "]"),
#             "Seedling_myc",
#             paste0("hetECM_15m[", mid_hetECM, "]"))
# )
# 
# ## Build dataframe cleanly
# df_D_AM <- data.frame(
#   x         = pred_D_AM$x,
#   predicted = pred_D_AM$predicted,
#   conf.low  = pred_D_AM$conf.low,
#   conf.high = pred_D_AM$conf.high,
#   seed_type = factor(pred_D_AM$facet, levels = c("AM","ECM")),
#   hetAM_raw = pred_D_AM$group
# )
# 
# ## Assign human-readable labels for spillover levels
# df_D_AM$hetAM <- factor(
#   df_D_AM$hetAM_raw,
#   levels = unique(pred_D_AM$group),
#   labels = c("Low","Medium","High")
# )
# 
# ## Build combo variable (AM_Low, ECM_Low, etc.)
# df_D_AM$combo <- interaction(df_D_AM$seed_type, df_D_AM$hetAM, sep="_")
# df_D_AM$combo <- factor(df_D_AM$combo, levels = desired_order)
# 
# ## Plot
# plot_D_AM <- ggplot(df_D_AM,
#                     aes(x = x, y = predicted, color = combo, fill = combo)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
#               alpha = 0.18, color = NA) +
#   geom_line(size = 0.8) +
#   facet_wrap(~ seed_type, ncol = 2,
#              labeller = labeller(seed_type = c(
#                AM  = "AM seedlings",
#                ECM = "ECM seedlings"))) +
#   scale_color_manual(values = color_map_combo,
#                      name = "AM BA Level") +
#   scale_fill_manual(values = color_map_combo, guide = "none") +
#   labs(title = "Panel D1 — Full model: AM spillover × Conspecific BA",
#        x = "Conspecific basal area (m²)",
#        y = "Predicted survival") +
#   base_theme
# 
# plot_D_AM
# 
# # =========================================================
# ## -------- PANEL D2 -------- Consp × hetECM (Panel D logic)
# ## Predict 
# pred_D_ECM <- ggpredict(
#   all_sp_15mBCD,
#   terms = c("consp_15m",
#             paste0("hetECM_15m[", paste(hetECM_quart, collapse = ","), "]"),
#             "Seedling_myc",
#             paste0("hetAM_15m[", mid_hetAM, "]")))
# 
# ## Build dataframe cleanly
# df_D_ECM <- data.frame(
#   x         = pred_D_ECM$x,
#   predicted = pred_D_ECM$predicted,
#   conf.low  = pred_D_ECM$conf.low,
#   conf.high = pred_D_ECM$conf.high,
#   seed_type = factor(pred_D_ECM$facet, levels = c("AM","ECM")),
#   hetECM_raw = pred_D_ECM$group)
# 
# ## Assign human-readable labels
# df_D_ECM$hetECM <- factor(
#   df_D_ECM$hetECM_raw,
#   levels = unique(pred_D_ECM$group),
#   labels = c("Low","Medium","High"))
# 
# ## Color-group label: AM_Low, AM_Medium, …
# df_D_ECM$combo <- interaction(df_D_ECM$seed_type, df_D_ECM$hetECM, sep="_")
# df_D_ECM$combo <- factor(df_D_ECM$combo, levels = desired_order)
# 
# ## Plot
# plot_D_ECM <- ggplot(df_D_ECM,
#                      aes(x = x, y = predicted, color = combo, fill = combo)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
#               alpha = 0.18, color = NA) +
#   geom_line(size = 0.8) +
#   facet_wrap(~ seed_type, ncol = 2,
#              labeller = labeller(seed_type = c(
#                AM="AM seedlings",
#                ECM="ECM seedlings"))) +
#   scale_color_manual(values = color_map_combo,
#                      name = "ECM BA Level") +
#   scale_fill_manual(values = color_map_combo, guide="none") +
#   labs(title = "Panel D2 — Additive & Interactive Spillover (ECM × Consp BA)",
#        x = "Conspecific basal area (m²)",
#        y = "Predicted survival") +
#   base_theme
# 
# plot_D_ECM


#########################################################
#########################################################
################# All Species Models - Proportional Basal Area

################# ################# ################# 
## Null model
all_sp_prop_null <- glmmTMB(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = nbinom2, data = allSpDat) 
summary(all_sp_prop_null)
# R2
r2_null_prop <- performance::r2(all_sp_prop_null)
r2_null_prop

################# ################# ################# 
## Panel A (proportional BA)
all_sp_propA <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = nbinom2, data = allSpDat) 
summary(all_sp_propA)
## Significant negative effect of propConsp
## propConsp:Seedling_myc interaction is negative and marginally significant
# VIF and R2
vif_propA <- compute_vif( ~ propConsp*Seedling_myc +
                        logPrior + log0yr, allSpDat)
vif_propA # all <4
r2_propA <- performance::r2(all_sp_propA)
r2_propA

################# ################# ################# 
## Panel B (proportional BA)
all_sp_propB <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propHetAM + 
                       propConsp:propHetAM + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propB)
## Significant negative effects of propConsp and propAM
## propConsp:propAM is significant and positive
# VIF and R2
vif_propB <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            logPrior + log0yr, allSpDat)
vif_propB # All <10
r2_propB <- performance::r2(all_sp_propB)
r2_propB

################# ################# ################# 
## Panel C (proportional BA)
all_sp_propC <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propHetAM + 
                        Seedling_myc:propHetAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propC)
## Significant positive and negative effects of propConsp and propAM, respectively
## propConsp:Seedling_myc and Seedling_myc:AM_15m interactions not significant 
# VIF and R2
vif_propC <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + Seedling_myc:propHetAM + 
                            logPrior + log0yr, allSpDat)
vif_propC # All <5
r2_propC <- performance::r2(all_sp_propC)
r2_propC

################# ################# ################# 
## Panel D (B+C) (proportional BA)
all_sp_propD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propHetAM + 
                        propConsp:propHetAM + # Panel B
                        Seedling_myc:propHetAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propD)
## Significant negative effects of propConsp and propAM
## propConsp:propAM interaction is significant and positive
# VIF and R2
vif_propD <- compute_vif( ~ propConsp*Seedling_myc + 
                            propHetAM + propConsp:propHetAM + 
                            Seedling_myc:propHetAM + 
                            logPrior + log0yr, allSpDat)
vif_propD # All <7
r2_propD <- performance::r2(all_sp_propD)
r2_propD

################# ################# ################# 
## Super interactive model (proportional BA)
all_sp_propBCD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propHetAM + 
                       propConsp:propHetAM + # Panel B
                       Seedling_myc:propHetAM + # Panel C
                       propConsp:Seedling_myc:propHetAM + # 3-way interaction
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propBCD)
## Significant negative effects of propConsp and propAM 
## propConsp:Seedling_myc and propConsp:propAMm interactions are significant and positive!
## Seedling_myc:propAM interaction not significant
## propConsp:Seedling_myc:propAM interaction is significant and negative
# VIF and R2
vif_propBCD <- compute_vif( ~ propConsp*Seedling_myc + 
                              propHetAM + propConsp:propHetAM + 
                              Seedling_myc:propHetAM + 
                              propConsp:Seedling_myc:propHetAM + 
                              logPrior + log0yr, allSpDat)
vif_propBCD # four values >10 due to interactions
r2_propBCD <- performance::r2(all_sp_propBCD)
r2_propBCD

## Model comparison
### Compare AICs 
AIC(all_sp_prop_null, all_sp_propA, all_sp_propB, all_sp_propC,
    all_sp_propD, all_sp_propBCD) # lowest is all_sp_propB 
anova(all_sp_propB, all_sp_propD) # Panel D is NOT significantly better!

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Diagnostics and Post-hoc tests for Panel B (best supported)
simBprop <- DHARMa::simulateResiduals(all_sp_propB, n = 1000)
DHARMa::testDispersion(simBprop) # No overdispersion
DHARMa::testZeroInflation(simBprop) # No zero-inflation
DHARMa::testUniformity(simBprop) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBprop) # ~0.32% of observations are outliers
# Model diagnostics indicated excellent fit, with no evidence of overdispersion or zero inflation. Although DHARMa’s uniformity and outlier tests showed significant deviations, the magnitude was trivial (<0.4% outliers) and expected for large-sample GLMMs, and residual plots showed no structured patterns.
# Post-hoc comparisons revealed no reliable differences in seedling abundance between AM and ECM species across proportional conspecific or heterospecific basal area gradients. All nine AM–ECM contrasts yielded ratios near 1.0 with non-significant p-values except for a single marginal comparison at high conspecific and high AM basal area, which we interpret cautiously given the lack of a consistent pattern.

########## ########## ########## ########## ########## 
########## Plotting effect sizes - Panel B
### Extract fixed effects 
fixef_df <- broom.mixed::tidy(
  all_sp_propB, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)

### Map raw term names to clean labels you want shown on the plot
fixef_df <- fixef_df %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propHetAM" ~ "Prop. AM BA",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "propConsp:propHetAM" ~ "Prop. Consp. BA:Prop. AM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))

# 3. Order components on the y-axis in the order you prefer
desired_order <- c("Prop. Conspecific BA","Prop. AM BA",
                   "Mycorrhizal type",  "Prop. Consp. BA:Myc. type ECM",
                   "Prop. Consp. BA:Prop. AM BA", 
                   "Log(Recruit Density)","Log(Prior Density)")
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = desired_order))

# 4. Significance flag (alpha = 0.05)
fixef_df <- fixef_df %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))

# Make your factor for y-axis
fixef_df$component <- factor(fixef_df$component, levels = rev(desired_order))

### Plot
ggplot(fixef_df, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error), 
                orientation = "y",width = 0.1, color = "black")+
  geom_point(aes(shape = Significant, fill = Significant), 
             size = 3.5, color = "black", stroke = 1) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

############ Estimated marginal means at specific predictor combinations
# Gives predicted survival for AM vs ECM seedlings at low/medium/high levels of the proportion of conspecific and AM BA (i.e., predictions under specific values of confounders)
emm_Bprop <- emmeans(all_sp_propB,
  ~ Seedling_myc | propConsp * propHetAM,
  type = "response",
  at = list(propConsp  = quantile(allSpDat$propConsp, c(.1, .5, .9)),
            propAM  = quantile(allSpDat$propHetAM, c(.1, .5, .9))))
emm_Bprop
pairs(emm_Bprop, by = c("propConsp", "propHetAM"))

############ marginal effects / average partial effects
# Computes average marginal effect (AME) of each continuous predictor
# Gives the closest analog to “the effect of X on survival,” even in the presence of interactions

ame_continuous <- function(model, var){
  emtrends(model, specs = ~1, var = var)
}
# Marginal effect of proportion of conspecific BA
ameProp_consp <- ame_continuous(all_sp_propB, "propConsp")
ameProp_consp %>%
  summary(infer = TRUE)
# Marginal effect of proportion of heterospecific AM BA
ameProp_hetAM <- ame_continuous(all_sp_propB, "propHetAM")
ameProp_hetAM %>%
  summary(infer = TRUE)
# Marginal effect of seedling mycorrhizal guild
ameProp_myc <- emmeans(all_sp_propB, pairwise ~ Seedling_myc)
ameProp_myc


################################################## 
##### Plot predictions for Panel B

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A and B
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
mid_propConsp  <- median(allSpDat$propConsp,  na.rm = TRUE)
mid_propHetAM  <- median(allSpDat$propHetAM,  na.rm = TRUE)

propHetAM_quart  <- quantile(allSpDat$propHetAM,  probs = c(.1,.5,.9), na.rm = TRUE)

### Colors
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
        plot.title = element_text(size=15, face="bold"))

# =========================================================
## -------- PANEL A — Classic feedback (propConsp × seedling type)
## Predict 
pred_A <- ggpredict(
  all_sp_propB,
  terms = c("propConsp", "Seedling_myc",
            paste0("propHetAM[",  mid_propHetAM,  "]")),
  bias_correction = TRUE)

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
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year seedling counts") +
  base_theme

plot_A

# =========================================================
## ----------------PANEL B1 ------ propConsp × propHetAM 
## Predict 
pred_B_AM <- ggpredict(all_sp_propB,
  terms = c("propConsp",
            paste0("propHetAM[", paste(propHetAM_quart, collapse = ","), "]"),
            "Seedling_myc"))

## Build dataframe cleanly
df_B_AM <- data.frame(x         = pred_B_AM$x,
                      predicted = pred_B_AM$predicted,
                      conf.low  = pred_B_AM$conf.low,
                      conf.high = pred_B_AM$conf.high,
                      seed_type = factor(pred_B_AM$facet, levels = c("AM","ECM")),
                      hetAM_raw = pred_B_AM$group)

## Assign human-readable labels
df_B_AM$hetAM <- factor(
  df_B_AM$hetAM_raw,
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
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_B_AM


#########################################################
#########################################################
#########################################################
################# Mycorrhizal Type Models - Raw Basal Area

## Subset out just AM and just ECM species
allSpDat_AM = subset(allSpDat, Seedling_myc == "AM")
allSpDat_ECM = subset(allSpDat, Seedling_myc == "ECM")

# For ECM seedlings, negative-binomial GLMMs consistently showed non-identifiable dispersion parameters (approaching 0 or extremely large), indicating that after accounting for random effects there was no meaningful overdispersion to estimate. Poisson GLMMs provided stable convergence, lower AIC values, and acceptable residual diagnostics. Therefore, all ECM models were fit using Poisson errors.
# Because AM and ECM seedlings differed strongly in their variance structure, we fit each guild using the distribution that best matched its data-generating process. AM seedling counts were overdispersed even after accounting for random effects, and negative binomial GLMMs gave stable fits and lower AIC values. In contrast, ECM seedling counts showed no estimable overdispersion—negative binomial models produced boundary estimates or convergence warnings—while Poisson GLMMs yielded stable fits and acceptable residual diagnostics. Thus, Poisson error was used for ECM seedlings and negative binomial error for AM seedlings. Importantly, both model sets used the same fixed and random effects, the same link function, and the same biological predictors, ensuring comparability of ecological interpretations across guilds.

########## ########## ########## ########## ########## 
##### Null model - AM
AM_sp_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15m_null)
# R2
r2_nullAM <- performance::r2(AM_sp_15m_null)
r2_nullAM

##### Null model - ECM
ECM_sp_15m_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15m_null)
# R2
r2_nullECM <- performance::r2(ECM_sp_15m_null)
r2_nullECM

########## ########## ########## ########## ########## 
##### Panel A - CDD - AM
AM_sp_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15mA)
## Significant negative effect of ConspBA 
# VIF and R2
vif_AMa <- compute_vif(~ consp_15m +
                        logPrior + log0yr, allSpDat_AM)
vif_AMa # All <2
r2_A_AM <- performance::r2(AM_sp_15mA)
r2_A_AM

##### Panel A - CDD - ECM
ECM_sp_15mA <- glmer(Year_1_counts ~ consp_15m +
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                       family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15mA)
## Significant negative effect of ConspBA 
# VIF and R2
vif_ECMa <- compute_vif(~ consp_15m +
                         logPrior + log0yr, allSpDat_ECM)
vif_ECMa # All <2
r2_A_ECM <- performance::r2(ECM_sp_15mA)
r2_A_ECM

########## ########## ########## ########## ########## 
######## Panel B - AM
AM_sp_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                       hetAM_15m + hetECM_15m + 
                       consp_15m:hetAM_15m + # Panel B
                       consp_15m:hetECM_15m + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat_AM)
summary(AM_sp_15mB)
## Significant positive main effects of conspBA, AM BA, and ECM BA
## consp_15m:AM_15m and consp_15m:ECM_15m are negative and significant
# VIF and R2
vif_AMb <- compute_vif(~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr, allSpDat_AM)
vif_AMb # All <6
r2_B_AM <- performance::r2(AM_sp_15mB)
r2_B_AM

######## Panel B - ECM
# For ECM seedlings, negative-binomial GLMMs failed to converge and estimated a dispersion parameter effectively equal to zero (α < 1×10⁻⁸), indicating no detectable overdispersion after accounting for random effects. We therefore used a Poisson GLMM for this group, which showed stable convergence and acceptable residual diagnostics.
ECM_sp_15mB <- glmer(Year_1_counts ~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr +
                         (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15mB)
## Non-significant main effects of conspBA, AM BA, and ECM BA
## consp_15m:AM_15m is negative and significant
## consp_15m:ECM_15m is positive and significant
# VIF and R2
vif_ECMb <- compute_vif(~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr, allSpDat_ECM)
vif_ECMb # All <7
r2_B_ECM <- performance::r2(ECM_sp_15mB)
r2_B_ECM

########## ########## ########## ########## ########## 
##### Model comparison - AM
### Compare AICs  
AIC(AM_sp_15m_null, AM_sp_15mA, AM_sp_15mB) # Lowest is AM_sp_15mB
## Model comparison - ECM
### Compare AICs  
AIC(ECM_sp_15m_null, ECM_sp_15mA, ECM_sp_15mB) # Lowest is ECM_sp_15mB
anova(ECM_sp_15mA, ECM_sp_15mB) # B is significantly better

### Diagnostics - AM
sim_AM_B <- DHARMa::simulateResiduals(AM_sp_15mB, n = 1000)
DHARMa::testDispersion(sim_AM_B) # No overdispersion
DHARMa::testZeroInflation(sim_AM_B) # No zero-inflation
DHARMa::testUniformity(sim_AM_B) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_AM_B) # ~0.30% of observations are outliers

### Diagnostics - ECM
sim_ECM_B <- DHARMa::simulateResiduals(ECM_sp_15mB, n = 1000)
DHARMa::testDispersion(sim_ECM_B) # No overdispersion
DHARMa::testZeroInflation(sim_ECM_B) # No zero-inflation
DHARMa::testUniformity(sim_ECM_B) # Non-uniform, but normal with large datasets 
DHARMa::testOutliers(sim_ECM_B) # ~0.09% of observations are outliers

########## Plotting effect sizes - for Panel B (most complex mycorrhizal version)
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(component = case_when(
        term == "consp_15m" ~ "Conspecific BA",
        term == "hetAM_15m" ~ "AM BA",
        term == "hetECM_15m" ~ "ECM BA",
        term == "consp_15m:hetAM_15m" ~ "Conspecific BA:AM BA",
        term == "consp_15m:hetECM_15m" ~ "Conspecific BA:ECM BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(component = factor(component, levels = c(
        "Conspecific BA", "AM BA", "ECM BA", 
        "Conspecific BA:AM BA","Conspecific BA:ECM BA",
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for both models
fixef_AM  <- extract_fixef(AM_sp_15mB, "AM")
fixef_ECM <- extract_fixef(ECM_sp_15mB, "ECM")

## Combine
fixef_df <- bind_rows(fixef_AM, fixef_ECM)

# Define desired vertical order
desired_order <- c("Conspecific BA", "AM BA", 
                   "ECM BA", "Conspecific BA:AM BA","Conspecific BA:ECM BA",
                   "Log(Recruit Density)", "Log(Prior Density)")

# Prepare data
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = rev(desired_order)),
    Model = factor(Model, levels = c("AM", "ECM")),
     Significant = factor(Significant, levels = c("yes", "no")))

######## Plot
offset <- 0.09
ggplot(fixef_df, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes( xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error),
    width = 0.1, color = "black",
    position = position_nudge(y = ifelse(fixef_df$Model == "AM", offset, -offset))) +
  geom_point(aes(shape = Model, fill = Significant),
    size = 3, color = "black", stroke = 1,
    position = position_nudge(y = ifelse(fixef_df$Model == "AM", offset, -offset))) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes =
              list(shape = 21, color = "black", size = 3.5)))

############ Estimated marginal means at specific predictor combinations
# Gives predicted survival for AM seedlings at low/medium/high levels of all three BA predictors (i.e., predictions under specific values of confounders)
# Across levels of hetAM BA, holding hetECM constant
emm_consp_hetAM <- emmeans(AM_sp_15mB,
  ~ consp_15m * hetAM_15m,
  type = "response",
  at = list(consp_15m = quantile(allSpDat_AM$consp_15m, c(.1, .5, .9)),
    hetAM_15m = quantile(allSpDat_AM$hetAM_15m, c(.1, .5, .9)),
    hetECM_15m = median(allSpDat_AM$hetECM_15m)))
pairs(emm_consp_hetAM, by = "consp_15m")
pairs(emm_consp_hetAM, by = "hetAM_15m")
# Across levels of hetECM BA, holding hetAM constant
emm_consp_hetECM <- emmeans(AM_sp_15mB,
  ~ consp_15m * hetECM_15m,
  type = "response",
  at = list(consp_15m  = quantile(allSpDat_AM$consp_15m, c(.1, .5, .9)),
    hetECM_15m = quantile(allSpDat_AM$hetECM_15m, c(.1, .5, .9)),
    hetAM_15m  = median(allSpDat_AM$hetAM_15m)))
pairs(emm_consp_hetECM, by = "consp_15m")
pairs(emm_consp_hetECM, by = "hetECM_15m")

# Gives predicted survival for ECM seedlings at low/medium/high levels of all three BA predictors (i.e., predictions under specific values of confounders)
# Across levels of hetAM BA, holding hetECM constant
emm_consp_hetECM <- emmeans(ECM_sp_15mB,
                           ~ consp_15m * hetAM_15m,
                           type = "response",
             at = list(consp_15m = quantile(allSpDat_ECM$consp_15m, c(.1, .5, .9)),
                      hetAM_15m = quantile(allSpDat_ECM$hetAM_15m, c(.1, .5, .9)),
                      hetECM_15m = median(allSpDat_ECM$hetECM_15m)))
pairs(emm_consp_hetECM, by = "consp_15m")
pairs(emm_consp_hetECM, by = "hetAM_15m")
# Across levels of hetECM BA, holding hetAM constant
emm_consp_hetECM <- emmeans(ECM_sp_15mB,
                            ~ consp_15m * hetECM_15m,
                            type = "response",
                  at = list(consp_15m  = quantile(allSpDat_ECM$consp_15m, c(.1, .5, .9)),
                            hetECM_15m = quantile(allSpDat_ECM$hetECM_15m, c(.1, .5, .9)),
                            hetAM_15m  = median(allSpDat_ECM$hetAM_15m)))
pairs(emm_consp_hetECM, by = "consp_15m")
pairs(emm_consp_hetECM, by = "hetECM_15m")


############ Marginal effects / average partial effects
# Computes average marginal effect (AME) of each continuous predictor
# Gives the closest analog to “the effect of X on survival,” even in the presence of interactions
ame_continuous <- function(model, var){
  emtrends(model, specs = ~1, var = var)
}
##### AM Seedlings
# Marginal effect of conspecific BA
ame_consp <- ame_continuous(AM_sp_15mB, "consp_15m")
ame_consp %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific AM BA
ame_hetAM <- ame_continuous(AM_sp_15mB, "hetAM_15m")
ame_hetAM %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific ECM BA
ame_hetECM <- ame_continuous(AM_sp_15mB, "hetECM_15m")
ame_hetECM %>%
  summary(infer = TRUE)
##### ECM Seedlings
# Marginal effect of conspecific BA
ame_consp <- ame_continuous(ECM_sp_15mB, "consp_15m")
ame_consp %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific AM BA
ame_hetAM <- ame_continuous(ECM_sp_15mB, "hetAM_15m")
ame_hetAM %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific ECM BA
ame_hetECM <- ame_continuous(ECM_sp_15mB, "hetECM_15m")
ame_hetECM %>%
  summary(infer = TRUE)


################################################## 
##### Plot predictions for Interactive spillover

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A–D
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
mid_consp  <- median(allSpDat_AM$consp_15m,  na.rm = TRUE)
mid_hetAM  <- median(allSpDat_AM$hetAM_15m,  na.rm = TRUE)
mid_hetECM <- median(allSpDat_AM$hetECM_15m, na.rm = TRUE)

hetAM_quart  <- quantile(allSpDat_AM$hetAM_15m,  probs = c(.1,.5,.9), na.rm = TRUE)
hetECM_quart <- quantile(allSpDat_AM$hetECM_15m, probs = c(.1,.5,.9), na.rm = TRUE)

### Colors
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
        plot.title = element_text(size=15, face="bold"))

########################## AM Seedlings
# =========================================================
## -------- PANEL A — Classic feedback (consp × seedling type)
## Predict 
pred_A <- ggpredict(AM_sp_15mB,
  terms = c("consp_15m",
            paste0("hetAM_15m[",  mid_hetAM,  "]"),
            paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = TRUE)

## Build dataframe cleanly
df_A <- data.frame(
  x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high)

## Single AM line
df_A$combo <- factor("AM_Medium", levels = c("AM_Medium"))

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted,
                     color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(size = 0.9) +
  scale_color_manual(values = color_map_combo,
                     guide = "none") +
  scale_fill_manual(values = color_map_combo,
                    guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year AM seedling counts") +
  base_theme

plot_A

# =========================================================
## ---------------------PANEL B1 ------ Consp × hetAM 
## Predict 
pred_B1_AM <- ggpredict(AM_sp_15mB,
  terms = c("consp_15m",
    paste0("hetAM_15m[", paste(hetAM_quart, collapse=","), "]"),
    paste0("hetECM_15m[", mid_hetECM, "]")),bias_correction = TRUE)

## Build dataframe cleanly
df_B1_AM <- data.frame(
  x         = pred_B1_AM$x,
  predicted = pred_B1_AM$predicted,
  conf.low  = pred_B1_AM$conf.low,
  conf.high = pred_B1_AM$conf.high,
  hetAM_raw = pred_B1_AM$group)

df_B1_AM$hetAM <- factor(
  df_B1_AM$hetAM_raw,
  levels = unique(pred_B1_AM$group),
  labels = c("Low","Medium","High"))

df_B1_AM$combo <- factor(
  paste0("AM_", df_B1_AM$hetAM),
  levels = c("AM_Low","AM_Medium","AM_High"))


## Plot
plot_B1_AM <- ggplot(df_B1_AM,
                     aes(x = x, y = predicted,
                         color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_combo,
                     name = "Heterospecific AM BA", guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year AM seedling counts") +
  base_theme
plot_B1_AM

# =========================================================
## -------- PANEL B2 ----------- Consp × hetECM
## Predict 
pred_B2_AM <- ggpredict(AM_sp_15mB,
  terms = c( "consp_15m",
    paste0("hetECM_15m[", paste(hetECM_quart, collapse=","), "]"),
    paste0("hetAM_15m[", mid_hetAM, "]")),bias_correction = TRUE)

## Build dataframe cleanly
df_B2_AM <- data.frame(
  x         = pred_B2_AM$x,
  predicted = pred_B2_AM$predicted,
  conf.low  = pred_B2_AM$conf.low,
  conf.high = pred_B2_AM$conf.high,
  hetECM_raw = pred_B2_AM$group)

df_B2_AM$hetECM <- factor(
  df_B2_AM$hetECM_raw,
  levels = unique(pred_B2_AM$group),
  labels = c("Low","Medium","High"))

df_B2_AM$combo <- factor(
  paste0("ECM_", df_B2_AM$hetECM),
  levels = c("ECM_Low","ECM_Medium","ECM_High"))

color_map_AM_ECMlevels <- c(
  "ECM_Low"    = AM_grad[1],
  "ECM_Medium" = AM_grad[2],
  "ECM_High"   = AM_grad[3])

## Plot
plot_B2_AM <- ggplot(df_B2_AM,
                     aes(x = x, y = predicted,
                         color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_AM_ECMlevels,
                     name = "Heterospecific ECM BA", guide = "none") +
  scale_fill_manual(values = color_map_AM_ECMlevels, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year AM seedling counts") +
  base_theme
plot_B2_AM


########################## ECM Seedlings
# =========================================================
## -------- PANEL A — Classic feedback (consp × seedling type)
## Predict 
pred_A <- ggpredict(ECM_sp_15mB,
  terms = c("consp_15m",
            paste0("hetAM_15m[", mid_hetAM, "]"),
            paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = FALSE)

## Build dataframe cleanly
df_A <- data.frame(
  x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high)

## Single AM line
df_A$combo <- factor("ECM_Medium", levels = c("ECM_Medium"))

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted,
                     color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_combo,
                     guide = "none") +
  scale_fill_manual(values = color_map_combo,
                    guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year ECM seedling counts") +
  base_theme

plot_A

# =========================================================
## --------------------- PANEL B1 — ECM seedlings: Consp × hetAM
pred_B1_ECM <- ggpredict(ECM_sp_15mB,
  terms = c("consp_15m",
    paste0("hetAM_15m[", paste(hetAM_quart, collapse = ","), "]"),
    paste0("hetECM_15m[", mid_hetECM, "]")),
  bias_correction = FALSE)

## Build dataframe
df_B1_ECM <- data.frame(
  x         = pred_B1_ECM$x,
  predicted = pred_B1_ECM$predicted,
  conf.low  = pred_B1_ECM$conf.low,
  conf.high = pred_B1_ECM$conf.high,
  hetAM_raw = pred_B1_ECM$group)

df_B1_ECM$hetAM <- factor(
  df_B1_ECM$hetAM_raw,
  levels = unique(pred_B1_ECM$group),
  labels = c("Low", "Medium", "High"))

df_B1_ECM$combo <- factor(
  paste0("AM_", df_B1_ECM$hetAM),
  levels = c("AM_Low", "AM_Medium", "AM_High"))

## ECM = blue gradient
color_map_ECM_AMlevels <- c(
  "AM_Low"    = ECM_grad[1],
  "AM_Medium" = ECM_grad[2],
  "AM_High"   = ECM_grad[3])

## Plot
plot_B1_ECM <- ggplot(df_B1_ECM,
  aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_ECM_AMlevels,
                     name = "Heterospecific AM BA", guide = "none") +
  scale_fill_manual(values = color_map_ECM_AMlevels, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
    y = "Predicted one-year ECM seedling counts") +
  base_theme

plot_B1_ECM


# =========================================================
## -------- PANEL B2 — ECM seedlings: Consp × hetECM
pred_B2_ECM <- ggpredict(ECM_sp_15mB,
  terms = c("consp_15m",
    paste0("hetECM_15m[", paste(hetECM_quart, collapse = ","), "]"),
    paste0("hetAM_15m[", mid_hetAM, "]")),
  bias_correction = FALSE)

## Build dataframe
df_B2_ECM <- data.frame(
  x          = pred_B2_ECM$x,
  predicted  = pred_B2_ECM$predicted,
  conf.low   = pred_B2_ECM$conf.low,
  conf.high  = pred_B2_ECM$conf.high,
  hetECM_raw = pred_B2_ECM$group)

df_B2_ECM$hetECM <- factor(
  df_B2_ECM$hetECM_raw,
  levels = unique(pred_B2_ECM$group),
  labels = c("Low", "Medium", "High"))

df_B2_ECM$combo <- factor(
  paste0("ECM_", df_B2_ECM$hetECM),
  levels = c("ECM_Low", "ECM_Medium", "ECM_High"))

## ECM blue gradient
color_map_ECM_ECMlevels <- c(
  "ECM_Low"    = ECM_grad[1],
  "ECM_Medium" = ECM_grad[2],
  "ECM_High"   = ECM_grad[3])

## Plot
plot_B2_ECM <- ggplot(df_B2_ECM,
  aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_ECM_ECMlevels,
                     name = "Heterospecific ECM BA", guide = "none") +
  scale_fill_manual(values = color_map_ECM_ECMlevels, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
    y = "Predicted one-year ECM seedling counts") +
  base_theme

plot_B2_ECM


#########################################################
#########################################################
################# Mycorrhizal Type Models - Proportional Basal Area

########## ########## ########## ########## ########## 
##### Null model - AM (Proportional BA)
AM_sp_15m_nullprop <- glmmTMB(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear) + (1|Sp),
                          family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15m_nullprop)
# R2
r2_nullAMprop <- performance::r2(AM_sp_15m_nullprop)
r2_nullAMprop

##### Null model - ECM (Proportional BA)
ECM_sp_15m_nullprop <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15m_nullprop)
# R2
r2_nullECMprop <- performance::r2(ECM_sp_15m_nullprop)
r2_nullECMprop

########## ########## ########## ########## ########## 
#### Panel A - AM (Proportional BA)
AM_sp_15mpropA <- glmmTMB(Year_1_counts ~ propConsp +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15mpropA)
## Significant negative effect of ConspBA 
# VIF and R2
vif_AMaProp <- compute_vif(~ propConsp +
                         logPrior + log0yr, allSpDat_AM)
vif_AMaProp # All <2
r2_A_AMprop <- performance::r2(AM_sp_15mpropA)
r2_A_AMprop

#### Panel A - ECM (Proportional BA)
ECM_sp_15mpropA <- glmer(Year_1_counts ~ propConsp +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15mpropA)
## Main effect of ConspBA is not significant
# VIF and R2
vif_ECMaProp <- compute_vif(~ propConsp +
                             logPrior + log0yr, allSpDat_ECM)
vif_ECMaProp # All <2
r2_A_ECMprop <- performance::r2(ECM_sp_15mpropA)
r2_A_ECMprop

########## ########## ########## ########## ########## 
#### Panel B - AM (Proportional BA)
AM_sp_15mpropB <- glmmTMB(Year_1_counts ~ propConsp + 
                      propHetAM + propConsp:propHetAM + 
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat_AM)
summary(AM_sp_15mpropB)
## Significant negative main effects of propConsp and propHetAM
## propConsp:propAM is positive and significant
# VIF and R2
vif_AMbProp <- compute_vif(~ propConsp + propHetAM + propConsp:propHetAM +
                             logPrior + log0yr, allSpDat_AM)
vif_AMbProp # All <2
r2_B_AMprop <- performance::r2(AM_sp_15mpropB)
r2_B_AMprop

#### Panel B - ECM (Proportional BA)
ECM_sp_15mpropB <- glmer(Year_1_counts ~ propConsp + 
                       propHetAM + propConsp:propHetAM + 
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15mpropB)
## Non-significant main effects of propConsp and propHetAM
## propConsp:propHetAM is negative and significant
# VIF and R2
vif_ECMbProp <- compute_vif(~ propConsp + 
                             propHetAM + propConsp:propHetAM + 
                             logPrior + log0yr, allSpDat_ECM)
vif_ECMbProp # All <2
r2_B_ECMprop <- performance::r2(ECM_sp_15mpropB)
r2_B_ECMprop

########## ########## ########## ########## ########## 
## Model comparison - AM
AIC(AM_sp_15m_nullprop, AM_sp_15mpropA, AM_sp_15mpropB) # AM_sp_15mpropB is lowest
## Model comparison - ECM
AIC(ECM_sp_15m_nullprop, ECM_sp_15mpropA, ECM_sp_15mpropB) # ECM_sp_15mpropB is lowest

### Diagnostics - AM
sim_AMprop_B <- DHARMa::simulateResiduals(AM_sp_15mpropB, n = 1000)
DHARMa::testDispersion(sim_AMprop_B) # No overdispersion
DHARMa::testZeroInflation(sim_AMprop_B) # No zero-inflation
DHARMa::testUniformity(sim_AMprop_B) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_AMprop_B) # ~0.32% of observations are outliers

### Diagnostics - ECM
sim_ECMprop_B <- DHARMa::simulateResiduals(ECM_sp_15mpropB, n = 1000)
DHARMa::testDispersion(sim_ECMprop_B) # No overdispersion
DHARMa::testZeroInflation(sim_ECMprop_B) # No zero-inflation
DHARMa::testUniformity(sim_ECMprop_B) # Non-uniform, but normal with large datasets 
DHARMa::testOutliers(sim_ECMprop_B) # ~0.12% of observations are outliers

########## Plotting effect sizes - for Panel B - proportional BA
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate( component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propHetAM" ~ "Prop. AM BA",
        term == "propConsp:propHetAM" ~ "Prop. Conspecific BA:Prop AM BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(component = factor(component, levels = c(
        "Prop. Conspecific BA", "Prop. AM BA", 
        "Prop. Conspecific BA:Prop AM BA",
        "Log(Recruit Density)","Log(Prior Density)")))
}

## Extract for both models
fixef_AM  <- extract_fixef(AM_sp_15mpropB, "AM")
fixef_ECM <- extract_fixef(ECM_sp_15mpropB, "ECM")

## Combine
fixef_df <- bind_rows(fixef_AM, fixef_ECM)

# Define desired vertical order
desired_order <- c("Prop. Conspecific BA", "Prop. AM BA", 
                   "Prop. Conspecific BA:Prop AM BA",
                   "Log(Recruit Density)","Log(Prior Density)")

# Prepare data
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = rev(desired_order)),
    Model = factor(Model, levels = c("AM", "ECM")),
    Significant = factor(Significant, levels = c("yes", "no")))

######## Plot
# Small vertical offset for AM vs ECM
offset <- 0.09
ggplot(fixef_df, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes( xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                width = 0.1, color = "black",
    position = position_nudge(y = ifelse(fixef_df$Model == "AM", offset, -offset))) +
  geom_point(aes(shape = Model, fill = Significant),
             size = 3, color = "black", stroke = 1,
    position = position_nudge(y = ifelse(fixef_df$Model == "AM", offset, -offset))) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes =
                               list(shape = 21, color = "black", size = 3.5)))

############ Marginal effects / average partial effects
# Computes average marginal effect (AME) of each continuous predictor
# Gives the closest analog to “the effect of X on survival,” even in the presence of interactions
ame_continuous <- function(model, var){
  emtrends(model, specs = ~1, var = var)
}
##### AM Seedlings
# Marginal effect of proportional conspecific BA
ame_consp <- ame_continuous(AM_sp_15mpropB, "propConsp")
ame_consp %>%
  summary(infer = TRUE)
# Marginal effect of proportional AM BA
ame_hetAM <- ame_continuous(AM_sp_15mpropB, "propHetAM")
ame_hetAM %>%
  summary(infer = TRUE)
##### ECM Seedlings
# Marginal effect of conspecific BA
ame_consp <- ame_continuous(ECM_sp_15mpropB, "propConsp")
ame_consp %>%
  summary(infer = TRUE)
# Marginal effect of heterospecific AM BA
ame_hetAM <- ame_continuous(ECM_sp_15mpropB, "propHetAM")
ame_hetAM %>%
  summary(infer = TRUE)


################################################## 
##### Plot predictions for Interactive spillover

# These plots are not showing the marginal effects, but rather are conditional predictions, not derivatives. They are the appropriate method for: 
# showing how effects vary across the interaction space
# illustrating ecological mechanisms tested in Panels A–D
# visualizing the implications of our best-supported model 

# =========================================================
# 1. Central values and quartiles
# =========================================================
mid_propConsp <- median(allSpDat_AM$propConsp,  na.rm = TRUE)
mid_propHetAM <- median(allSpDat_AM$propHetAM,  na.rm = TRUE)

propHetAM_quart <- quantile(allSpDat_AM$propHetAM,
  probs = c(.1, .5, .9), na.rm = TRUE)

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
        plot.title = element_text(size=15, face="bold"))

########################## AM Seedlings
# =========================================================
## -------- PANEL A — Classic feedback (consp × seedling type)
pred_A <- ggpredict(AM_sp_15mpropB,
  terms = c("propConsp",
    paste0("propHetAM[", mid_propHetAM, "]")),
  bias_correction = TRUE)

df_A <- data.frame(x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high)

df_A$combo <- factor("AM_Medium", levels = "AM_Medium")

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year AM seedling counts") +
  base_theme
plot_A

# =========================================================
## ---------------------PANEL B1 ------ Consp × hetAM 
pred_B1_AM <- ggpredict(AM_sp_15mpropB,
  terms = c( "propConsp",
    paste0("propHetAM[", paste(propHetAM_quart, collapse = ","), "]")),
  bias_correction = TRUE)

df_B1_AM <- data.frame(
  x          = pred_B1_AM$x,
  predicted  = pred_B1_AM$predicted,
  conf.low   = pred_B1_AM$conf.low,
  conf.high  = pred_B1_AM$conf.high,
  propHetAM_raw = pred_B1_AM$group)

df_B1_AM$propHetAM <- factor(
  df_B1_AM$propHetAM_raw,
  levels = unique(pred_B1_AM$group),
  labels = c("Low","Medium","High"))

df_B1_AM$combo <- factor(
  paste0("AM_", df_B1_AM$propHetAM),
  levels = c("AM_Low","AM_Medium","AM_High"))

plot_B1_AM <- ggplot(df_B1_AM,
                     aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_combo,
                     name = "Prop. AM BA") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year AM seedling counts") +
  base_theme

plot_B1_AM


########################## ECM Seedlings
# =========================================================
## -------- PANEL A — Classic feedback (consp × seedling type)
## Predict 
pred_A <- ggpredict(ECM_sp_15mpropB,
  terms = c("propConsp",
    paste0("propHetAM[", mid_propHetAM, "]")),
  bias_correction = FALSE)

df_A <- data.frame(
  x         = pred_A$x,
  predicted = pred_A$predicted,
  conf.low  = pred_A$conf.low,
  conf.high = pred_A$conf.high)

df_A$combo <- factor("ECM_Medium", levels = "ECM_Medium")

plot_A <- ggplot(df_A,
                 aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year ECM seedling counts") +
  base_theme

plot_A


# =========================================================
## --------------------- PANEL B1 — ECM seedlings: Consp × hetAM
pred_B1_ECM <- ggpredict(ECM_sp_15mpropB,
  terms = c("propConsp",
    paste0("propHetAM[", paste(propHetAM_quart, collapse = ","), "]")),
  bias_correction = FALSE)

df_B1_ECM <- data.frame(x          = pred_B1_ECM$x,
  predicted  = pred_B1_ECM$predicted,
  conf.low   = pred_B1_ECM$conf.low,
  conf.high  = pred_B1_ECM$conf.high,
  propHetAM_raw = pred_B1_ECM$group)

df_B1_ECM$propHetAM <- factor(
  df_B1_ECM$propHetAM_raw,
  levels = unique(pred_B1_ECM$group),
  labels = c("Low","Medium","High"))

df_B1_ECM$combo <- factor(
  paste0("AM_", df_B1_ECM$propHetAM),
  levels = c("AM_Low","AM_Medium","AM_High"))

color_map_ECM_AMlevels <- c(
  "AM_Low"    = ECM_grad[1],
  "AM_Medium" = ECM_grad[2],
  "AM_High"   = ECM_grad[3])

plot_B1_ECM <- ggplot(df_B1_ECM,
                      aes(x = x, y = predicted, color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = color_map_ECM_AMlevels,
                     name = "Prop. AM BA") +
  scale_fill_manual(values = color_map_ECM_AMlevels, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year ECM seedling counts") +
  base_theme

plot_B1_ECM


################################################
################################################
##### Simple figure illustrating one-year survivorship

# Add survivorship column, dropping cases where Year_0_counts == 0
allSpDat <- allSpDat %>%
  mutate(survivorship = ifelse(Year_0_counts > 0,
                               Year_1_counts / Year_0_counts,
                               NA))

# Example lookup tables (fill in with your actual names)
species_names <- c("bonro" = "Quercus velutina",
                   "wo" = "Quercus alba",
                   "nro" = "Quercus rubra",
  "rm" = "Acer rubrum",
  "bc" = "Prunus serotina",
  "ab" = "Fagus grandifolia",
  "sm" = "Acer saccharum",
  "iw" = "Ostrya virginiana",
  "wa" = "Fraxinus americana",
  "bs" = "Tilia americana")

site_names <- c("1" = "49.01",
  "2" = "59.12",
  "3" = "8.21",
  "4" = "39.21",
  "5" = "16.35",
  "6" = "60.35",
  "7" = "7.40",
  "8" = "40.40",
  "9" = "21.43",
  "10" = "25.43",
  "11" = "6.45",
  "12" = "24.45")

# Compute site means
site_means <- allSpDat %>%
  group_by(SiteID, Sp) %>%
  summarise(mean_surv = mean(survivorship, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Sp = recode(Sp, !!!species_names),
    SiteID = recode(as.character(SiteID), !!!site_names),
    # Set factor levels to preserve the order in your lookup vectors
    Sp = factor(Sp, levels = species_names),
    SiteID = factor(SiteID, levels = site_names)
  )

# Create italic labels for species
species_labels <- setNames(
  lapply(species_names, function(x) bquote(italic(.(x)))),
  species_names
)

# Plot heatmap with italic species names
ggplot(site_means, aes(x = SiteID, y = Sp, fill = mean_surv)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "cividis", na.value = "grey90") +
  scale_y_discrete(labels = species_labels) +
  labs(x = "Site", y = "", fill = "Mean\nsurvivorship") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                  size = 11, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12),
        legend.position = "right")


