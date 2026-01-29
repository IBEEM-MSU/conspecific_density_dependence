
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 27 January 2026
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

# Marginal prediction R² (log-scale)
pred_marg <- predict(all_sp_15mA, type = "response", re.form = NA)
obs <- allSpDat$Year_1_counts
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mA, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond


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

# Marginal prediction R² (log-scale)
pred_marg <- predict(all_sp_15mB, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mB, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond

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

# Marginal prediction R² (log-scale)
pred_marg <- predict(all_sp_15mC, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mC, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond


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

# Marginal prediction R² (log-scale)
pred_marg <- predict(all_sp_15mD, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mD, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond

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

# Marginal prediction R² (log-scale)
pred_marg <- predict(all_sp_15mBCD, type = "response", re.form = NA)
r_marg <- cor(log1p(obs), log1p(pred_marg), use = "complete.obs")
R2_pred_marg <- r_marg^2
R2_pred_marg
# Conditional prediction r² (log-scale)
pred_cond <- predict(all_sp_15mBCD, type = "response")
r_cond <- cor(log1p(obs), log1p(pred_cond), use = "complete.obs")
R2_pred_cond <- r_cond^2
R2_pred_cond

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
    plot.title = element_text(size=15, face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

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
  scale_y_continuous(limits = c(0.5, 3.5), breaks = c(1,2,3)) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
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
  scale_y_continuous(limits = c(0.2, 4.9), breaks = c(1,2,3,4)) +
  scale_color_manual(values = color_map_combo, name = "AM BA Level", guide = "none") +
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
  scale_y_continuous(limits = c(0.2, 4.9), breaks = c(1,2,3,4)) +
  scale_color_manual(values = color_map_combo,
                     name = "ECM BA Level", guide="none") +
  scale_fill_manual(values = color_map_combo, guide="none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_B_ECM

#############################
##### Combined plot
df_all <- bind_rows(
  df_B_AM  %>% mutate(neighbor_type = "AM neighbors"),
  df_B_ECM %>% mutate(neighbor_type = "ECM neighbors")
)

df_all$seed_type <- factor(df_all$seed_type, levels = c("AM", "ECM"))
df_all$neighbor_type <- factor(df_all$neighbor_type,
                               levels = c("AM neighbors", "ECM neighbors"))


ggplot(df_all,
       aes(x = x, y = predicted,
           color = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(seed_type ~ neighbor_type,
             labeller = labeller(
               seed_type = c(AM = "AM seedlings", ECM = "ECM seedlings")
             )) +
  scale_y_continuous(limits = c(0.2, 4.9), breaks = c(1,2,3,4)) +
  scale_color_manual(values = color_map_combo, guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Conspecific basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme


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
  scale_y_continuous(limits = c(0.8, 2.5), breaks = c(1,2)) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
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
  scale_y_continuous(limits = c(0.8, 2.5), breaks = c(1,2)) +
  scale_color_manual(values = color_map_combo,
                     name = "Seedling type",
                     labels = c("AM","ECM"), guide = "none") +
  scale_fill_manual(values = color_map_combo, guide="none") +
  labs(x = "Heterospecific ECM basal area (m²)",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_C_ECM

######### Combined plot

(plot_C_AM | plot_C_ECM) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(axis.title.y = element_text(),
    axis.text.y  = element_text(),
    axis.ticks.y = element_line())


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
  scale_y_continuous(limits = c(0.6, 2.1), breaks = c(1,2)) +
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
  scale_y_continuous(limits = c(0.6, 2.1), breaks = c(1,2)) +
  scale_color_manual(values = color_map_combo, name = "Prop. AM BA",
                     guide = "none") +
  scale_fill_manual(values = color_map_combo, guide = "none") +
  labs(x = "Proportion of conspecific basal area",
       y = "Predicted one-year seedling counts") +
  base_theme
plot_B_AM





