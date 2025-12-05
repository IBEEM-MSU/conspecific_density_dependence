
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### *REVISED VERSIONS BASED ON CHRIS' FEEDBACK*
####### Date: 17 November 2025
#################################################

## Read in packages
library(tidyverse) 
library(lme4) 
library(emmeans)
library(dplyr)
library(ggplot2)
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
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1

## Create two new columns, hetAM_15m and hetECM_15m, that exclude conspecific BA from calculations of AM and ECM associating basal area

# Define species as AM or ECM associating
AM_species  <- c("bc", "rm", "sm", "wa")
ECM_species <- c("ab", "bonro", "bs", "iw", "nro", "wo")

# Create placeholder columns
allSpDat$hetAM_15m  <- NA
allSpDat$hetECM_15m <- NA

# Match and substract
allSpDat <- allSpDat %>%
  rowwise() %>%
  mutate(
    # Identify the column containing the focal species' basal area
    focal_col = paste0(Sp, "_15m"),
    focal_value = get(focal_col),
    # Calculate heterospecific AM and ECM basal area
    hetAM_15m = ifelse(Sp %in% AM_species,
                       AM_15m - focal_value,   # subtract conspecific
                       AM_15m),                 # all AM are heterospecific for ECM species
    hetECM_15m = ifelse(Sp %in% ECM_species,
                        ECM_15m - focal_value,  # subtract conspecific
                        ECM_15m)) %>%        # all ECM are heterospecific for AM species
  
  ungroup() %>%
  select(-focal_col, -focal_value)   # cleanup temporary columns


### glmmTMB still has log-link built in with family = nbinom2

#########################################################
#########################################################
#########################################################
################# All Species Models - Raw Basal Area

### Compared AIC values with and without logPrior - significantly greater support for models including logPrior, so we should leave in
### Also compared best negative binomial model in each case to a Poisson - considerably better fit and greater support for use of the negative binomial distribution

##  Helper function: VIF for glmmTMB fixed effects
compute_vif <- function(formula, data) {
  mm <- model.matrix(formula, data = data)
  mm_df <- as.data.frame(mm[, -1])  # Remove intercept
  fake_y <- rnorm(nrow(mm_df))
  lm_fake <- lm(fake_y ~ ., data = mm_df)
  car::vif(lm_fake)
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

########## Plotting effect sizes - for Panel A
### Extract fixed effects 
fixef_df_A <- broom.mixed::tidy(
  all_sp_15mA,
  effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Relabel terms
fixef_df_A <- fixef_df_A %>%
  mutate(component = case_when(
    term == "consp_15m" ~ "Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "consp_15m:Seedling_mycECM" ~ "Consp. BA:Mycorrhizal type",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
### Desired order (only terms in Model A)
desired_order_A <- c( "Conspecific BA",
  "Mycorrhizal type","Consp. BA:Mycorrhizal type",
  "Log(Recruit Density)","Log(Prior Density)")
fixef_df_A <- fixef_df_A %>%
  mutate(component = factor(component, levels = desired_order_A))
### Add significance flag
fixef_df_A <- fixef_df_A %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
### Reverse order for plotting
fixef_df_A$component <- factor(fixef_df_A$component, levels = rev(desired_order_A))
### Plot
ggplot(fixef_df_A, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error
    ),width = 0.1, orientation = "y", color = "black") +
  geom_point(aes(shape = Significant, fill = Significant),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"))

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

########## Plotting effect sizes - for Panel B
### Extract fixed effects 
fixef_df_B <- broom.mixed::tidy(
  all_sp_15mB,
  effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Relabel terms
fixef_df_B <- fixef_df_B %>%
  mutate(component = case_when(
    term == "consp_15m" ~ "Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "hetAM_15m" ~ "AM BA",
    term == "hetECM_15m" ~ "ECM BA",
    term == "consp_15m:Seedling_mycECM" ~ "Consp. BA:Mycorrhizal type",
    term == "consp_15m:hetAM_15m" ~ "Consp. BA:AM BA",
    term == "consp_15m:hetECM_15m" ~ "Consp. BA:ECM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
### Desired order (only terms in Model B)
desired_order_B <- c( "Conspecific BA","AM BA","ECM BA",
                      "Mycorrhizal type","Consp. BA:AM BA",
                      "Consp. BA:ECM BA","Consp. BA:Mycorrhizal type",
                      "Log(Recruit Density)","Log(Prior Density)")
fixef_df_B <- fixef_df_B %>%
  mutate(component = factor(component, levels = desired_order_B))
### Add significance flag
fixef_df_B <- fixef_df_B %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
### Reverse order for plotting
fixef_df_B$component <- factor(fixef_df_B$component, levels = rev(desired_order_B))
### Plot
ggplot(fixef_df_B, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error
  ),width = 0.1, orientation = "y", color = "black") +
  geom_point(aes(shape = Significant, fill = Significant),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

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

########## Plotting effect sizes - for Panel C
### Extract fixed effects 
fixef_df_C <- broom.mixed::tidy(
  all_sp_15mC,
  effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Relabel terms
fixef_df_C <- fixef_df_C %>%
  mutate(component = case_when(
    term == "consp_15m" ~ "Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "hetAM_15m" ~ "AM BA",
    term == "hetECM_15m" ~ "ECM BA",
    term == "consp_15m:Seedling_mycECM" ~ "Consp. BA:Mycorrhizal type",
    term == "Seedling_mycECM:hetAM_15m" ~ "Mycorrhizal type:AM BA",
    term == "Seedling_mycECM:hetECM_15m" ~ "Mycorrhizal type:ECM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
### Desired order (only terms in Model B)
desired_order_C <- c( "Conspecific BA","AM BA","ECM BA",
                      "Mycorrhizal type","Consp. BA:Mycorrhizal type",
                      "Mycorrhizal type:AM BA","Mycorrhizal type:ECM BA",
                      "Log(Recruit Density)","Log(Prior Density)")
fixef_df_C <- fixef_df_C %>%
  mutate(component = factor(component, levels = desired_order_C))
### Add significance flag
fixef_df_C <- fixef_df_C %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
### Reverse order for plotting
fixef_df_C$component <- factor(fixef_df_C$component, levels = rev(desired_order_C))
### Plot
ggplot(fixef_df_C, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error
  ),width = 0.1, orientation = "y", color = "black") +
  geom_point(aes(shape = Significant, fill = Significant),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

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

########## Plotting effect sizes - for Panel D
### Extract fixed effects 
fixef_df_D <- broom.mixed::tidy(
  all_sp_15mD,
  effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Relabel terms
fixef_df_D <- fixef_df_D %>%
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
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
### Reorder
desired_order_D <- c(
  "Conspecific BA","AM BA","ECM BA","Mycorrhizal type",
  "Consp. BA:AM BA","Consp. BA:ECM BA","Consp. BA:Mycorrhizal type",
  "Mycorrhizal type:AM BA","Mycorrhizal type:ECM BA",
  "Log(Recruit Density)","Log(Prior Density)")
fixef_df_D <- fixef_df_D %>%
  mutate(component = factor(component, levels = desired_order_D))
### Add in significance
fixef_df_D <- fixef_df_D %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
### Reorder for plotting
fixef_df_D$component <- factor(fixef_df_D$component, levels = rev(desired_order_D))
### Plot
ggplot(fixef_df_D, aes(x = estimate, y = component)) +
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
vif_BCD <- compute_vif(
  ~ consp_15m * Seedling_myc +
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
anova(all_sp_15mB, all_sp_15mBCD) # And is significantly better than B

#### Diagnostics and Post-hoc tests for super-interactive model (best supported)
### Diagnostics
simBCD <- DHARMa::simulateResiduals(all_sp_15mBCD, n = 1000)
DHARMa::testDispersion(simBCD) # No overdispersion
DHARMa::testZeroInflation(simBCD) # No zero-inflation
DHARMa::testUniformity(simBCD) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCD) # ~0.28% of observations are outliers - not concerning
### Post-hoc tests
# Define a grid of ecologically meaningful predictor values
emm_BCD <- emmeans(
  all_sp_15mBCD,
  ~ Seedling_myc | consp_15m * hetAM_15m * hetECM_15m,
  type = "response",
  at = list(consp_15m   = quantile(allSpDat$consp_15m, c(.1, .5, .9)),
    hetAM_15m   = quantile(allSpDat$hetAM_15m, c(.1, .5, .9)),
    hetECM_15m  = quantile(allSpDat$hetECM_15m, c(.1, .5, .9))))
emm_BCD
pairs(emm_BCD, by = c("consp_15m", "hetAM_15m", "hetECM_15m"))
# AM vs ECM seedling mycorrhizal type has no meaningful main effect or interactions with conspecific density or heterospecific AM/ECM densities.
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


################### Exploratory figure making

# Panel A predictions
predA <- ggpredict(all_sp_15mA,
  terms = c("consp_15m [all]", "Seedling_myc"))

plotA <- ggplot(predA, aes(x = x, y = predicted, color = group)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.15, color = NA) +
  labs(x = "Conspecific BA", y = "Predicted recruitment",
       color = "Mycorrhizal type", fill = "Mycorrhizal type",
       title = "Panel A – Classic Feedback") +
  theme_bw(base_size = 12)

predB <- ggpredict(all_sp_15mB,
  terms = c("consp_15m [all]", "hetECM_15m [quartiles]", "Seedling_myc"))

plotB <- plot(predB) +
  labs(title = "Panel B – Interactive Spillover")

predB2 <- ggpredict(all_sp_15mB,
                   terms = c("consp_15m [all]", "hetAM_15m [quartiles]", "Seedling_myc"))

plotB2 <- plot(predB2) +
  labs(title = "Panel B – Interactive Spillover")

predC <- ggpredict(all_sp_15mC,
  terms = c("hetECM_15m [all]", "Seedling_myc"))

plotC <- plot(predC) +
  labs(title = "Panel C – Additive Spillover")

predC2 <- ggpredict(all_sp_15mC,
                   terms = c("hetAM_15m [all]", "Seedling_myc"))

plotC2 <- plot(predC2) +
  labs(title = "Panel C – Additive Spillover")

predD <- ggpredict(all_sp_15mD,
  terms = c("consp_15m [all]", "hetECM_15m [quartiles]", "Seedling_myc"))

plotD <- plot(predD) +
  labs(title = "Panel D – Interactive + Additive Spillover")

predBCD <- ggpredict(
  all_sp_15mBCD,
  terms = c("consp_15m [all]", "hetECM_15m [quartiles]", "Seedling_myc")
)

plotBCD <- plot(predBCD) +
  labs(title = "Super Interactive Model (BCD+)")

predBCD <- ggpredict(
  all_sp_15mBCD,
  terms = c("consp_15m [all]", "Seedling_myc")
)

plotBCD <- plot(predBCD) +
  labs(title = "Super Interactive Model (BCD+)")

predBCD <- ggpredict(
  all_sp_15mBCD,
  terms = c("hetECM_15m [all]", "Seedling_myc")
)

plotBCD <- plot(predBCD) +
  labs(title = "Super Interactive Model (BCD+)")

predBCD <- ggpredict(
  all_sp_15mBCD,
  terms = c("hetAM_15m [all]", "Seedling_myc")
)

plotBCD <- plot(predBCD) +
  labs(title = "Super Interactive Model (BCD+)")

############ Marginal effects

ame_continuous <- function(model, var){
  emtrends(model, specs = ~1, var = var)
}

ame_consp <- ame_continuous(all_sp_15mBCD, "consp_15m")
ame_consp

ame_hetAM <- ame_continuous(all_sp_15mBCD, "hetAM_15m")
ame_hetAM

ame_hetECM <- ame_continuous(all_sp_15mBCD, "hetECM_15m")
ame_hetECM

ame_myc <- emmeans(all_sp_15mBCD, pairwise ~ Seedling_myc)
ame_myc

ame_table <- bind_rows(
  tidy(ame_consp)  %>% mutate(term = "consp_15m"),
  tidy(ame_hetAM)  %>% mutate(term = "hetAM_15m"),
  tidy(ame_hetECM) %>% mutate(term = "hetECM_15m"),
  tidy(ame_myc$contrasts) %>% mutate(term = "Seedling_myc (ECM – AM)")
)



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

########## Plotting effect sizes - Panel A (proportional)
### Extract fixed effects 
fixef_df_aProp <- broom.mixed::tidy(
  all_sp_propA, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Map raw term names to clean labels you want shown on the plot
fixef_df_aProp <- fixef_df_aProp %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
# Order components on the y-axis in the order you prefer
desired_order_aProp <- c("Prop. Conspecific BA","Mycorrhizal type", 
                   "Prop. Consp. BA:Myc. type ECM",
                   "Log(Recruit Density)","Log(Prior Density)")
fixef_df_aProp <- fixef_df_aProp %>%
  mutate(component = factor(component, levels = desired_order_aProp))
# Significance flag (alpha = 0.05)
fixef_df_aProp <- fixef_df_aProp %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
# Make your factor for y-axis
fixef_df_aProp$component <- factor(fixef_df_aProp$component, levels = rev(desired_order_aProp))
### Plot
ggplot(fixef_df_aProp, aes(x = estimate, y = component)) +
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

################# ################# ################# 
## Panel B (proportional BA)
all_sp_propB <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propAM + 
                       propConsp:propAM + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propB)
## Significant negative effects of propConsp and propAM
## propConsp:propAM is significant and positive
# VIF and R2
vif_propB <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + propConsp:propAM + 
                            logPrior + log0yr, allSpDat)
vif_propB # All <10
r2_propB <- performance::r2(all_sp_propB)
r2_propB

########## Plotting effect sizes - Panel B (proportional)
### Extract fixed effects 
fixef_df_bProp <- broom.mixed::tidy(
  all_sp_propB, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Map raw term names to clean labels you want shown on the plot
fixef_df_bProp <- fixef_df_bProp %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propAM" ~ "Prop. AM BA",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "propConsp:propAM" ~ "Prop. Consp. BA:Prop. AM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
# Order components on the y-axis in the order you prefer
desired_order_bProp <- c("Prop. Conspecific BA","Prop. AM BA", "Mycorrhizal type", 
                         "Prop. Consp. BA:Prop. AM BA", "Prop. Consp. BA:Myc. type ECM",
                         "Log(Recruit Density)","Log(Prior Density)")
fixef_df_bProp <- fixef_df_bProp %>%
  mutate(component = factor(component, levels = desired_order_bProp))
# Significance flag (alpha = 0.05)
fixef_df_bProp <- fixef_df_bProp %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
# Make your factor for y-axis
fixef_df_bProp$component <- factor(fixef_df_bProp$component, levels = rev(desired_order_bProp))
### Plot
ggplot(fixef_df_bProp, aes(x = estimate, y = component)) +
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

################# ################# ################# 
## Panel C (proportional BA)
all_sp_propC <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propAM + 
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propC)
## Significant positive and negative effects of propConsp and propAM, respectively
## propConsp:Seedling_myc interaction is significant and negative
## Seedling_myc:AM_15m interaction is not significant 
# VIF and R2
vif_propC <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + Seedling_myc:propAM + 
                            logPrior + log0yr, allSpDat)
vif_propC # All <9
r2_propC <- performance::r2(all_sp_propC)
r2_propC

########## Plotting effect sizes - Panel C (proportional)
### Extract fixed effects 
fixef_df_cProp <- broom.mixed::tidy(
  all_sp_propC, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Map raw term names to clean labels you want shown on the plot
fixef_df_cProp <- fixef_df_cProp %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propAM" ~ "Prop. AM BA",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "Seedling_mycECM:propAM" ~ "Prop. AM BA:Myc. type ECM",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
# Order components on the y-axis in the order you prefer
desired_order_cProp <- c("Prop. Conspecific BA","Prop. AM BA", "Mycorrhizal type", 
                         "Prop. Consp. BA:Myc. type ECM","Prop. AM BA:Myc. type ECM",
                         "Log(Recruit Density)","Log(Prior Density)")
fixef_df_cProp <- fixef_df_cProp %>%
  mutate(component = factor(component, levels = desired_order_cProp))
# Significance flag (alpha = 0.05)
fixef_df_cProp <- fixef_df_cProp %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
# Make your factor for y-axis
fixef_df_cProp$component <- factor(fixef_df_cProp$component, levels = rev(desired_order_cProp))
### Plot
ggplot(fixef_df_cProp, aes(x = estimate, y = component)) +
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

################# ################# ################# 
## Panel D (B+C) (proportional BA)
all_sp_propD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propAM + 
                        propConsp:propAM + # Panel B
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propD)
## Significant negative effects of propConsp and propAM
## propConsp:Seedling_myc interaction is marginally significant and positive
## propConsp:propAM interaction is significant and positive
# VIF and R2
vif_propD <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + propConsp:propAM + 
                            Seedling_myc:propAM + 
                            logPrior + log0yr, allSpDat)
vif_propD # All <10
r2_propD <- performance::r2(all_sp_propD)
r2_propD

########## Plotting effect sizes - Panel C (proportional)
### Extract fixed effects 
fixef_df_dProp <- broom.mixed::tidy(
  all_sp_propD, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)
### Map raw term names to clean labels you want shown on the plot
fixef_df_dProp <- fixef_df_dProp %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propAM" ~ "Prop. AM BA",
    term == "propConsp:propAM" ~ "Prop. Consp. BA:Prop. AM BA",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "Seedling_mycECM:propAM" ~ "Prop. AM BA:Myc. type ECM",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))
# Order components on the y-axis in the order you prefer
desired_order_dProp <- c("Prop. Conspecific BA","Prop. AM BA", "Mycorrhizal type", 
                         "Prop. Consp. BA:Prop. AM BA",
                         "Prop. Consp. BA:Myc. type ECM","Prop. AM BA:Myc. type ECM",
                         "Log(Recruit Density)","Log(Prior Density)")
fixef_df_dProp <- fixef_df_dProp %>%
  mutate(component = factor(component, levels = desired_order_dProp))
# Significance flag (alpha = 0.05)
fixef_df_dProp <- fixef_df_dProp %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))
# Make your factor for y-axis
fixef_df_dProp$component <- factor(fixef_df_dProp$component, levels = rev(desired_order_dProp))
### Plot
ggplot(fixef_df_dProp, aes(x = estimate, y = component)) +
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

################# ################# ################# 
## Super interactive model (proportional BA)
all_sp_propBCD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propAM + 
                       propConsp:propAM + # Panel B
                       Seedling_myc:propAM + # Panel C
                       propConsp:Seedling_myc:propAM + # 3-way interaction
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propBCD)
## Significant negative effects of propConsp and propAM 
## propConsp:Seedling_myc, propConsp:propAMm, and Seedling_myc:propAM interactions are all significant and positive!
## propConsp:Seedling_myc:propAM interaction is significant and negative
# VIF and R2
vif_propBCD <- compute_vif( ~ propConsp*Seedling_myc + 
                              propAM + propConsp:propAM + 
                              Seedling_myc:propAM + 
                              propConsp:Seedling_myc:propAM + 
                              logPrior + log0yr, allSpDat)
vif_propBCD # two values >10 due to interactions
r2_propBCD <- performance::r2(all_sp_propBCD)
r2_propBCD

## Model comparison
### Compare AICs 
AIC(all_sp_prop_null, all_sp_propA, all_sp_propB, all_sp_propC,
    all_sp_propD, all_sp_propBCD) # lowest is all_sp_propBCD 
anova(all_sp_propB, all_sp_propBCD) # Super-interactive model is significantly better!

#### Diagnostics and Post-hoc tests for super-interactive model (best supported)
### Diagnostics
simBCDprop <- DHARMa::simulateResiduals(all_sp_propBCD, n = 1000)
DHARMa::testDispersion(simBCDprop) # No overdispersion
DHARMa::testZeroInflation(simBCDprop) # No zero-inflation
DHARMa::testUniformity(simBCDprop) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCDprop) # ~0.25% of observations are outliers
### Post-hoc tests
# Define a grid of ecologically meaningful predictor values
emm_BCDprop <- emmeans(
  all_sp_propBCD,
  ~ Seedling_myc | propConsp * propAM,
  type = "response",
  at = list(propConsp  = quantile(allSpDat$propConsp, c(.1, .5, .9)),
            propAM  = quantile(allSpDat$propAM, c(.1, .5, .9))))
emm_BCDprop
pairs(emm_BCDprop, by = c("propConsp", "propAM"))
# Model diagnostics indicated excellent fit, with no evidence of overdispersion (dispersion = 0.85, p = 0.79) or zero inflation (ratioObsSim = 0.88, p = 0.21). Although DHARMa’s uniformity and outlier tests showed significant deviations, the magnitude was trivial (<0.3% outliers) and expected for large-sample GLMMs, and residual plots showed no structured patterns.
# Post-hoc comparisons revealed no reliable differences in seedling abundance between AM and ECM species across proportional conspecific or heterospecific basal area gradients. All nine AM–ECM contrasts yielded ratios near 1.0 with non-significant p-values except for a single marginal comparison at high conspecific and high AM basal area, which we interpret cautiously given the lack of a consistent pattern.

########## ########## ########## ########## ########## 
########## Plotting effect sizes - Superinteractive model
### Extract fixed effects 
fixef_df <- broom.mixed::tidy(
  all_sp_propBCD, effects = "fixed",
  component = "cond") %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, std.error, p.value)

### Map raw term names to clean labels you want shown on the plot
fixef_df <- fixef_df %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Mycorrhizal type",
    term == "propAM" ~ "Prop. AM BA",
    term == "propConsp:Seedling_mycECM" ~ "Prop. Consp. BA:Myc. type ECM",
    term == "propConsp:propAM" ~ "Prop. Consp. BA:Prop. AM BA",
    term == "Seedling_mycECM:propAM" ~  "Myc. type ECM:Prop. AM BA",
    term == "propConsp:Seedling_mycECM:propAM" ~ 
    "Prop. Consp. BA:Myc. type ECM:\nProp. AM BA",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))

# 3. Order components on the y-axis in the order you prefer
desired_order <- c("Prop. Conspecific BA","Prop. AM BA",
                   "Mycorrhizal type",  "Prop. Consp. BA:Myc. type ECM",
                   "Prop. Consp. BA:Prop. AM BA", "Myc. type ECM:Prop. AM BA",
                   "Prop. Consp. BA:Myc. type ECM:\nProp. AM BA",
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

########## Plotting effect sizes - Panel A
extract_fixef_A <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(component = case_when(
        term == "consp_15m" ~ "Conspecific BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(component = factor(component, levels = c(
        "Conspecific BA",  "Log(Recruit Density)", "Log(Prior Density)")))
}
## Extract for both models
fixef_AM_A  <- extract_fixef_A(AM_sp_15mA, "AM")
fixef_ECM_A <- extract_fixef_A(ECM_sp_15mA, "ECM")
## Combine
fixef_df_A_myc <- bind_rows(fixef_AM_A, fixef_ECM_A)
# Define desired vertical order
desired_order_A_myc <- c("Conspecific BA", 
                   "Log(Recruit Density)", "Log(Prior Density)")
# Prepare data
fixef_df_A_myc <- fixef_df_A_myc %>%
  mutate(component = factor(component, levels = rev(desired_order_A_myc)),
         Model = factor(Model, levels = c("AM", "ECM")),
         Significant = factor(Significant, levels = c("yes", "no")))
######## Plot
offset <- 0.05
ggplot(fixef_df_A_myc, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes( xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                width = 0.1, color = "black",
      position = position_nudge(y = ifelse(fixef_df_A_myc$Model == "AM", offset, -offset))) +
  geom_point(aes(shape = Model, fill = Significant),
             size = 3, color = "black", stroke = 1,
  position = position_nudge(y = ifelse(fixef_df_A_myc$Model == "AM", offset, -offset))) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes =
                               list(shape = 21, color = "black", size = 3.5)))


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

predB <- ggpredict(AM_sp_15mB,
                  terms = c("consp_15m [all]", "hetAM_15m [quartiles]"))

plotB <- plot(predB) +
  labs(title = "Panel B – Interactive Spillover")

predBECM <- ggpredict(ECM_sp_15mB,
                   terms = c("consp_15m [all]", "hetAM_15m [quartiles]"))

plotBECM <- plot(predBECM) +
  labs(title = "Panel B – Interactive Spillover")

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
DHARMa::testOutliers(sim_ECM_B) # ~0.16% of observations are outliers

########## Plotting effect sizes - for Panel B (most complex mycorrhizal version)
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
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
    mutate(
      component = factor(component, levels = c(
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
offset <- 0.05
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

########## Plotting effect sizes - Panel A - proportional BA
extract_fixef_propA <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate( component = case_when(
      term == "propConsp" ~ "Prop. Conspecific BA",
      term == "logPrior" ~ "Log(Prior Density)",
      term == "log0yr" ~ "Log(Recruit Density)",
      TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(component = factor(component, levels = c(
      "Prop. Conspecific BA", 
      "Log(Recruit Density)","Log(Prior Density)")))
}

## Extract for both models
fixef_AM_propA  <- extract_fixef_propA(AM_sp_15mpropA, "AM")
fixef_ECM_propA <- extract_fixef_propA(ECM_sp_15mpropA, "ECM")

## Combine
fixef_df_propA <- bind_rows(fixef_AM_propA, fixef_ECM_propA)

# Define desired vertical order
desired_order_propA <- c("Prop. Conspecific BA", 
                   "Log(Recruit Density)","Log(Prior Density)")

# Prepare data
fixef_df_propA <- fixef_df_propA %>%
  mutate(component = factor(component, levels = rev(desired_order_propA)),
         Model = factor(Model, levels = c("AM", "ECM")),
         Significant = factor(Significant, levels = c("yes", "no")))

######## Plot
# Small vertical offset for AM vs ECM
offset <- 0.05
ggplot(fixef_df_propA, aes(x = estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes( xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                width = 0.1, color = "black",
                position = position_nudge(y = ifelse(fixef_df_propA$Model == "AM", offset, -offset))) +
  geom_point(aes(shape = Model, fill = Significant),
             size = 3, color = "black", stroke = 1,
             position = position_nudge(y = ifelse(fixef_df_propA$Model == "AM", offset, -offset))) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes =
                               list(shape = 21, color = "black", size = 3.5)))

########## ########## ########## ########## ########## 
#### Panel B - AM (Proportional BA)
AM_sp_15mpropB <- glmmTMB(Year_1_counts ~ propConsp + 
                      propAM + propConsp:propAM + 
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat_AM)
summary(AM_sp_15mpropB)
## Significant negative main effects of propConsp and propAM
## propConsp:propAM is positive and significant
# VIF and R2
vif_AMbProp <- compute_vif(~ propConsp +
                             logPrior + log0yr, allSpDat_AM)
vif_AMbProp # All <2
r2_B_AMprop <- performance::r2(AM_sp_15mpropB)
r2_B_AMprop

#### Panel B - ECM (Proportional BA)
ECM_sp_15mpropB <- glmer(Year_1_counts ~ propConsp + 
                       propAM + propConsp:propAM + 
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15mpropB)
## Non-significant main effects of propConsp and propAM
## propConsp:propAM is negative and significant
# VIF and R2
vif_ECMbProp <- compute_vif(~ propConsp +
                             logPrior + log0yr, allSpDat_ECM)
vif_ECMbProp # All <2
r2_B_ECMprop <- performance::r2(ECM_sp_15mpropB)
r2_B_ECMprop


########## ########## ########## ########## ########## 
## Model comparison - AM
AIC(AM_sp_15m_nullprop, AM_sp_15mpropA, AM_sp_15mpropB) # AM_sp_15mpropB is lowest
## Model comparison - ECM
AIC(ECM_sp_15m_nullprop, ECM_sp_15mpropA, ECM_sp_15mpropB) # ECM_sp_15mpropB is lowest
anova(ECM_sp_15mpropA, ECM_sp_15mpropB) # And is significantly more supported

### Diagnostics - AM
sim_AMprop_B <- DHARMa::simulateResiduals(AM_sp_15mpropB, n = 1000)
DHARMa::testDispersion(sim_AMprop_B) # No overdispersion
DHARMa::testZeroInflation(sim_AMprop_B) # No zero-inflation
DHARMa::testUniformity(sim_AMprop_B) # Non-uniform, but normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(sim_AMprop_B) # ~0.37% of observations are outliers

### Diagnostics - ECM
sim_ECMprop_B <- DHARMa::simulateResiduals(ECM_sp_15mpropB, n = 1000)
DHARMa::testDispersion(sim_ECMprop_B) # No overdispersion
DHARMa::testZeroInflation(sim_ECMprop_B) # No zero-inflation
DHARMa::testUniformity(sim_ECMprop_B) # Non-uniform, but normal with large datasets 
DHARMa::testOutliers(sim_ECMprop_B) # ~0.09% of observations are outliers

########## Plotting effect sizes - for Panel B - proportional BA
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate( component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propAM" ~ "Prop. AM BA",
        term == "propConsp:propAM" ~ "Prop. Conspecific BA:Prop AM BA",
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
offset <- 0.05
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


#########################################################
#########################################################
################# Species-Level Model - RM - Raw Basal Area

##### Subset RM data
rm_surv <- subset(allSpDat, Sp == "rm")

########## ########## ########## ########## ########## 
##### Null model - RM
rm_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                     (1|SiteID) + (1|survYear),
                     family = nbinom2, data = rm_surv) 
summary(rm_15m_null)
# R2
rm_15m_null_r2 <- performance::r2(rm_15m_null)
rm_15m_null_r2

########## ########## ########## ########## ########## 
##### Partial Panel A - CNDD - RM
rm_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = rm_surv)
rm_15mA2 <- glmer(Year_1_counts ~ consp_15m +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = rm_surv)
summary(rm_15mA)
# Significant negative effect of conspBA 
# VIF and R2
vif_rm_15mA <- compute_vif(~ consp_15m +
                             logPrior + log0yr, rm_surv)
vif_rm_15mA # All <2
rm_15mA_r2 <- performance::r2(rm_15mA)
rm_15mA_r2

########## ########## ########## ########## ########## 
###### Partial Panel B - RM
rm_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = rm_surv)
summary(rm_15mB)
## Significant negative and positive effects of AM_15m and ECM_15m, respectively
## Consp_15m:AM_15m is significant and negative
# VIF and R2
vif_rm_15mB <- compute_vif(~ consp_15m + 
                             hetAM_15m + hetECM_15m + 
                             consp_15m:hetAM_15m + 
                             consp_15m:hetECM_15m + 
                             logPrior + log0yr, rm_surv)
vif_rm_15mB # Two are >10 (but interactions, so okay)
rm_15mB_r2 <- performance::r2(rm_15mB)
rm_15mB_r2

########## ########## ########## ########## ########## 
## Model comparison - RM
AIC(rm_15m_null, rm_15mA, rm_15mB) # B is significantly better

### Diagnostics - RM
sim_rm_B <- DHARMa::simulateResiduals(rm_15mB, n = 1000)
DHARMa::testDispersion(sim_rm_B) # No overdispersion
DHARMa::testZeroInflation(sim_rm_B) # Mild under-inflation
DHARMa::testUniformity(sim_rm_B) # Non-uniform 
DHARMa::testOutliers(sim_rm_B) # ~0.12% of observations are outliers


#########################################################
################# Species-Level Test - SM - Raw Basal Area

##### Subset sm data
sm_surv <- subset(allSpDat, Sp == "sm")

########## ########## ########## ########## ########## 
##### Null model - SM
sm_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                         (1|SiteID) + (1|survYear),
                       family = nbinom2, data = sm_surv) 
summary(sm_15m_null)
# R2
sm_15m_null_r2 <- performance::r2(sm_15m_null)
sm_15m_null_r2

########## ########## ########## ########## ########## 
##### Partial Panel A - CNDD - sm
sm_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = nbinom2, data = sm_surv)
summary(sm_15mA)
# Significant positive effect of conspBA 
# VIF and R2
vif_sm_15mA <- compute_vif(~ consp_15m +
                             logPrior + log0yr, sm_surv)
vif_sm_15mA # All <2
sm_15mA_r2 <- performance::r2(sm_15mA)
sm_15mA_r2

########## ########## ########## ########## ########## 
###### Partial Panel B - sm
sm_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                     hetAM_15m + hetECM_15m + 
                     consp_15m:hetAM_15m + # Panel B
                     consp_15m:hetECM_15m + # Panel B
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = nbinom2, data = sm_surv)
summary(sm_15mB)
## Significant positive effects of consp, AM_15m, and ECM_15m
## Consp_15m:AM_15m is significant and negative
# VIF and R2
vif_sm_15mB <- compute_vif(~ consp_15m + 
                             hetAM_15m + hetECM_15m + 
                             consp_15m:hetAM_15m + 
                             consp_15m:hetECM_15m + 
                             logPrior + log0yr, sm_surv)
vif_sm_15mB # All <10
sm_15mB_r2 <- performance::r2(sm_15mB)
sm_15mB_r2

########## ########## ########## ########## ########## 
## Model comparison - sm
AIC(sm_15m_null, sm_15mA, sm_15mB) # B is significantly better
anova(sm_15mA, sm_15mB)

### Diagnostics - sm
sim_sm_B <- DHARMa::simulateResiduals(sm_15mB, n = 1000)
DHARMa::testDispersion(sim_sm_B) # No overdispersion
DHARMa::testZeroInflation(sim_sm_B) # Mild under-inflation
DHARMa::testUniformity(sim_sm_B) # Non-uniform 
DHARMa::testOutliers(sim_sm_B) # ~0.47% of observations are outliers

#########################################################
################# Species-Level Test - bc - Raw Basal Area

##### Subset bc data
bc_surv <- subset(allSpDat, Sp == "bc")

########## ########## ########## ########## ########## 
##### Null model - BC
bc_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_null)
bc_15m_null_r2 <- performance::r2(bc_15m_null)
bc_15m_null_r2

########## ########## ########## ########## ########## 
##### Partial Panel A - CNDD - BC
bc_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15mA) 
## Significant negative effect of conspBA 
# VIF and R2
car::vif(bc_15mA) # All <2
bc_15mA_r2 <- performance::r2(bc_15mA)
bc_15mA_r2

########## ########## ########## ########## ########## 
##### Partial Panel B - BC
bc_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15mB)
## Significant negative effects of conspBA and AM_15m
## consp_15m:AM_15m is positive and significant
# VIF and R2
car::vif(bc_15mB) # One >10
bc_15mB_r2 <- performance::r2(bc_15mB)
bc_15mB_r2

########## ########## ########## ########## ########## 
## Model comparison - bc
AIC(bc_15m_null, bc_15mA, bc_15mB) 
anova(bc_15mA, bc_15mB) # B is not significantly better than A

### Diagnostics - bc
sim_bc_B <- DHARMa::simulateResiduals(bc_15mB, n = 1000)
DHARMa::testDispersion(sim_bc_B) # No overdispersion
DHARMa::testZeroInflation(sim_bc_B) # No zero-inflation
DHARMa::testUniformity(sim_bc_B) # Non-uniform 
DHARMa::testOutliers(sim_bc_B) # ~0.08% of observations are outliers

#########################################################
################# Species-Level Test - wa - Raw Basal Area

##### Subset wa data
wa_surv <- subset(allSpDat, Sp == "wa")

########## ########## ########## ########## ########## 
##### Null model - wa
wa_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                         (1|SiteID) + (1|survYear),
                       family = nbinom2, data = wa_surv) 
summary(wa_15m_null)
# R2
wa_15m_null_r2 <- performance::r2(wa_15m_null)
wa_15m_null_r2

########## ########## ########## ########## ########## 
##### Partial Panel A - CNDD - wa
wa_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = nbinom2, data = wa_surv)
summary(wa_15mA)
# Significant negative effect of conspBA 
# VIF and R2
vif_wa_15mA <- compute_vif(~ consp_15m +
                             logPrior + log0yr, wa_surv)
vif_wa_15mA # All <2
wa_15mA_r2 <- performance::r2(wa_15mA)
wa_15mA_r2

########## ########## ########## ########## ########## 
###### Partial Panel B - wa
wa_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                     hetAM_15m + hetECM_15m + 
                     consp_15m:hetAM_15m + # Panel B
                     consp_15m:hetECM_15m + # Panel B
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = nbinom2, data = wa_surv)
summary(wa_15mB)
## Significant negative effects of conspBA
## Consp_15m:AM_15m is marginally significant and positive
# VIF and R2
vif_wa_15mB <- compute_vif(~ consp_15m + 
                             hetAM_15m + hetECM_15m + 
                             consp_15m:hetAM_15m + 
                             consp_15m:hetECM_15m + 
                             logPrior + log0yr, wa_surv)
vif_wa_15mB # Two are >10
wa_15mB_r2 <- performance::r2(wa_15mB)
wa_15mB_r2

########## ########## ########## ########## ########## 
## Model comparison - wa
AIC(wa_15m_null, wa_15mA, wa_15mB) 
anova(wa_15mA, wa_15mB) # B is not significantly better

### Diagnostics - wa
sim_wa_B <- DHARMa::simulateResiduals(wa_15mB, n = 1000)
DHARMa::testDispersion(sim_wa_B) # No overdispersion
DHARMa::testZeroInflation(sim_wa_B) # No zero-inflation
DHARMa::testUniformity(sim_wa_B) # Non-uniform 
DHARMa::testOutliers(sim_wa_B) # ~0.44% of observations are outliers

#########################################################
################# Species-Level Test - ab - Raw Basal Area - Excludes RE of SiteID

##### Subset ab data
ab_surv <- subset(allSpDat, Sp == "ab")

########## ########## ########## ########## ########## 
##### Null model - AB
ab_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_null)
# R2
ab_15m_null_r2 <- performance::r2(ab_15m_null)
ab_15m_null_r2

########## ########## ########## ########## ########## 
##### Partial Panel A - CNDD - AB
ab_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_15mA) 
# conspBA is not significant
# VIF and R2
car::vif(ab_15mA) # All <2
ab_15m_A_r2 <- performance::r2(ab_15mA)
ab_15m_A_r2

########## ########## ########## ########## ########## 
##### Partial Panel B - AB
ab_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                  (1|survYear),
                 family = poisson(link = "log"), data = ab_surv,
                 control = glmerControl(optimizer = "bobyqa"))
summary(ab_15mB)
# VIF and R2
car::vif(ab_15mB) # A few >15 - interactions
ab_15m_B_r2 <- performance::r2(ab_15mB)
ab_15m_B_r2

########## ########## ########## ########## ########## 
## Model comparison - AB
AIC(ab_15m_null, ab_15mA, ab_15mB) 
anova(ab_15mA, ab_15mB) # B is not significantly better

### Diagnostics - wa
sim_ab_B <- DHARMa::simulateResiduals(ab_15mB, n = 1000)
DHARMa::testDispersion(sim_ab_B) # No overdispersion
DHARMa::testZeroInflation(sim_ab_B) # No zero-inflation
DHARMa::testUniformity(sim_ab_B) # Non-uniform 
DHARMa::testOutliers(sim_ab_B) # ~0.40% of observations are outliers

#########################################################
################# Species-Level Test - nro - Raw Basal Area - removed RE of site

##### Subset nro data
nro_surv <- subset(allSpDat, Sp == "nro")

## Null model - species-level
nro_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_null)
performance::r2(nro_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
nro_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_15mA) # conspBA is significant
car::vif(nro_15mA)
performance::r2(nro_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
nro_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   #consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_15mB)
## main effect of ECM_15m is significant
car::vif(nro_15mB) 
performance::r2(nro_15mB) 

## Model comparison - species-level
### Compare AICs 
AIC(nro_15m_null, nro_15mA, nro_15mB)
anova(nro_15mA, nro_15mB) # B is significantly better - but excludes consp_15m:hetAM_15m

#########################################################
################# Species-Level Test - bonro - Raw Basal Area - No RE of site

##### Subset bonro data
bonro_surv <- subset(allSpDat, Sp == "bonro")

## Null model - species-level
bonro_15m_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = bonro_surv) 
summary(bonro_15m_null)
performance::r2(bonro_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
bonro_15mA <- glmer(Year_1_counts ~ consp_15m +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15mA) # conspBA is not significant
car::vif(bonro_15mA)
performance::r2(bonro_15mA)

## Testing Partial Panel B conceptual figure - Interactive spillover
bonro_15mB <- glmer(Year_1_counts ~ consp_15m + 
                    hetAM_15m + hetECM_15m + 
                    #consp_15m:hetAM_15m + # Panel B
                    consp_15m:hetECM_15m + # Panel B
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15mB)
## none of main effects or interactions significant
car::vif(bonro_15mB) # One VIF >15
performance::r2(bonro_15mB) 

## Model comparison - species-level
### Compare AICs 
AIC(bonro_15m_null, bonro_15mA, bonro_15mB)
anova(bonro_15mA, bonro_15mB) # B is not significantly better

#########################################################
################# Species-Level Test - wo - Raw Basal Area - Excludes RE of site

##### Subset wo data
wo_surv <- subset(allSpDat, Sp == "wo")

## Null model - species-level
wo_15m_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = wo_surv) 
summary(wo_15m_null)
performance::r2(wo_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
wo_15mA <- glmer(Year_1_counts ~ consp_15m +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_15mA) # conspBA is not significant
car::vif(wo_15mA)
performance::r2(wo_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
wo_15mB <- glmer(Year_1_counts ~ consp_15m + 
                    hetAM_15m + hetECM_15m + 
                    consp_15m:hetAM_15m + # Panel B
                    consp_15m:hetECM_15m + # Panel B
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_15mB)
## none of main effects or their interactions are significant
car::vif(wo_15mB) # One VIF >15
performance::r2(wo_15mB) 

## Model comparison - species-level
### Compare AICs 
AIC(wo_15m_null, wo_15mA, wo_15mB)
anova(wo_15mA, wo_15mB) # B is NOT significantly better

#########################################################
################# Species-Level Test - iw - Raw Basal Area - Excludes RE of site

##### Subset iw data
iw_surv <- subset(allSpDat, Sp == "iw")

## Null model - species-level
iw_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_null)
performance::r2(iw_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
iw_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15mA) # conspBA is not significant
car::vif(iw_15mA)
performance::r2(iw_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
iw_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   #consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15mB)
## none of main effects or their interactions are significant
car::vif(iw_15mB) # One VIF >30
performance::r2(iw_15mB) 

## Model comparison - species-level
### Compare AICs 
AIC(iw_15m_null, iw_15mA, iw_15mB)
anova(iw_15mA, iw_15mB) # B is NOT significantly better

##################################################
########## Plotting effect sizes - Partial Panel B - species-level
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
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
    mutate(
      component = factor(component, levels = c(
        "Conspecific BA", "AM BA", "ECM BA",
        "Conspecific BA:AM BA","Conspecific BA:ECM BA",
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM  <- extract_fixef(rm_15mB, "rm")
fixef_SM  <- extract_fixef(sm_15mB, "sm")
fixef_BC  <- extract_fixef(bc_15mB, "bc")
fixef_WA  <- extract_fixef(wa_15mB, "wa")
fixef_AB  <- extract_fixef(ab_15mB, "ab")
fixef_NRO  <- extract_fixef(nro_15mB, "nro")
fixef_BONRO  <- extract_fixef(bonro_15mB, "bonro")
fixef_WO  <- extract_fixef(wo_15mB, "wo")
fixef_IW  <- extract_fixef(iw_15mB, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

# Define desired vertical order
desired_order <- c( "Conspecific BA", "AM BA", "ECM BA",
                    "Conspecific BA:AM BA","Conspecific BA:ECM BA",
                    "Log(Recruit Density)", 
                    "Log(Prior Density)")

## Species labels (italics, ordered alphabetically within AM/ECM)
species_labels <- c(
  ab    = expression(italic("Fagus grandifolia")),
  iw    = expression(italic("Ostrya virginiana")),
  wo    = expression(italic("Quercus alba")),
  nro   = expression(italic("Quercus rubra")),
  bonro = expression(italic("Quercus velutina")),
  rm    = expression(italic("Acer rubrum")),
  sm    = expression(italic("Acer saccharum")),
  wa    = expression(italic("Fraxinus americana")),
  bc    = expression(italic("Prunus serotina")))

## Colors matched to alphabetical order
species_cols <- c(
  ab    = "#C6E0EB",   # Fagus grandifolia
  iw    = "#4A90C2",   # Ostrya virginiana
  wo    = "#143D7A",   # Quercus alba
  nro   = "#1F5BA1",   # Quercus rubra
  bonro = "#7FB8D0",   # Quercus velutina
  rm    = "#FB8A74",   # Acer rubrum
  sm    = "#E63E2D",   # Acer saccharum
  wa    = "#991F1F",   # Fraxinus americana
  bc    = "#FDCBC0")    # Prunus serotina)

## Sizes for percent BA
species_sizes <- data.frame(species = c("ab", "bonro", "iw", "nro", 
                                        "wo", "bc", "rm", "sm", "wa"), 
                            percBA = c(2.43, 15.20, 0.11, 29.7, 10.30, 
                                       1.03, 9.40, 15.80, 1.82))

# Representative % BA values for legend
legend_percBA <- c(1, 10, 20)

###### Both AM and ECM
fixef_df_plot <- fixef_df %>%
  filter(estimate > -4 & std.error < 5) %>%  # exclude extreme outliers
  mutate(component = factor(component, levels = desired_order),
         # species legend order (alphabetical within ECM/AM)
         Model = factor(Model, levels = c("ab","iw","wo","nro","bonro",
                                          "rm","sm","wa","bc")),
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
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"))

#### Plot - All species
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error,
                     color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot$component)),
                     labels = levels(fixef_df_plot$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

##### Plot - AM
ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

##### Plot - ECM
ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))


##################################################
########## Plotting effect sizes - Partial Panel A - species-level
extract_fixef_A <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "consp_15m" ~ "Conspecific BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Conspecific BA", 
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM_A  <- extract_fixef_A(rm_15mA, "rm")
fixef_SM_A  <- extract_fixef_A(sm_15mA, "sm")
fixef_BC_A  <- extract_fixef_A(bc_15mA, "bc")
fixef_WA_A  <- extract_fixef_A(wa_15mA, "wa")
fixef_AB_A  <- extract_fixef_A(ab_15mA, "ab")
fixef_NRO_A  <- extract_fixef_A(nro_15mA, "nro")
fixef_BONRO_A  <- extract_fixef_A(bonro_15mA, "bonro")
fixef_WO_A  <- extract_fixef_A(wo_15mA, "wo")
fixef_IW_A  <- extract_fixef_A(iw_15mA, "iw")

## Combine
fixef_df_A <- bind_rows(fixef_RM_A, fixef_SM_A, fixef_BC_A,
                      fixef_WA_A, fixef_AB_A, fixef_NRO_A,
                      fixef_BONRO_A, fixef_WO_A, fixef_IW_A)

# Define desired vertical order
desired_order_A <- c("Conspecific BA", 
                    "Log(Recruit Density)", 
                    "Log(Prior Density)")

###### Both AM and ECM
fixef_df_plot_A <- fixef_df_A %>%
  filter(estimate > -4) %>%  # exclude extreme outliers
  mutate(component = factor(component, levels = desired_order_A),
         # species legend order (alphabetical within ECM/AM)
         Model = factor(Model, levels = c("ab","iw","wo","nro","bonro",
                                          "rm","sm","wa","bc")),
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
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"))

#### Plot - All species
ggplot(fixef_df_plot_A, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot_A$component)),
                     labels = levels(fixef_df_plot_A$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include AM
fixef_df_AM_A <- filter(fixef_df_plot_A, Myc_type == "AM")

##### Plot - AM
ggplot(fixef_df_AM_A, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include ECM
fixef_df_ECM_A <- filter(fixef_df_plot_A, Myc_type == "ECM")

##### Plot - ECM
ggplot(fixef_df_ECM_A, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType", 
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(
    name = "Species", values = species_cols, labels = species_labels) +
  scale_size(range = c(2, 6),
             breaks = log1p(legend_percBA),
             labels = paste0(legend_percBA, "%"),
             name = "% BA" ) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

#################################################
#################################################
################# Species-level Proportional BA Models
################# Species-Level Test - RM - Proportional Basal Area

## Null model - species-level (proportional BA)
rm_15m_prop_null <- glmmTMB(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = nbinom2, data = rm_surv) 
summary(rm_15m_prop_null)
performance::r2(rm_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
rm_propA <- glmmTMB(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = nbinom2, data = rm_surv)
summary(rm_propA) # significant negative effect of propConsp  
car::vif(rm_propA)
performance::r2(rm_propA) 

## Testing Partial Panel B (with AM) - species-level (proportional BA)
rm_propB <- glmmTMB(Year_1_counts ~ propConsp*propAM + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = nbinom2, data = rm_surv)
summary(rm_propB) # significant negative effect of propAM  
car::vif(rm_propB)
performance::r2(rm_propB) 

## Model comparison - RM (proportional BA)
### Compare AICs
AIC(rm_15m_prop_null, rm_propA, rm_propB) # Model B significantly better

################# 
################# Species-Level Test - sm - Proportional Basal Area

## Null model - species-level (proportional BA)
sm_15m_prop_null <- glmmTMB(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = nbinom2, data = sm_surv) 
summary(sm_15m_prop_null)
performance::r2(sm_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
sm_propA <- glmmTMB(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = nbinom2, data = sm_surv)
summary(sm_propA) # Significant positive effect of propConsp
car::vif(sm_propA)
performance::r2(sm_propA) 

## Testing Partial Panel B - species-level (proportional BA)
sm_propB <- glmmTMB(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = nbinom2, data = sm_surv)
summary(sm_propB) # Significant negative effect of propAM 
# Significant positive propConsp:propAM interaction
car::vif(sm_propB)
performance::r2(sm_propB) 

## Model comparison - SM (proportional BA)
### Compare AICs
AIC(sm_15m_prop_null, sm_propA, sm_propB)
anova(sm_propA, sm_propB) # More complex model is better!

################# 
################# Species-Level Test - bc - Proportional Basal Area

## Null model - species-level (proportional BA)
bc_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_prop_null)
performance::r2(bc_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
bc_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_propA) # Significant negative effect of propConsp
car::vif(bc_propA)
performance::r2(bc_propA) 

## Testing Partial Panel B - species-level (proportional BA)
bc_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_propB) # Only negative effect of propConsp significant
car::vif(bc_propB)
performance::r2(bc_propB) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(bc_15m_prop_null, bc_propA, bc_propB)
anova(bc_propA, bc_propB) # More complex model is NOT better

#################
################# Species-Level Test - wa - Proportional Basal Area

## Null model - species-level (proportional BA)
wa_15m_prop_null <- glmmTMB(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = nbinom2, data = wa_surv) 
summary(wa_15m_prop_null)
performance::r2(wa_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
wa_propA <- glmmTMB(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = nbinom2, data = wa_surv)
summary(wa_propA) # Significant negative effect of propConsp
car::vif(wa_propA)
performance::r2(wa_propA) 

## Testing Partial Panel B - species-level (proportional BA)
wa_propB <- glmmTMB(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = nbinom2, data = wa_surv)
summary(wa_propB) # Only negative main effect of propConsp significant
car::vif(wa_propB)
performance::r2(wa_propB) 

## Model comparison - WA (proportional BA)
### Compare AICs
AIC(wa_15m_prop_null, wa_propA, wa_propB)
anova(wa_propA, wa_propB) # More complex model is NOT better

################# 
################# Species-Level Test - ab - Proportional Basal Area

## Null model - species-level (proportional BA)
ab_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_prop_null)
performance::r2(ab_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
ab_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                     (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_propA) # Non-significant effect of propConsp
car::vif(ab_propA)
performance::r2(ab_propA) 

## Testing Partial Panel B - species-level (proportional BA)
ab_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                     (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_propB) # No terms significant
car::vif(ab_propB)
performance::r2(ab_propB) 

## Model comparison - AB (proportional BA)
### Compare AICs
AIC(ab_15m_prop_null, ab_propA, ab_propB)
anova(ab_propA, ab_propB) # More complex model is NOT better

################# 
################# Species-Level Test - nro - Proportional Basal Area

## Null model - species-level (proportional BA)
nro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_prop_null)
performance::r2(nro_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
nro_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_propA) # Non-significant effect of propConsp
car::vif(nro_propA)
performance::r2(nro_propA) 

## Testing Partial Panel B - species-level (proportional BA)
nro_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_propB) # Significant negative effects of propConsp and propAM 
car::vif(nro_propB)
performance::r2(nro_propB) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(nro_15m_prop_null, nro_propA, nro_propB)
anova(nro_propA, nro_propB) # More complex model IS better

################# 
################# Species-Level Test - bonro - Proportional Basal Area

## Null model - species-level (proportional BA)
bonro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear),
                           family = poisson(link = "log"), data = bonro_surv) 
summary(bonro_15m_prop_null)
performance::r2(bonro_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
bonro_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                      (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_propA) # Non-significant effect of propConsp
car::vif(bonro_propA)
performance::r2(bonro_propA) 

## Testing Partial Panel B - species-level (proportional BA)
bonro_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_propB) # No terms significant
car::vif(bonro_propB)
performance::r2(bonro_propB) 

## Model comparison - BONRO (proportional BA)
### Compare AICs
AIC(bonro_15m_prop_null, bonro_propA, bonro_propB)
anova(bonro_propA, bonro_propB) # More complex model is NOT better

################# 
################# Species-Level Test - wo - Proportional Basal Area

## Null model - species-level (proportional BA)
wo_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                             (1|survYear),
                           family = poisson(link = "log"), data = wo_surv)
summary(wo_15m_prop_null)
performance::r2(wo_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
wo_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_propA) # Non-significant effect of propConsp
car::vif(wo_propA)
performance::r2(wo_propA) 

## Testing Partial Panel B - species-level (proportional BA)
wo_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_propB) # No terms significant
car::vif(wo_propB)
performance::r2(wo_propB) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(wo_15m_prop_null, wo_propA, wo_propB)
anova(wo_propA, wo_propB) 

################# 
################# Species-Level Test - iw - Proportional Basal Area

## Null model - species-level (proportional BA)
iw_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_prop_null)
performance::r2(iw_15m_prop_null) 

## Testing Panel A - CNDD - species-level (proportional BA)
iw_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_propA) # Non-significant effect of propConsp
car::vif(iw_propA)
performance::r2(iw_propA) 

## Testing Partial Panel B - species-level (proportional BA)
iw_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_propB) # No terms significant
car::vif(iw_propB)
performance::r2(iw_propB) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(iw_15m_prop_null, iw_propA, iw_propB)
anova(iw_propB, iw_propB2) # More complex model is NOT better


##################################################
########## Plotting effect sizes - Partial Panel B - species-level
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propAM" ~ "Prop. AM BA",
        term == "propConsp:propAM" ~ "Prop. Consp. BA:Prop. AM BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Prop. Conspecific BA", "Prop. AM BA",
        "Prop. Consp. BA:Prop. AM BA",
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM  <- extract_fixef(rm_propB, "rm")
fixef_SM  <- extract_fixef(sm_propB, "sm")
fixef_BC  <- extract_fixef(bc_propB, "bc")
fixef_WA  <- extract_fixef(wa_propB, "wa")
fixef_AB  <- extract_fixef(ab_propB, "ab")
fixef_NRO <- extract_fixef(nro_propB, "nro")
fixef_BONRO <- extract_fixef(bonro_propB, "bonro")
fixef_WO  <- extract_fixef(wo_propB, "wo")
fixef_IW  <- extract_fixef(iw_propB, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

# Define desired vertical order
desired_order <- c("Prop. Conspecific BA", 
                   "Prop. AM BA", 
                   "Prop. Consp. BA:Prop. AM BA",
                   "Log(Recruit Density)", 
                   "Log(Prior Density)")

## Species labels (italics, alphabetized within AM & ECM)
species_labels <- c(
  ab    = expression(italic("Fagus grandifolia")),
  iw    = expression(italic("Ostrya virginiana")),
  wo    = expression(italic("Quercus alba")),
  nro   = expression(italic("Quercus rubra")),
  bonro = expression(italic("Quercus velutina")),
  rm    = expression(italic("Acer rubrum")),
  sm    = expression(italic("Acer saccharum")),
  wa    = expression(italic("Fraxinus americana")),
  bc    = expression(italic("Prunus serotina")))

## Species colors (alphabetical within AM & ECM)
species_cols <- c(
  ab    = "#C6E0EB",   # Fagus grandifolia
  iw    = "#4A90C2",   # Ostrya virginiana
  wo    = "#143D7A",   # Quercus alba
  nro   = "#1F5BA1",   # Quercus rubra
  bonro = "#7FB8D0",   # Quercus velutina
  rm    = "#FB8A74",   # Acer rubrum
  sm    = "#E63E2D",   # Acer saccharum
  wa    = "#991F1F",   # Fraxinus americana
  bc    = "#FDCBC0"    # Prunus serotina
)

## Mean % basal area
species_sizes <- data.frame(species = c("ab", "bonro", "iw", "nro", 
                                        "wo", "bc", "rm", "sm", "wa"), 
                            percBA = c(2.43, 15.20, 0.11, 29.7, 10.30, 
                                       1.03, 9.40, 15.80, 1.82))

# Representative % BA values for legend
legend_percBA <- c(1, 10, 20)

###### Both AM and ECM
fixef_df_plot <- fixef_df %>%
  filter(estimate > -10 & estimate < 10) %>%
  mutate(component = factor(component, levels = desired_order),
         Model = factor(Model, levels = c("ab","iw","wo","nro","bonro",
                                          "rm","sm","wa","bc")),
         Significant = factor(Significant, levels = c("yes", "no")),
         Myc_type = case_when(Model %in% c("rm","sm","bc","wa") ~ "AM",
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

#### Plot all species
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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

##### AM only
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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

##### ECM only
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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


##################################################
########## Plotting effect sizes - Partial Panel A - species-level
extract_fixef_Aprop <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Prop. Conspecific BA", 
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM_propA  <- extract_fixef_Aprop(rm_propA, "rm")
fixef_SM_propA  <- extract_fixef_Aprop(sm_propA, "sm")
fixef_BC_propA  <- extract_fixef_Aprop(bc_propA, "bc")
fixef_WA_propA  <- extract_fixef_Aprop(wa_propA, "wa")
fixef_AB_propA  <- extract_fixef_Aprop(ab_propA, "ab")
fixef_NRO_propA <- extract_fixef_Aprop(nro_propA, "nro")
fixef_BONRO_propA <- extract_fixef_Aprop(bonro_propA, "bonro")
fixef_WO_propA  <- extract_fixef_Aprop(wo_propA, "wo")
fixef_IW_propA  <- extract_fixef_Aprop(iw_propA, "iw")

## Combine
fixef_df_propA <- bind_rows(fixef_RM_propA, fixef_SM_propA, fixef_BC_propA,
                      fixef_WA_propA, fixef_AB_propA, fixef_NRO_propA,
                      fixef_BONRO_propA, fixef_WO_propA, fixef_IW_propA)

# Define desired vertical order
desired_order_propA <- c("Prop. Conspecific BA", 
                   "Log(Recruit Density)", 
                   "Log(Prior Density)")

###### Both AM and ECM
fixef_df_plot_propA <- fixef_df_propA %>%
  filter(estimate > -10) %>%
  mutate(component = factor(component, levels = desired_order_propA),
         Model = factor(Model, levels = c("ab","iw","wo","nro","bonro",
                                          "rm","sm","wa","bc")),
         Significant = factor(Significant, levels = c("yes", "no")),
         Myc_type = case_when(Model %in% c("rm","sm","bc","wa") ~ "AM",
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

#### Plot all species
ggplot(fixef_df_plot_propA, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df_plot_propA$component)),
                     labels = levels(fixef_df_plot_propA$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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

##### AM only
fixef_df_AM_propA <- filter(fixef_df_plot_propA, Myc_type == "AM")

ggplot(fixef_df_AM_propA, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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

##### ECM only
fixef_df_ECM_propA <- filter(fixef_df_plot_propA, Myc_type == "ECM")

ggplot(fixef_df_ECM_propA, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbar(aes(xmin = estimate - 1.96 * std.error,
                    xmax = estimate + 1.96 * std.error,
                    color = Model), width = 0.05) +
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)), stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse") +
  scale_shape_manual(name = "Mycorrhizal\nType",
                     values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(name = "Species", values = species_cols,
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


