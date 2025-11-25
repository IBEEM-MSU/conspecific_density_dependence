
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
library(DHARMa)
library(broom.mixed)
library(tibble)
library(forcats)
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
# Significant negative effect of consp
# VIF
vif_A <- compute_vif(~ consp_15m * Seedling_myc +
                        logPrior + log0yr, allSpDat)
vif_A # All <=10
# R2
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
# VIF
vif_B <- compute_vif(~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    logPrior + log0yr,  allSpDat)
vif_B # All <7
# R2
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
# VIF
vif_C <- compute_vif( ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    logPrior + log0yr, allSpDat)
vif_C # One >16
# R2
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
# VIF
vif_D <- compute_vif( ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    logPrior + log0yr, allSpDat)
vif_D
# R2
r2_D <- performance::r2(all_sp_15mD)
r2_D

# Diagnostics
simD <- DHARMa::simulateResiduals(all_sp_15mD, n = 1000)
plot(simD)              
DHARMa::testDispersion(simD)
DHARMa::testZeroInflation(simD)
DHARMa::testUniformity(simD)
DHARMa::testOutliers(simD)

######### ######### ######### ######### ######### 
######### PANEL B+C+D — SUPERINTERACTIVE MODEL
all_sp_15mBCD <- glmmTMB(
  Year_1_counts ~ consp_15m * Seedling_myc +
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
## Negative interactions between consp:hetAM and consp:hetECM 
## Non-significant Seedling_mycECM:hetAM_15m and Seedling_mycECM:hetECM_15m 
## consp:Seedling_mycECM:hetECM_15m Significant positive; ECM seedlings in ECM-heterogeneous microsites experience less negative or even reversed competitive effects.
# VIF
vif_BCD <- compute_vif(
  ~ consp_15m * Seedling_myc +
    hetAM_15m + hetECM_15m +
    consp_15m:hetAM_15m + consp_15m:hetECM_15m +
    Seedling_myc:hetAM_15m + Seedling_myc:hetECM_15m +
    consp_15m:Seedling_myc:hetAM_15m +
    consp_15m:Seedling_myc:hetECM_15m +
    logPrior + log0yr,allSpDat)
vif_BCD # None of these high VIFs on interactions threaten the validity of the model; they simply reflect interaction-induced multicollinearity.
# R2
r2_BCD <- performance::r2(all_sp_15mBCD)
r2_BCD # 47% of variance explained by fixed effects, 13% by REs
# Diagnostics
simBCD <- DHARMa::simulateResiduals(all_sp_15mBCD, n = 1000)
plot(simBCD)              
DHARMa::testDispersion(simBCD) # No overdispersion
DHARMa::testZeroInflation(simBCD) # No zero-inflation
DHARMa::testUniformity(simBCD) # Non-uniform, but very normal with large datasets - not an indicator of poor fit
DHARMa::testOutliers(simBCD) # ~0.28% of obs are outliers - not concerning


### Compare AICs
### All sites combined
AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, 
    all_sp_15mC, all_sp_15mD, all_sp_15mBCD) # D is lowest
anova(all_sp_15mB2, all_sp_15mD2)

#### Post-hoc tests

# 1. Define a grid of ecologically meaningful predictor values
emm_BCD <- emmeans(
  all_sp_15mBCD,
  ~ Seedling_myc | consp_15m * hetAM_15m * hetECM_15m,
  type = "response",
  at = list(
    consp_15m   = quantile(allSpDat$consp_15m, c(.1, .5, .9)),
    hetAM_15m   = quantile(allSpDat$hetAM_15m, c(.1, .5, .9)),
    hetECM_15m  = quantile(allSpDat$hetECM_15m, c(.1, .5, .9))))
emm_BCD
pairs(emm_BCD, by = c("consp_15m", "hetAM_15m", "hetECM_15m"))
# Across realistic ranges of conspecific density and heterospecific neighbor composition, AM and ECM seedlings have indistinguishable recruitment rates.

## Model comparison
### Compare AICs - lowest is all_sp_15mB but all_sp_15mD close 2nd
AIC(all_sp_15m_null, all_sp_15mA, all_sp_15mB, all_sp_15mC,
    all_sp_15mD, all_sp_15mBCD)
# Global model is significantly better!!
anova(all_sp_15mB, all_sp_15mBCD) 


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
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error), height = 0.1, color = "black")+
  geom_point(aes(shape = Significant, fill = Significant),
             size = 3.2, color = "black", stroke = 0.9) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± 95% CI (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"))


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
performance::r2(all_sp_prop_null)

################# ################# ################# 
## Testing Panel A conceptual figure (proportional BA)
all_sp_propA <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = nbinom2, data = allSpDat) 
summary(all_sp_propA)
## main effect of propConsp is marginally significant
## main effect of Seedling_myc is not signficant 
## propConsp:Seedling_myc interaction is negative and significant
# VIF
vif_propA <- compute_vif( ~ propConsp*Seedling_myc +
                        logPrior + log0yr, allSpDat)
vif_propA
# R2
r2_propA <- performance::r2(all_sp_propA)
r2_propA

## Testing Panel B conceptual figure (proportional BA)
all_sp_propB <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propAM + 
                       propConsp:propAM + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propB)
## main effects of propConsp and propAM are significant
## propConsp:Seedling_myc is marginally significant and positive (small)
## propConsp:propAM is significant and positive
# VIF
vif_propB <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + propConsp:propAM + 
                            logPrior + log0yr, allSpDat)
vif_propB
# R2
r2_propB <- performance::r2(all_sp_propB)
r2_propB

## Testing Panel C conceptual figure (proportional BA)
all_sp_propC <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propAM + 
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propC)
## main effects of propConsp and propAM are significant, but not of Seedling_myc
## propConsp:Seedling_myc interaction is significant and negative
## Seedling_myc:AM_15m interaction is not significant 
# VIF
vif_propC <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + Seedling_myc:propAM + 
                            logPrior + log0yr +, allSpDat)
vif_propC
# R2
r2_propC <- performance::r2(all_sp_propC)
r2_propC

## Testing Panel D (B+C) conceptual figure (proportional BA)
all_sp_propD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                        propAM + 
                        propConsp:propAM + # Panel B
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat)
summary(all_sp_propD)
## main effects of propConsp and propAM are significant, but not of Seedling_myc
## propConsp:Seedling_myc interaction is significant and positive
## propConsp:propAM interaction is significant and positive
## Seedling_myc:propAM interaction is not significant 
# VIF
vif_propD <- compute_vif( ~ propConsp*Seedling_myc + 
                            propAM + propConsp:propAM + 
                            Seedling_myc:propAM + 
                            logPrior + log0yr , allSpDat)
vif_propD
# R2
r2_propD <- performance::r2(all_sp_propD)
r2_propD

## Testing Panels B, C, and D - Global model (proportional BA)
## Super interactive model (testing everything)
all_sp_propBCD <- glmmTMB(Year_1_counts ~ propConsp*Seedling_myc + 
                       propAM + 
                       propConsp:propAM + # Panel B
                       Seedling_myc:propAM + # Panel C
                       propConsp:Seedling_myc:propAM + # Panel D
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat)
summary(all_sp_propBCD)
## main effects of propConsp and propAM are significant, but not of Seedling_myc
## propConsp:Seedling_myc interaction is significant and positive
## propConsp:propAM interaction is significant and positive
## Seedling_myc:propAM interaction is now significant and positive
## propConsp:Seedling_myc:propAM interaction is significant and negative
# VIF
vif_propBCD <- compute_vif( ~ propConsp*Seedling_myc + 
                              propAM + propConsp:propAM + 
                              Seedling_myc:propAM + 
                              propConsp:Seedling_myc:propAM + 
                              logPrior + log0yr, allSpDat)
vif_propBCD
# R2
r2_propBCD <- performance::r2(all_sp_propBCD)
r2_propBCD

## Model comparison
### Compare AICs - lowest is all_sp_propBCD 
AIC(all_sp_prop_null, all_sp_propA, all_sp_propB, all_sp_propC,
    all_sp_propD, all_sp_propBCD)
anova(all_sp_propB, all_sp_propBCD) # Global model is significantly better!


########## Plotting effect sizes 
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
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error), 
                 height = 0.1, color = "black")+
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

##### Null model - AM
AM_sp_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15m_null)
# R2
r2_null <- performance::r2(AM_sp_15m_null)
r2_null

##### Null model - ECM
ECM_sp_15m_null <- glmmTMB(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear) + (1|Sp),
                          family = nbinom2, data = allSpDat_ECM) 
summary(ECM_sp_15m_null)
# R2
r2_null <- performance::r2(ECM_sp_15m_null)
r2_null

##### Testing Panel A conceptual figure - CDD - AM
AM_sp_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = nbinom2, data = allSpDat_AM) 
summary(AM_sp_15mA)
## Main effect of ConspBA is negative and significant
# VIF
vif_AMa <- compute_vif(~ consp_15m +
                        logPrior + log0yr, allSpDat_AM)
vif_AMa # All <2
# R2
performance::r2(AM_sp_15mA)

##### Testing Panel A conceptual figure - CDD - ECM
ECM_sp_15mA <- glmmTMB(Year_1_counts ~ consp_15m +
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat_ECM) 
summary(ECM_sp_15mA)
## Main effect of ConspBA is negative and significant
# VIF
vif_ECMa <- compute_vif(~ consp_15m +
                         logPrior + log0yr, allSpDat_ECM)
vif_ECMa # All <2
# R2
performance::r2(ECM_sp_15mA)

######## Testing Panel B conceptual figure - AM
AM_sp_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                       hetAM_15m + hetECM_15m + 
                       consp_15m:hetAM_15m + # Panel B
                       consp_15m:hetECM_15m + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = nbinom2, data = allSpDat_AM)
summary(AM_sp_15mB)
## Significant positive main effects of conspBA, AM BA, and ECM BA
## consp_15m:AM_15m is positive and significant
## consp_15m:ECM_15m is negative and significant
# VIF
vif_AMb <- compute_vif(~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr, allSpDat_AM)
vif_AMb # All <2
# R2
performance::r2(AM_sp_15mB)

## Testing Panel B conceptual figure - ECM - HAS TO BE POISSON
ECM_sp_15mB <- glmmTMB(Year_1_counts ~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr +
                         (1|SiteID) + 
                         (1|survYear) + (1|Sp),
                  family = poisson, data = allSpDat_ECM,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(ECM_sp_15mB)


# Compare AIC
AIC(m_pois, m_nb1, m_nb2)
anova(m_pois, m_nb1)

# basic checks
table(allSpDat_ECM$Sp)             # observations per species
table(allSpDat_ECM$SiteID)         # obs per site (should be 4 sites * plots)
length(unique(allSpDat_ECM$Sp))
length(unique(allSpDat_ECM$SiteID))

# check random-effect variances from a simpler model
mod_simple <- glmmTMB(Year_1_counts ~ 1 + (1|SiteID) + (1|survYear) + (1|Sp),
                      family = nbinom2, data = allSpDat_ECM)
summary(mod_simple)
VarCorr(mod_simple)

## Main effects of conspBA, AM BA, and ECM BA all not significant
## consp_15m:AM_15m is negative and significant
## consp_15m:ECM_15m is positive and significant
# VIF
vif_ECMb <- compute_vif(~ consp_15m + 
                         hetAM_15m + hetECM_15m + 
                         consp_15m:hetAM_15m + # Panel B
                         consp_15m:hetECM_15m + # Panel B
                         logPrior + log0yr, allSpDat_ECM)
vif_ECMb # All <2
# R2
performance::r2(ECM_sp_15mB)

## Model comparison - AM
### Compare AICs - lowest is most complex - AM_sp_15mB
AIC(AM_sp_15m_null, AM_sp_15mA, AM_sp_15mA2, AM_sp_15mB, AM_sp_15mB2)
## Model comparison - ECM
AIC(ECM_sp_15m_null, ECM_sp_15mA, ECM_sp_15mA2, ECM_sp_15mB, ECM_sp_15mB2)
anova(ECM_sp_15mA,ECM_sp_15mB) # B is better

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

# --- prepare data
fixef_df <- fixef_df %>%
  mutate(# enforce consistent factor order
    component = factor(component, levels = rev(desired_order)),
    # tidy up Model and Significant
    Model = factor(Model, levels = c("AM", "ECM")),
    Significant = factor(Significant, levels = c("yes", "no")),
    # shift large estimates for plotting
    plot_est = ifelse(estimate > 1, estimate - 1, estimate),
    plot_xmin = ifelse(estimate > 1, estimate - 1, estimate) - std.error,
    plot_xmax = ifelse(estimate > 1, estimate - 1, estimate) + std.error)

######## Plot
offset <- 0.05
ggplot(fixef_df, aes(x = plot_est, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = plot_xmin, xmax = plot_xmax,
                     y = as.numeric(component) + 
                       ifelse(Model == "AM", offset, -offset)),
                 height = 0.05, color = "black") +
  geom_point(aes(shape = Model, fill = Significant,
        y = as.numeric(component) + 
           ifelse(Model == "AM", offset, -offset)),
             size = 3, color = "black", stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component)) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
    labels = c("-0.5", "0", "0.5", "1.5", "2.0", "2.5")) +
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes = 
          list(shape = 21, color = "black", size = 3.5)))

#########################################################
#########################################################
################# Mycorrhizal Type Models - Proportional Basal Area

## Null model - AM - Proportional BA
AM_sp_15m_nullprop <- glmer(Year_1_counts ~ 1 +
                          (1|SiteID) + (1|survYear) + (1|Sp),
                        family = poisson(link = "log"), data = allSpDat_AM) 
summary(AM_sp_15m_nullprop)
r.squaredGLMM(AM_sp_15m_nullprop) 

## Null model - ECM - Proportional BA
ECM_sp_15m_nullprop <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15m_nullprop)
r.squaredGLMM(ECM_sp_15m_nullprop) 

## Testing Panel A conceptual figure - AM - Proportional BA
AM_sp_15mpropA <- glmer(Year_1_counts ~ propConsp +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = poisson(link = "log"), data = allSpDat_AM) 
summary(AM_sp_15mpropA)
## Main effect of ConspBA is negative and significant
car::vif(AM_sp_15mpropA)
r.squaredGLMM(AM_sp_15mpropA) 

## Testing Panel A conceptual figure - ECM - Proportional BA
ECM_sp_15mpropA <- glmer(Year_1_counts ~ propConsp +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15mpropA)
## Main effect of ConspBA is not significant
car::vif(ECM_sp_15mpropA)
r.squaredGLMM(ECM_sp_15mpropA) 

## Testing Panel B conceptual figure - AM - Proportional BA
AM_sp_15mpropB <- glmer(Year_1_counts ~ propConsp + 
                      propAM + propConsp:propAM + 
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_AM)
summary(AM_sp_15mpropB)
## Significant negative main effects of propConsp and propAM
## propConsp:propAM is positive and significant
car::vif(AM_sp_15mpropB)
r.squaredGLMM(AM_sp_15mpropB) 

## Testing Panel B conceptual figure - ECM - Proportional BA
ECM_sp_15mpropB <- glmer(Year_1_counts ~ propConsp + 
                       propAM + propConsp:propAM + 
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15mpropB)
## Non-significant main effects of propConsp and propAM
## propConsp:propAM is negative and significant
car::vif(ECM_sp_15mpropB)
r.squaredGLMM(ECM_sp_15mpropB) 


## Model comparison - AM
### Compare AICs - lowest is most complex - AM_sp_15mpropB
AIC(AM_sp_15m_nullprop, AM_sp_15mpropA, AM_sp_15mpropB)
## Model comparison - ECM
### Compare AICs - lowest is most complex - ECM_sp_15mpropB
AIC(ECM_sp_15m_nullprop, ECM_sp_15mpropA, ECM_sp_15mpropB)
anova(ECM_sp_15mpropA, ECM_sp_15mpropB)

########## Plotting effect sizes - for Panel B - proportional BA
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
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
fixef_AM  <- extract_fixef(AM_sp_15m1prop, "AM")
fixef_ECM <- extract_fixef(ECM_sp_15m1prop, "ECM")

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
  # Error bars with same offset as points
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error,
          y = as.numeric(component) + ifelse(Model == "AM", offset, -offset)),
                 height = 0.05, color = "black") +
  # Points with the same vertical offset
  geom_point(aes(shape = Model, fill = Significant,
        y = as.numeric(component) + ifelse(Model == "AM", offset, -offset)),
             size = 3, color = "black", stroke = 1) +
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component)) +
  scale_shape_manual(values = c("AM" = 24, "ECM" = 22)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(override.aes = 
                    list(shape = 21, color = "black", size = 3.5)))

#########################################################
#########################################################
#########################################################
################# Species-Level Model - RM - Raw Basal Area

##### Subset RM data
rm_surv <- subset(allSpDat, Sp == "rm")

## Null model - species-level
rm_15m_null <- glmer(Year_1_counts ~ 1 +
                     (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv) 
summary(rm_15m_null)
r.squaredGLMM(rm_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
rm_15mA <- glmer(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_15mA) # conspBA highly significant
car::vif(rm_15mA)
r.squaredGLMM(rm_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
rm_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_15mB)
## main effects of AM_15m and ECM_15m are significant
## main effect of conspBA is not significant
## Consp_15m:AM_15m is significant and negative
## Consp_15m:ECM_15m is marginally significant and negative
car::vif(rm_15mB) # One VIF >12
r.squaredGLMM(rm_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# rm_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                    hetAM_15m + hetECM_15m + 
#                    logPrior + log0yr +
#                    (1|SiteID) + (1|survYear),
#                  family = poisson(link = "log"), data = rm_surv)
# summary(rm_15mB2)
# ## Main effects of conspBA and AM_15m are significant
# ## Main effects of ECM_15m is not significant
# car::vif(rm_15mB2) 
# r.squaredGLMM(rm_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(rm_15m_null, rm_15mA, rm_15mB)
#AIC(rm_15m_null, rm_15mA, rm_15mB, rm_15mB2)
# More complex model is better
anova(rm_15mA, rm_15mB)


#########################################################
################# Species-Level Test - SM - Raw Basal Area

##### Subset sm data
sm_surv <- subset(allSpDat, Sp == "sm")

## Null model - species-level
sm_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = sm_surv) 
summary(sm_15m_null)
r.squaredGLMM(sm_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
sm_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_15mA) # conspBA highly significant
car::vif(sm_15mA)
r.squaredGLMM(sm_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
sm_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_15mB)
## main effects of AM_15m is not significant
## main effects of conspBA and ECM_15m are significant
## Consp_15m:AM_15m is significant and negative
## Consp_15m:ECM_15m is not significant
car::vif(sm_15mB) 
r.squaredGLMM(sm_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# sm_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                     (1|SiteID) + (1|survYear),
#                   family = poisson(link = "log"), data = sm_surv)
# summary(sm_15mB2)
# ## main effects of AM_15m, conspBA, and ECM_15m are all significant
# car::vif(sm_15mB2) 
# r.squaredGLMM(sm_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(sm_15m_null, sm_15mA, sm_15mB)
#AIC(sm_15m_null, sm_15mA, sm_15mB, sm_15mB2)
# More complex model is better
anova(sm_15mA, sm_15mB)

#########################################################
################# Species-Level Test - bc - Raw Basal Area

##### Subset bc data
bc_surv <- subset(allSpDat, Sp == "bc")

## Null model - species-level
bc_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_null)
r.squaredGLMM(bc_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
bc_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15mA) # conspBA highly significant
car::vif(bc_15mA)
r.squaredGLMM(bc_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
bc_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15mB)
## main effects of conspBA and AM_15m are significant
## main effect ECM_15m is not significant
## consp_15m:AM_15m is postitive and significant
## consp_15m:ECM_15m is not significant
car::vif(bc_15mB) # One VIF >14
r.squaredGLMM(bc_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# bc_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                     (1|SiteID) + (1|survYear),
#                   family = poisson(link = "log"), data = bc_surv)
# summary(bc_15mB2)
# ## main effect of conspBA is significant
# ## main effects of AM_15m and ECM_15m are not significant
# car::vif(bc_15mB2) 
# r.squaredGLMM(bc_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(bc_15m_null, bc_15mA, bc_15mB)
# AIC(bc_15m_null, bc_15mA, bc_15mB, bc_15mB2)
# More complex model is NOT better
anova(bc_15mA, bc_15mB)

#########################################################
################# Species-Level Test - wa - Raw Basal Area

##### Subset wa data
wa_surv <- subset(allSpDat, Sp == "wa")

## Null model - species-level
wa_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = wa_surv) 
summary(wa_15m_null)
r.squaredGLMM(wa_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
wa_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_15mA) # conspBA highly significant
car::vif(wa_15mA)
r.squaredGLMM(wa_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
wa_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_15mB)
## Main effect of conspBA is significant 
## None of other main effects or their interactions are significant 
car::vif(wa_15mB) # One VIF >30
r.squaredGLMM(wa_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# wa_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                     (1|SiteID) + (1|survYear),
#                   family = poisson(link = "log"), data = wa_surv)
# summary(wa_15mB2)
# ## main effect of conspBA is significant
# ## main effects of AM_15m and ECM_15m are not significant
# car::vif(wa_15mB2) 
# r.squaredGLMM(wa_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(wa_15m_null, wa_15mA, wa_15mB)
# AIC(wa_15m_null, wa_15mA, wa_15mB, wa_15mB2)
# More complex model is NOT better
anova(wa_15mA, wa_15mB)

#########################################################
################# Species-Level Test - ab - Raw Basal Area - Excludes RE of SiteID

##### Subset ab data
ab_surv <- subset(allSpDat, Sp == "ab")

## Null model - species-level
ab_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_null)
r.squaredGLMM(ab_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
ab_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_15mA) # conspBA is not significant
car::vif(ab_15mA)
r.squaredGLMM(ab_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
ab_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   #consp_15m:hetAM_15m + # Panel B
                   consp_15m:hetECM_15m + # Panel B
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_15mB)
## none of the main effects or their interactions are significant 
car::vif(ab_15mB) 
r.squaredGLMM(ab_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# ab_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                      (1|survYear),
#                   family = poisson(link = "log"), data = ab_surv)
# summary(ab_15mB2)
# ## main effects of conspBA, AM_15m, and ECM_15m are not significant
# car::vif(ab_15mB2) 
# r.squaredGLMM(ab_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(ab_15m_null, ab_15mA, ab_15mB)
# AIC(ab_15m_null, ab_15mA, ab_15mB, ab_15mB2)
# More complex model is NOT better
anova(ab_15mA, ab_15mB)

#########################################################
################# Species-Level Test - nro - Raw Basal Area - removed RE of site

##### Subset nro data
nro_surv <- subset(allSpDat, Sp == "nro")

## Null model - species-level
nro_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_null)
r.squaredGLMM(nro_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
nro_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_15mA) # conspBA is significant
car::vif(nro_15mA)
r.squaredGLMM(nro_15mA) 

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
r.squaredGLMM(nro_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# nro_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                     (1|survYear),
#                   family = poisson(link = "log"), data = nro_surv)
# summary(nro_15mB2)
# ## main effects of conspBA, AM_15m, and ECM_15m are all significant
# car::vif(nro_15mB2) 
# r.squaredGLMM(nro_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(nro_15m_null, nro_15mA, nro_15mB)
# AIC(nro_15m_null, nro_15mA, nro_15mB, nro_15mB2)
# More complex model is better
anova(nro_15mA, nro_15mB)

#########################################################
################# Species-Level Test - bonro - Raw Basal Area - No RE of site

##### Subset bonro data
bonro_surv <- subset(allSpDat, Sp == "bonro")

## Null model - species-level
bonro_15m_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = bonro_surv) 
summary(bonro_15m_null)
r.squaredGLMM(bonro_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
bonro_15mA <- glmer(Year_1_counts ~ consp_15m +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15mA) # conspBA is not significant
car::vif(bonro_15mA)
r.squaredGLMM(bonro_15mA) 

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
r.squaredGLMM(bonro_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# bonro_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                      hetAM_15m + hetECM_15m + 
#                      logPrior + log0yr +
#                      (1|survYear),
#                    family = poisson(link = "log"), data = bonro_surv)
# summary(bonro_15mB2)
# ## none of terms are significant
# car::vif(bonro_15mB2) 
# r.squaredGLMM(bonro_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(bonro_15m_null, bonro_15mA, bonro_15mB)
# AIC(bonro_15m_null, bonro_15mA, bonro_15mB, bonro_15mB2)
# More complex model is NOT better
anova(bonro_15mA, bonro_15mB)

#########################################################
################# Species-Level Test - wo - Raw Basal Area - Excludes RE of site

##### Subset wo data
wo_surv <- subset(allSpDat, Sp == "wo")

## Null model - species-level
wo_15m_null <- glmer(Year_1_counts ~ 1 +
                        (1|survYear),
                      family = poisson(link = "log"), data = wo_surv) 
summary(wo_15m_null)
r.squaredGLMM(wo_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
wo_15mA <- glmer(Year_1_counts ~ consp_15m +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_15mA) # conspBA is not significant
car::vif(wo_15mA)
r.squaredGLMM(wo_15mA) 

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
r.squaredGLMM(wo_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# wo_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                      hetAM_15m + hetECM_15m + 
#                      logPrior + log0yr +
#                      (1|survYear),
#                    family = poisson(link = "log"), data = wo_surv)
# summary(wo_15mB2)
# ## none of main effects or their interactions are significant
# car::vif(wo_15mB2) 
# r.squaredGLMM(wo_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(wo_15m_null, wo_15mA, wo_15mB)
# AIC(wo_15m_null, wo_15mA, wo_15mB, wo_15mB2)
# More complex model is NOT better
anova(wo_15mA, wo_15mB)

#########################################################
################# Species-Level Test - iw - Raw Basal Area - Excludes RE of site

##### Subset iw data
iw_surv <- subset(allSpDat, Sp == "iw")

## Null model - species-level
iw_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|survYear),
                     family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_null)
r.squaredGLMM(iw_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
iw_15mA <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15mA) # conspBA is not significant
car::vif(iw_15mA)
r.squaredGLMM(iw_15mA) 

## Testing Partial Panel B conceptual figure - Interactive spillover
## ***I THINK THIS IS ONE SPECIFIC JUSTIFIABLE CASE FOR NOT INCLUDING THE INTERACTIONS DUE TO HIGHLY INFLATED PARAMETER ESTIMATES (RESULTING FROM VERY LOW IW BA)
iw_15mB <- glmer(Year_1_counts ~ consp_15m + 
                   hetAM_15m + hetECM_15m + 
                   #consp_15m:hetAM_15m + # Panel B
                   #consp_15m:hetECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15mB)
## none of main effects or their interactions are significant
car::vif(iw_15mB) # One VIF >30
r.squaredGLMM(iw_15mB) 

# ## Testing Partial Panel B conceptual figure - Interactive spillover
# ## Removed consp_15m:AM_15m and consp_15m:ECM_15m interactions
# iw_15mB2 <- glmer(Year_1_counts ~ consp_15m + 
#                     hetAM_15m + hetECM_15m + 
#                     logPrior + log0yr +
#                     (1|survYear),
#                   family = poisson(link = "log"), data = iw_surv)
# summary(iw_15mB2)
# ## none of main effects or their interactions are significant
# car::vif(iw_15mB2) 
# r.squaredGLMM(iw_15mB2)

## Model comparison - species-level
### Compare AICs 
AIC(iw_15m_null, iw_15mA, iw_15mB)
# AIC(iw_15m_null, iw_15mA, iw_15mB, iw_15mB2)
# More complex model is NOT better
anova(iw_15mA, iw_15mB)

##################################################
########## Plotting effect sizes - Partial Panel C - species-level
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
  # ECM
  ab    = expression(italic("Fagus grandifolia")),
  iw    = expression(italic("Ostrya virginiana")),
  wo    = expression(italic("Quercus alba")),
  nro   = expression(italic("Quercus rubra")),
  bonro = expression(italic("Quercus velutina")),
  # AM
  rm    = expression(italic("Acer rubrum")),
  sm    = expression(italic("Acer saccharum")),
  wa    = expression(italic("Fraxinus americana")),
  bc    = expression(italic("Prunus serotina"))
)

## Colors matched to alphabetical order
species_cols <- c(
  # ECM (alphabetical by Latin name)
  ab    = "#C6E0EB",   # Fagus grandifolia
  iw    = "#4A90C2",   # Ostrya virginiana
  wo    = "#143D7A",   # Quercus alba
  nro   = "#1F5BA1",   # Quercus rubra
  bonro = "#7FB8D0",   # Quercus velutina
  # AM (alphabetical by Latin name)
  rm    = "#FB8A74",   # Acer rubrum
  sm    = "#E63E2D",   # Acer saccharum
  wa    = "#991F1F",   # Fraxinus americana
  bc    = "#FDCBC0"    # Prunus serotina
)

## Sizes for percent BA
species_sizes <- data.frame(species = c("ab", "bonro", "iw", "nro", 
                                        "wo", "bc", "rm", "sm", "wa"), 
                            percBA = c(2.43, 15.20, 0.11, 29.7, 10.30, 
                                       1.03, 9.40, 15.80, 1.82))

# Representative % BA values for legend
legend_percBA <- c(1, 10, 20)

###### Both AM and ECM
# --- Prepare plotting data
fixef_df_plot <- fixef_df %>%
  filter(estimate > -5) %>%  # exclude extreme outliers
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


fixef_df_plot %>%
  filter(is.na(estimate) | is.na(std.error) |
           is.na(estimate - std.error) |
           is.na(estimate + std.error) |
           is.na(y_pos) |
           !is.finite(estimate) |
           !is.finite(std.error) |
           !is.finite(estimate - std.error) |
           !is.finite(estimate + std.error) |
           !is.finite(log1p(percBA)))

#### Plot - All species
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

##### Plot - AM
ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

##### Plot - ECM
ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

#################################################
#################################################
################# Species-level Proportional BA Models
### PICK UP HERE!!!
################# Species-Level Test - RM - Proportional Basal Area

## Null model - species-level (proportional BA)
rm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv) 
summary(rm_15m_prop_null)
r.squaredGLMM(rm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
rm_propA <- glmer(Year_1_counts ~ propConsp +
                  logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = rm_surv)
summary(rm_propA)
car::vif(rm_propA)
r.squaredGLMM(rm_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
rm_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_propB) # Only main effect of propAM is significant
car::vif(rm_propB)
r.squaredGLMM(rm_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
rm_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = rm_surv)
summary(rm_propB2) # Both main effects of propConsp and propAM are significant
car::vif(rm_propB2)
r.squaredGLMM(rm_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(rm_15m_prop_null, rm_propA, rm_propB, rm_propB2)
# More complex model is NOT better
anova(rm_propB, rm_propB2)

################# 
################# Species-Level Test - sm - Proportional Basal Area

## Null model - species-level (proportional BA)
sm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = sm_surv) 
summary(sm_15m_prop_null)
r.squaredGLMM(sm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
sm_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(sm_propA) # Significant effect of propConsp
car::vif(sm_propA)
r.squaredGLMM(sm_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
sm_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(sm_propB) # Both main effects of propConsp and propAM are significant
car::vif(sm_propB)
r.squaredGLMM(sm_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
sm_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = sm_surv)
summary(sm_propB2) # Both main effects of propConsp and propAM are significant
car::vif(sm_propB2)
r.squaredGLMM(sm_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(sm_15m_prop_null, sm_propA, sm_propB, sm_propB2)
# More complex model is NOT better
anova(sm_propB, sm_propB2)

################# 
################# Species-Level Test - bc - Proportional Basal Area

## Null model - species-level (proportional BA)
bc_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_prop_null)
r.squaredGLMM(bc_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
bc_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_propA) # Significant effect of propConsp
car::vif(bc_propA)
r.squaredGLMM(bc_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
bc_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_propB) # Only main effect of propConsp significant
car::vif(bc_propB)
r.squaredGLMM(bc_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
bc_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = bc_surv)
summary(bc_propB2) # Only main effect of propConsp significant
car::vif(bc_propB2)
r.squaredGLMM(bc_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(bc_15m_prop_null, bc_propA, bc_propB, bc_propB2)
# More complex model is NOT better
anova(bc_propA, bc_propB)

#################
################# Species-Level Test - wa - Proportional Basal Area

## Null model - species-level (proportional BA)
wa_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = wa_surv) 
summary(wa_15m_prop_null)
r.squaredGLMM(wa_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
wa_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(wa_propA) # Significant effect of propConsp
car::vif(wa_propA)
r.squaredGLMM(wa_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
wa_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(wa_propB) # Only main effect of propConsp significant
car::vif(wa_propB)
r.squaredGLMM(wa_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
wa_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = wa_surv)
summary(wa_propB2) # Only main effect of propConsp significant
car::vif(wa_propB2)
r.squaredGLMM(wa_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(wa_15m_prop_null, wa_propA, wa_propB, wa_propB2)
# More complex model is NOT better
anova(wa_propA, wa_propB2)

################# 
################# Species-Level Test - ab - Proportional Basal Area

## Null model - species-level (proportional BA)
ab_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_prop_null)
r.squaredGLMM(ab_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
ab_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                     (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_propA) # Non-significant effect of propConsp
car::vif(ab_propA)
r.squaredGLMM(ab_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
ab_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                     (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_propB) # No terms significant
car::vif(ab_propB)
r.squaredGLMM(ab_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
ab_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                      (1|survYear),
                   family = poisson(link = "log"), data = ab_surv)
summary(ab_propB2) # Only main effect of propConsp significant
car::vif(ab_propB2)
r.squaredGLMM(ab_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(ab_15m_prop_null, ab_propA, ab_propB, ab_propB2)
# More complex model is NOT better
anova(ab_propA, ab_propB)

################# 
################# Species-Level Test - nro - Proportional Basal Area

## Null model - species-level (proportional BA)
nro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_prop_null)
r.squaredGLMM(nro_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
nro_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_propA) # Non-significant effect of propConsp
car::vif(nro_propA)
r.squaredGLMM(nro_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
nro_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_propB) # Main effects of propConsp and propAM are significant
car::vif(nro_propB)
r.squaredGLMM(nro_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
nro_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear),
                   family = poisson(link = "log"), data = nro_surv)
summary(nro_propB2) # Main effects of propConsp and propAM are significant
car::vif(nro_propB2)
r.squaredGLMM(nro_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(nro_15m_prop_null, nro_propA, nro_propB, nro_propB2)
# More complex model is NOT better
anova(nro_propB, nro_propB2)

################# 
################# Species-Level Test - bonro - Proportional Basal Area

## Null model - species-level (proportional BA)
bonro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                           (1|survYear),
                           family = poisson(link = "log"), data = bonro_surv) 
summary(bonro_15m_prop_null)
r.squaredGLMM(bonro_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
bonro_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_propA) # Non-significant effect of propConsp
car::vif(bonro_propA)
r.squaredGLMM(bonro_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
bonro_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_propB) # No terms significant
car::vif(bonro_propB)
r.squaredGLMM(bonro_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
bonro_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|survYear),
                   family = poisson(link = "log"), data = bonro_surv)
summary(bonro_propB2) # Only main effect of propConsp significant
car::vif(bonro_propB2)
r.squaredGLMM(bonro_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(bonro_15m_prop_null, bonro_propA, bonro_propB, bonro_propB2)
# More complex model is NOT better
anova(bonro_propA, bonro_propB2)

################# 
################# Species-Level Test - wo - Proportional Basal Area

## Null model - species-level (proportional BA)
wo_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                             (1|survYear),
                           family = poisson(link = "log"), data = wo_surv)
summary(wo_15m_prop_null)
r.squaredGLMM(wo_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
wo_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_propA) # Non-significant effect of propConsp
car::vif(wo_propA)
r.squaredGLMM(wo_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
wo_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_propB) # No terms significant
car::vif(wo_propB)
r.squaredGLMM(wo_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
wo_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|survYear),
                   family = poisson(link = "log"), data = wo_surv)
summary(wo_propB2) # Only main effect of propConsp significant
car::vif(wo_propB2)
r.squaredGLMM(wo_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(wo_15m_prop_null, wo_propA, wo_propB, wo_propB2)
# More complex model is NOT better
anova(wo_propA, wo_propB2)

################# 
################# Species-Level Test - iw - Proportional Basal Area

## Null model - species-level (proportional BA)
iw_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_prop_null)
r.squaredGLMM(iw_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
iw_propA <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_propA) # Non-significant effect of propConsp
car::vif(iw_propA)
r.squaredGLMM(iw_propA) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
iw_propB <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_propB) # No terms significant
car::vif(iw_propB)
r.squaredGLMM(iw_propB) 

## Testing Partial Panel B conceptual figure (with AM) - species-level (proportional BA)
# No interaction
iw_propB2 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                     (1|survYear),
                   family = poisson(link = "log"), data = iw_surv)
summary(iw_propB2) # Only main effect of propAM significant
car::vif(iw_propB2)
r.squaredGLMM(iw_propB2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(iw_15m_prop_null, iw_propA, iw_propB, iw_propB2)
# More complex model is NOT better
anova(iw_propB, iw_propB2)


##################################################
########## Plotting effect sizes - Partial Panel C - species-level
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propAM" ~ "Prop. AM BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Prop. Conspecific BA", "Prop. AM BA",
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM  <- extract_fixef(rm_propB2, "rm")
fixef_SM  <- extract_fixef(sm_propB2, "sm")
fixef_BC  <- extract_fixef(bc_propB2, "bc")
fixef_WA  <- extract_fixef(wa_propB2, "wa")
fixef_AB  <- extract_fixef(ab_propB2, "ab")
fixef_NRO <- extract_fixef(nro_propB2, "nro")
fixef_BONRO <- extract_fixef(bonro_propB2, "bonro")
fixef_WO  <- extract_fixef(wo_propB2, "wo")
fixef_IW  <- extract_fixef(iw_propB2, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

# Define desired vertical order
desired_order <- c("Prop. Conspecific BA", 
                   "Prop. AM BA", 
                   "Log(Recruit Density)", 
                   "Log(Prior Density)")

## Species labels (italics, alphabetized within AM & ECM)
species_labels <- c(
  # ECM
  ab    = expression(italic("Fagus grandifolia")),
  iw    = expression(italic("Ostrya virginiana")),
  wo    = expression(italic("Quercus alba")),
  nro   = expression(italic("Quercus rubra")),
  bonro = expression(italic("Quercus velutina")),
  # AM
  rm    = expression(italic("Acer rubrum")),
  sm    = expression(italic("Acer saccharum")),
  wa    = expression(italic("Fraxinus americana")),
  bc    = expression(italic("Prunus serotina"))
)

## Species colors (alphabetical within AM & ECM)
species_cols <- c(
  # ECM
  ab    = "#C6E0EB",   # Fagus grandifolia
  iw    = "#4A90C2",   # Ostrya virginiana
  wo    = "#143D7A",   # Quercus alba
  nro   = "#1F5BA1",   # Quercus rubra
  bonro = "#7FB8D0",   # Quercus velutina
  # AM
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
  filter(estimate > -5) %>%
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
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### AM only
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### ECM only
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error,
                     color = Model), height = 0.05) +
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
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
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


