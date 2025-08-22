
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
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
library(forcats)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1

#########################################################
#########################################################
#########################################################
################# All Species Models - Raw Basal Area

## Null model
all_sp_15m_null <- glmer(Year_1_counts ~ 1 +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_15m_null)
r.squaredGLMM(all_sp_15m_null) 

## Testing Panel A conceptual figure
all_sp_15m <- glmer(Year_1_counts ~ consp_15m*Seedling_myc +
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_15m)
car::vif(all_sp_15m)
r.squaredGLMM(all_sp_15m) 

## Testing Panel B conceptual figure
all_sp_15m1 <- glmer(Year_1_counts ~ consp_15m + Seedling_myc +
                       AM_15m + ECM_15m + 
                       consp_15m:AM_15m + # Panel B
                       consp_15m:ECM_15m + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m1)
car::vif(all_sp_15m1)
r.squaredGLMM(all_sp_15m1) 

## Testing Panel C conceptual figure
all_sp_15m2 <- glmer(Year_1_counts ~ consp_15m + Seedling_myc + 
                       AM_15m + ECM_15m + 
                       Seedling_myc:AM_15m + # Panel C
                       Seedling_myc:ECM_15m + # Panel C
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m2)
car::vif(all_sp_15m2)
r.squaredGLMM(all_sp_15m2) 

## Testing Panel D conceptual figure
all_sp_15m3 <- glmer(Year_1_counts ~ consp_15m + Seedling_myc + 
                       AM_15m + ECM_15m + 
                       consp_15m:Seedling_myc:AM_15m + # Panel D
                       consp_15m:Seedling_myc:ECM_15m + # Panel D
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m3) 
car::vif(all_sp_15m3) ## Some VIFs >10
r.squaredGLMM(all_sp_15m3) 

## Testing Panels B, C, and D of conceptual figure - Global model
all_sp_15m4 <- glmer(Year_1_counts ~ consp_15m + Seedling_myc + 
                    AM_15m + ECM_15m + 
                    consp_15m:AM_15m + # Panel B
                    consp_15m:ECM_15m + # Panel B
                    Seedling_myc:AM_15m + # Panel C
                    Seedling_myc:ECM_15m + # Panel C
                    consp_15m:Seedling_myc:AM_15m + # Panel D
                    consp_15m:Seedling_myc:ECM_15m + # Panel D
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear) + (1|Sp),
                  family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m4)
car::vif(all_sp_15m4)
r.squaredGLMM(all_sp_15m4) 

## Model comparison
### Compare AICs - lowest is for all_sp_15m4 but all_sp_15m3 close 2nd
AIC(all_sp_15m_null, all_sp_15m, all_sp_15m1, all_sp_15m2, all_sp_15m3, all_sp_15m4)
# More complex model fits better
anova(all_sp_15m3, all_sp_15m4) # Global model significantly better than Panel D

########## Plotting effect sizes - for Panel D
# 1. Extract fixed-effect summary from the model
fixef_df <- broom.mixed::tidy(all_sp_15m3, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%            # drop intercept
  select(term, estimate, std.error, p.value)

# 2. Map raw term names to clean labels you want shown on the plot
fixef_df <- fixef_df %>%
  mutate(component = case_when(
    term == "consp_15m" ~ "Conspecific BA",
    term == "Seedling_mycECM" ~ "Seedling mycorrhizal type",
    term == "AM_15m" ~ "Arbuscular mycorrhizal BA",
    term == "ECM_15m" ~ "Ectomycorrhizal BA",
    term == "consp_15m:Seedling_mycAM:AM_15m" ~ "Consp:Seedling myc AM:AM",
    term == "consp_15m:Seedling_mycECM:AM_15m" ~ "Consp:Seedling myc ECM:AM",
    term == "consp_15m:Seedling_mycAM:ECM_15m" ~ "Consp:Seedling myc AM:ECM",
    term == "consp_15m:Seedling_mycECM:ECM_15m" ~ "Consp:Seedling myc ECM:ECM",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))

# 3. Order components on the y-axis in the order you prefer
desired_order <- c("Conspecific BA","Arbuscular mycorrhizal BA",
                   "Ectomycorrhizal BA", "Seedling mycorrhizal type",
                   "Consp:Seedling myc AM:AM", "Consp:Seedling myc ECM:AM",
                   "Consp:Seedling myc AM:ECM", "Consp:Seedling myc ECM:ECM",
                   "Log(Recruit Density)","Log(Prior Density)")
fixef_df <- fixef_df %>%
  mutate(component = factor(component, levels = desired_order))

# 4. Significance flag (alpha = 0.05)
fixef_df <- fixef_df %>%
  mutate(Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"))

# Identify large values and shift them
fixef_df <- fixef_df %>%
  mutate(plot_est = ifelse(estimate > 1, estimate - 1, estimate))

# Make your factor for y-axis
fixef_df$component <- factor(fixef_df$component, levels = rev(desired_order))

### Plot
ggplot(fixef_df, aes(x = plot_est, y = component)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_errorbarh(aes(xmin = plot_est - std.error, xmax = plot_est + std.error),
                 height = 0.1, color = "black") +
  geom_point(aes(shape = Significant, fill = Significant), 
             size = 3.5, color = "black", stroke = 1) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  # Custom axis labels with blank for dummy tick
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.25, 1.5),
    labels = c("-0.5", "0", "0.5", "2.0", "", "2.5")) +
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))
#########################################################

#########################################################
################# All Species Models - Proportional Basal Area

## Null model
all_sp_prop_null <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_prop_null)
r.squaredGLMM(all_sp_prop_null) 

## Testing Panel A conceptual figure (proportional BA)
all_sp_prop <- glmer(Year_1_counts ~ propConsp*Seedling_myc +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_prop)
car::vif(all_sp_prop)
r.squaredGLMM(all_sp_prop) 

## Testing Panel B conceptual figure (proportional BA)
all_sp_prop1 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                       propAM + # Panel B
                       propConsp:propAM + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop1)
car::vif(all_sp_prop1)
r.squaredGLMM(all_sp_prop1) 

## Testing Panel C conceptual figure (proportional BA)
all_sp_prop2 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                        propAM + 
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop2)
car::vif(all_sp_prop2)
r.squaredGLMM(all_sp_prop2) 

## Testing Panel D conceptual figure (proportional BA)
all_sp_prop3 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                        propAM + 
                        propConsp:Seedling_myc:propAM + # Panel D
                        logPrior + log0yr +
                        (1|SiteID) + (1|survYear) + (1|Sp),
                      family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop3)
car::vif(all_sp_prop3)
r.squaredGLMM(all_sp_prop3) 

## Testing Panels B, C, and D of conceptual figure - Global model (proportional BA)
all_sp_prop4 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                       propAM + 
                       propConsp:propAM + # Panel B
                       Seedling_myc:propAM + # Panel C
                       propConsp:Seedling_myc:propAM + # Panel D
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop4)
car::vif(all_sp_prop4)
r.squaredGLMM(all_sp_prop4) 

### Compare AICs - 
AIC(all_sp_prop_null, all_sp_prop, all_sp_prop1, 
    all_sp_prop2, all_sp_prop3, all_sp_prop4)
# More complex model fits better
anova(all_sp_prop3, all_sp_prop4) # Global model not better than Panel D
anova(all_sp_prop1, all_sp_prop3) # Panel D significantly better than Panel B
# We lose sensitivity when we use the proportional data - consider AICs

## ALL Model AIC comparison (raw BA and proportional)
# Named list of models
models <- list(
  all_sp_15m_null = all_sp_15m_null,
  all_sp_15m = all_sp_15m,
  all_sp_15m1 = all_sp_15m1,
  all_sp_15m2 = all_sp_15m2,
  all_sp_15m3 = all_sp_15m3,
  all_sp_15m4 = all_sp_15m4,
  all_sp_prop_null = all_sp_prop_null,
  all_sp_prop = all_sp_prop,
  all_sp_prop1 = all_sp_prop1,
  all_sp_prop2 = all_sp_prop2,
  all_sp_prop3 = all_sp_prop3,
  all_sp_prop4 = all_sp_prop4
)

# Extract and order AICs
aic_vals <- sapply(models, AIC)
aic_vals[order(aic_vals)]
## The raw BA models are better!

########## Plotting effect sizes - for Panel D
# 1. Extract fixed-effect summary from the model
fixef_df <- broom.mixed::tidy(all_sp_prop3, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%            # drop intercept
  select(term, estimate, std.error, p.value)

# 2. Map raw term names to clean labels you want shown on the plot
fixef_df <- fixef_df %>%
  mutate(component = case_when(
    term == "propConsp" ~ "Prop. Conspecific BA",
    term == "Seedling_mycECM" ~ "Seedling mycorrhizal type",
    term == "propAM" ~ "Prop. Arbuscular mycorrhizal BA",
    term == "propConsp:Seedling_mycAM:propAM" ~ "propConsp:Seedling myc AM:propAM",
    term == "propConsp:Seedling_mycECM:propAM" ~ "propConsp:Seedling myc ECM:propAM",
    term == "logPrior" ~ "Log(Prior Density)",
    term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ term))

# 3. Order components on the y-axis in the order you prefer
desired_order <- c("Prop. Conspecific BA","Prop. Arbuscular mycorrhizal BA",
                   "Seedling mycorrhizal type",
                   "propConsp:Seedling myc AM:propAM",
                   "propConsp:Seedling myc ECM:propAM",
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
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error),
                 height = 0.1, color = "black") +
  geom_point(aes(shape = Significant, fill = Significant), 
             size = 3.5, color = "black", stroke = 1) +
  scale_shape_manual(values = c("yes" = 16, "no" = 21)) +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))
#########################################################
#########################################################
#########################################################

#########################################################
#########################################################
################# Mycorrhizal Type Models - Raw Basal Area

## Subset out just AM and just ECM species
allSpDat_AM = subset(allSpDat, Seedling_myc == "AM")
allSpDat_ECM = subset(allSpDat, Seedling_myc == "ECM")

## Null model - AM
AM_sp_15m_null <- glmer(Year_1_counts ~ 1 +
                           (1|SiteID) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat_AM) 
summary(AM_sp_15m_null)
r.squaredGLMM(AM_sp_15m_null) 

## Null model - ECM
ECM_sp_15m_null <- glmer(Year_1_counts ~ 1 +
                          (1|SiteID) + (1|survYear) + (1|Sp),
                        family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15m_null)
r.squaredGLMM(ECM_sp_15m_null) 

## Testing Panel A conceptual figure - AM
AM_sp_15m <- glmer(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_AM) 
summary(AM_sp_15m)
car::vif(AM_sp_15m)
r.squaredGLMM(AM_sp_15m) 

## Testing Panel A conceptual figure - ECM
ECM_sp_15m <- glmer(Year_1_counts ~ consp_15m +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15m)
car::vif(ECM_sp_15m)
r.squaredGLMM(ECM_sp_15m) 

## Testing Panel B conceptual figure - AM
AM_sp_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                       AM_15m + ECM_15m + 
                       consp_15m:AM_15m + # Panel B
                       consp_15m:ECM_15m + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_AM)
summary(AM_sp_15m1)
car::vif(AM_sp_15m1)
r.squaredGLMM(AM_sp_15m1) 

## Testing Panel B conceptual figure - ECM
ECM_sp_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                      AM_15m + ECM_15m + 
                      consp_15m:AM_15m + # Panel B
                      consp_15m:ECM_15m + # Panel B
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15m1)
car::vif(ECM_sp_15m1)
r.squaredGLMM(ECM_sp_15m1) 

## Model comparison - AM
### Compare AICs - lowest is most complex - AM_sp_15m1
AIC(AM_sp_15m_null, AM_sp_15m, AM_sp_15m1)
## Model comparison - ECM
### Compare AICs - lowest is most complex (just barely) - ECM_sp_15m1
AIC(ECM_sp_15m_null, ECM_sp_15m, ECM_sp_15m1)
anova(ECM_sp_15m,ECM_sp_15m1)

########## Plotting effect sizes - for Panel B (most complex mycorrhizal version)
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "consp_15m" ~ "Conspecific BA",
        term == "AM_15m" ~ "Arbuscular mycorrhizal BA",
        term == "ECM_15m" ~ "Ectomycorrhizal BA",
        term == "consp_15m:AM_15m" ~ "Consp:AM",
        term == "consp_15m:ECM_15m" ~ "Consp:ECM",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Conspecific BA", "Arbuscular mycorrhizal BA", 
        "Ectomycorrhizal BA", "Consp:AM", "Consp:ECM",
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for both models
fixef_AM  <- extract_fixef(AM_sp_15m1, "AM")
fixef_ECM <- extract_fixef(ECM_sp_15m1, "ECM")

## Combine
fixef_df <- bind_rows(fixef_AM, fixef_ECM)

# Define desired vertical order
desired_order <- c("Conspecific BA",
  "Arbuscular mycorrhizal BA","Ectomycorrhizal BA",
  "Consp:AM","Consp:ECM",
  "Log(Recruit Density)","Log(Prior Density)")

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
AM_sp_15mprop <- glmer(Year_1_counts ~ propConsp +
                     logPrior + log0yr +
                     (1|SiteID) + (1|survYear) + (1|Sp),
                   family = poisson(link = "log"), data = allSpDat_AM) 
summary(AM_sp_15mprop)
car::vif(AM_sp_15mprop)
r.squaredGLMM(AM_sp_15mprop) 

## Testing Panel A conceptual figure - ECM - Proportional BA
ECM_sp_15mprop <- glmer(Year_1_counts ~ propConsp +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_ECM) 
summary(ECM_sp_15mprop)
car::vif(ECM_sp_15mprop)
r.squaredGLMM(ECM_sp_15mprop) 

## Testing Panel B conceptual figure - AM - Proportional BA
AM_sp_15m1prop <- glmer(Year_1_counts ~ propConsp + 
                      propAM + propConsp:propAM + 
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_AM)
summary(AM_sp_15m1prop)
car::vif(AM_sp_15m1prop)
r.squaredGLMM(AM_sp_15m1prop) 

## Testing Panel B conceptual figure - ECM - Proportional BA
ECM_sp_15m1prop <- glmer(Year_1_counts ~ propConsp + 
                       propAM + propConsp:propAM + 
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_ECM)
summary(ECM_sp_15m1prop)
car::vif(ECM_sp_15m1prop)
r.squaredGLMM(ECM_sp_15m1prop) 

## Model comparison - AM
### Compare AICs - lowest is most complex - AM_sp_15m1prop
AIC(AM_sp_15m_nullprop, AM_sp_15mprop, AM_sp_15m1prop)
## Model comparison - ECM
### Compare AICs - lowest is most complex - ECM_sp_15m1prop
AIC(ECM_sp_15m_nullprop, ECM_sp_15mprop, ECM_sp_15m1prop)
anova(ECM_sp_15mprop, ECM_sp_15m1prop)

########## Plotting effect sizes - for Panel B - proportional BA
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propAM" ~ "Prop. Arbuscular mycorrhizal BA",
        term == "propConsp:propAM" ~ "propConsp:propAM",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(component = factor(component, levels = c(
        "Prop. Conspecific BA", "Prop. Arbuscular mycorrhizal BA", 
        "propConsp:propAM","Log(Recruit Density)","Log(Prior Density)")))
}

## Extract for both models
fixef_AM  <- extract_fixef(AM_sp_15m1prop, "AM")
fixef_ECM <- extract_fixef(ECM_sp_15m1prop, "ECM")

## Combine
fixef_df <- bind_rows(fixef_AM, fixef_ECM)

# Define desired vertical order
desired_order <- c("Prop. Conspecific BA", "Prop. Arbuscular mycorrhizal BA", 
                   "propConsp:propAM","Log(Recruit Density)","Log(Prior Density)")

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
rm_15m <- glmer(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|SiteID) + (1|survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_15m)
car::vif(rm_15m)
r.squaredGLMM(rm_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
rm_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_15m1)
car::vif(rm_15m1)
r.squaredGLMM(rm_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
rm_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                       AM_15m + ECM_15m + 
                       consp_15m:AM_15m + # Panel B
                       consp_15m:ECM_15m + # Panel B
                       logPrior + log0yr +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv)
summary(rm_15m2)
car::vif(rm_15m2)
r.squaredGLMM(rm_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for rm_15m4 but rm_15m3 close 2nd
AIC(rm_15m_null, rm_15m, rm_15m1, rm_15m2)
# More complex model is NOT better
anova(rm_15m1, rm_15m2)

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
sm_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = sm_surv)
summary(sm_15m)
car::vif(sm_15m)
r.squaredGLMM(sm_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
sm_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_15m1)
car::vif(sm_15m1)
r.squaredGLMM(sm_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
sm_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_15m2)
car::vif(sm_15m2)
r.squaredGLMM(sm_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for sm_15m4 but sm_15m3 close 2nd
AIC(sm_15m_null, sm_15m, sm_15m1, sm_15m2)
# More complex model is NOT better
anova(sm_15m1, sm_15m2)

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
bc_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = bc_surv)
summary(bc_15m)
car::vif(bc_15m)
r.squaredGLMM(bc_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
bc_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15m1)
car::vif(bc_15m1)
r.squaredGLMM(bc_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
bc_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15m2)
car::vif(bc_15m2)
r.squaredGLMM(bc_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for bc_15m4 but bc_15m3 close 2nd
AIC(bc_15m_null, bc_15m, bc_15m1, bc_15m2)
# More complex model is NOT better than bc_15m
anova(bc_15m, bc_15m2)

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
wa_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = wa_surv)
summary(wa_15m)
car::vif(wa_15m)
r.squaredGLMM(wa_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
wa_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_15m1)
car::vif(wa_15m1)
r.squaredGLMM(wa_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
wa_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_15m2)
car::vif(wa_15m2)
r.squaredGLMM(wa_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for wa_15m4 but wa_15m3 close 2nd
AIC(wa_15m_null, wa_15m, wa_15m1, wa_15m2)
# More complex model is NOT better
anova(wa_15m, wa_15m2)

#########################################################
################# Species-Level Test - ab - Raw Basal Area - Excludes RE of SiteID

##### Subset ab data
ab_surv <- subset(allSpDat, Sp == "ab")

## Null model - species-level
ab_15m_null <- glmer(Year_1_counts ~ 1 + (1|survYear),
                     family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_null)
r.squaredGLMM(ab_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
ab_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr + (1|survYear),
                family = poisson(link = "log"), data = ab_surv)
summary(ab_15m)
car::vif(ab_15m)
r.squaredGLMM(ab_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
ab_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_15m1)
car::vif(ab_15m1)
r.squaredGLMM(ab_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
ab_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_15m2)
car::vif(ab_15m2)
r.squaredGLMM(ab_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for ab_15m4 but ab_15m3 close 2nd
AIC(ab_15m_null, ab_15m, ab_15m1, ab_15m2)
# More complex model is NOT better
anova(ab_15m, ab_15m1)

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
nro_15m <- glmer(Year_1_counts ~ consp_15m +
                   logPrior + log0yr + (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_15m)
car::vif(nro_15m)
r.squaredGLMM(nro_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
nro_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                    AM_15m + ECM_15m + 
                    logPrior + log0yr + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_15m1)
car::vif(nro_15m1)
r.squaredGLMM(nro_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
nro_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                    AM_15m + ECM_15m + 
                    consp_15m:AM_15m + # Panel B
                    consp_15m:ECM_15m + # Panel B
                    logPrior + log0yr +
                     (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_15m2)
car::vif(nro_15m2)
r.squaredGLMM(nro_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for nro_15m4 but nro_15m3 close 2nd
AIC(nro_15m_null, nro_15m, nro_15m1, nro_15m2)
# More complex model is NOT better
anova(nro_15m1, nro_15m2)

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
bonro_15m <- glmer(Year_1_counts ~ consp_15m +
                     logPrior + log0yr +
                     (1|survYear),
                   family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15m)
car::vif(bonro_15m)
r.squaredGLMM(bonro_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
bonro_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                      AM_15m + ECM_15m + 
                      logPrior + log0yr +
                      (1|survYear),
                    family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15m1)
car::vif(bonro_15m1)
r.squaredGLMM(bonro_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
bonro_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                      AM_15m + ECM_15m + 
                      consp_15m:AM_15m + # Panel B
                      consp_15m:ECM_15m + # Panel B
                      logPrior + log0yr +
                      (1|survYear),
                    family = poisson(link = "log"), data = bonro_surv)
summary(bonro_15m2)
car::vif(bonro_15m2)
r.squaredGLMM(bonro_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for bonro_15m4 but bonro_15m3 close 2nd
AIC(bonro_15m_null, bonro_15m, bonro_15m1, bonro_15m2)
# More complex model is NOT better
anova(bonro_15m, bonro_15m1)

#########################################################
################# Species-Level Test - wo - Raw Basal Area - Excludes RE of site

##### Subset wo data
wo_surv <- subset(allSpDat, Sp == "wo")

## Null model - species-level
wo_15m_null <- glmer(Year_1_counts ~ 1 + (1|survYear),
                     family = poisson(link = "log"), data = wo_surv) 
summary(wo_15m_null)
r.squaredGLMM(wo_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
wo_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr + 
                  (1|survYear),
                family = poisson(link = "log"), data = wo_surv)
summary(wo_15m)
car::vif(wo_15m)
r.squaredGLMM(wo_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
wo_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_15m1)
car::vif(wo_15m1)
r.squaredGLMM(wo_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
wo_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = wo_surv)
summary(wo_15m2)
car::vif(wo_15m2)
r.squaredGLMM(wo_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for wo_15m4 but wo_15m3 close 2nd
AIC(wo_15m_null, wo_15m, wo_15m1, wo_15m2)
# More complex model is NOT better
anova(wo_15m, wo_15m1)

#########################################################
################# Species-Level Test - iw - Raw Basal Area - Excludes RE of site

##### Subset iw data
iw_surv <- subset(allSpDat, Sp == "iw")

## Null model - species-level
iw_15m_null <- glmer(Year_1_counts ~ 1 + (1|survYear),
                     family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_null)
r.squaredGLMM(iw_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
iw_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr + 
                  (1|survYear),
                family = poisson(link = "log"), data = iw_surv)
summary(iw_15m)
car::vif(iw_15m)
r.squaredGLMM(iw_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
iw_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15m1)
car::vif(iw_15m1)
r.squaredGLMM(iw_15m1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level
iw_15m2 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   consp_15m:AM_15m + # Panel B
                   consp_15m:ECM_15m + # Panel B
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_15m2)
car::vif(iw_15m2)
r.squaredGLMM(iw_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for iw_15m4 but iw_15m3 close 2nd
AIC(iw_15m_null, iw_15m, iw_15m1, iw_15m2)
# More complex model is NOT better
anova(iw_15m, iw_15m1)

##################################################
########## Plotting effect sizes - Partial Panel C - species-level
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "consp_15m" ~ "Conspecific BA",
        term == "AM_15m" ~ "Arbuscular mycorrhizal BA",
        term == "ECM_15m" ~ "Ectomycorrhizal BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Conspecific BA", "Arbuscular mycorrhizal BA", 
        "Ectomycorrhizal BA", "Log(Recruit Density)", 
        "Log(Prior Density)")))
}

## Extract for each species
fixef_RM  <- extract_fixef(rm_15m1, "rm")
fixef_SM  <- extract_fixef(sm_15m1, "sm")
fixef_BC  <- extract_fixef(bc_15m1, "bc")
fixef_WA  <- extract_fixef(wa_15m1, "wa")
fixef_AB  <- extract_fixef(ab_15m1, "ab")
fixef_NRO  <- extract_fixef(nro_15m1, "nro")
fixef_BONRO  <- extract_fixef(bonro_15m1, "bonro")
fixef_WO  <- extract_fixef(wo_15m1, "wo")
fixef_IW  <- extract_fixef(iw_15m1, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

# Define desired vertical order
desired_order <- c("Conspecific BA", "Arbuscular mycorrhizal BA", 
                   "Ectomycorrhizal BA", "Log(Recruit Density)", 
                   "Log(Prior Density)")

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

###### Both AM and ECM
# --- Prepare plotting data
fixef_df_plot <- fixef_df %>%
  filter(estimate > -5) %>%  # exclude extreme outliers
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
      Model == "rm" ~ 0.12,
      Model == "bc" ~ 0.09,
      Model == "sm" ~ 0.06,
      Model == "wa" ~ 0.09,
      Model == "ab" ~ -0,
      Model == "nro" ~ -0.03,
      Model == "bonro" ~ -0.06,
      Model == "wo" ~ -0.09,
      Model == "iw" ~ -0.12,
      TRUE ~ 0),
    # final y-position for plotting
    y_pos = as.numeric(component) + y_offset,
    # fill color: species color if significant, white otherwise
    fill_col = ifelse(Significant == "yes", 
                      species_cols[Model], "white"))

#### Plot
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
    labels = levels(fixef_df$component),
    trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"))

###### AM vs EM only - Prepare plotting data
fixef_df_plot <- fixef_df %>%
  filter(estimate > -5) %>%  # exclude extreme outliers
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
         percBA = species_sizes$percBA[match(Model, 
                                    species_sizes$species)],
         # color for species
         color_col = species_cols[Model],
         # species-specific vertical offset
         y_offset = case_when(
           Model == "rm" ~ 0.06,
           Model == "bc" ~ 0.03,
           Model == "sm" ~ 0,
           Model == "wa" ~ -0.03,
           Model == "ab" ~ 0.06,
           Model == "nro" ~ 0.03,
           Model == "bonro" ~ 0,
           Model == "wo" ~ -0.03,
           Model == "iw" ~ -0.06,
           TRUE ~ 0),
         # final y-position for plotting
         y_pos = as.numeric(component) + y_offset,
         # fill color: species color if significant, white otherwise
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"))

##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

##### Plot - AM
ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

##### Plot - ECM
ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

#################################################
#################################################
################# Species-level Proportional BA Models

################# 
################# 
################# Species-Level Test - RM - Proportional Basal Area

## Null model - species-level (proportional BA)
rm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                       (1|SiteID) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv) 
summary(rm_15m_prop_null)
r.squaredGLMM(rm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
rm_prop <- glmer(Year_1_counts ~ propConsp +
                  logPrior + log0yr +
                  (1|SiteID) + (1|survYear),
                family = poisson(link = "log"), data = rm_surv)
summary(rm_prop)
car::vif(rm_prop)
r.squaredGLMM(rm_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
rm_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_prop1)
car::vif(rm_prop1)
r.squaredGLMM(rm_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
rm_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_prop2)
car::vif(rm_prop2)
r.squaredGLMM(rm_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(rm_15m_prop_null, rm_prop, rm_prop1, rm_prop2)
# More complex model is NOT better
anova(rm_prop1, rm_prop2)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  rm_15m_null = rm_15m_null,
  rm_15m = rm_15m,
  rm_15m1 = rm_15m1,
  rm_15m2 = rm_15m2,
  rm_15m_prop_null = rm_15m_prop_null,
  rm_prop = rm_prop,
  rm_prop1 = rm_prop1,
  rm_prop2 = rm_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are better!

################# 
################# 
################# Species-Level Test - sm - Proportional Basal Area

## Null model - species-level (proportional BA)
sm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = sm_surv) 
summary(sm_15m_prop_null)
r.squaredGLMM(sm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
sm_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_prop)
car::vif(sm_prop)
r.squaredGLMM(sm_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
sm_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(sm_prop1)
car::vif(sm_prop1)
r.squaredGLMM(sm_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
sm_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(sm_prop2)
car::vif(sm_prop2)
r.squaredGLMM(sm_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(sm_15m_prop_null, sm_prop, sm_prop1, sm_prop2)
# More complex model is NOT better
anova(sm_prop1, sm_prop2)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  sm_15m_null = sm_15m_null,
  sm_15m = sm_15m,
  sm_15m1 = sm_15m1,
  sm_15m2 = sm_15m2,
  sm_15m_prop_null = sm_15m_prop_null,
  sm_prop = sm_prop,
  sm_prop1 = sm_prop1,
  sm_prop2 = sm_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are better!

################# 
################# 
################# Species-Level Test - bc - Proportional Basal Area

## Null model - species-level (proportional BA)
bc_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_prop_null)
r.squaredGLMM(bc_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
bc_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_prop)
car::vif(bc_prop)
r.squaredGLMM(bc_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
bc_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_prop1)
car::vif(bc_prop1)
r.squaredGLMM(bc_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
bc_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_prop2)
car::vif(bc_prop2)
r.squaredGLMM(bc_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(bc_15m_prop_null, bc_prop, bc_prop1, bc_prop2)
# More complex model is NOT better
anova(bc_prop, bc_prop2)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  bc_15m_null = bc_15m_null,
  bc_15m = bc_15m,
  bc_15m1 = bc_15m1,
  bc_15m2 = bc_15m2,
  bc_15m_prop_null = bc_15m_prop_null,
  bc_prop = bc_prop,
  bc_prop1 = bc_prop1,
  bc_prop2 = bc_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are better, but not by much

################# 
#################
################# Species-Level Test - wa - Proportional Basal Area

## Null model - species-level (proportional BA)
wa_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = wa_surv) 
summary(wa_15m_prop_null)
r.squaredGLMM(wa_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
wa_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_prop)
car::vif(wa_prop)
r.squaredGLMM(wa_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
wa_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(wa_prop1)
car::vif(wa_prop1)
r.squaredGLMM(wa_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
wa_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(wa_prop2)
car::vif(wa_prop2)
r.squaredGLMM(wa_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(wa_15m_prop_null, wa_prop, wa_prop1, wa_prop2)
# More complex model is NOT better
anova(wa_prop, wa_prop1)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  wa_15m_null = wa_15m_null,
  wa_15m = wa_15m,
  wa_15m1 = wa_15m1,
  wa_15m2 = wa_15m2,
  wa_15m_prop_null = wa_15m_prop_null,
  wa_prop = wa_prop,
  wa_prop1 = wa_prop1,
  wa_prop2 = wa_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are much better!

################# 
################# 
################# Species-Level Test - ab - Proportional Basal Area

## Null model - species-level (proportional BA)
ab_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = ab_surv) 
summary(ab_15m_prop_null)
r.squaredGLMM(ab_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
ab_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = ab_surv)
summary(ab_prop)
car::vif(ab_prop)
r.squaredGLMM(ab_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
ab_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_prop1)
car::vif(ab_prop1)
r.squaredGLMM(ab_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
ab_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(ab_prop2)
car::vif(ab_prop2)
r.squaredGLMM(ab_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(ab_15m_prop_null, ab_prop, ab_prop1, ab_prop2)
# More complex model is NOT better
anova(ab_prop, ab_prop1)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  ab_15m_null = ab_15m_null,
  ab_15m = ab_15m,
  ab_15m1 = ab_15m1,
  ab_15m2 = ab_15m2,
  ab_15m_prop_null = ab_15m_prop_null,
  ab_prop = ab_prop,
  ab_prop1 = ab_prop1,
  ab_prop2 = ab_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are better, but just barely

################# 
################# 
################# Species-Level Test - nro - Proportional Basal Area

## Null model - species-level (proportional BA)
nro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|SiteID) + (1|survYear),
                          family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_prop_null)
r.squaredGLMM(nro_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
nro_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|SiteID) + (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_prop)
car::vif(nro_prop)
r.squaredGLMM(nro_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
nro_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_prop1)
car::vif(nro_prop1)
r.squaredGLMM(nro_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
nro_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|SiteID) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_prop2)
car::vif(nro_prop2)
r.squaredGLMM(nro_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(nro_15m_prop_null, nro_prop, nro_prop1, nro_prop2)
# More complex model is NOT better
anova(nro_prop1, nro_prop2)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  nro_15m_null = nro_15m_null,
  nro_15m = nro_15m,
  nro_15m1 = nro_15m1,
  nro_15m2 = nro_15m2,
  nro_15m_prop_null = nro_15m_prop_null,
  nro_prop = nro_prop,
  nro_prop1 = nro_prop1,
  nro_prop2 = nro_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The raw BA models are better, but not by much
anova(nro_15m1, nro_prop1)

################# 
################# 
################# Species-Level Test - bonro - Proportional Basal Area

## Null model - species-level (proportional BA)
bonro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                           (1|survYear),
                           family = poisson(link = "log"), data = bonro_surv) 
summary(bonro_15m_prop_null)
r.squaredGLMM(bonro_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
bonro_prop <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                   (1|survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonro_prop)
car::vif(bonro_prop)
r.squaredGLMM(bonro_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
bonro_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                   (1|survYear),
                   family = poisson(link = "log"), data = bonro_surv)
summary(bonro_prop1)
car::vif(bonro_prop1)
r.squaredGLMM(bonro_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
bonro_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                     logPrior + log0yr +
                   (1|survYear),
                   family = poisson(link = "log"), data = bonro_surv)
summary(bonro_prop2)
car::vif(bonro_prop2)
r.squaredGLMM(bonro_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(bonro_15m_prop_null, bonro_prop, bonro_prop1, bonro_prop2)
# More complex model is NOT better
anova(bonro_prop, bonro_prop1)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  bonro_15m_null = bonro_15m_null,
  bonro_15m = bonro_15m,
  bonro_15m1 = bonro_15m1,
  bonro_15m2 = bonro_15m2,
  bonro_15m_prop_null = bonro_15m_prop_null,
  bonro_prop = bonro_prop,
  bonro_prop1 = bonro_prop1,
  bonro_prop2 = bonro_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The proportional BA model is slightly better than raw BA, but not by much
anova(bonro_prop, bonro_prop1)

################# 
################# 
################# Species-Level Test - wo - Proportional Basal Area

## Null model - species-level (proportional BA)
wo_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                             (1|survYear),
                           family = poisson(link = "log"), data = wo_surv)
summary(wo_15m_prop_null)
r.squaredGLMM(wo_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
wo_prop <- glmer(Year_1_counts ~ propConsp +
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(wo_prop)
car::vif(wo_prop)
r.squaredGLMM(wo_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
wo_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                     logPrior + log0yr +
                      (1|survYear),
                   family = poisson(link = "log"), data = wo_surv)
summary(wo_prop1)
car::vif(wo_prop1)
r.squaredGLMM(wo_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
wo_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                     logPrior + log0yr +
                     (1|survYear),
                   family = poisson(link = "log"), data = wo_surv)
summary(wo_prop2)
car::vif(wo_prop2)
r.squaredGLMM(wo_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(wo_15m_prop_null, wo_prop, wo_prop1, wo_prop2)
# More complex model is NOT better
anova(wo_prop, wo_prop1)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  wo_15m_null = wo_15m_null,
  wo_15m = wo_15m,
  wo_15m1 = wo_15m1,
  wo_15m2 = wo_15m2,
  wo_15m_prop_null = wo_15m_prop_null,
  wo_prop = wo_prop,
  wo_prop1 = wo_prop1,
  wo_prop2 = wo_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The proportional BA models are better, but not by much
anova(wo_prop, wo_prop1)


################# 
################# 
################# Species-Level Test - iw - Proportional Basal Area

## Null model - species-level (proportional BA)
iw_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|survYear),
                          family = poisson(link = "log"), data = iw_surv) 
summary(iw_15m_prop_null)
r.squaredGLMM(iw_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
iw_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|survYear),
                 family = poisson(link = "log"), data = iw_surv)
summary(iw_prop)
car::vif(iw_prop)
r.squaredGLMM(iw_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
iw_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_prop1)
car::vif(iw_prop1)
r.squaredGLMM(iw_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
iw_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iw_prop2)
car::vif(iw_prop2)
r.squaredGLMM(iw_prop2) 

## Model comparison - species-level (proportional BA)
### Compare AICs
AIC(iw_15m_prop_null, iw_prop, iw_prop1, iw_prop2)
# More complex model is NOT better
anova(iw_prop1, iw_prop2)

## ALL Model AIC comparison - species-level (raw BA and proportional)
# Named list of models
models_sp <- list(
  iw_15m_null = iw_15m_null,
  iw_15m = iw_15m,
  iw_15m1 = iw_15m1,
  iw_15m2 = iw_15m2,
  iw_15m_prop_null = iw_15m_prop_null,
  iw_prop = iw_prop,
  iw_prop1 = iw_prop1,
  iw_prop2 = iw_prop2)
# Extract and order AICs
aic_vals_sp <- sapply(models_sp, AIC)
aic_vals_sp[order(aic_vals_sp)]
## The proportional BA models are better, but not by much
anova(iw_prop1, iw_prop2)


##################################################
########## Plotting effect sizes - Partial Panel C - species-level
extract_fixef <- function(model, model_type) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(
      component = case_when(
        term == "propConsp" ~ "Prop. Conspecific BA",
        term == "propAM" ~ "Prop. Arbuscular mycorrhizal BA",
        term == "logPrior" ~ "Log(Prior Density)",
        term == "log0yr" ~ "Log(Recruit Density)",
        TRUE ~ term),
      Significant = ifelse(!is.na(p.value) & 
                             p.value <= 0.05, "yes", "no"),
      Model = model_type) %>%
    mutate(
      component = factor(component, levels = c(
        "Prop. Conspecific BA", "Prop. Arbuscular mycorrhizal BA", 
        "Log(Recruit Density)", "Log(Prior Density)")))
}

## Extract for each species
fixef_RM  <- extract_fixef(rm_prop1, "rm")
fixef_SM  <- extract_fixef(sm_prop1, "sm")
fixef_BC  <- extract_fixef(bc_prop1, "bc")
fixef_WA  <- extract_fixef(wa_prop1, "wa")
fixef_AB  <- extract_fixef(ab_prop1, "ab")
fixef_NRO  <- extract_fixef(nro_prop1, "nro")
fixef_BONRO  <- extract_fixef(bonro_prop1, "bonro")
fixef_WO  <- extract_fixef(wo_prop1, "wo")
fixef_IW  <- extract_fixef(iw_prop1, "iw")

## Combine
fixef_df <- bind_rows(fixef_RM, fixef_SM, fixef_BC,
                      fixef_WA, fixef_AB, fixef_NRO,
                      fixef_BONRO, fixef_WO, fixef_IW)

# Define desired vertical order
desired_order <- c("Prop. Conspecific BA", 
                   "Prop. Arbuscular mycorrhizal BA", 
                   "Log(Recruit Density)", "Log(Prior Density)")

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

###### Both AM and ECM
# --- Prepare plotting data
fixef_df_plot <- fixef_df %>%
  filter(estimate > -5) %>%  # exclude extreme outliers
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
           Model == "rm" ~ 0.12,
           Model == "bc" ~ 0.09,
           Model == "sm" ~ 0.06,
           Model == "wa" ~ 0.09,
           Model == "ab" ~ -0,
           Model == "nro" ~ -0.03,
           Model == "bonro" ~ -0.06,
           Model == "wo" ~ -0.09,
           Model == "iw" ~ -0.12,
           TRUE ~ 0),
         # final y-position for plotting
         y_pos = as.numeric(component) + y_offset,
         # fill color: species color if significant, white otherwise
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"))

#### Plot
ggplot(fixef_df_plot, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

###### AM vs EM only - Prepare plotting data
fixef_df_plot <- fixef_df %>%
  filter(estimate > -5) %>%  # exclude extreme outliers
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
         percBA = species_sizes$percBA[match(Model, 
                                             species_sizes$species)],
         # color for species
         color_col = species_cols[Model],
         # species-specific vertical offset
         y_offset = case_when(
           Model == "rm" ~ 0.06,
           Model == "bc" ~ 0.03,
           Model == "sm" ~ 0,
           Model == "wa" ~ -0.03,
           Model == "ab" ~ 0.06,
           Model == "nro" ~ 0.03,
           Model == "bonro" ~ 0,
           Model == "wo" ~ -0.03,
           Model == "iw" ~ -0.06,
           TRUE ~ 0),
         # final y-position for plotting
         y_pos = as.numeric(component) + y_offset,
         # fill color: species color if significant, white otherwise
         fill_col = ifelse(Significant == "yes", 
                           species_cols[Model], "white"))

##### Filter to just include AM
fixef_df_AM <- filter(fixef_df_plot, Myc_type == "AM")

##### Plot - AM
ggplot(fixef_df_AM, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

##### Filter to just include ECM
fixef_df_ECM <- filter(fixef_df_plot, Myc_type == "ECM")

##### Plot - ECM
ggplot(fixef_df_ECM, aes(x = estimate, y = y_pos)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Error bars colored by species
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error,
                     color = Model), height = 0.05) +
  # Points: shape by Myc_type, fill by significance, outline by species
  geom_point(aes(shape = Myc_type, fill = fill_col,
                 color = Model, size = log1p(percBA)),
             stroke = 1) +
  # Y-axis with desired top-to-bottom order
  scale_y_continuous(breaks = 1:length(levels(fixef_df$component)),
                     labels = levels(fixef_df$component),
                     trans = "reverse"   # ensures top-to-bottom ordering visually
  ) +
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
  # Labels & theme
  labs(x = "Effect size ± SE (log count scale)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))





