
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 16 April 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)
library(ggplot2)

## Most complex model to fit (15-m radius)
# glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conSp_BA + 
# conMycNS * hetMyc + log10(conSp_Prior3_counts + 1) + log10(hetSp_Prior3_counts + 1) + log10(total_count_yr0 + 1) + (1 | Site) + (1 | survYear), family = binomial(link = "logit"), data = data)

# For species with sufficient observations, we included a random effect of site. For species with limited observations across sites, we excluded this term due to singular fits and lack of variance attributable to site.

##### Read in file that has the counts by species for every TransectID-ObsYear combination
allSp_TI_ObsYear <- read.csv("./data/L1/allSp_prior_cohort_densities.csv")

#################################################
#################################################
#################################################
##### Red Maple

### Read in data file from the R/L1 folder 
rm_surv_temp <- read.csv("./data/L1/rm_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
rm_surv_temp <- dplyr::filter(rm_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
rm_surv_temp$TI_FY = paste(rm_surv_temp$TransectID,rm_surv_temp$First.Year)
# All species 3-year count data
allSp_TI_ObsYear$TI_FY = paste(allSp_TI_ObsYear$TransectID,allSp_TI_ObsYear$Obs.Year)
# Calculate sum of heterospecific species previous 3-year counts for RM
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,6:14], na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
rm_surv_temp$rm_Prior3_counts = 
              allSp_TI_ObsYear$rm_Prior3_counts[match(rm_surv_temp$TI_FY,
                                                      allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
rm_surv_temp$rm_Prior3_counts[is.na(rm_surv_temp$rm_Prior3_counts)]<- 0
# Heterospecific
rm_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(rm_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

rm_surv <- rm_surv_temp

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#################### 15-m radius
rm_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     rm_15m + conMycNS_15m * hetMyc_15m + 
                     log10(rm_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = rm_surv)
summary(rm_mod_15m)
# Variance inflation factors - all < 3
car::vif(rm_mod_15m)
r.squaredGLMM(rm_mod_15m) 
# Fixed effects: R2m = 0.05563300 
# Fixed + random effects: R2c = 0.5629413

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(rm_surv$rm_15m, rm_surv$conMycNS_15m, 
              rm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- rm_surv %>%
  summarize(across(c(rm_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
          rm_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    rm_15m = if (focal == "rm_15m") ba_seq else median_vals$rm_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    rm_Prior3_counts = median_vals$rm_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = rm_surv$Site[1],survYear = rm_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("rm_15m", "Conspecific BA"),
  make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
  make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(rm_mod_15m, newdata = new_data, 
                      type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
              ifelse(focal_label == "Conspecific BA", rm_15m,
              ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
    y = "Predicted Survivorship", 
    color = "",title = "Red Maple") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
    axis.title = element_text(size = 12), legend.key = element_blank(), 
    plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 10), 
    panel.grid.minor = element_blank(),
    legend.position = c(0.05, 0.05),  # Bottom left inside plot
    legend.justification = c(0, 0))



#################################################
#################################################
##### Sugar Maple

### Read in data file from the R/L1 folder 
sm_surv_temp <- read.csv("./data/L1/sm_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
sm_surv_temp <- dplyr::filter(sm_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
sm_surv_temp$TI_FY = paste(sm_surv_temp$TransectID,sm_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for sm
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5,7:14)], na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
sm_surv_temp$sm_Prior3_counts = 
  allSp_TI_ObsYear$sm_Prior3_counts[match(sm_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
sm_surv_temp$sm_Prior3_counts[is.na(sm_surv_temp$sm_Prior3_counts)]<- 0
# Heterospecific
sm_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(sm_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

sm_surv <- sm_surv_temp

#### Make site and year into factors
sm_surv$Site = as.factor(sm_surv$Site)
sm_surv$survYear = as.factor(sm_surv$survYear)

#################### 15-m radius
sm_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      sm_15m + conMycNS_15m * hetMyc_15m + 
                      log10(sm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = sm_surv)
summary(sm_mod_15m)
# Variance inflation factors - all < 3
car::vif(sm_mod_15m)
r.squaredGLMM(sm_mod_15m)
# Fixed effects: R2m = 0.03176506 
# Fixed + random effects: R2c = 0.8307705

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(sm_surv$sm_15m, sm_surv$conMycNS_15m, 
                     sm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- sm_surv %>%
  summarize(across(c(sm_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     sm_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    sm_15m = if (focal == "sm_15m") ba_seq else median_vals$sm_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    sm_Prior3_counts = median_vals$sm_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = sm_surv$Site[1],survYear = sm_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("sm_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(sm_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", sm_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Sugar Maple") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### American Beech

### Read in data file from the R/L1 folder 
ab_surv_temp <- read.csv("./data/L1/ab_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
ab_surv_temp <- dplyr::filter(ab_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
ab_surv_temp$TI_FY = paste(ab_surv_temp$TransectID,ab_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for ab
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5,6,8:14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
ab_surv_temp$ab_Prior3_counts = 
  allSp_TI_ObsYear$ab_Prior3_counts[match(ab_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
ab_surv_temp$ab_Prior3_counts[is.na(ab_surv_temp$ab_Prior3_counts)]<- 0
# Heterospecific
ab_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(ab_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

ab_surv <- ab_surv_temp

#### Make site and year into factors
ab_surv$Site = as.factor(ab_surv$Site)
ab_surv$survYear = as.factor(ab_surv$survYear)


#################### 15-m radius - no RE of site
### Most complex model - not converging properly with all teabs
ab_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      ab_15m + conMycNS_15m * hetMyc_15m + 
                      log10(ab_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = ab_surv)
summary(ab_mod_15m)
# Variance inflation factors - <10
car::vif(ab_mod_15m)
r.squaredGLMM(ab_mod_15m)
# Fixed effects: R2m = 0.03955228 
# Fixed + random effects (survYear only): R2c = 0.4981347

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(ab_surv$ab_15m, ab_surv$conMycNS_15m, 
                     ab_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- ab_surv %>%
  summarize(across(c(ab_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     ab_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    ab_15m = if (focal == "ab_15m") ba_seq else median_vals$ab_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    ab_Prior3_counts = median_vals$ab_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    survYear = ab_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("ab_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(ab_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects 

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", ab_15m,
               ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "American Beech") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### White ash

### Read in data file from the R/L1 folder 
wa_surv_temp <- read.csv("./data/L1/wa_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
wa_surv_temp <- dplyr::filter(wa_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
wa_surv_temp$TI_FY = paste(wa_surv_temp$TransectID,wa_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for wa
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:7,9:14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
wa_surv_temp$wa_Prior3_counts = 
  allSp_TI_ObsYear$wa_Prior3_counts[match(wa_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
wa_surv_temp$wa_Prior3_counts[is.na(wa_surv_temp$wa_Prior3_counts)]<- 0
# Heterospecific
wa_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(wa_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

wa_surv <- wa_surv_temp

#### Make site and year into factors
wa_surv$Site = as.factor(wa_surv$Site)
wa_surv$survYear = as.factor(wa_surv$survYear)

#################### 15-m radius
wa_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      wa_15m + conMycNS_15m * hetMyc_15m + 
                      log10(wa_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wa_surv)
summary(wa_mod_15m) # we may need to account for the mortality more closely
# Do an individual set of analyses from 2016 to present
# Variance inflation factors - all < 5
car::vif(wa_mod_15m)
r.squaredGLMM(wa_mod_15m)
# Fixed effects: R2m = 0.02632910 
# Fixed + random effects: R2c = 0.8731862

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(wa_surv$wa_15m, wa_surv$conMycNS_15m, 
                     wa_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- wa_surv %>%
  summarize(across(c(wa_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     wa_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    wa_15m = if (focal == "wa_15m") ba_seq else median_vals$wa_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    wa_Prior3_counts = median_vals$wa_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = wa_surv$Site[1],survYear = wa_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("wa_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(wa_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", wa_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "White Ash") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.05),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### Black Cherry

### Read in data file from the R/L1 folder 
bc_surv_temp <- read.csv("./data/L1/bc_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
bc_surv_temp <- dplyr::filter(bc_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
bc_surv_temp$TI_FY = paste(bc_surv_temp$TransectID,bc_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for bc
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:8,10:14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
bc_surv_temp$bc_Prior3_counts = 
  allSp_TI_ObsYear$bc_Prior3_counts[match(bc_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
bc_surv_temp$bc_Prior3_counts[is.na(bc_surv_temp$bc_Prior3_counts)]<- 0
# Heterospecific
bc_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(bc_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

bc_surv <- bc_surv_temp

#### Make site and year into factors
bc_surv$Site = as.factor(bc_surv$Site)
bc_surv$survYear = as.factor(bc_surv$survYear)

#################### 15-m radius
bc_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      bc_15m + conMycNS_15m * hetMyc_15m + 
                      log10(bc_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bc_surv)
summary(bc_mod_15m)
# Variance inflation factors - all < 5
car::vif(bc_mod_15m)
r.squaredGLMM(bc_mod_15m)
# Fixed effects: R2m = 0.02266342 
# Fixed + random effects: R2c = 0.4276095

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(bc_surv$bc_15m, bc_surv$conMycNS_15m, 
                     bc_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- bc_surv %>%
  summarize(across(c(bc_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     bc_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    bc_15m = if (focal == "bc_15m") ba_seq else median_vals$bc_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    bc_Prior3_counts = median_vals$bc_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = bc_surv$Site[1],survYear = bc_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("bc_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(bc_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", bc_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Black Cherry") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.05),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### Northern Red Oak

### Read in data file from the R/L1 folder 
nro_surv_temp <- read.csv("./data/L1/nro_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
nro_surv_temp <- dplyr::filter(nro_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
nro_surv_temp$TI_FY = paste(nro_surv_temp$TransectID,nro_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for nro
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:9,11:14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
nro_surv_temp$nro_Prior3_counts = 
  allSp_TI_ObsYear$nro_Prior3_counts[match(nro_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
nro_surv_temp$nro_Prior3_counts[is.na(nro_surv_temp$nro_Prior3_counts)]<- 0
# Heterospecific
nro_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(nro_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

nro_surv <- nro_surv_temp

#### Make site and year into factors
nro_surv$Site = as.factor(nro_surv$Site)
nro_surv$survYear = as.factor(nro_surv$survYear)

#################### 15-m radius
nro_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      nro_15m + conMycNS_15m * hetMyc_15m + 
                      log10(nro_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = nro_surv)
summary(nro_mod_15m)
# Variance inflation factors - all < 7
car::vif(nro_mod_15m)
r.squaredGLMM(nro_mod_15m)
# Fixed effects: R2m = 0.02860598 
# Fixed + random effects: R2c = 0.4900101

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(nro_surv$nro_15m, nro_surv$conMycNS_15m, 
                     nro_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- nro_surv %>%
  summarize(across(c(nro_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     nro_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    nro_15m = if (focal == "nro_15m") ba_seq else median_vals$nro_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    nro_Prior3_counts = median_vals$nro_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = nro_surv$Site[1],survYear = nro_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("nro_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(nro_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", nro_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Northern Red Oak") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### Black oak

### Read in data file from the R/L1 folder 
bonro_surv_temp <- read.csv("./data/L1/bonro_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
bonro_surv_temp <- dplyr::filter(bonro_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
bonro_surv_temp$TI_FY = paste(bonro_surv_temp$TransectID,bonro_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for bonro
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:10,12:14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
bonro_surv_temp$bonro_Prior3_counts = 
  allSp_TI_ObsYear$bonro_Prior3_counts[match(bonro_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
bonro_surv_temp$bonro_Prior3_counts[is.na(bonro_surv_temp$bonro_Prior3_counts)]<- 0
# Heterospecific
bonro_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(bonro_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

bonro_surv <- bonro_surv_temp

#### Make site and year into factors
bonro_surv$Site = as.factor(bonro_surv$Site)
bonro_surv$survYear = as.factor(bonro_surv$survYear)

#################### 15-m radius
bonro_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      bonro_15m + conMycNS_15m * hetMyc_15m + 
                      log10(bonro_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bonro_surv)
summary(bonro_mod_15m)
# Variance inflation factors - all < 7
car::vif(bonro_mod_15m)
r.squaredGLMM(bonro_mod_15m)
# Fixed effects: R2m = 0.03257658 
# Fixed + random effects: R2c = 0.4265718

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(bonro_surv$bonro_15m, bonro_surv$conMycNS_15m, 
                     bonro_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- bonro_surv %>%
  summarize(across(c(bonro_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     bonro_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    bonro_15m = if (focal == "bonro_15m") ba_seq else median_vals$bonro_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    bonro_Prior3_counts = median_vals$bonro_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = bonro_surv$Site[1],survYear = bonro_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("bonro_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(bonro_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", bonro_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Black Oak") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### White oak

### Read in data file from the R/L1 folder 
wo_surv_temp <- read.csv("./data/L1/wo_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
wo_surv_temp <- dplyr::filter(wo_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
wo_surv_temp$TI_FY = paste(wo_surv_temp$TransectID,wo_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for wo
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:11,13,14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
wo_surv_temp$wo_Prior3_counts = 
  allSp_TI_ObsYear$wo_Prior3_counts[match(wo_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
wo_surv_temp$wo_Prior3_counts[is.na(wo_surv_temp$wo_Prior3_counts)]<- 0
# Heterospecific
wo_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(wo_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

wo_surv <- wo_surv_temp

#### Make site and year into factors
wo_surv$Site = as.factor(wo_surv$Site)
wo_surv$survYear = as.factor(wo_surv$survYear)

#################### 15-m radius
wo_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      wo_15m + conMycNS_15m * hetMyc_15m + 
                      log10(wo_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wo_surv)
summary(wo_mod_15m)
# Variance inflation factors - all < 7
car::vif(wo_mod_15m)
r.squaredGLMM(wo_mod_15m)
# Fixed effects: R2m = 0.09389812 
# Fixed + random effects: R2c = 0.4235249

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(wo_surv$wo_15m, wo_surv$conMycNS_15m, 
                     wo_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- wo_surv %>%
  summarize(across(c(wo_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     wo_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    wo_15m = if (focal == "wo_15m") ba_seq else median_vals$wo_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    wo_Prior3_counts = median_vals$wo_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = wo_surv$Site[1],survYear = wo_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("wo_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(wo_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", wo_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "White Oak") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.05),  # Bottom left inside plot
        legend.justification = c(0, 0))


#################################################
#################################################
##### Ironwood

### Read in data file from the R/L1 folder 
iw_surv_temp <- read.csv("./data/L1/iw_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
iw_surv_temp <- dplyr::filter(iw_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
iw_surv_temp$TI_FY = paste(iw_surv_temp$TransectID,iw_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for iw
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:12,14)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
iw_surv_temp$iw_Prior3_counts = 
  allSp_TI_ObsYear$iw_Prior3_counts[match(iw_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
iw_surv_temp$iw_Prior3_counts[is.na(iw_surv_temp$iw_Prior3_counts)]<- 0
# Heterospecific
iw_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(iw_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

iw_surv <- iw_surv_temp

#### Make site and year into factors
iw_surv$Site = as.factor(iw_surv$Site)
iw_surv$survYear = as.factor(iw_surv$survYear)

#################### 15-m radius
iw_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      iw_15m + conMycNS_15m * hetMyc_15m + 
                      log10(iw_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = iw_surv)
summary(iw_mod_15m)
# Variance inflation factors - all < 5
car::vif(iw_mod_15m)
r.squaredGLMM(iw_mod_15m)
# Fixed effects: R2m = 0.01951896 
# Fixed + random effects: R2c = 0.2874703

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(iw_surv$iw_15m, iw_surv$conMycNS_15m, 
                     iw_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- iw_surv %>%
  summarize(across(c(iw_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     iw_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal variable
make_pred_data <- function(focal, label) {
  df <- data.frame(
    iw_15m = if (focal == "iw_15m") ba_seq else median_vals$iw_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    iw_Prior3_counts = median_vals$iw_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    Site = iw_surv$Site[1],survYear = iw_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("iw_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(iw_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", iw_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Ironwood") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))



#################################################
#################################################
##### American Basswood

### Read in data file from the R/L1 folder 
bs_surv_temp <- read.csv("./data/L1/bs_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
bs_surv_temp <- dplyr::filter(bs_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
bs_surv_temp$TI_FY = paste(bs_surv_temp$TransectID,bs_surv_temp$First.Year)

# Calculate sum of heterospecific species previous 3-year counts for bs
allSp_TI_ObsYear$hetSp_Prior3_counts = rowSums(allSp_TI_ObsYear[,c(5:13)],
                                               na.rm = TRUE)

## Add conspecific and heterospecific species previous 3-year counts columns with match
# Conspecific
bs_surv_temp$bs_Prior3_counts = 
  allSp_TI_ObsYear$bs_Prior3_counts[match(bs_surv_temp$TI_FY,
                                          allSp_TI_ObsYear$TI_FY)]
# Convert NAs to 0s
bs_surv_temp$bs_Prior3_counts[is.na(bs_surv_temp$bs_Prior3_counts)]<- 0
# Heterospecific
bs_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(bs_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

bs_surv <- bs_surv_temp

#### Make site and year into factors
bs_surv$Site = as.factor(bs_surv$Site)
bs_surv$survYear = as.factor(bs_surv$survYear)

#################### 15-m radius
bs_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      bs_15m + conMycNS_15m + hetMyc_15m + 
                      log10(bs_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bs_surv)
summary(bs_mod_15m)
# Variance inflation factors - all < 3
car::vif(bs_mod_15m)
r.squaredGLMM(bs_mod_15m)
# Fixed effects: R2m = 0.01815265 
# Fixed + random effects (survYear only): R2c = 0.6103411

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across years

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(bs_surv$bs_15m, bs_surv$conMycNS_15m, 
                     bs_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- bs_surv %>%
  summarize(across(c(bs_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     bs_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Create prediction datasets for each focal varibsle
make_pred_data <- function(focal, label) {
  df <- data.frame(
    bs_15m = if (focal == "bs_15m") ba_seq else median_vals$bs_15m,
    conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
    hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
    bs_Prior3_counts = median_vals$bs_Prior3_counts,
    hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
    total_count_yr0 = median_vals$total_count_yr0,
    survYear = bs_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(make_pred_data("bs_15m", "Conspecific BA"),
                      make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
                      make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(bs_mod_15m, newdata = new_data, 
                              type = "response", re.form = NA)
# re.form = NA ignores RE, uses only the fixed effects 

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", bs_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", 
       color = "",title = "Basswood") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.75),  # Bottom left inside plot
        legend.justification = c(0, 0))

###############################
###############################
###############################
###############################

# ## Model selection/evaluation using AIC
# # Create a data frame with model names and AIC values
# aic_rm_5m <- tibble(Model = c("cs_5m", "cm_5m", "cmNs_5m", "hm_5m", 
#                               "hmcm_5m", "hmcmNS_5m"),
#                     AIC = sapply(list(cs_5m, cm_5m, cmNs_5m, hm_5m, 
#                                       hmcm_5m, hmcmNS_5m), AIC)) %>%
#   arrange(AIC)  # Sort by AIC
# # Print the sorted AIC values
# print(aic_rm_5m)
# # Pseudo R2 for best models
# r.squaredGLMM(hmcm_5m)
# r.squaredGLMM(cm_5m)



##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods
##### Panel for each site!
# All this really does is change the intercept by site - not super informative

# Step 1: Create a sequence of basal area values
ba_seq <- seq(0, max(rm_surv$rm_15m, rm_surv$conMycNS_15m, 
                     rm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)

# Step 2: Set median values for the non-focal covariates
median_vals <- rm_surv %>%
  summarize(across(c(rm_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
                     rm_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))

# Step 3: Get all site levels and choose a valid survYear
site_levels <- unique(rm_surv$Site)
surv_year_level <- rm_surv$survYear[1]  # just to ensure model structure is respected

# Step 4: Function to generate prediction data across sites
make_pred_data <- function(focal, label) {
  expand.grid(
    x_val = ba_seq,
    Site = site_levels
  ) %>%
    mutate(
      rm_15m = if (focal == "rm_15m") x_val else median_vals$rm_15m,
      conMycNS_15m = if (focal == "conMycNS_15m") x_val else median_vals$conMycNS_15m,
      hetMyc_15m = if (focal == "hetMyc_15m") x_val else median_vals$hetMyc_15m,
      rm_Prior3_counts = median_vals$rm_Prior3_counts,
      hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
      total_count_yr0 = median_vals$total_count_yr0,
      survYear = surv_year_level,
      focal_label = label
    )
}

# Step 5: Combine prediction data for all three components
new_data <- bind_rows(
  make_pred_data("rm_15m", "Conspecific BA"),
  make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
  make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA")
)

# Step 6: Predict survivorship (average effect, no REs)
new_data$pred_surv <- predict(rm_mod_15m, newdata = new_data, 
                              type = "response", re.form = ~(1 | Site))

# Step 7: Add focal x-axis
new_data$x_focal <- with(new_data, 
                         ifelse(focal_label == "Conspecific BA", rm_15m,
                                ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Optional: Make Site a factor for consistent facet ordering
new_data$Site <- factor(new_data$Site, levels = sort(unique(new_data$Site)))

# Step 8: Plot with facets per site
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Site, ncol = 4) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", color = "") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = 'black', size = 10), 
    axis.title = element_text(size = 12), 
    legend.key = element_blank(), 
    plot.title = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

