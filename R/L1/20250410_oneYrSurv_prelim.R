
#################################################
####### Analyses of one-year survivorship of Red Maple WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol
####### Date: 04 April 2025
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

# For species with sufficient observations per site (n ≥ XX), we included a random effect of site to account for spatial clustering. For species with limited replication across sites, we excluded this term due to singular fits and lack of variance attributable to site.

### Summary by species
## RM - 15-m
## SM - 10-m
## AB - 15-m (site excluded)
## WA - 15-m
## BC - 15-m 
## NRO - 15-m
## BO - 15-m 
## WO - 10-m
## IW - 10-m
## BS - 15-m  (site excluded)
#### Good justification to go with the 15-m radius

#################################################
#################################################
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

###############################
###############################

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


##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods
##### Averaged across sites
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
    # Choose any one site and year for valid factor levels
    Site = rm_surv$Site[1],survYear = rm_surv$survYear[1],
    focal_label = label)
  df
}

# Step 4: Combine prediction data for all three components
new_data <- bind_rows(
  make_pred_data("rm_15m", "Conspecific BA"),
  make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
  make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA")
)

# Step 5: Predict survivorship without random effects (average effect)
new_data$pred_surv <- predict(rm_mod_15m, newdata = new_data, 
                      type = "response", re.form = NA)
# re.form = NA ignores the RE and uses only the fixed effects (i.e., an “average” prediction across sites and years).

# Step 6: Add x-axis for the focal BA value
new_data$x_focal <- with(new_data, 
              ifelse(focal_label == "Conspecific BA", rm_15m,
              ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))

# Step 7: Plot
ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
  geom_line(linewidth = 1) +
  labs(x = "Basal Area in 15-m Radius (m²)", 
       y = "Predicted Survivorship", color = "") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")


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

###############################
###############################

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


#################################################
#################################################
##### American Beech

### Read in data file from the R/L1 folder 
ab_surv_temp <- read.csv("./data/L1/ab_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
ab_surv_temp <- dplyr::filter(ab_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
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

###############################
###############################

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

###############################
###############################

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

###############################
###############################

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

###############################
###############################

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


#################################################
#################################################
##### Black oak

### Read in data file from the R/L1 folder 
bonro_surv_temp <- read.csv("./data/L1/bonro_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
bonro_surv_temp <- dplyr::filter(bonro_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
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

###############################
###############################

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

###############################
###############################

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

###############################
###############################

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


#################################################
#################################################
##### American Basswood

### Read in data file from the R/L1 folder 
bs_surv_temp <- read.csv("./data/L1/bs_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
bs_surv_temp <- dplyr::filter(bs_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data, using a concatenate and match 
# Survivorship data
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

###############################
###############################

#################### 15-m radius
bs_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      bs_15m + conMycNS_15m * hetMyc_15m + 
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

###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### 
###### Plot survivorship versus TOTAL conmycorrhizal BA with fit lines from:
# 1) cm_5m ; 2) cmNs_5m; 3) cs_5m
# Extract model coefficients (logit scale)
models <- data.frame(
  Model = c("Con-Mycorrhizal - RM", "Con-Mycorrhizal - no RM", "RM"),
  Intercept = c(fixef(cm_5m)[1], fixef(cmNs_5m)[1], fixef(cs_5m)[1]),
  Slope = c(fixef(cm_5m)[2], fixef(cmNs_5m)[2], fixef(cs_5m)[2])
)

# Generate predictor values (common range for x-axis)
x_values <- seq(min(rm_surv$conMyc_5m), max(rm_surv$conMyc_5m), length.out = 100)

# Compute back-transformed survival probabilities
fit_data <- do.call(rbind, lapply(1:nrow(models), function(i) {
  logit_pred <- models$Intercept[i] + models$Slope[i] * x_values
  prob_pred <- plogis(logit_pred)  # Transform logit to probability
  data.frame(conMyc_5m = x_values, survivorship = prob_pred, Model = models$Model[i])
}))

## Re-order levels
fit_data$Model <- factor(fit_data$Model, 
                levels = c("Con-Mycorrhizal - RM", "Con-Mycorrhizal - no RM", "RM"))

# Scatterplot with fitted lines from all models
ggplot(rm_surv, aes(x = conMyc_5m, y = total_count_yr1 / total_count_yr0)) +
  geom_point(alpha = 0.5) +  # Observed survivorship
  geom_line(data = fit_data, aes(x = conMyc_5m, y = survivorship, color = Model), linewidth = 1) +
  scale_color_manual(values = c("blue", "red", "gold")) +  # Custom colors
  labs(x = "Con-Mycorrhizal Basal Area in 5-m Neighborhood", 
       y = "Seedling Survivorship", 
       title = "RM Seedling Survivorship") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 12), legend.key = element_blank(), 
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods
##### Panel for each site!

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

