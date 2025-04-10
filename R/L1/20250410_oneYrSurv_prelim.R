
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

## Most complex model to fit (at 5-m, 10-m, and 15-m scales): 
# glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conSp_BA + 
# conMycNS * hetMyc + log10(conSp_Prior3_counts + 1) + log10(hetSp_Prior3_counts + 1) + log10(total_count_yr0 + 1) + (1 | Site) + (1 | survYear), family = binomial(link = "logit"), data = data)

#################################################
#################################################
##### Read in file that has the counts by species for every TransectID-ObsYear combination
allSp_TI_ObsYear <- read.csv("./data/L1/allSp_prior_cohort_densities.csv")

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

#################### 5-m radius
# cor(rm_surv[, c("rm_5m", "conMycNS_5m", "hetMyc_5m",
#                 "rm_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")

### Most complex model
# Scaling the count covariates using log10() to account for difference in magnitude 
rm_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_5m + conMycNS_5m * hetMyc_5m + 
                      log10(rm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                      family = binomial(link = "logit"), 
                      data = rm_surv)
summary(rm_mod_5m)
# Variance inflation factors - all <2
car::vif(rm_mod_5m)
r.squaredGLMM(rm_mod_5m)

#################### 10-m radius
### Most complex model
# Scaling the count covariates using log10() to account for difference in magnitude 
rm_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     rm_10m + conMycNS_10m * hetMyc_10m + 
                     log10(rm_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = rm_surv)
summary(rm_mod_10m)
# Variance inflation factors - all <2
car::vif(rm_mod_10m)
r.squaredGLMM(rm_mod_10m)

#################### 15-m radius
### Most complex model - not converging properly with all terms
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

############## PICK UP HERE

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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
sm_surv$percEMF_5m = (sm_surv$hetMyc_5m / (sm_surv$hetMyc_5m + sm_surv$conMyc_5m))
# %EMF - 10m
sm_surv$percEMF_10m = (sm_surv$hetMyc_10m / (sm_surv$hetMyc_10m + sm_surv$conMyc_10m))
# %EMF - 15m
sm_surv$percEMF_15m = (sm_surv$hetMyc_15m / (sm_surv$hetMyc_15m + sm_surv$conMyc_15m))

#################### 5-m radius
# cor(sm_surv[, c("sm_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "sm_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model
sm_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     sm_5m * percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(sm_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = sm_surv)
summary(sm_mod_5m)
# Variance inflation factors - all < 4
car::vif(sm_mod_5m)

#################### 10-m radius
### Most complex model
sm_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      sm_10m * percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(sm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = sm_surv)
summary(sm_mod_10m)
# Variance inflation factors - all < 7
car::vif(sm_mod_10m)

#################### 15-m radius
### Most complex model - convergence issues
sm_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #sm_15m * percEMF_15m + 
                      sm_15m + percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(sm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = sm_surv)
summary(sm_mod_15m)
# Variance inflation factors - all < 8, but questionable
car::vif(sm_mod_15m)

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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
ab_surv$percEMF_5m = (ab_surv$hetMyc_5m / (ab_surv$hetMyc_5m + ab_surv$conMyc_5m))
# %EMF - 10m
ab_surv$percEMF_10m = (ab_surv$hetMyc_10m / (ab_surv$hetMyc_10m + ab_surv$conMyc_10m))
# %EMF - 15m
ab_surv$percEMF_15m = (ab_surv$hetMyc_15m / (ab_surv$hetMyc_15m + ab_surv$conMyc_15m))

#################### 5-m radius
# cor(ab_surv[, c("ab_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "ab_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model - cannot include all terms
ab_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     #ab_5m * percEMF_5m + 
                     ab_5m + percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(ab_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     #(1 | Site) + 
                     (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = ab_surv)
summary(ab_mod_5m)
# Variance inflation factors - all < 6
car::vif(ab_mod_5m)

#################### 10-m radius
### Most complex model - not working... 
ab_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #ab_10m * percEMF_10m + 
                      ab_10m + #percEMF_10m + 
                      conMycNS_10m + 
                      hetMyc_10m + 
                      log10(ab_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = ab_surv)
summary(ab_mod_10m)
# Variance inflation factors - some >10
car::vif(ab_mod_10m)

#################### 15-m radius
### Most complex model - not working... 
ab_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #ab_15m * percEMF_15m + 
                      ab_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(ab_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = ab_surv)
summary(ab_mod_15m)
# Variance inflation factors - some >10
car::vif(ab_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
wa_surv$percEMF_5m = (wa_surv$hetMyc_5m / (wa_surv$hetMyc_5m + wa_surv$conMyc_5m))
# %EMF - 10m
wa_surv$percEMF_10m = (wa_surv$hetMyc_10m / (wa_surv$hetMyc_10m + wa_surv$conMyc_10m))
# %EMF - 15m
wa_surv$percEMF_15m = (wa_surv$hetMyc_15m / (wa_surv$hetMyc_15m + wa_surv$conMyc_15m))

#################### 5-m radius
cor(wa_surv[, c("wa_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
                 "wa_Prior3_counts", "hetSp_Prior3_counts",
                "total_count_yr0")], use = "complete.obs")
### Most complex model 
wa_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     wa_5m * percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(wa_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = wa_surv)
summary(wa_mod_5m)
# Variance inflation factors - all < 4
car::vif(wa_mod_5m)

#################### 10-m radius
### Most complex model
wa_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      wa_10m * percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(wa_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wa_surv)
summary(wa_mod_10m)
# Variance inflation factors - all < 9 - questionable
car::vif(wa_mod_10m)

#################### 15-m radius
### Most complex model
wa_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #wa_15m * percEMF_15m + 
                      wa_15m + #percEMF_15m + 
                      conMycNS_15m * hetMyc_15m + 
                      log10(wa_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wa_surv)
summary(wa_mod_15m)
# Variance inflation factors - >10
car::vif(wa_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
bc_surv$percEMF_5m = (bc_surv$hetMyc_5m / (bc_surv$hetMyc_5m + bc_surv$conMyc_5m))
# %EMF - 10m
bc_surv$percEMF_10m = (bc_surv$hetMyc_10m / (bc_surv$hetMyc_10m + bc_surv$conMyc_10m))
# %EMF - 15m
bc_surv$percEMF_15m = (bc_surv$hetMyc_15m / (bc_surv$hetMyc_15m + bc_surv$conMyc_15m))

#################### 5-m radius
# cor(bc_surv[, c("bc_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "bc_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model 
bc_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     bc_5m * percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(bc_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = bc_surv)
summary(bc_mod_5m)
# Variance inflation factors - all < 4
car::vif(bc_mod_5m)

#################### 10-m radius
### Most complex model 
bc_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #bc_10m * percEMF_10m + 
                      bc_10m + #percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(bc_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bc_surv)
summary(bc_mod_10m)
# Variance inflation factors - all < 10 but questionable
car::vif(bc_mod_10m)

#################### 15-m radius
### Most complex model
bc_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #bc_15m * percEMF_15m + 
                      bc_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(bc_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bc_surv)
summary(bc_mod_15m)
# Variance inflation factors - all < 10, but questionable
car::vif(bc_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
nro_surv$percEMF_5m = (nro_surv$hetMyc_5m / (nro_surv$hetMyc_5m + nro_surv$conMyc_5m))
# %EMF - 10m
nro_surv$percEMF_10m = (nro_surv$hetMyc_10m / (nro_surv$hetMyc_10m + nro_surv$conMyc_10m))
# %EMF - 15m
nro_surv$percEMF_15m = (nro_surv$hetMyc_15m / (nro_surv$hetMyc_15m + nro_surv$conMyc_15m))

#################### 5-m radius
# cor(nro_surv[, c("nro_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "nro_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model 
nro_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     #nro_5m * percEMF_5m +
                      nro_5m + percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(nro_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = nro_surv)
summary(nro_mod_5m)
# Variance inflation factors - all < 4
car::vif(nro_mod_5m)

#################### 10-m radius
### Most complex model - not working... 
nro_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #nro_10m * percEMF_10m + 
                      nro_10m + percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(nro_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = nro_surv)
summary(nro_mod_10m)
# Variance inflation factors - all < 9, but questionable
car::vif(nro_mod_10m)

#################### 15-m radius
### Most complex model
nro_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #nro_15m * percEMF_15m + 
                      nro_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(nro_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = nro_surv)
summary(nro_mod_15m)
# Variance inflation factors - some >10
car::vif(nro_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
bonro_surv$percEMF_5m = (bonro_surv$hetMyc_5m / (bonro_surv$hetMyc_5m + bonro_surv$conMyc_5m))
# %EMF - 10m
bonro_surv$percEMF_10m = (bonro_surv$hetMyc_10m / (bonro_surv$hetMyc_10m + bonro_surv$conMyc_10m))
# %EMF - 15m
bonro_surv$percEMF_15m = (bonro_surv$hetMyc_15m / (bonro_surv$hetMyc_15m + bonro_surv$conMyc_15m))

#################### 5-m radius
# cor(bonro_surv[, c("bonro_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "bonro_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model
bonro_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     #bonro_5m * percEMF_5m + 
                     bonro_5m + percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(bonro_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = bonro_surv)
summary(bonro_mod_5m)
# Variance inflation factors - all < 3
car::vif(bonro_mod_5m)

#################### 10-m radius
### Most complex model - not working... 
bonro_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #bonro_10m * percEMF_10m + 
                      bonro_10m + #percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(bonro_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bonro_surv)
summary(bonro_mod_10m)
# Variance inflation factors - all < 9, but questionable
car::vif(bonro_mod_10m)

#################### 15-m radius
### Most complex model
bonro_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #bonro_15m * percEMF_15m + 
                      bonro_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(bonro_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bonro_surv)
summary(bonro_mod_15m)
# Variance inflation factors - all < 10, but questionbonrole
car::vif(bonro_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
wo_surv$percEMF_5m = (wo_surv$hetMyc_5m / (wo_surv$hetMyc_5m + wo_surv$conMyc_5m))
# %EMF - 10m
wo_surv$percEMF_10m = (wo_surv$hetMyc_10m / (wo_surv$hetMyc_10m + wo_surv$conMyc_10m))
# %EMF - 15m
wo_surv$percEMF_15m = (wo_surv$hetMyc_15m / (wo_surv$hetMyc_15m + wo_surv$conMyc_15m))

#################### 5-m radius
# cor(wo_surv[, c("wo_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "wo_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model 
wo_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     wo_5m * percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(wo_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = wo_surv)
summary(wo_mod_5m)
# Variance inflation factors - all < 3
car::vif(wo_mod_5m)

#################### 10-m radius
### Most complex model 
wo_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      wo_10m * percEMF_10m + 
                      #wo_10m + percEMF_10m + 
                      conMycNS_10m + #hetMyc_10m + 
                      log10(wo_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wo_surv)
summary(wo_mod_10m)
# Variance inflation factors - all < 3 if hetMyc excluded
car::vif(wo_mod_10m)

#################### 15-m radius
### Most complex model
wo_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      wo_15m * percEMF_15m + 
                      #wo_15m + percEMF_15m + 
                      conMycNS_15m + #hetMyc_15m + 
                      log10(wo_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = wo_surv)
summary(wo_mod_15m)
# Variance inflation factors - all < 4 if hetMyc excluded
car::vif(wo_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
iw_surv$percEMF_5m = (iw_surv$hetMyc_5m / (iw_surv$hetMyc_5m + iw_surv$conMyc_5m))
# %EMF - 10m
iw_surv$percEMF_10m = (iw_surv$hetMyc_10m / (iw_surv$hetMyc_10m + iw_surv$conMyc_10m))
# %EMF - 15m
iw_surv$percEMF_15m = (iw_surv$hetMyc_15m / (iw_surv$hetMyc_15m + iw_surv$conMyc_15m))

#################### 5-m radius
# cor(iw_surv[, c("iw_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
#                  "iw_Prior3_counts", "hetSp_Prior3_counts",
#                 "total_count_yr0")], use = "complete.obs")
### Most complex model - cannot include all terms
iw_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     #iw_5m * percEMF_5m + 
                     iw_5m + percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(iw_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = iw_surv)
summary(iw_mod_5m)
# Variance inflation factors - all < 4
car::vif(iw_mod_5m)

#################### 10-m radius
### Most complex model - not working... 
iw_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #iw_10m * percEMF_10m + 
                      iw_10m + #percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(iw_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = iw_surv)
summary(iw_mod_10m)
# Variance inflation factors - all < 2 in simpler model
car::vif(iw_mod_10m)

#################### 15-m radius
### Most complex model
iw_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      #iw_15m * percEMF_15m + 
                      iw_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(iw_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = iw_surv)
summary(iw_mod_15m)
# Variance inflation factors - all < 2 in simpler model
car::vif(iw_mod_15m)


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

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
bs_surv$percEMF_5m = (bs_surv$hetMyc_5m / (bs_surv$hetMyc_5m + bs_surv$conMyc_5m))
# %EMF - 10m
bs_surv$percEMF_10m = (bs_surv$hetMyc_10m / (bs_surv$hetMyc_10m + bs_surv$conMyc_10m))
# %EMF - 15m
bs_surv$percEMF_15m = (bs_surv$hetMyc_15m / (bs_surv$hetMyc_15m + bs_surv$conMyc_15m))

#################### 5-m radius
cor(bs_surv[, c("bs_5m", "percEMF_5m", "conMycNS_5m", "hetMyc_5m",
                 "bs_Prior3_counts", "hetSp_Prior3_counts",
                "total_count_yr0")], use = "complete.obs")
### Most complex model - cannot include all terms
bs_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     #bs_5m * percEMF_5m + 
                     bs_5m + percEMF_5m + 
                     conMycNS_5m + hetMyc_5m + 
                     log10(bs_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     #(1 | Site) + 
                     (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = bs_surv)
summary(bs_mod_5m)
# Variance inflation factors - all < 4
car::vif(bs_mod_5m)

#################### 10-m radius
### Most complex model - not working... 
bs_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      #bs_10m * percEMF_10m + 
                      bs_10m + #percEMF_10m + 
                      conMycNS_10m + hetMyc_10m + 
                      log10(bs_Prior3_counts + 1) +
                      log10(hetSp_Prior3_counts + 1) +
                      log10(total_count_yr0 + 1) +
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bs_surv)
summary(bs_mod_10m)
# Variance inflation factors - all < 3
car::vif(bs_mod_10m)

#################### 15-m radius
### Most complex model
bs_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                     # bs_15m * percEMF_15m + 
                      bs_15m + #percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(bs_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      #(1 | Site) + 
                      (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = bs_surv)
summary(bs_mod_15m)
# Variance inflation factors - all < 3 with simplified structure
car::vif(bs_mod_15m)



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



