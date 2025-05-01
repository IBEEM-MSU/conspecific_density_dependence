#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 1 May 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)
library(ggplot2)
library(DHARMa)

### Most complex model to fit (15-m radius)
###  Biotic model of one-year survivorship - with adult densities (BA)
# (#1 yr) ~ Consp BA:log(#0 yr) + Con-myc BA:log(#0 yr) + Het-myc BA:log(#0 yr) + log(hetSp_Current_Prior3 + 1):log(#0 yr) + Consp BA:log(#0 yr):Con-myc BA + Consp BA:log(#0 yr):Het-myc BA + Het-myc BA:log(#0 yr):Con-myc BA + (1 | Site) + (1 | Year)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
# (#1 yr) ~ ConspRelativeAbund:log(#0 yr) + %ECM:log(#0 yr) + log(hetSp_Current_Prior3 + 1):log(#0 yr) + ConspRelativeAbund:log(#0 yr):%ECM + (1 | Site) + (1 | Year)

# For species with sufficient observations, we included a random effect of site. For species with limited observations across sites, we excluded this term due to singular fits and lack of variance attributable to site.

##### Read in file that has the counts by species for every TransectID-ObsYear combination
allSp_TI_ObsYear <- read.csv("./data/L1/allSp_prior_cohort_densities.csv")

## Create TI_FY concatenated variable
allSp_TI_ObsYear$TI_FY = paste(allSp_TI_ObsYear$TransectID,
                               allSp_TI_ObsYear$Obs.Year)

##### Going to bring in the total_count_yr0 counts for each species, matching by Obs.Year (First.Year) and Transect ID
###### rm
## Read in data file from the R/L1 folder 
rm_surv_temp <- read.csv("./data/L1/rm_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
rm_surv_temp$TI_FY = paste(rm_surv_temp$TransectID,rm_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$rm_0yr_count = 
  rm_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     rm_surv_temp$TI_FY)]
###### sm
## Read in data file from the R/L1 folder 
sm_surv_temp <- read.csv("./data/L1/sm_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
sm_surv_temp$TI_FY = paste(sm_surv_temp$TransectID,sm_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$sm_0yr_count = 
  sm_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     sm_surv_temp$TI_FY)]
###### ab
## Read in data file from the R/L1 folder 
ab_surv_temp <- read.csv("./data/L1/ab_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
ab_surv_temp$TI_FY = paste(ab_surv_temp$TransectID,ab_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$ab_0yr_count = 
  ab_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     ab_surv_temp$TI_FY)]
###### wa
## Read in data file from the R/L1 folder 
wa_surv_temp <- read.csv("./data/L1/wa_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
wa_surv_temp$TI_FY = paste(wa_surv_temp$TransectID,wa_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$wa_0yr_count = 
  wa_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     wa_surv_temp$TI_FY)]
###### bc
## Read in data file from the R/L1 folder 
bc_surv_temp <- read.csv("./data/L1/bc_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
bc_surv_temp$TI_FY = paste(bc_surv_temp$TransectID,bc_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$bc_0yr_count = 
  bc_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     bc_surv_temp$TI_FY)]
###### nro
## Read in data file from the R/L1 folder 
nro_surv_temp <- read.csv("./data/L1/nro_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
nro_surv_temp$TI_FY = paste(nro_surv_temp$TransectID,nro_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$nro_0yr_count = 
  nro_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                      nro_surv_temp$TI_FY)]
###### bonro
## Read in data file from the R/L1 folder 
bonro_surv_temp <- read.csv("./data/L1/bonro_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
bonro_surv_temp$TI_FY = paste(bonro_surv_temp$TransectID,bonro_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$bonro_0yr_count = 
  bonro_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                        bonro_surv_temp$TI_FY)]
###### wo
## Read in data file from the R/L1 folder 
wo_surv_temp <- read.csv("./data/L1/wo_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
wo_surv_temp$TI_FY = paste(wo_surv_temp$TransectID,wo_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$wo_0yr_count = 
  wo_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     wo_surv_temp$TI_FY)]
###### iw
## Read in data file from the R/L1 folder 
iw_surv_temp <- read.csv("./data/L1/iw_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
iw_surv_temp$TI_FY = paste(iw_surv_temp$TransectID,iw_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$iw_0yr_count = 
  iw_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     iw_surv_temp$TI_FY)]
###### bs
## Read in data file from the R/L1 folder 
bs_surv_temp <- read.csv("./data/L1/bs_oneYrSurv_matureBA.csv", row.names = 1)
# Create concatenated Transect ID - First Year identifier
bs_surv_temp$TI_FY = paste(bs_surv_temp$TransectID,bs_surv_temp$First.Year)
# Bring over 0-yr counts from species' data using match
allSp_TI_ObsYear$bs_0yr_count = 
  bs_surv_temp$total_count_yr0[match(allSp_TI_ObsYear$TI_FY,
                                     bs_surv_temp$TI_FY)]

#################################################
#################################################
#################################################
##### Red Maple

### Filter to include 2001 as first First.Year
rm_surv_temp <- dplyr::filter(rm_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,17:25)], na.rm = TRUE)
# Heterospecific
rm_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(rm_surv_temp$TI_FY,allSp_TI_ObsYear$TI_FY)]

# Calculate sum of heterospecific species previous 3-year and current year counts for RM
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(6:14,17:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
rm_surv$propEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
rm_surv$propRM_15m = (rm_surv$rm_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))
rm_surv$percEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Log transform the count of 0-yr old seedlings
rm_surv$log0yr <- log(rm_surv$total_count_yr0)  

#################### 15-m radius
###  Biotic model of one-year survivorship - with adult densities (BA)
# (#1 yr) ~ Consp BA:log(#0 yr) + Con-myc BA:log(#0 yr) + Het-myc BA:log(#0 yr) + log(Prior3_counts + 1):log(#0 yr) + Consp BA:log(#0 yr):Con-myc BA + Consp BA:log(#0 yr):Het-myc BA + Het-myc BA:log(#0 yr):Con-myc BA + (1 | Site) + (1 | Year)

rmBA_15m <- glmer(total_count_yr1 ~ rm_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    rm_15m:log0yr:conMycNS_15m + 
                    rm_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = rm_surv)
summary(rmBA_15m)

## Simplified to exclude 3-way interactions
rmBA_15m2 <- glmer(total_count_yr1 ~ rm_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = rm_surv)
summary(rmBA_15m2)

## VIFs and pseudo r-squared values
car::vif(rmBA_15m) # all < 10
r.squaredGLMM(rmBA_15m) 
car::vif(rmBA_15m2) # all < 2
r.squaredGLMM(rmBA_15m2) 
## AIC
aic_values <- AIC(rmBA_15m, rmBA_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # rm_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = rmBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
# (#1 yr) ~ ConspRelativeAbund:log(#0 yr) + %ECM:log(#0 yr) + log(Prior3_counts + 1):log(#0 yr) + ConspRelativeAbund:log(#0 yr):%ECM + (1 | Site) + (1 | Year)

rmProp_15m <- glmer(total_count_yr1 ~ propRM_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propRM_15m:log0yr:propEMF_15m +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rmProp_15m)

## Simplified to exclude 3-way interactions
rmProp_15m2 <- glmer(total_count_yr1 ~ propRM_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = rm_surv)
summary(rmProp_15m2)

## VIFs and pseudo r-squared values
car::vif(rmProp_15m) 
r.squaredGLMM(rmProp_15m) 
car::vif(rmProp_15m2) 
r.squaredGLMM(rmProp_15m2) 
## AIC
aic_values <- AIC(rmProp_15m, rmProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # rm_mod_15m1 is best


#################################################
#################################################
##### Sugar Maple

### Filter to include 2001 as first First.Year
sm_surv_temp <- dplyr::filter(sm_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16,18:25)], na.rm = TRUE)
# Heterospecific
sm_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(sm_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

sm_surv <- sm_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
sm_surv$propEMF_15m = (sm_surv$hetMyc_15m / (sm_surv$hetMyc_15m + sm_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
sm_surv$propsm_15m = (sm_surv$sm_15m / (sm_surv$hetMyc_15m + sm_surv$conMyc_15m))

#### Make site and year into factors
sm_surv$Site = as.factor(sm_surv$Site)
sm_surv$survYear = as.factor(sm_surv$survYear)

#### Log transfosm the count of 0-yr old seedlings
sm_surv$log0yr <- log(sm_surv$total_count_yr0)  

#################### 15-m radius
## Full model
smBA_15m <- glmer(total_count_yr1 ~ sm_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    sm_15m:log0yr:conMycNS_15m + 
                    sm_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(smBA_15m)

## Simplified to exclude 3-way interactions
smBA_15m2 <- glmer(total_count_yr1 ~ sm_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = sm_surv)
summary(smBA_15m2)

## VIFs and pseudo r-squared values
car::vif(smBA_15m) # some >10
r.squaredGLMM(smBA_15m) 
car::vif(smBA_15m2) # all < 2
r.squaredGLMM(smBA_15m2) 
## AIC
aic_values <- AIC(smBA_15m, smBA_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # sm_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = smBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
smProp_15m <- glmer(total_count_yr1 ~ propsm_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propsm_15m:log0yr:propEMF_15m +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = sm_surv)
summary(smProp_15m)

## Simplified to exclude 3-way interactions
smProp_15m2 <- glmer(total_count_yr1 ~ propsm_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = sm_surv)
summary(smProp_15m2)

## VIFs and pseudo r-squared values
car::vif(smProp_15m) # all < 3
r.squaredGLMM(smProp_15m) 
car::vif(smProp_15m2) # all < 2
r.squaredGLMM(smProp_15m2) 
## AIC
aic_values <- AIC(smProp_15m, smProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])

#################################################
#################################################
##### American Beech

### Filter to include 2001 as first First.Year
ab_surv_temp <- dplyr::filter(ab_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:17,19:25)], na.rm = TRUE)
# Heterospecific
ab_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(ab_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

ab_surv <- ab_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
ab_surv$propEMF_15m = (ab_surv$conMyc_15m / (ab_surv$hetMyc_15m + ab_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
ab_surv$propab_15m = (ab_surv$ab_15m / (ab_surv$hetMyc_15m + ab_surv$conMyc_15m))

#### Make site and year into factors
ab_surv$Site = as.factor(ab_surv$Site)
ab_surv$survYear = as.factor(ab_surv$survYear)

#### Log transfoab the count of 0-yr old seedlings
ab_surv$log0yr <- log(ab_surv$total_count_yr0)  

#################### 15-m radius
## Full model - no RE of Site
abBA_15m <- glmer(total_count_yr1 ~ ab_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    ab_15m:log0yr:conMycNS_15m + 
                    ab_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(abBA_15m)

## Simplified to exclude 3-way interactions
abBA_15m2 <- glmer(total_count_yr1 ~ ab_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | survYear),
                   family = poisson(link = "log"), data = ab_surv)
summary(abBA_15m2)

## VIFs and pseudo r-squared values
car::vif(abBA_15m) # some >10
r.squaredGLMM(abBA_15m) 
car::vif(abBA_15m2) # all < 5
r.squaredGLMM(abBA_15m2) 
## AIC
aic_values <- AIC(abBA_15m, abBA_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # ab_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = abBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
abProp_15m <- glmer(total_count_yr1 ~ propab_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propab_15m:log0yr:propEMF_15m +
                      (1 | survYear),
                    family = poisson(link = "log"), data = ab_surv)
summary(abProp_15m)

## Simplified to exclude 3-way interactions
abProp_15m2 <- glmer(total_count_yr1 ~ propab_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = ab_surv)
summary(abProp_15m2)

## VIFs and pseudo r-squared values
car::vif(abProp_15m) # some >10
r.squaredGLMM(abProp_15m) 
car::vif(abProp_15m2) # all < 4
r.squaredGLMM(abProp_15m2) 
## AIC
aic_values <- AIC(abProp_15m, abProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) 


#################################################
#################################################
##### White ash

### Filter to include 2001 as first First.Year
wa_surv_temp <- dplyr::filter(wa_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:18,20:25)], na.rm = TRUE)
# Heterospecific
wa_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(wa_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

wa_surv <- wa_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
wa_surv$propEMF_15m = (wa_surv$hetMyc_15m / (wa_surv$hetMyc_15m + wa_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
wa_surv$propwa_15m = (wa_surv$wa_15m / (wa_surv$hetMyc_15m + wa_surv$conMyc_15m))

#### Make site and year into factors
wa_surv$Site = as.factor(wa_surv$Site)
wa_surv$survYear = as.factor(wa_surv$survYear)

#### Log transfowa the count of 0-yr old seedlings
wa_surv$log0yr <- log(wa_surv$total_count_yr0)  

#################### 15-m radius
## Full model
waBA_15m <- glmer(total_count_yr1 ~ wa_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    wa_15m:log0yr:conMycNS_15m + 
                    wa_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(waBA_15m)

## Simplified to exclude 3-way interactions
waBA_15m2 <- glmer(total_count_yr1 ~ wa_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = wa_surv)
summary(waBA_15m2)

## VIFs and pseudo r-squared values
car::vif(waBA_15m) # some >10
r.squaredGLMM(waBA_15m) 
car::vif(waBA_15m2) # all < 2
r.squaredGLMM(waBA_15m2) 
## AIC
aic_values <- AIC(waBA_15m, waBA_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # wa_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = waBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
waProp_15m <- glmer(total_count_yr1 ~ propwa_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propwa_15m:log0yr:propEMF_15m +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = wa_surv)
summary(waProp_15m)

## Simplified to exclude 3-way interactions
waProp_15m2 <- glmer(total_count_yr1 ~ propwa_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wa_surv)
summary(waProp_15m2)

## VIFs and pseudo r-squared values
car::vif(waProp_15m) # all < 7
r.squaredGLMM(waProp_15m) 
car::vif(waProp_15m2) # all < 2
r.squaredGLMM(waProp_15m2) 
## AIC
aic_values <- AIC(waProp_15m, waProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### Black Cherry

### Filter to include 2001 as first First.Year
bc_surv_temp <- dplyr::filter(bc_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:19,21:25)], na.rm = TRUE)
# Heterospecific
bc_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(bc_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

bc_surv <- bc_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
bc_surv$propEMF_15m = (bc_surv$hetMyc_15m / (bc_surv$hetMyc_15m + bc_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
bc_surv$propbc_15m = (bc_surv$bc_15m / (bc_surv$hetMyc_15m + bc_surv$conMyc_15m))

#### Make site and year into factors
bc_surv$Site = as.factor(bc_surv$Site)
bc_surv$survYear = as.factor(bc_surv$survYear)

#### Log transfobc the count of 0-yr old seedlings
bc_surv$log0yr <- log(bc_surv$total_count_yr0)  

#################### 15-m radius
## Full model - DOES NOT CONVERGE
bcBA_15m <- glmer(total_count_yr1 ~ bc_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    bc_15m:log0yr:conMycNS_15m + 
                    bc_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bcBA_15m)

## Simplified to exclude 3-way interactions
bcBA_15m2 <- glmer(total_count_yr1 ~ bc_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = bc_surv)
summary(bcBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(bcBA_15m) # some >10
# r.squaredGLMM(bcBA_15m) 
car::vif(bcBA_15m2) # all < 5
r.squaredGLMM(bcBA_15m2) 
## AIC
# aic_values <- AIC(bcBA_15m, bcBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) # bc_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = bcBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
bcProp_15m <- glmer(total_count_yr1 ~ propbc_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propbc_15m:log0yr:propEMF_15m +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = bc_surv)
summary(bcProp_15m)

## Simplified to exclude 3-bcy interactions
bcProp_15m2 <- glmer(total_count_yr1 ~ propbc_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = bc_surv)
summary(bcProp_15m2)

## VIFs and pseudo r-squared values
car::vif(bcProp_15m) # one > 10
r.squaredGLMM(bcProp_15m) 
car::vif(bcProp_15m2) # all < 4
r.squaredGLMM(bcProp_15m2) 
## AIC
aic_values <- AIC(bcProp_15m, bcProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### Northern Red Oak

### Filter to include 2001 as first First.Year
nro_surv_temp <- dplyr::filter(nro_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:20,22:25)], na.rm = TRUE)
# Heterospecific
nro_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(nro_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

nro_surv <- nro_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
nro_surv$propEMF_15m = (nro_surv$conMyc_15m / (nro_surv$hetMyc_15m + nro_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
nro_surv$propnro_15m = (nro_surv$nro_15m / (nro_surv$hetMyc_15m + nro_surv$conMyc_15m))

#### Make site and year into factors
nro_surv$Site = as.factor(nro_surv$Site)
nro_surv$survYear = as.factor(nro_surv$survYear)

#### Log transfonro the count of 0-yr old seedlings
nro_surv$log0yr <- log(nro_surv$total_count_yr0)  

#################### 15-m radius
## Full model - DOES NOT CONVERGE
nroBA_15m <- glmer(total_count_yr1 ~ nro_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     nro_15m:log0yr:conMycNS_15m + 
                     nro_15m:log0yr:hetMyc_15m + 
                     hetMyc_15m:log0yr:conMycNS_15m +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = nro_surv)
summary(nroBA_15m)

## Simplified to exclude 3-way interactions
nroBA_15m2 <- glmer(total_count_yr1 ~ nro_15m:log0yr + conMycNS_15m:log0yr + 
                      hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = nro_surv)
summary(nroBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(nroBA_15m) # some >10
# r.squaredGLMM(nroBA_15m) 
car::vif(nroBA_15m2) # all < 3
r.squaredGLMM(nroBA_15m2) 
## AIC
# aic_values <- AIC(nroBA_15m, nroBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) # nro_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = nroBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
nroProp_15m <- glmer(total_count_yr1 ~ propnro_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       propnro_15m:log0yr:propEMF_15m +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = nro_surv)
summary(nroProp_15m)

## Simplified to exclude 3-way interactions - strong correlation btwn propnro_15m and propEMF_15m 
nroProp_15m2 <- glmer(total_count_yr1 ~ propnro_15m:log0yr + propEMF_15m:log0yr + 
                        log(Prior3_counts + 1):log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = nro_surv)
summary(nroProp_15m2)

## VIFs and pseudo r-squared values
car::vif(nroProp_15m) # many > 15
r.squaredGLMM(nroProp_15m) 
car::vif(nroProp_15m2) # some >10
r.squaredGLMM(nroProp_15m2) 
## AIC
aic_values <- AIC(nroProp_15m, nroProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### Black oak

### Filter to include 2001 as first First.Year
bonro_surv_temp <- dplyr::filter(bonro_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:21,23:25)], na.rm = TRUE)
# Heterospecific
bonro_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(bonro_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

bonro_surv <- bonro_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
bonro_surv$propEMF_15m = (bonro_surv$conMyc_15m / (bonro_surv$hetMyc_15m + bonro_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
bonro_surv$propbonro_15m = (bonro_surv$bonro_15m / (bonro_surv$hetMyc_15m + bonro_surv$conMyc_15m))

#### Make site and year into factors
bonro_surv$Site = as.factor(bonro_surv$Site)
bonro_surv$survYear = as.factor(bonro_surv$survYear)

#### Log transfobonro the count of 0-yr old seedlings
bonro_surv$log0yr <- log(bonro_surv$total_count_yr0)  

#################### 15-m radius
## Full model - DOES NOT CONVERGE
bonroBA_15m <- glmer(total_count_yr1 ~ bonro_15m:log0yr + conMycNS_15m:log0yr + 
                       hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                       bonro_15m:log0yr:conMycNS_15m + 
                       bonro_15m:log0yr:hetMyc_15m + 
                       hetMyc_15m:log0yr:conMycNS_15m +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = bonro_surv)
summary(bonroBA_15m)

## Simplified to exclude 3-way interactions
bonroBA_15m2 <- glmer(total_count_yr1 ~ bonro_15m:log0yr + conMycNS_15m:log0yr + 
                        hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = bonro_surv)
summary(bonroBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(bonroBA_15m) # some >10
# r.squaredGLMM(bonroBA_15m) 
car::vif(bonroBA_15m2) # all < 6
r.squaredGLMM(bonroBA_15m2) 
## AIC
# aic_values <- AIC(bonroBA_15m, bonroBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) # bonro_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = bonroBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
bonroProp_15m <- glmer(total_count_yr1 ~ propbonro_15m:log0yr + propEMF_15m:log0yr + 
                         log(Prior3_counts + 1):log0yr +
                         propbonro_15m:log0yr:propEMF_15m +
                         (1 | Site) + (1 | survYear),
                       family = poisson(link = "log"), data = bonro_surv)
summary(bonroProp_15m)

## Simplified to exclude 3-way interactions 
bonroProp_15m2 <- glmer(total_count_yr1 ~ propbonro_15m:log0yr + propEMF_15m:log0yr + 
                          log(Prior3_counts + 1):log0yr +
                          (1 | Site) + (1 | survYear),
                        family = poisson(link = "log"), data = bonro_surv)
summary(bonroProp_15m2)

## VIFs and pseudo r-squared values
car::vif(bonroProp_15m) # many > 20
r.squaredGLMM(bonroProp_15m) 
car::vif(bonroProp_15m2) # some >15
r.squaredGLMM(bonroProp_15m2) 
## AIC
aic_values <- AIC(bonroProp_15m, bonroProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### White oak

### Filter to include 2001 as first First.Year
wo_surv_temp <- dplyr::filter(wo_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:22,24:25)], na.rm = TRUE)
# Heterospecific
wo_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(wo_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

wo_surv <- wo_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
wo_surv$propEMF_15m = (wo_surv$conMyc_15m / (wo_surv$hetMyc_15m + wo_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
wo_surv$propwo_15m = (wo_surv$wo_15m / (wo_surv$hetMyc_15m + wo_surv$conMyc_15m))

#### Make site and year into factors
wo_surv$Site = as.factor(wo_surv$Site)
wo_surv$survYear = as.factor(wo_surv$survYear)

#### Log transfowo the count of 0-yr old seedlings
wo_surv$log0yr <- log(wo_surv$total_count_yr0)  

#################### 15-m radius
## Full model - DOES NOT CONVERGE
woBA_15m <- glmer(total_count_yr1 ~ wo_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    wo_15m:log0yr:conMycNS_15m + 
                    wo_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = wo_surv)
summary(woBA_15m)

## Simplified to exclude 3-way interactions
woBA_15m2 <- glmer(total_count_yr1 ~ wo_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = wo_surv)
summary(woBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(woBA_15m) # some >10
# r.squaredGLMM(woBA_15m) 
car::vif(woBA_15m2) # all < 3
r.squaredGLMM(woBA_15m2) 
## AIC
# aic_values <- AIC(woBA_15m, woBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) # wo_mod_15m1 is best

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = woBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
woProp_15m <- glmer(total_count_yr1 ~ propwo_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propwo_15m:log0yr:propEMF_15m +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = wo_surv)
summary(woProp_15m)

## Simplified to exclude 3-way interactions - correlation btwn propwo_15m and propEMF_15m
woProp_15m2 <- glmer(total_count_yr1 ~ propwo_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wo_surv)
summary(woProp_15m2)

## VIFs and pseudo r-squared values
car::vif(woProp_15m) # many > 19
r.squaredGLMM(woProp_15m) 
car::vif(woProp_15m2) # some >19
r.squaredGLMM(woProp_15m2) 
## AIC
aic_values <- AIC(woProp_15m, woProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### Ironwood

### Filter to include 2001 as first First.Year
iw_surv_temp <- dplyr::filter(iw_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:23,25)], na.rm = TRUE)
# Heterospecific
iw_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(iw_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

iw_surv <- iw_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
iw_surv$propEMF_15m = (iw_surv$conMyc_15m / (iw_surv$hetMyc_15m + iw_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
iw_surv$propiw_15m = (iw_surv$iw_15m / (iw_surv$hetMyc_15m + iw_surv$conMyc_15m))

#### Make site and year into factors
iw_surv$Site = as.factor(iw_surv$Site)
iw_surv$survYear = as.factor(iw_surv$survYear)

#### Log transfoiw the count of 0-yr old seedlings
iw_surv$log0yr <- log(iw_surv$total_count_yr0)  

#################### 15-m radius
## Full model - no RE of Site - DOES NOT CONVERGE
iwBA_15m <- glmer(total_count_yr1 ~ iw_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    iw_15m:log0yr:conMycNS_15m + 
                    iw_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | survYear),
                  family = poisson(link = "log"), data = iw_surv)
summary(iwBA_15m)

## Simplified to exclude 3-way interactions
iwBA_15m2 <- glmer(total_count_yr1 ~ iw_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | survYear),
                   family = poisson(link = "log"), data = iw_surv)
summary(iwBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(iwBA_15m) 
# r.squaredGLMM(iwBA_15m) 
car::vif(iwBA_15m2) # all < 9
r.squaredGLMM(iwBA_15m2) 
## AIC
# aic_values <- AIC(iwBA_15m, iwBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) 

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = iwBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM)
iwProp_15m <- glmer(total_count_yr1 ~ propiw_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propiw_15m:log0yr:propEMF_15m +
                      (1 | survYear),
                    family = poisson(link = "log"), data = iw_surv)
summary(iwProp_15m)

## Simplified to exclude 3-way interactions
iwProp_15m2 <- glmer(total_count_yr1 ~ propiw_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = iw_surv)
summary(iwProp_15m2)

## VIFs and pseudo r-squared values
car::vif(iwProp_15m) # all < 2
r.squaredGLMM(iwProp_15m) 
car::vif(iwProp_15m2) # all < 2
r.squaredGLMM(iwProp_15m2) 
## AIC
aic_values <- AIC(iwProp_15m, iwProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])


#################################################
#################################################
##### American Basswood

### Filter to include 2001 as first First.Year
bs_surv_temp <- dplyr::filter(bs_surv_temp, First.Year > 2000)

# Calculate sum of previous 3-year counts of hetero- and conspecific seedlings 
allSp_TI_ObsYear$Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:14,16:24)], na.rm = TRUE)
# Heterospecific
bs_surv_temp$Prior3_counts = 
  allSp_TI_ObsYear$Prior3_counts[match(bs_surv_temp$TI_FY,
                                       allSp_TI_ObsYear$TI_FY)]

bs_surv <- bs_surv_temp

#### Calculate the proportion of EMF hosts in 15-m radius
bs_surv$propEMF_15m = (bs_surv$conMyc_15m / (bs_surv$hetMyc_15m + bs_surv$conMyc_15m))

#### Calculate the relative proportion of conspecific adults in 15-m radius
bs_surv$propbs_15m = (bs_surv$bs_15m / (bs_surv$hetMyc_15m + bs_surv$conMyc_15m))

#### Make site and year into factors
bs_surv$Site = as.factor(bs_surv$Site)
bs_surv$survYear = as.factor(bs_surv$survYear)

#### Log transfobs the count of 0-yr old seedlings
bs_surv$log0yr <- log(bs_surv$total_count_yr0)  

#################### 15-m radius
## Full model - no RE of Site - DOES NOT CONVERGE
bsBA_15m <- glmer(total_count_yr1 ~ bs_15m:log0yr + conMycNS_15m:log0yr + 
                    hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                    bs_15m:log0yr:conMycNS_15m + 
                    bs_15m:log0yr:hetMyc_15m + 
                    hetMyc_15m:log0yr:conMycNS_15m +
                    (1 | survYear),
                  family = poisson(link = "log"), data = bs_surv)
summary(bsBA_15m)

## Simplified to exclude 3-way interactions - also does not converge...
bsBA_15m2 <- glmer(total_count_yr1 ~ bs_15m:log0yr + conMycNS_15m:log0yr + 
                     hetMyc_15m:log0yr + log(Prior3_counts + 1):log0yr +
                     (1 | survYear),
                   family = poisson(link = "log"), data = bs_surv)
summary(bsBA_15m2)

## VIFs and pseudo r-squared values
# car::vif(bsBA_15m) 
# r.squaredGLMM(bsBA_15m) 
car::vif(bsBA_15m2) # all < 9
r.squaredGLMM(bsBA_15m2) 
## AIC
# aic_values <- AIC(bsBA_15m, bsBA_15m2)
# # Sort AIC values from lowest to highest and show
# print(aic_values[order(aic_values$AIC), ]) 

# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = bsBA_15m, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)

### Biotic model of one-year survivorship - with adult proportions (%ECM) - does not converge
bsProp_15m <- glmer(total_count_yr1 ~ propbs_15m:log0yr + propEMF_15m:log0yr + 
                      log(Prior3_counts + 1):log0yr +
                      propbs_15m:log0yr:propEMF_15m +
                      (1 | survYear),
                    family = poisson(link = "log"), data = bs_surv)
summary(bsProp_15m)

## Simplified to exclude 3-way interactions - does not converge
bsProp_15m2 <- glmer(total_count_yr1 ~ propbs_15m:log0yr + propEMF_15m:log0yr + 
                       log(Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = bs_surv)
summary(bsProp_15m2)

## VIFs and pseudo r-squared values
car::vif(bsProp_15m) # all < 2
r.squaredGLMM(bsProp_15m) 
car::vif(bsProp_15m2) # all < 2
r.squaredGLMM(bsProp_15m2) 
## AIC
aic_values <- AIC(bsProp_15m, bsProp_15m2)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ])



# ##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years
# 
# # Step 1: Create a sequence of basal area values
# ba_seq <- seq(0, max(rm_surv$rm_15m, rm_surv$conMycNS_15m, 
#               rm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)
# 
# # Step 2: Set median values for the non-focal covariates
# median_vals <- rm_surv %>%
#   summarize(across(c(rm_Prior3_counts, hetSp_Prior3_counts, total_count_yr0,
#           rm_15m, conMycNS_15m, hetMyc_15m), ~ median(.x, na.rm = TRUE)))
# 
# # Step 3: Create prediction datasets for each focal variable
# make_pred_data <- function(focal, label) {
#   df <- data.frame(
#     rm_15m = if (focal == "rm_15m") ba_seq else median_vals$rm_15m,
#     conMycNS_15m = if (focal == "conMycNS_15m") ba_seq else median_vals$conMycNS_15m,
#     hetMyc_15m = if (focal == "hetMyc_15m") ba_seq else median_vals$hetMyc_15m,
#     rm_Prior3_counts = median_vals$rm_Prior3_counts,
#     hetSp_Prior3_counts = median_vals$hetSp_Prior3_counts,
#     total_count_yr0 = median_vals$total_count_yr0,
#     Site = rm_surv$Site[1],survYear = rm_surv$survYear[1],
#     focal_label = label)
#   df
# }
# 
# # Step 4: Combine prediction data for all three components
# new_data <- bind_rows(make_pred_data("rm_15m", "Conspecific BA"),
#   make_pred_data("conMycNS_15m", "Con-mycorrhizal BA"),
#   make_pred_data("hetMyc_15m", "Hetero-mycorrhizal BA"))
# 
# # Step 5: Predict survivorship without random effects (average effect)
# new_data$pred_surv <- predict(rm_mod_15m, newdata = new_data, 
#                       type = "response", re.form = NA)
# # re.form = NA ignores RE, uses only the fixed effects (i.e., an “average” prediction across sites and years).
# 
# # Step 6: Add x-axis for the focal BA value
# new_data$x_focal <- with(new_data, 
#               ifelse(focal_label == "Conspecific BA", rm_15m,
#               ifelse(focal_label == "Con-mycorrhizal BA", conMycNS_15m, hetMyc_15m)))
# 
# # Step 7: Plot
# ggplot(new_data, aes(x = x_focal, y = pred_surv, color = focal_label)) +
#   geom_line(linewidth = 1) +
#   labs(x = "Basal Area in 15-m Radius (m²)", 
#     y = "Predicted Survivorship", 
#     color = "",title = "Red Maple") +
#   theme_bw() +
#   theme(axis.text = element_text(colour = 'black', size = 10), 
#     axis.title = element_text(size = 12), legend.key = element_blank(), 
#     plot.title = element_text(size = 14, hjust = 0.5),  # Slightly larger title
#     legend.text = element_text(size = 10), 
#     legend.title = element_text(size = 10), 
#     panel.grid.minor = element_blank(),
#     legend.position = c(0.05, 0.05),  # Bottom left inside plot
#     legend.justification = c(0, 0))