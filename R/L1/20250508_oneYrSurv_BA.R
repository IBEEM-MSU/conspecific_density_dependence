
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 8 May 2025
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
library(broom.mixed)
library(tibble)

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

##### Read in file with species' average cumulative BA per site
speciesSite_sum = read.csv("./data/L1/mean_totalSpSiteBA.csv")
speciesSite_sum$Site = factor(speciesSite_sum$Site)
speciesSite_sum = speciesSite_sum[,c(1,2,4)]

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

# Rename data
rm_surv <- rm_surv_temp

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Log transform the count of 0-yr old seedlings
rm_surv$log0yr <- log(rm_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
rm_surv$logPrior = log(rm_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
rm_surv$propEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))
rm_surv$propEMF_15m[rm_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "rm"
# Step 2: Add Sp column to match join key
rm_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
rm_surv <- rm_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
rm_surv$total_yr0_std = rm_surv$total_count_yr0/rm_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(rm_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Red maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(rm_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Red maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(rm_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Red maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(rm_surv, aes(x = rm_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Red maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius
rmBA_15m <- glmer(total_count_yr1 ~ rm_15m + conMycNS_15m + 
                     hetMyc_15m + rm_15m:hetMyc_15m +
                     conMycNS_15m:hetMyc_15m +
                     logPrior + log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = rm_surv)
summary(rmBA_15m)

## VIFs and pseudo r-squared values
car::vif(rmBA_15m) # all < 7
r.squaredGLMM(rmBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = rmBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
sm_surv <- sm_surv_temp

#### Make site and year into factors
sm_surv$Site = as.factor(sm_surv$Site)
sm_surv$survYear = as.factor(sm_surv$survYear)

#### Log transform the count of 0-yr old seedlings
sm_surv$log0yr <- log(sm_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
sm_surv$logPrior = log(sm_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
sm_surv$propEMF_15m = (sm_surv$hetMyc_15m / (sm_surv$hetMyc_15m + sm_surv$conMyc_15m))
sm_surv$propEMF_15m[sm_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "sm"
# Step 2: Add Sp column to match join key
sm_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
sm_surv <- sm_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
sm_surv$total_yr0_std = sm_surv$total_count_yr0/sm_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(sm_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Sugar maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(sm_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Sugar maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(sm_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Sugar maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(sm_surv, aes(x = sm_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Sugar maple") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius
smBA_15m <- glmer(total_count_yr1 ~ sm_15m + conMycNS_15m + 
                    hetMyc_15m + sm_15m:hetMyc_15m +
                    conMycNS_15m:hetMyc_15m +
                    logPrior + log0yr +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(smBA_15m)

## VIFs and pseudo r-squared values
car::vif(smBA_15m) # all < 3
r.squaredGLMM(smBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = smBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
ab_surv <- ab_surv_temp

#### Make site and year into factors
ab_surv$Site = as.factor(ab_surv$Site)
ab_surv$survYear = as.factor(ab_surv$survYear)

#### Log transform the count of 0-yr old seedlings
ab_surv$log0yr <- log(ab_surv$total_count_yr0) 
## Calculate the log-transformed prior 3 year count variable + 1
ab_surv$logPrior = log(ab_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
ab_surv$propEMF_15m = (ab_surv$conMyc_15m / (ab_surv$hetMyc_15m + ab_surv$conMyc_15m))
ab_surv$propEMF_15m[ab_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "ab"
# Step 2: Add Sp column to match join key
ab_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
ab_surv <- ab_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
ab_surv$total_yr0_std = ab_surv$total_count_yr0/ab_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(ab_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "American beech") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(ab_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "American beech") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(ab_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "American beech") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(ab_surv, aes(x = ab_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "American beech") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius - no RE of site
abBA_15m <- glmer(total_count_yr1 ~ ab_15m + conMycNS_15m + 
                    hetMyc_15m + ab_15m:hetMyc_15m +
                    conMycNS_15m:hetMyc_15m +
                    logPrior + log0yr +
                    (1 | survYear),
                  family = poisson(link = "log"), data = ab_surv)
summary(abBA_15m)

## VIFs and pseudo r-squared values
car::vif(abBA_15m) # all < 10
r.squaredGLMM(abBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = abBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
wa_surv <- wa_surv_temp

#### Make site and year into factors
wa_surv$Site = as.factor(wa_surv$Site)
wa_surv$survYear = as.factor(wa_surv$survYear)

#### Log transform the count of 0-yr old seedlings
wa_surv$log0yr <- log(wa_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
wa_surv$logPrior = log(wa_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
wa_surv$propEMF_15m = (wa_surv$hetMyc_15m / (wa_surv$hetMyc_15m + wa_surv$conMyc_15m))
wa_surv$propEMF_15m[wa_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "wa"
# Step 2: Add Sp column to match join key
wa_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
wa_surv <- wa_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
wa_surv$total_yr0_std = wa_surv$total_count_yr0/wa_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(wa_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White ash") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(wa_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White ash") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(wa_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White ash") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(wa_surv, aes(x = wa_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White ash") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius
waBA_15m <- glmer(total_count_yr1 ~ wa_15m + conMycNS_15m + 
                    hetMyc_15m + wa_15m:hetMyc_15m +
                    conMycNS_15m:hetMyc_15m +
                    logPrior + log0yr +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(waBA_15m)

## VIFs and pseudo r-squared values
car::vif(waBA_15m) # all < 9
r.squaredGLMM(waBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = waBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
bc_surv <- bc_surv_temp

#### Make site and year into factors
bc_surv$Site = as.factor(bc_surv$Site)
bc_surv$survYear = as.factor(bc_surv$survYear)

#### Log transform the count of 0-yr old seedlings
bc_surv$log0yr <- log(bc_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
bc_surv$logPrior = log(bc_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
bc_surv$propEMF_15m = (bc_surv$hetMyc_15m / (bc_surv$hetMyc_15m + bc_surv$conMyc_15m))
bc_surv$propEMF_15m[bc_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "bc"
# Step 2: Add Sp column to match join key
bc_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
bc_surv <- bc_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
bc_surv$total_yr0_std = bc_surv$total_count_yr0/bc_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(bc_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black cherry") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(bc_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black cherry") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(bc_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black cherry") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(bc_surv, aes(x = bc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black cherry") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius
bcBA_15m <- glmer(total_count_yr1 ~ bc_15m + conMycNS_15m + 
                    hetMyc_15m + bc_15m:hetMyc_15m +
                    conMycNS_15m:hetMyc_15m +
                    logPrior + log0yr +
                    (1 | Site) + (1 | survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bcBA_15m)

## VIFs and pseudo r-squared values
car::vif(bcBA_15m) # all < 6
r.squaredGLMM(bcBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = bcBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
nro_surv <- nro_surv_temp

#### Make site and year into factors
nro_surv$Site = as.factor(nro_surv$Site)
nro_surv$survYear = as.factor(nro_surv$survYear)

#### Log transform the count of 0-yr old seedlings
nro_surv$log0yr <- log(nro_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
nro_surv$logPrior = log(nro_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
nro_surv$propEMF_15m = (nro_surv$conMyc_15m / (nro_surv$hetMyc_15m + nro_surv$conMyc_15m))
nro_surv$propEMF_15m[nro_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "nro"
# Step 2: Add Sp column to match join key
nro_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
nro_surv <- nro_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
nro_surv$total_yr0_std = nro_surv$total_count_yr0/nro_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(nro_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Northern red oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(nro_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Northern red oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(nro_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Northern red oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(nro_surv, aes(x = nro_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Northern red oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


#################### 15-m radius
nroBA_15m <- glmer(total_count_yr1 ~ nro_15m + conMycNS_15m + 
                     hetMyc_15m + nro_15m:hetMyc_15m +
                     conMycNS_15m:hetMyc_15m +
                     logPrior + log0yr +
                     (1 | Site) + (1 | survYear),
                   family = poisson(link = "log"), data = nro_surv)
summary(nroBA_15m)

## VIFs and pseudo r-squared values
car::vif(nroBA_15m) # all < 10
r.squaredGLMM(nroBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = nroBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
bonro_surv <- bonro_surv_temp

#### Make site and year into factors
bonro_surv$Site = as.factor(bonro_surv$Site)
bonro_surv$survYear = as.factor(bonro_surv$survYear)

#### Log transform the count of 0-yr old seedlings
bonro_surv$log0yr <- log(bonro_surv$total_count_yr0) 
## Calculate the log-transformed prior 3 year count variable + 1
bonro_surv$logPrior = log(bonro_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
bonro_surv$propEMF_15m = (bonro_surv$conMyc_15m / (bonro_surv$hetMyc_15m + bonro_surv$conMyc_15m))
bonro_surv$propEMF_15m[bonro_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "bonro"
# Step 2: Add Sp column to match join key
bonro_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
bonro_surv <- bonro_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
bonro_surv$total_yr0_std = bonro_surv$total_count_yr0/bonro_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(bonro_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(bonro_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black oak") +
  scale_x_continuous(limits = c(0, 1.25)) +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(bonro_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(bonro_surv, aes(x = bonro_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Black oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius 
bonroBA_15m <- glmer(total_count_yr1 ~ bonro_15m + conMycNS_15m + 
                    hetMyc_15m + #bonro_15m:hetMyc_15m +
                    conMycNS_15m:hetMyc_15m +
                    logPrior + log0yr +
                    #(1 | Site) + 
                    (1 | survYear),
                  family = poisson(link = "log"), data = bonro_surv)
summary(bonroBA_15m)

## VIFs and pseudo r-squared values
car::vif(bonroBA_15m) # all < 8
r.squaredGLMM(bonroBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = bonroBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
wo_surv <- wo_surv_temp

#### Make site and year into factors
wo_surv$Site = as.factor(wo_surv$Site)
wo_surv$survYear = as.factor(wo_surv$survYear)

#### Log transform the count of 0-yr old seedlings
wo_surv$log0yr <- log(wo_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
wo_surv$logPrior = log(wo_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
wo_surv$propEMF_15m = (wo_surv$conMyc_15m / (wo_surv$hetMyc_15m + wo_surv$conMyc_15m))
wo_surv$propEMF_15m[wo_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "wo"
# Step 2: Add Sp column to match join key
wo_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
wo_surv <- wo_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
wo_surv$total_yr0_std = wo_surv$total_count_yr0/wo_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(wo_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(wo_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(wo_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(wo_surv, aes(x = wo_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "White oak") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius 
woBA_15m <- glmer(total_count_yr1 ~ wo_15m + conMycNS_15m + 
                       hetMyc_15m + #wo_15m:hetMyc_15m +
                       conMycNS_15m:hetMyc_15m +
                       logPrior + log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wo_surv)
summary(woBA_15m)

## VIFs and pseudo r-squared values
car::vif(woBA_15m) # all < 7
r.squaredGLMM(woBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = woBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)

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

# Rename
iw_surv <- iw_surv_temp

#### Make site and year into factors
iw_surv$Site = as.factor(iw_surv$Site)
iw_surv$survYear = as.factor(iw_surv$survYear)

#### Log transform the count of 0-yr old seedlings
iw_surv$log0yr <- log(iw_surv$total_count_yr0)  
## Calculate the log-transformed prior 3 year count variable + 1
iw_surv$logPrior = log(iw_surv$Prior3_counts + 1)

#### Calculate the proportion of EMF hosts in 15-m radius
iw_surv$propEMF_15m = (iw_surv$conMyc_15m / (iw_surv$hetMyc_15m + iw_surv$conMyc_15m))
iw_surv$propEMF_15m[iw_surv$propEMF_15m == "NaN"] <- 0

# Step 1: Define species code from the name
sp_code <- "iw"
# Step 2: Add Sp column to match join key
iw_surv$Sp <- sp_code
# Step 3: Join with the species BA dataframe
iw_surv <- iw_surv %>%
  left_join(speciesSite_sum, by = c("Site", "Sp")) 
# Make a column for the standardized 0-yr counts
iw_surv$total_yr0_std = iw_surv$total_count_yr0/iw_surv$spBA_corr

##### Plot - Standardized 0 yr versus Proportion ECM Basal Area
ggplot(iw_surv, aes(x = propEMF_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Proportion ECM Basal Area", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Ironwood") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Heteromycorrhizal Basal Area
ggplot(iw_surv, aes(x = hetMyc_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Heteromycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Ironwood") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conmycorrhizal Basal Area
ggplot(iw_surv, aes(x = conMycNS_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conmycorrhizal Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Ironwood") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))
##### Plot - Standardized 0 yr versus Total Conspecific Basal Area
ggplot(iw_surv, aes(x = iw_15m, y = total_yr0_std)) +
  geom_point() +  
  theme_bw() +
  labs(x = "Total Conspecific Basal Area (m²)", 
       y = "Standardized Count (No. 0-Year Seedlings/Adult Basal Area, m²)",
       title = "Ironwood") +
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

#################### 15-m radius
iwBA_15m <- glmer(total_count_yr1 ~ iw_15m + conMycNS_15m + 
                       hetMyc_15m + #iw_15m:hetMyc_15m +
                       conMycNS_15m:hetMyc_15m +
                       logPrior + log0yr +
                       #(1 | Site) + 
                       (1 | survYear),
                     family = poisson(link = "log"), data = iw_surv)
summary(iwBA_15m)

## VIFs and pseudo r-squared values
car::vif(iwBA_15m) # one = 11
r.squaredGLMM(iwBA_15m) 

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = iwBA_15m, n = 1000)
# QQ plot of residuals
plotQQunif(sim_res)


################################################ 
######## Plotting model coefficients +/- SE
# Function to pull out relevant terms
extract_model_summary <- function(model, species_abbrev) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    mutate(species = species_abbrev) %>%
    select(species, term, estimate, std.error, statistic, p.value)
}

# Create a list of our models
models <- list(rm = rmBA_15m, sm = smBA_15m, ab = abBA_15m, wa = waBA_15m,
               bc = bcBA_15m, nro = nroBA_15m, bonro = bonroBA_15m,
               wo = woBA_15m, iw = iwBA_15m)

# Map terms by model
allSp_results <- purrr::imap_dfr(models, extract_model_summary)

# Add uniform names for covariates across species
allSp_results <- allSp_results %>%
  rowwise() %>%
  mutate(component = case_when(
    term == "conMycNS_15m" ~ "Conmycorrhizal_BA",
    term == "hetMyc_15m" ~ "Heteromycorrhizal_BA",
    term == "conMycNS_15m:hetMyc_15m" ~ "Conmyc_x_Heteromycorrhizal_BA",
    term == "logPrior" ~ "Log_Prior_Density",
    term == "log0yr" ~ "Log_Recruit_Density",
    grepl(paste0("^", species, "_15m:hetMyc_15m$"), term) ~
      "Consp_x_Heteromycorrhizal_BA",
    grepl(paste0("^", species, "_15m$"), term) ~ "Conspecific_BA",
    TRUE ~ NA_character_
  )) %>%
  ungroup()
# Make into an ordered factor
allSp_results$component = factor(allSp_results$component,
                        levels = c( "Log_Prior_Density", 
                        "Log_Recruit_Density",
                        "Conmyc_x_Heteromycorrhizal_BA",
                        "Consp_x_Heteromycorrhizal_BA",
                        "Heteromycorrhizal_BA",
                        "Conmycorrhizal_BA",
                        "Conspecific_BA"))

# Add AM or ECM to each row by species
allSp_results$Myc_Type = ifelse(allSp_results$species == "rm", "AM",
                         ifelse(allSp_results$species == "sm", "AM",
                         ifelse(allSp_results$species == "wa", "AM",
                         ifelse(allSp_results$species == "bc", "AM","ECM"))))

# Add column for significance of parameter estimate
allSp_results$Significant = ifelse(allSp_results$p.value <= 0.05, "yes", "no")

View(allSp_results)

######################
## Create a color vector with red shades for AM and blue shades for ECM
species_cols <- c(ab = "#B4D4DA", bonro = "#6EB6D9", iw = "#3F8BBA",
  nro = "#1F74B1", wo = "#244D7F",
  bc = "#FEB1A3", rm = "#F77964",  
  sm = "#D2241F", wa = "#AE123A")

## Create a size vector for species to scale points by mean percent ba
species_sizes <- data.frame(species = c("ab", "bonro", "iw", "nro", "wo",
                                        "bc", "rm", "sm", "wa"),
                            percBA = c(2.43, 15.20, 0.11, 29.7, 10.30,
                                    1.03, 9.40, 15.80, 1.82))
# Add data to results df
allSp_results$percBA = species_sizes$percBA[match(allSp_results$species,
                                                  species_sizes$species)]

# Create ordered factor for species by Myc_Type group and alphabetical within each
allSp_results <- allSp_results %>%
  mutate(species_label = factor(
      species,levels = c(
        sort(unique(species[Myc_Type == "AM"])),
        sort(unique(species[Myc_Type == "ECM"])))),
    component = factor(component, levels = c(
      "Log_Prior_Density", 
      "Log_Recruit_Density",
      "Conmyc_x_Heteromycorrhizal_BA",
      "Consp_x_Heteromycorrhizal_BA",
      "Heteromycorrhizal_BA",
      "Conmycorrhizal_BA",
      "Conspecific_BA")),
    Significant = factor(Significant, levels = c("yes", "no")))

## Calculate mean parameter estimates for AM versus ECM
group_means <- allSp_results %>%
  group_by(Myc_Type, component) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    std.error = sd(estimate, na.rm = TRUE) / sqrt(n()),  
    .groups = "drop")

# Filter to exclude one outlier
allSp_results2 = dplyr::filter(allSp_results, std.error <= 5)

# Plot - All Species (ECM and AM)
ggplot(allSp_results2, aes(x = estimate, y = component, 
                           color = species_label, shape = Significant)) +
  geom_point(aes(size = log1p(percBA)), 
             position = position_dodge(width = 0.6)) +   # size only here
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error), 
                 height = 0, position = position_dodge(width = 0.6)) +
  ## Add ECM and AM group means
  geom_point(data = group_means, 
             aes(x = mean_estimate, y = component, shape = Myc_Type), 
             color = "black", size = 3, inherit.aes = FALSE) +
  geom_errorbarh(data = group_means,
                 aes(xmin = mean_estimate - std.error, 
                     xmax = mean_estimate + std.error, y = component),
                 color = "black", height = 0, inherit.aes = FALSE) +
  scale_color_manual(values = species_cols) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1, "AM" = 17, "ECM" = 15)) +
  scale_size(range = c(1.5, 5), 
             name = "log(% BA + 1)") +  # nicely scaled
  theme_bw() +
  labs(x = "Effect Size ± SE (Log Count Scale)", 
       y = "", color = "Species", shape = "") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

# Plot - ECM only
allSp_results2 %>% filter(., Myc_Type == "ECM")  %>% 
ggplot(., aes(x = estimate,  y = component, 
                           color = species_label,  shape = Significant)) +
  geom_point(aes(size = log1p(percBA)), 
             position = position_dodge(width = 0.6)) +   # size only here
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error), 
                 height = 0, position = position_dodge(width = 0.6)) +
  ## Add ECM and AM group means
  geom_point(data = group_means %>% filter(Myc_Type == "ECM"), 
             aes(x = mean_estimate, y = component, shape = Myc_Type), 
             color = "black", size = 3, inherit.aes = FALSE) +
  geom_errorbarh(data = group_means %>% filter(Myc_Type == "ECM"),
                 aes(xmin = mean_estimate - std.error, 
                     xmax = mean_estimate + std.error, y = component),
                 color = "black", height = 0, inherit.aes = FALSE) +
  scale_color_manual(values = species_cols) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1, 
                                "AM" = 17, "ECM" = 15)) +
  scale_size(range = c(1.5, 5), 
             name = "log(% BA + 1)") +  # nicely scaled
  theme_bw() +
  labs(x = "Effect Size ± SE (Log Count Scale)", 
       y = "", color = "Species", shape = "") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))

# Plot - AM only
allSp_results2 %>% filter(., Myc_Type == "AM")  %>% 
  ggplot(., aes(x = estimate,  y = component, 
                color = species_label,  shape = Significant)) +
  geom_point(aes(size = log1p(percBA)), 
             position = position_dodge(width = 0.6)) +   # size only here
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error), 
                 height = 0, position = position_dodge(width = 0.6)) +
  ## Add ECM and AM group means
  geom_point(data = group_means %>% filter(Myc_Type == "AM"), 
             aes(x = mean_estimate, y = component, shape = Myc_Type), 
             color = "black", size = 3, inherit.aes = FALSE) +
  geom_errorbarh(data = group_means %>% filter(Myc_Type == "AM"),
                 aes(xmin = mean_estimate - std.error, 
                     xmax = mean_estimate + std.error, y = component),
                 color = "black", height = 0, inherit.aes = FALSE) +
  scale_color_manual(values = species_cols) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1, 
                                "AM" = 17, "ECM" = 15)) +
  scale_size(range = c(1.5, 5), 
             name = "log(% BA + 1)") +  # nicely scaled
  theme_bw() +
  labs(x = "Effect Size ± SE (Log Count Scale)", 
       y = "", color = "Species", shape = "") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))


########################################################################
##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods - Averaged across sites and years
########################################################################
########################################################################
# ######### RED MAPLE EXAMPLE
# ######### Median only
# # Step 1: Create a sequence of basal area values
# ba_seq_list <- list(
#   rm_15m = seq(0, max(rm_surv$rm_15m, na.rm = TRUE), length.out = 100),
#   conMycNS_15m = seq(0, max(rm_surv$conMycNS_15m, na.rm = TRUE), length.out = 100),
#   hetMyc_15m = seq(0, max(rm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100)
# )
# 
# # Step 2: Get median values of other covariates
# med_vals <- list(
#   log0yr = median(rm_surv$log0yr, na.rm = TRUE),
#   logPrior = median(rm_surv$logPrior, na.rm = TRUE),
#   rm_15m = median(rm_surv$rm_15m, na.rm = TRUE),
#   conMycNS_15m = median(rm_surv$conMycNS_15m, na.rm = TRUE),
#   hetMyc_15m = median(rm_surv$hetMyc_15m, na.rm = TRUE)
# )
# 
# # Step 3: Prediction data function
# make_pred_data_median <- function(focal_var, ba_seq) {
#   n <- length(ba_seq)
#   
#   data.frame(
#     rm_15m = if (focal_var == "rm_15m") ba_seq else 
#       rep(med_vals$rm_15m, n),
#     conMycNS_15m = if (focal_var == "conMycNS_15m") ba_seq else 
#       rep(med_vals$conMycNS_15m, n),
#     hetMyc_15m = if (focal_var == "hetMyc_15m") ba_seq else 
#       rep(med_vals$hetMyc_15m, n),
#     log0yr = rep(med_vals$log0yr, n),
#     logPrior = rep(med_vals$logPrior, n),
#     Site = rep(rm_surv$Site[1], n),
#     survYear = rep(rm_surv$survYear[1], n),
#     focal_label = rep(focal_var, n)
#   )
# }
# 
# # Step 4: Generate prediction data
# pred_data <- bind_rows(
#   make_pred_data_median("rm_15m", ba_seq_list$rm_15m),
#   make_pred_data_median("conMycNS_15m", ba_seq_list$conMycNS_15m),
#   make_pred_data_median("hetMyc_15m", ba_seq_list$hetMyc_15m)
# )
# 
# # Step 5: Predict using only fixed effects (re.form = NA)
# pred_data$pred_count <- predict(rmBA_15m, newdata = pred_data, type = "response", re.form = NA)
# 
# # Add x-axis variable
# pred_data$x_focal <- with(pred_data,
#                           ifelse(focal_label == "rm_15m", rm_15m,
#                                  ifelse(focal_label == "conMycNS_15m", conMycNS_15m, hetMyc_15m)))
# 
# # Make sure focal_label is a factor with levels matching the color manual
# pred_data$focal_label <- factor(pred_data$focal_label,
#                  levels = c("rm_15m", "conMycNS_15m", "hetMyc_15m"),
#                 labels = c("Conspecific BA", "Conmycorrhizal BA", "Heteromycorrhizal BA"))
# 
# # Step 6: Plot
# ggplot(pred_data, aes(x = x_focal, y = pred_count, color = focal_label)) +
#   geom_line(size = 1) +
#   labs(x = "Basal Area in 15-m Radius (m²)",
#        y = "Predicted 1-Year Seedling Count",
#        color = "Neighborhood Type") +
#   scale_color_manual(values = c("Conspecific BA" = "black", 
#                                 "Conmycorrhizal BA" = "#DA3742", 
#                                 "Heteromycorrhizal BA" = "#3885B6")) +
#   theme_bw() +
#   theme(axis.text = element_text(color = "black"),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 11))
# 
# ########################################################################
# ######### 25% and 75% quantiles
# # Step 1: Create sequence of focal BA values
# ba_seq_list <- list(
#   rm_15m = seq(0, max(rm_surv$rm_15m, na.rm = TRUE), length.out = 100),
#   conMycNS_15m = seq(0, max(rm_surv$conMycNS_15m, na.rm = TRUE), length.out = 100),
#   hetMyc_15m = seq(0, max(rm_surv$hetMyc_15m, na.rm = TRUE), length.out = 100))
# 
# # Step 2: Get quantile values of all covariates
# q_vals <- list(
#   q25 = list(
#     log0yr = quantile(rm_surv$log0yr, 0.25, na.rm = TRUE),
#     logPrior = quantile(rm_surv$logPrior, 0.25, na.rm = TRUE),
#     rm_15m = quantile(rm_surv$rm_15m, 0.25, na.rm = TRUE),
#     conMycNS_15m = quantile(rm_surv$conMycNS_15m, 0.25, na.rm = TRUE),
#     hetMyc_15m = quantile(rm_surv$hetMyc_15m, 0.25, na.rm = TRUE)
#   ),
#   q75 = list(
#     log0yr = quantile(rm_surv$log0yr, 0.75, na.rm = TRUE),
#     logPrior = quantile(rm_surv$logPrior, 0.75, na.rm = TRUE),
#     rm_15m = quantile(rm_surv$rm_15m, 0.75, na.rm = TRUE),
#     conMycNS_15m = quantile(rm_surv$conMycNS_15m, 0.75, na.rm = TRUE),
#     hetMyc_15m = quantile(rm_surv$hetMyc_15m, 0.75, na.rm = TRUE)
#   )
# )
# 
# # Step 3: Prediction data function
# make_pred_data_conditional <- function(focal_var, ba_seq, qval) {
#   n <- length(ba_seq)
#   data.frame(
#     rm_15m = if (focal_var == "rm_15m") ba_seq else 
#       rep(qval$rm_15m, n),
#     conMycNS_15m = if (focal_var == "conMycNS_15m") ba_seq else 
#       rep(qval$conMycNS_15m, n),
#     hetMyc_15m = if (focal_var == "hetMyc_15m") ba_seq else 
#       rep(qval$hetMyc_15m, n),
#     log0yr = rep(qval$log0yr, n),
#     logPrior = rep(qval$logPrior, n),
#     Site = rep(rm_surv$Site[1], n),
#     survYear = rep(rm_surv$survYear[1], n),
#     focal_label = rep(focal_var, n))
# }
# 
# # Step 4: Generate prediction datasets
# pred_data_25 <- bind_rows(
#   make_pred_data_conditional("rm_15m", ba_seq_list$rm_15m, q_vals$q25),
#   make_pred_data_conditional("conMycNS_15m", ba_seq_list$conMycNS_15m, q_vals$q25),
#   make_pred_data_conditional("hetMyc_15m", ba_seq_list$hetMyc_15m, q_vals$q25)
# )
# 
# pred_data_75 <- bind_rows(
#   make_pred_data_conditional("rm_15m", ba_seq_list$rm_15m, q_vals$q75),
#   make_pred_data_conditional("conMycNS_15m", ba_seq_list$conMycNS_15m, q_vals$q75),
#   make_pred_data_conditional("hetMyc_15m", ba_seq_list$hetMyc_15m, q_vals$q75)
# )
# 
# # Step 5: Predict using model (fixed effects only)
# pred_data_25$pred_count <- predict(rmBA_15m, newdata = pred_data_25, type = "response", re.form = NA)
# pred_data_75$pred_count <- predict(rmBA_15m, newdata = pred_data_75, type = "response", re.form = NA)
# 
# # Step 6: Format for plotting
# for (df_name in c("pred_data_25", "pred_data_75")) {
#   df <- get(df_name)
#   
#   df$x_focal <- with(df,
#           ifelse(focal_label == "rm_15m", rm_15m,
#           ifelse(focal_label == "conMycNS_15m", conMycNS_15m, hetMyc_15m)))
#   
#   df$focal_label <- factor(df$focal_label,
#           levels = c("rm_15m", "conMycNS_15m", "hetMyc_15m"),
#     labels = c("Conspecific BA", "Conmycorrhizal BA", "Heteromycorrhizal BA"))
#   
#   assign(df_name, df)  # Re-assign updated data frame
# }
# 
# # Step 7: Plot for 25th percentile
# ggplot(pred_data_25, aes(x = x_focal, y = pred_count, color = focal_label)) +
#   geom_line(size = 1) +
#   labs(x = "Basal Area in 15-m Radius (m²)",
#        y = "Predicted 1-Year Seedling Count",
#        color = "Neighborhood Type",
#        title = "Predictions at 25th Percentile of All Other Covariates") +
#   scale_color_manual(values = c("Conspecific BA" = "black", 
#                                 "Conmycorrhizal BA" = "#DA3742", 
# "Heteromycorrhizal BA" = "#3885B6")) +
#   theme_bw()
# 
# # Step 8: Plot for 75th percentile
# ggplot(pred_data_75, aes(x = x_focal, y = pred_count, color = focal_label)) +
#   geom_line(size = 1) +
#   labs(x = "Basal Area in 15-m Radius (m²)",
#        y = "Predicted 1-Year Seedling Count",
#        color = "Neighborhood Type",
#        title = "Predictions at 75th Percentile of All Other Covariates") +
#   scale_color_manual(values = c("Conspecific BA" = "black", 
#                                 "Conmycorrhizal BA" = "#DA3742", 
#                                 "Heteromycorrhizal BA" = "#3885B6")) +
#   theme_bw()
# 
# 

########################################################################
########################################################################
##### Plot to visualize effects of conspecific, con-myco, and hetero-myco neighborhoods on the 1-yr seedling counts for each species together!
# Your model list
models <- list(rm = rmBA_15m, sm = smBA_15m, ab = abBA_15m, wa = waBA_15m,
               bc = bcBA_15m, nro = nroBA_15m, bonro = bonroBA_15m,
               wo = woBA_15m, iw = iwBA_15m)

# Your datasets (must match model order)
data_list <- list(rm = rm_surv, sm = sm_surv, ab = ab_surv, wa = wa_surv,
                  bc = bc_surv, nro = nro_surv, bonro = bonro_surv,
                  wo = wo_surv, iw = iw_surv)

# Your color vector
species_cols <- c(ab = "#B4D4DA", bonro = "#6EB6D9", iw = "#3F8BBA",
                  nro = "#1F74B1", wo = "#244D7F", bc = "#FEB1A3",
                  rm = "#F77964", sm = "#D2241F", wa = "#AE123A")

# Define conspecific BA variables per species
conspecific_vars <- c(rm = "rm_15m", sm = "sm_15m", ab = "ab_15m", wa = "wa_15m",
                      bc = "bc_15m", nro = "nro_15m", bonro = "bonro_15m",
                      wo = "wo_15m", iw = "iw_15m")


make_species_predictions <- function(focal_var) {
  pred_list <- list()
  
  for (sp in names(models)) {
    mod <- models[[sp]]
    dat <- data_list[[sp]]
    conspecific_var <- conspecific_vars[[sp]]
    
    # Define x-axis variable and interaction handling
    interacting_var <- NULL
    if (focal_var == "conspecific") {
      x_var <- conspecific_var
      x_label <- x_var
    } else if (focal_var == "conspecific:hetMyc_15m") {
      interaction_term <- paste0(conspecific_var, ":hetMyc_15m")
      if (!grepl(interaction_term, paste(deparse(formula(mod)), collapse = ""))) next
      x_var <- "hetMyc_15m"
      interacting_var <- conspecific_var
      x_label <- interaction_term
    } else if (focal_var == "conMycNS_15m:hetMyc_15m") {
      interaction_term <- "conMycNS_15m:hetMyc_15m"
      if (!grepl(interaction_term, paste(deparse(formula(mod)), collapse = ""))) next
      x_var <- "conMycNS_15m"
      interacting_var <- "hetMyc_15m"
      x_label <- interaction_term
    } else {
      x_var <- focal_var
      x_label <- focal_var
    }
    
    # Skip if focal or interacting variable missing
    if (!x_var %in% names(dat) || all(is.na(dat[[x_var]]))) next
    if (!is.null(interacting_var) && (!interacting_var %in% names(dat) || all(is.na(dat[[interacting_var]])))) next
    
    # Sequence of values for prediction
    max_val <- max(dat[[x_var]], na.rm = TRUE)
    if (!is.finite(max_val)) next
    seq_vals <- seq(0, max_val, length.out = 100)
    
    # Median covariates
    med_vals <- list(
      log0yr = median(dat$log0yr, na.rm = TRUE),
      logPrior = median(dat$logPrior, na.rm = TRUE),
      conMycNS_15m = median(dat$conMycNS_15m, na.rm = TRUE),
      hetMyc_15m = median(dat$hetMyc_15m, na.rm = TRUE),
      conspecific = median(dat[[conspecific_var]], na.rm = TRUE)
    )
    
    # Create base row for prediction data
    pred_data <- data.frame(
      log0yr = med_vals$log0yr,
      logPrior = med_vals$logPrior,
      conMycNS_15m = med_vals$conMycNS_15m,
      hetMyc_15m = med_vals$hetMyc_15m,
      stringsAsFactors = FALSE
    )
    pred_data[[conspecific_var]] <- med_vals$conspecific
    
    # Expand to prediction sequence
    newdata <- pred_data[rep(1, length(seq_vals)), ]
    
    # Modify variables for current prediction type
    if (focal_var == "conspecific:hetMyc_15m") {
      newdata[[x_var]] <- seq_vals
      newdata[[interacting_var]] <- med_vals$conspecific
    } else if (focal_var == "conMycNS_15m:hetMyc_15m") {
      newdata[[x_var]] <- seq_vals
      newdata[[interacting_var]] <- seq_vals
    } else {
      newdata[[x_var]] <- seq_vals
    }
    
    # Add dummy random effects values
    newdata$Site <- dat$Site[1]
    newdata$survYear <- dat$survYear[1]
    
    # Safe prediction
    newdata$pred_count <- suppressWarnings(predict(mod, newdata = newdata, type = "response", re.form = NA))
    
    # Store x value and species for plotting
    newdata$x_focal <- newdata[[x_var]]
    newdata$species <- sp
    pred_list[[sp]] <- newdata
  }
  
  all_preds <- bind_rows(pred_list)
  all_preds$species <- factor(all_preds$species, levels = names(species_cols))
  return(all_preds)
}

# Titles and focal variable mapping
plot_titles <- c(
  "Conspecific Basal Area", 
  "Conmycorrhizal Basal Area", 
  "Heteromycorrhizal Basal Area", 
  "Log(Prior Seedling Count + 1)",
  "Log(Recruit Seedling Count)", 
  "Conmyc × Heteromyc Basal Area", 
  "Conspecific × Heteromyc Basal Area"
)

focal_vars <- c(
  "conspecific", 
  "conMycNS_15m", 
  "hetMyc_15m", 
  "logPrior", 
  "log0yr", 
  "conMycNS_15m:hetMyc_15m", 
  "conspecific:hetMyc_15m"
)

# Plot each focal variable
for (i in seq_along(focal_vars)) {
  focal_var <- focal_vars[i]
  pred_df <- make_species_predictions(focal_var)
  
  if (nrow(pred_df) == 0) next  # skip empty predictions
  
  p <- ggplot(pred_df, aes(x = x_focal, y = pred_count, color = species)) +
    geom_line(size = 1) +
    labs(
      x = plot_titles[i],
      y = "Predicted 1-Year Seedling Count",
      color = "Species"
    ) +
    scale_color_manual(values = species_cols) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.position = "right"
    )
  
  print(p)
}

##########

# 1. Mycorrhizal type classification
am_species <- c("bc", "rm", "sm", "wa")  # example AM species codes
ecm_species <- c("ab", "bonro", "iw", "nro", "wo")  # example ECM species codes

# 2. Prediction function (unchanged)
make_species_predictions <- function(focal_var) {
  pred_list <- list()
  
  for (sp in names(models)) {
    mod <- models[[sp]]
    dat <- data_list[[sp]]
    conspecific_var <- conspecific_vars[[sp]]
    
    # Define x-axis variable and interaction handling
    interacting_var <- NULL
    if (focal_var == "conspecific") {
      x_var <- conspecific_var
      x_label <- x_var
    } else if (focal_var == "conspecific:hetMyc_15m") {
      interaction_term <- paste0(conspecific_var, ":hetMyc_15m")
      if (!grepl(interaction_term, paste(deparse(formula(mod)), collapse = ""))) next
      x_var <- "hetMyc_15m"
      interacting_var <- conspecific_var
      x_label <- interaction_term
    } else if (focal_var == "conMycNS_15m:hetMyc_15m") {
      interaction_term <- "conMycNS_15m:hetMyc_15m"
      if (!grepl(interaction_term, paste(deparse(formula(mod)), collapse = ""))) next
      x_var <- "conMycNS_15m"
      interacting_var <- "hetMyc_15m"
      x_label <- interaction_term
    } else {
      x_var <- focal_var
      x_label <- focal_var
    }
    
    if (!x_var %in% names(dat) || all(is.na(dat[[x_var]]))) next
    if (!is.null(interacting_var) && (!interacting_var %in% names(dat) || all(is.na(dat[[interacting_var]])))) next
    
    max_val <- max(dat[[x_var]], na.rm = TRUE)
    if (!is.finite(max_val)) next
    seq_vals <- seq(0, max_val, length.out = 100)
    
    med_vals <- list(
      log0yr = median(dat$log0yr, na.rm = TRUE),
      logPrior = median(dat$logPrior, na.rm = TRUE),
      conMycNS_15m = median(dat$conMycNS_15m, na.rm = TRUE),
      hetMyc_15m = median(dat$hetMyc_15m, na.rm = TRUE),
      conspecific = median(dat[[conspecific_var]], na.rm = TRUE)
    )
    
    pred_data <- data.frame(
      log0yr = med_vals$log0yr,
      logPrior = med_vals$logPrior,
      conMycNS_15m = med_vals$conMycNS_15m,
      hetMyc_15m = med_vals$hetMyc_15m,
      stringsAsFactors = FALSE
    )
    pred_data[[conspecific_var]] <- med_vals$conspecific
    
    newdata <- pred_data[rep(1, length(seq_vals)), ]
    
    if (focal_var == "conspecific:hetMyc_15m") {
      newdata[[x_var]] <- seq_vals
      newdata[[interacting_var]] <- med_vals$conspecific
    } else if (focal_var == "conMycNS_15m:hetMyc_15m") {
      newdata[[x_var]] <- seq_vals
      newdata[[interacting_var]] <- seq_vals
    } else {
      newdata[[x_var]] <- seq_vals
    }
    
    newdata$Site <- dat$Site[1]
    newdata$survYear <- dat$survYear[1]
    newdata$pred_count <- suppressWarnings(predict(mod, newdata = newdata, type = "response", re.form = NA))
    newdata$x_focal <- newdata[[x_var]]
    newdata$species <- sp
    pred_list[[sp]] <- newdata
  }
  
  all_preds <- bind_rows(pred_list)
  all_preds$species <- factor(all_preds$species, levels = names(species_cols))
  return(all_preds)
}

# 3. Plot labels and variables
plot_titles <- c(
  "Conspecific Basal Area", 
  "Conmycorrhizal Basal Area", 
  "Heteromycorrhizal Basal Area", 
  "Log(Prior Seedling Count + 1)",
  "Log(Recruit Seedling Count)", 
  "Conmyc × Heteromyc Basal Area", 
  "Conspecific × Heteromyc Basal Area"
)

focal_vars <- c(
  "conspecific", 
  "conMycNS_15m", 
  "hetMyc_15m", 
  "logPrior", 
  "log0yr", 
  "conMycNS_15m:hetMyc_15m", 
  "conspecific:hetMyc_15m"
)

# 4. Plotting loop (AM and ECM separated)
for (i in seq_along(focal_vars)) {
  focal_var <- focal_vars[i]
  pred_df <- make_species_predictions(focal_var)
  
  if (nrow(pred_df) == 0) next
  
  # AM species plot
  am_df <- pred_df %>% filter(species %in% am_species)
  if (nrow(am_df) > 0) {
    p_am <- ggplot(am_df, aes(x = x_focal, y = pred_count, color = species)) +
      geom_line(size = 1) +
      labs(
        x = plot_titles[i],
        y = "Predicted 1-Year Seedling Count",
        title = paste0(plot_titles[i], " — AM Species"),
        color = "Species"
      ) +
      scale_color_manual(values = species_cols[am_species]) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.position = "right"
      )
    print(p_am)
  }
  
  # ECM species plot
  ecm_df <- pred_df %>% filter(species %in% ecm_species)
  if (nrow(ecm_df) > 0) {
    p_ecm <- ggplot(ecm_df, aes(x = x_focal, y = pred_count, color = species)) +
      geom_line(size = 1) +
      labs(
        x = plot_titles[i],
        y = "Predicted 1-Year Seedling Count",
        title = paste0(plot_titles[i], " — ECM Species"),
        color = "Species"
      ) +
      scale_color_manual(values = species_cols[ecm_species]) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.position = "right"
      )
    print(p_ecm)
  }
}


######################################################


# Main function to generate predictions per focal variable
make_species_predictions <- function(focal_var) {
  pred_list <- list()
  
  for (sp in names(models)) {
    mod <- models[[sp]]
    dat <- data_list[[sp]]
    conspecific_var <- conspecific_vars[[sp]]
    
    # Set actual variable name to use
    var_map <- list(
      conspecific = conspecific_var,
      conMycNS_15m = "conMycNS_15m",
      hetMyc_15m = "hetMyc_15m",
      logPrior = "logPrior"
    )
    var_name <- var_map[[focal_var]]
    
    # Skip if that variable isn't present or all NA
    if (!var_name %in% names(dat) || all(is.na(dat[[var_name]]))) next
    
    # Create sequence for focal variable
    max_val <- max(dat[[var_name]], na.rm = TRUE)
    if (!is.finite(max_val)) next  # skip species with no usable data for var
    
    seq_vals <- seq(0, max_val, length.out = 100)
    
    # Median values for covariates
    med_vals <- list(
      log0yr = median(dat$log0yr, na.rm = TRUE),
      logPrior = median(dat$logPrior, na.rm = TRUE),
      conMycNS_15m = median(dat$conMycNS_15m, na.rm = TRUE),
      hetMyc_15m = median(dat$hetMyc_15m, na.rm = TRUE),
      conspecific = median(dat[[conspecific_var]], na.rm = TRUE)
    )
    
    # Build base prediction data
    pred_data <- data.frame(
      log0yr = med_vals$log0yr,
      logPrior = med_vals$logPrior,
      conMycNS_15m = med_vals$conMycNS_15m,
      hetMyc_15m = med_vals$hetMyc_15m,
      stringsAsFactors = FALSE
    )
    pred_data[[conspecific_var]] <- med_vals$conspecific
    
    # Repeat rows for prediction sequence
    newdata <- pred_data[rep(1, length(seq_vals)), ]
    newdata[[var_name]] <- seq_vals
    newdata$Site <- dat$Site[1]
    newdata$survYear <- dat$survYear[1]
    
    # Predict
    newdata$pred_count <- predict(mod, newdata = newdata, type = "response", re.form = NA)
    newdata$x_focal <- seq_vals
    newdata$species <- sp
    pred_list[[sp]] <- newdata
  }
  
  bind_rows(pred_list)
  
  all_preds <- bind_rows(pred_list)
  all_preds$species <- factor(all_preds$species, levels = names(species_cols))
  return(all_preds)
  
}


# Labels for the x-axis
plot_titles <- c("Conspecific Basal Area", "Conmycorrhizal Basal Area",
                 "Heteromycorrhizal Basal Area", "Log(Prior Seedling Count + 1)")
focal_vars <- c("conspecific", "conMycNS_15m", "hetMyc_15m", "logPrior")

for (i in seq_along(focal_vars)) {
  focal_var <- focal_vars[i]
  pred_df <- make_species_predictions(focal_var)
  
  p <- ggplot(pred_df, aes(x = x_focal, y = pred_count, color = species)) +
    geom_line(size = 1) +
    labs(x = plot_titles[i],
         y = "Predicted 1-Year Seedling Count",
         color = "Species") +
    scale_color_manual(values = species_cols) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
  
  print(p)
}


