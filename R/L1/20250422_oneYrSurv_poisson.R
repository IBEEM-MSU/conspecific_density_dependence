
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
# (#1 yr) ~ Consp BA:log(#0 yr) + % EMF hosts:log(#0 yr) + Con-myc BA:log(#0 yr) + log(hetSp_Current_Prior3 + 1):log(#0 yr) + (1 | Site) + (1 | Year)

# For species with sufficient observations, we included a random effect of site. For species with limited observations across sites, we excluded this term due to singular fits and lack of variance attributable to site.

##### Read in file that has the counts by species for every TransectID-ObsYear combination
allSp_TI_ObsYear <- read.csv("./data/L1/allSp_prior_cohort_densities.csv")

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

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
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
rm_surv$percEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Log transform the count of 0-yr old seedlings
rm_surv$log0yr <- log(rm_surv$total_count_yr0)  
 
#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
rm_mod_15m <- glmer(total_count_yr1 ~ rm_15m:log0yr + 
              percEMF_15m:log0yr + 
              conMycNS_15m:log0yr + 
              log(hetSp_Prior3_counts + 1):log0yr +
              (1 | Site) + (1 | survYear),
    family = poisson(link = "log"), data = rm_surv)
summary(rm_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
rm_mod_15m1 <- glmer(total_count_yr1 ~ rm_15m + 
                       percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = rm_surv)
summary(rm_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
rm_mod_15m2 <- glmer(total_count_yr1 ~ 
                rm_15m + percEMF_15m + 
                 conMycNS_15m + 
                 log(hetSp_Prior3_counts + 1) +
                 #rm_15m:percEMF_15m + 
                 rm_15m:log0yr + 
                 percEMF_15m:log0yr + 
                 conMycNS_15m:log0yr + 
                 log(hetSp_Prior3_counts + 1):log0yr +
                 (1 | Site) + (1 | survYear),
  family = poisson(link = "log"), data = rm_surv)
summary(rm_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
rm_mod_15m3 <- glmer(total_count_yr1 ~ 
                       rm_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       rm_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = rm_surv)
summary(rm_mod_15m3)

# Variance inflation factors 
car::vif(rm_mod_15m) # all < 2
car::vif(rm_mod_15m1) # all < 4
car::vif(rm_mod_15m2) # all < 6
car::vif(rm_mod_15m3) # all < 6
r.squaredGLMM(rm_mod_15m) 
r.squaredGLMM(rm_mod_15m1) # very close 2nd
r.squaredGLMM(rm_mod_15m2) # highest r2 values
r.squaredGLMM(rm_mod_15m3) 
## AIC
aic_values <- AIC(rm_mod_15m, rm_mod_15m1, rm_mod_15m2, rm_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # rm_mod_15m1 is best

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


#################################################
#################################################
##### Sugar Maple

### Filter to include 2001 as first First.Year
sm_surv_temp <- dplyr::filter(sm_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for sm
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5,7:14,16,18:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
sm_surv$percEMF_15m = (sm_surv$hetMyc_15m / (sm_surv$hetMyc_15m + sm_surv$conMyc_15m))

#### Make site and year into factors
sm_surv$Site = as.factor(sm_surv$Site)
sm_surv$survYear = as.factor(sm_surv$survYear)

#### Log transform the count of 0-yr old seedlings
sm_surv$log0yr <- log(sm_surv$total_count_yr0)  

#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
sm_mod_15m <- glmer(total_count_yr1 ~ sm_15m:log0yr + 
                      percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = sm_surv)
summary(sm_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
sm_mod_15m1 <- glmer(total_count_yr1 ~ sm_15m + 
                       percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = sm_surv)
summary(sm_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
sm_mod_15m2 <- glmer(total_count_yr1 ~ 
                       sm_15m + percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       sm_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = sm_surv)
summary(sm_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
sm_mod_15m3 <- glmer(total_count_yr1 ~ 
                       sm_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       sm_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = sm_surv)
summary(sm_mod_15m3)

# Variance inflation factors 
car::vif(sm_mod_15m) # all < 2
car::vif(sm_mod_15m1) # all < 2
car::vif(sm_mod_15m2) # all < 8
car::vif(sm_mod_15m3) # all < 7
r.squaredGLMM(sm_mod_15m) 
r.squaredGLMM(sm_mod_15m1) # highest r2 value - fixed
r.squaredGLMM(sm_mod_15m2) # highest r2 value - total
r.squaredGLMM(sm_mod_15m3) 
## AIC
aic_values <- AIC(sm_mod_15m, sm_mod_15m1, sm_mod_15m2, sm_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # sm_mod_15m1 is best

#################################################
#################################################
##### American Beech

### Filter to include 2001 as first First.Year
ab_surv_temp <- dplyr::filter(ab_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for ab
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:6,8:14,16:17,19:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
ab_surv$percEMF_15m = (ab_surv$conMyc_15m / (ab_surv$hetMyc_15m + ab_surv$conMyc_15m))

#### Make site and year into factors
ab_surv$Site = as.factor(ab_surv$Site)
ab_surv$survYear = as.factor(ab_surv$survYear)

#### Log transform the count of 0-yr old seedlings
ab_surv$log0yr <- log(ab_surv$total_count_yr0)  

#################### 15-m radius - NO SITE EFFECT
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
ab_mod_15m <- glmer(total_count_yr1 ~ ab_15m:log0yr + 
                      percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | survYear),
                    family = poisson(link = "log"), data = ab_surv)
summary(ab_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
ab_mod_15m1 <- glmer(total_count_yr1 ~ ab_15m + 
                       percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = ab_surv)
summary(ab_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
ab_mod_15m2 <- glmer(total_count_yr1 ~ 
                       ab_15m + percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       ab_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = ab_surv)
summary(ab_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
ab_mod_15m3 <- glmer(total_count_yr1 ~ 
                       ab_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       ab_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = ab_surv)
summary(ab_mod_15m3)

# Variance inflation factors 
car::vif(ab_mod_15m) # all < 2
car::vif(ab_mod_15m1) # all < 9
car::vif(ab_mod_15m2) # all < 8
car::vif(ab_mod_15m3) # all < 7
r.squaredGLMM(ab_mod_15m) 
r.squaredGLMM(ab_mod_15m1) # highest r2 values
r.squaredGLMM(ab_mod_15m2) 
r.squaredGLMM(ab_mod_15m3) 
## AIC
aic_values <- AIC(ab_mod_15m, ab_mod_15m1, ab_mod_15m2, ab_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # ab_mod_15m1 is best


#################################################
#################################################
##### White ash

### Filter to include 2001 as first First.Year
wa_surv_temp <- dplyr::filter(wa_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for wa
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:7,9:14,16:18,20:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
wa_surv$percEMF_15m = (wa_surv$hetMyc_15m / (wa_surv$hetMyc_15m + wa_surv$conMyc_15m))

#### Make site and year into factors
wa_surv$Site = as.factor(wa_surv$Site)
wa_surv$survYear = as.factor(wa_surv$survYear)

#### Log transform the count of 0-yr old seedlings
wa_surv$log0yr <- log(wa_surv$total_count_yr0)  

#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
wa_mod_15m <- glmer(total_count_yr1 ~ wa_15m:log0yr + 
                      percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = wa_surv)
summary(wa_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
wa_mod_15m1 <- glmer(total_count_yr1 ~ wa_15m + 
                       percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wa_surv)
summary(wa_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
wa_mod_15m2 <- glmer(total_count_yr1 ~ 
                       wa_15m + percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       wa_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wa_surv)
summary(wa_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
wa_mod_15m3 <- glmer(total_count_yr1 ~ 
                       wa_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       wa_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wa_surv)
summary(wa_mod_15m3)

# Variance inflation factors 
car::vif(wa_mod_15m) # all < 2
car::vif(wa_mod_15m1) # all < 3
car::vif(wa_mod_15m2) # all < 11
car::vif(wa_mod_15m3) # all < 6
r.squaredGLMM(wa_mod_15m) 
r.squaredGLMM(wa_mod_15m1) # next highest, no issues
r.squaredGLMM(wa_mod_15m2) # highest r2 values, but convergence issues
r.squaredGLMM(wa_mod_15m3) 
## AIC
aic_values <- AIC(wa_mod_15m, wa_mod_15m1, wa_mod_15m2, wa_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # wa_mod_15m1 is best


#################################################
#################################################
##### Black Cherry

### Filter to include 2001 as first First.Year
bc_surv_temp <- dplyr::filter(bc_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for bc
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:8,10:14,16:19,21:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
bc_surv$percEMF_15m = (bc_surv$hetMyc_15m / (bc_surv$hetMyc_15m + bc_surv$conMyc_15m))

#### Make site and year into factors
bc_surv$Site = as.factor(bc_surv$Site)
bc_surv$survYear = as.factor(bc_surv$survYear)

#### Log transform the count of 0-yr old seedlings
bc_surv$log0yr <- log(bc_surv$total_count_yr0)  

#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
bc_mod_15m <- glmer(total_count_yr1 ~ bc_15m:log0yr + 
                      percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = bc_surv)
summary(bc_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
bc_mod_15m1 <- glmer(total_count_yr1 ~ bc_15m + 
                        percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = bc_surv)
summary(bc_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
bc_mod_15m2 <- glmer(total_count_yr1 ~ 
                       bc_15m + percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       bc_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = bc_surv)
summary(bc_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
bc_mod_15m3 <- glmer(total_count_yr1 ~ 
                       bc_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       bc_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = bc_surv)
summary(bc_mod_15m3)

# Variance inflation factors 
car::vif(bc_mod_15m) # all < 4
car::vif(bc_mod_15m1) # all < 3
car::vif(bc_mod_15m2) # all < 6
car::vif(bc_mod_15m3) # all < 6
r.squaredGLMM(bc_mod_15m) # highest r2 values
r.squaredGLMM(bc_mod_15m1)
r.squaredGLMM(bc_mod_15m2) 
r.squaredGLMM(bc_mod_15m3) 
## AIC
aic_values <- AIC(bc_mod_15m, bc_mod_15m1, bc_mod_15m2, bc_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # bc_mod_15m3 is best, but not major differences between the models


#################################################
#################################################
##### Northern Red Oak

### Filter to include 2001 as first First.Year
nro_surv_temp <- dplyr::filter(nro_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for nro
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:9,11:14,16:20,22:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
nro_surv$percEMF_15m = (nro_surv$conMyc_15m / (nro_surv$hetMyc_15m + nro_surv$conMyc_15m))

#### Make site and year into factors
nro_surv$Site = as.factor(nro_surv$Site)
nro_surv$survYear = as.factor(nro_surv$survYear)

#### Log transform the count of 0-yr old seedlings
nro_surv$log0yr <- log(nro_surv$total_count_yr0)  

#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
nro_mod_15m <- glmer(total_count_yr1 ~ nro_15m:log0yr + 
                      percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | Site) + (1 | survYear),
                    family = poisson(link = "log"), data = nro_surv)
summary(nro_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
nro_mod_15m1 <- glmer(total_count_yr1 ~ nro_15m + 
                       percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = nro_surv)
summary(nro_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
nro_mod_15m2 <- glmer(total_count_yr1 ~ 
                       nro_15m + percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       nro_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = nro_surv)
summary(nro_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
nro_mod_15m3 <- glmer(total_count_yr1 ~ 
                       nro_15m:percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       nro_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = nro_surv)
summary(nro_mod_15m3)

# Variance inflation factors 
car::vif(nro_mod_15m) # multiple >18
car::vif(nro_mod_15m1) # all <5
car::vif(nro_mod_15m2) # multiple >30
car::vif(nro_mod_15m3) # multiple >24
r.squaredGLMM(nro_mod_15m) 
r.squaredGLMM(nro_mod_15m1) # highest r2
r.squaredGLMM(nro_mod_15m2) 
r.squaredGLMM(nro_mod_15m3) 
## AIC
aic_values <- AIC(nro_mod_15m, nro_mod_15m1, nro_mod_15m2, nro_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # nro_mod_15m1 is best


#################################################
#################################################
##### Black oak

### Filter to include 2001 as first First.Year
bonro_surv_temp <- dplyr::filter(bonro_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for bonro
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:10,12:14,16:21,23:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
bonro_surv$percEMF_15m = (bonro_surv$conMyc_15m / (bonro_surv$hetMyc_15m + bonro_surv$conMyc_15m))

#### Make site and year into factors
bonro_surv$Site = as.factor(bonro_surv$Site)
bonro_surv$survYear = as.factor(bonro_surv$survYear)

#### Log transform the count of 0-yr old seedlings
bonro_surv$log0yr <- log(bonro_surv$total_count_yr0)  

#################### 15-m radius - removing Site RE
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
bonro_mod_15m <- glmer(total_count_yr1 ~ bonro_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = bonro_surv)
summary(bonro_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
bonro_mod_15m1 <- glmer(total_count_yr1 ~ bonro_15m + 
                        percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        log0yr +
                        (1 | survYear),
                      family = poisson(link = "log"), data = bonro_surv)
summary(bonro_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
bonro_mod_15m2 <- glmer(total_count_yr1 ~ 
                        bonro_15m + percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        bonro_15m:log0yr + 
                        percEMF_15m:log0yr + 
                        conMycNS_15m:log0yr + 
                        log(hetSp_Prior3_counts + 1):log0yr +
                        (1 | survYear),
                      family = poisson(link = "log"), data = bonro_surv)
summary(bonro_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
bonro_mod_15m3 <- glmer(total_count_yr1 ~ 
                        bonro_15m:percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        bonro_15m:log0yr + 
                        percEMF_15m:log0yr + 
                        conMycNS_15m:log0yr + 
                        log(hetSp_Prior3_counts + 1):log0yr +
                        (1 | survYear),
                      family = poisson(link = "log"), data = bonro_surv)
summary(bonro_mod_15m3)

# Variance inflation factors 
car::vif(bonro_mod_15m) # multiple >22
car::vif(bonro_mod_15m1) # all < 3
car::vif(bonro_mod_15m2) # multiple >24
car::vif(bonro_mod_15m3) # multiple >24
r.squaredGLMM(bonro_mod_15m) 
r.squaredGLMM(bonro_mod_15m1) 
r.squaredGLMM(bonro_mod_15m2) # highest r2, but convergence issues
r.squaredGLMM(bonro_mod_15m3) 
## AIC
aic_values <- AIC(bonro_mod_15m, bonro_mod_15m1, bonro_mod_15m2, bonro_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # bonro_mod_15m is best, but all close


#################################################
#################################################
##### White oak

### Filter to include 2001 as first First.Year
wo_surv_temp <- dplyr::filter(wo_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for wo
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:11,13:14,16:22,24:25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
wo_surv$percEMF_15m = (wo_surv$conMyc_15m / (wo_surv$hetMyc_15m + wo_surv$conMyc_15m))

#### Make site and year into factors
wo_surv$Site = as.factor(wo_surv$Site)
wo_surv$survYear = as.factor(wo_surv$survYear)

#### Log transform the count of 0-yr old seedlings
wo_surv$log0yr <- log(wo_surv$total_count_yr0)  

#################### 15-m radius
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
wo_mod_15m <- glmer(total_count_yr1 ~ wo_15m:log0yr + 
                       percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | Site) + (1 | survYear),
                     family = poisson(link = "log"), data = wo_surv)
summary(wo_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
wo_mod_15m1 <- glmer(total_count_yr1 ~ wo_15m + 
                        percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = wo_surv)
summary(wo_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
wo_mod_15m2 <- glmer(total_count_yr1 ~ 
                        wo_15m + percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        wo_15m:log0yr + 
                        percEMF_15m:log0yr + 
                        conMycNS_15m:log0yr + 
                        log(hetSp_Prior3_counts + 1):log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = wo_surv)
summary(wo_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
wo_mod_15m3 <- glmer(total_count_yr1 ~ 
                        wo_15m:percEMF_15m + 
                        conMycNS_15m + 
                        log(hetSp_Prior3_counts + 1) +
                        wo_15m:log0yr + 
                        percEMF_15m:log0yr + 
                        conMycNS_15m:log0yr + 
                        log(hetSp_Prior3_counts + 1):log0yr +
                        (1 | Site) + (1 | survYear),
                      family = poisson(link = "log"), data = wo_surv)
summary(wo_mod_15m3)

# Variance inflation factors 
car::vif(wo_mod_15m) # multiple >17
car::vif(wo_mod_15m1) # all <2
car::vif(wo_mod_15m2) # multiple >34
car::vif(wo_mod_15m3) # multiple >29
r.squaredGLMM(wo_mod_15m) 
r.squaredGLMM(wo_mod_15m1) 
r.squaredGLMM(wo_mod_15m2) 
r.squaredGLMM(wo_mod_15m3) # highest r2, but convergence issues and huge VIFs
## AIC
aic_values <- AIC(wo_mod_15m, wo_mod_15m1, wo_mod_15m2, wo_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # wo_mod_15m and wo_mod_15m1 are tied for best, but all are within 3 units


#################################################
#################################################
##### Ironwood

### Filter to include 2001 as first First.Year
iw_surv_temp <- dplyr::filter(iw_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for iw
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:12,14,16:23,25)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
iw_surv$percEMF_15m = (iw_surv$conMyc_15m / (iw_surv$hetMyc_15m + iw_surv$conMyc_15m))

#### Make site and year into factors
iw_surv$Site = as.factor(iw_surv$Site)
iw_surv$survYear = as.factor(iw_surv$survYear)

#### Log transform the count of 0-yr old seedlings
iw_surv$log0yr <- log(iw_surv$total_count_yr0)  

#################### 15-m radius - no RE of site, and excluding %emf
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
iw_mod_15m <- glmer(total_count_yr1 ~ iw_15m:log0yr + 
                      #percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | survYear),
                    family = poisson(link = "log"), data = iw_surv)
summary(iw_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
iw_mod_15m1 <- glmer(total_count_yr1 ~ iw_15m + 
                       #percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr + (1 | survYear),
                     family = poisson(link = "log"), data = iw_surv)
summary(iw_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
iw_mod_15m2 <- glmer(total_count_yr1 ~ 
                       iw_15m + #percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       iw_15m:log0yr + 
                       #percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = iw_surv)
summary(iw_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
iw_mod_15m3 <- glmer(total_count_yr1 ~ 
                       iw_15m:percEMF_15m + 
                       #conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       iw_15m:log0yr + 
                       #percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = iw_surv)
summary(iw_mod_15m3)

# Variance inflation factors 
car::vif(iw_mod_15m) # multiple >18
car::vif(iw_mod_15m1) # multiple >14
car::vif(iw_mod_15m2) # multiple >19
car::vif(iw_mod_15m3) # multiple >19
r.squaredGLMM(iw_mod_15m) 
r.squaredGLMM(iw_mod_15m1) # highest r2, but some large VIFs - go away if exclude %emf
r.squaredGLMM(iw_mod_15m2) 
r.squaredGLMM(iw_mod_15m3) 
## AIC
aic_values <- AIC(iw_mod_15m, iw_mod_15m1, iw_mod_15m2, iw_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # iw_mod_15m1 is best


#################################################
#################################################
##### American Basswood

### Filter to include 2001 as first First.Year
bs_surv_temp <- dplyr::filter(bs_surv_temp, First.Year > 2000)

# Add in conspecific and heterospecific species previous 3-year counts into the one-year survivorship data
# Calculate sum of heterospecific species previous 3-year and current year counts for bs
allSp_TI_ObsYear$hetSp_Prior3_counts = 
  rowSums(allSp_TI_ObsYear[,c(5:13,16:24)], na.rm = TRUE)
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

#### Calculate the proportion of EMF hosts in 15-m radius
bs_surv$percEMF_15m = (bs_surv$conMyc_15m / (bs_surv$hetMyc_15m + bs_surv$conMyc_15m))

#### Make site and year into factors
bs_surv$Site = as.factor(bs_surv$Site)
bs_surv$survYear = as.factor(bs_surv$survYear)

#### Log transform the count of 0-yr old seedlings
bs_surv$log0yr <- log(bs_surv$total_count_yr0)  

#################### 15-m radius - no RE of site, and excluding %emf
### NONE OF THE BASSWOOD MODELS ARE CONVERGING
## Just interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr
bs_mod_15m <- glmer(total_count_yr1 ~ bs_15m:log0yr + 
                      #percEMF_15m:log0yr + 
                      conMycNS_15m:log0yr + 
                      log(hetSp_Prior3_counts + 1):log0yr +
                      (1 | survYear),
                    family = poisson(link = "log"), data = bs_surv)
summary(bs_mod_15m)

## Just main effects - consp, %EMF, conmyc, hetspCounts, and log0yr
bs_mod_15m1 <- glmer(total_count_yr1 ~ bs_15m + 
                       #percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       log0yr + (1 | survYear),
                     family = poisson(link = "log"), data = bs_surv)
summary(bs_mod_15m1)

## Version 2 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of consp, %EMF, conmyc, and hetspCounts
bs_mod_15m2 <- glmer(total_count_yr1 ~ 
                       bs_15m + #percEMF_15m + 
                       conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       bs_15m:log0yr + 
                       #percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = bs_surv)
summary(bs_mod_15m2)

## Version 3 - interactions between consp:log0yr, %EMF:log0yr, conmyc:log0yr, and hetspCounts:log0yr, plus main effects of conmyc and hetspCounts, and consp:%EMF
bs_mod_15m3 <- glmer(total_count_yr1 ~ 
                       bs_15m:percEMF_15m + 
                       #conMycNS_15m + 
                       log(hetSp_Prior3_counts + 1) +
                       bs_15m:log0yr + 
                       #percEMF_15m:log0yr + 
                       conMycNS_15m:log0yr + 
                       log(hetSp_Prior3_counts + 1):log0yr +
                       (1 | survYear),
                     family = poisson(link = "log"), data = bs_surv)
summary(bs_mod_15m3)


# Variance inflation factors 
car::vif(bs_mod_15m) # multiple >18
car::vif(bs_mod_15m1) # multiple >14
car::vif(bs_mod_15m2) # multiple >19
car::vif(bs_mod_15m3) # multiple >19
r.squaredGLMM(bs_mod_15m) 
r.squaredGLMM(bs_mod_15m1) # highest r2, but some large VIFs - go away if exclude %emf
r.squaredGLMM(bs_mod_15m2) 
r.squaredGLMM(bs_mod_15m3) 
## AIC
aic_values <- AIC(bs_mod_15m, bs_mod_15m1, bs_mod_15m2, bs_mod_15m3)
# Sort AIC values from lowest to highest and show
print(aic_values[order(aic_values$AIC), ]) # bs_mod_15m1 is best

