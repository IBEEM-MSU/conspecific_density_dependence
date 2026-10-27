
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

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1

#########################################################
#########################################################
################# All Species Models - Raw Basal Area

## Null model
all_sp_15m_null <- glmer(Year_1_counts ~ 1 +
                      (1|Site) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_15m_null)
r.squaredGLMM(all_sp_15m_null) 

## Testing Panel A conceptual figure
all_sp_15m <- glmer(Year_1_counts ~ consp_15m*Seedling_myc +
                       logPrior + log0yr +
                       (1|Site) + (1|survYear) + (1|Sp),
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
                       (1|Site) + (1|survYear) + (1|Sp),
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
                       (1|Site) + (1|survYear) + (1|Sp),
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
                       (1|Site) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m3) 
car::vif(all_sp_15m3) ## Some VIFs >10
r.squaredGLMM(all_sp_15m3) 

# ## Testing Panel D conceptual figure - log0yr as RE
# all_sp_15m3b <- glmer(Year_1_counts ~ consp_15m + Seedling_myc + 
#                        AM_15m + ECM_15m + 
#                        consp_15m:Seedling_myc:AM_15m + # Panel D
#                        consp_15m:Seedling_myc:ECM_15m + # Panel D
#                        logPrior + #log0yr +
#                        (1|Site) + (1|survYear) + (1|Sp) + (1|log0yr),
#                      family = poisson(link = "log"), data = allSpDat)
# summary(all_sp_15m3b) 
# car::vif(all_sp_15m3b) ## Some VIFs >10
# r.squaredGLMM(all_sp_15m3b) 

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
                    (1|Site) + (1|survYear) + (1|Sp),
                  family = poisson(link = "log"), data = allSpDat)
summary(all_sp_15m4)
car::vif(all_sp_15m4)
r.squaredGLMM(all_sp_15m4) 

## Model comparison
### Compare AICs - lowest is for all_sp_15m4 but all_sp_15m3 close 2nd
AIC(all_sp_15m_null, all_sp_15m, all_sp_15m1, all_sp_15m2, all_sp_15m3, all_sp_15m4)
# More complex model fits better
anova(all_sp_15m3, all_sp_15m4) # Global model significantly better than Panel D

#########################################################
################# All Species Models - Proportional Basal Area

## Null model
all_sp_prop_null <- glmer(Year_1_counts ~ 1 +
                           (1|Site) + (1|survYear) + (1|Sp),
                         family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_prop_null)
r.squaredGLMM(all_sp_prop_null) 

## Testing Panel A conceptual figure (proportional BA)
all_sp_prop <- glmer(Year_1_counts ~ propConsp*Seedling_myc +
                      logPrior + log0yr +
                      (1|Site) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat) 
summary(all_sp_prop)
car::vif(all_sp_prop)
r.squaredGLMM(all_sp_prop) 

## Testing Panel B conceptual figure (proportional BA)
all_sp_prop1 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                       propAM + # Panel B
                       propConsp:propAM + # Panel B
                       logPrior + log0yr +
                       (1|Site) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop1)
car::vif(all_sp_prop1)
r.squaredGLMM(all_sp_prop1) 

## Testing Panel C conceptual figure (proportional BA)
all_sp_prop2 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                        propAM + 
                        Seedling_myc:propAM + # Panel C
                        logPrior + log0yr +
                        (1|Site) + (1|survYear) + (1|Sp),
                      family = poisson(link = "log"), data = allSpDat)
summary(all_sp_prop2)
car::vif(all_sp_prop2)
r.squaredGLMM(all_sp_prop2) 

## Testing Panel D conceptual figure (proportional BA)
all_sp_prop3 <- glmer(Year_1_counts ~ propConsp + Seedling_myc + 
                        propAM + 
                        propConsp:Seedling_myc:propAM + # Panel D
                        logPrior + log0yr +
                        (1|Site) + (1|survYear) + (1|Sp),
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
                       (1|Site) + (1|survYear) + (1|Sp),
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


#########################################################
#########################################################
################# Species-Level Test - RM - Raw Basal Area

##### Subset RM data
rm_surv <- subset(allSpDat, Sp == "rm")

## Null model - species-level
rm_15m_null <- glmer(Year_1_counts ~ 1 +
                     (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv) 
summary(rm_15m_null)
r.squaredGLMM(rm_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
rm_15m <- glmer(Year_1_counts ~ consp_15m +
                      logPrior + log0yr +
                      (1|Site) + (1|survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_15m)
car::vif(rm_15m)
r.squaredGLMM(rm_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
rm_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv)
summary(rm_15m2)
car::vif(rm_15m2)
r.squaredGLMM(rm_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for rm_15m4 but rm_15m3 close 2nd
AIC(rm_15m_null, rm_15m, rm_15m1, rm_15m2)
# More complex model is NOT better
anova(rm_15m1, rm_15m2)

################# Species-Level Test - RM - Proportional Basal Area

## Null model - species-level (proportional BA)
rm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = rm_surv) 
summary(rm_15m_prop_null)
r.squaredGLMM(rm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
rm_prop <- glmer(Year_1_counts ~ propConsp +
                  logPrior + log0yr +
                  (1|Site) + (1|survYear),
                family = poisson(link = "log"), data = rm_surv)
summary(rm_prop)
car::vif(rm_prop)
r.squaredGLMM(rm_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
rm_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = rm_surv)
summary(rm_prop1)
car::vif(rm_prop1)
r.squaredGLMM(rm_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
rm_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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

#########################################################
################# Species-Level Test - SM - Raw Basal Area

##### Subset sm data
sm_surv <- subset(allSpDat, Sp == "sm")

## Null model - species-level
sm_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = sm_surv) 
summary(sm_15m_null)
r.squaredGLMM(sm_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
sm_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|Site) + (1|survYear),
                family = poisson(link = "log"), data = sm_surv)
summary(sm_15m)
car::vif(sm_15m)
r.squaredGLMM(sm_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
sm_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_15m2)
car::vif(sm_15m2)
r.squaredGLMM(sm_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for sm_15m4 but sm_15m3 close 2nd
AIC(sm_15m_null, sm_15m, sm_15m1, sm_15m2)
# More complex model is NOT better
anova(sm_15m1, sm_15m2)

################# Species-Level Test - sm - Proportional Basal Area

## Null model - species-level (proportional BA)
sm_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|Site) + (1|survYear),
                          family = poisson(link = "log"), data = sm_surv) 
summary(sm_15m_prop_null)
r.squaredGLMM(sm_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
sm_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = sm_surv)
summary(sm_prop)
car::vif(sm_prop)
r.squaredGLMM(sm_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
sm_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
                  family = poisson(link = "log"), data = sm_surv)
summary(sm_prop1)
car::vif(sm_prop1)
r.squaredGLMM(sm_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
sm_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
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

#########################################################
################# Species-Level Test - bc - Raw Basal Area

##### Subset bc data
bc_surv <- subset(allSpDat, Sp == "bc")

## Null model - species-level
bc_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_null)
r.squaredGLMM(bc_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
bc_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|Site) + (1|survYear),
                family = poisson(link = "log"), data = bc_surv)
summary(bc_15m)
car::vif(bc_15m)
r.squaredGLMM(bc_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
bc_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_15m2)
car::vif(bc_15m2)
r.squaredGLMM(bc_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for bc_15m4 but bc_15m3 close 2nd
AIC(bc_15m_null, bc_15m, bc_15m1, bc_15m2)
# More complex model is NOT better than bc_15m
anova(bc_15m, bc_15m2)

################# Species-Level Test - bc - Proportional Basal Area

## Null model - species-level (proportional BA)
bc_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|Site) + (1|survYear),
                          family = poisson(link = "log"), data = bc_surv) 
summary(bc_15m_prop_null)
r.squaredGLMM(bc_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
bc_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = bc_surv)
summary(bc_prop)
car::vif(bc_prop)
r.squaredGLMM(bc_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
bc_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
                  family = poisson(link = "log"), data = bc_surv)
summary(bc_prop1)
car::vif(bc_prop1)
r.squaredGLMM(bc_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
bc_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
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

#########################################################
################# Species-Level Test - wa - Raw Basal Area

##### Subset wa data
wa_surv <- subset(allSpDat, Sp == "wa")

## Null model - species-level
wa_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = wa_surv) 
summary(wa_15m_null)
r.squaredGLMM(wa_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
wa_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|Site) + (1|survYear),
                family = poisson(link = "log"), data = wa_surv)
summary(wa_15m)
car::vif(wa_15m)
r.squaredGLMM(wa_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
wa_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_15m2)
car::vif(wa_15m2)
r.squaredGLMM(wa_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for wa_15m4 but wa_15m3 close 2nd
AIC(wa_15m_null, wa_15m, wa_15m1, wa_15m2)
# More complex model is NOT better
anova(wa_15m, wa_15m2)

################# Species-Level Test - wa - Proportional Basal Area

## Null model - species-level (proportional BA)
wa_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|Site) + (1|survYear),
                          family = poisson(link = "log"), data = wa_surv) 
summary(wa_15m_prop_null)
r.squaredGLMM(wa_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
wa_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = wa_surv)
summary(wa_prop)
car::vif(wa_prop)
r.squaredGLMM(wa_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
wa_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
                  family = poisson(link = "log"), data = wa_surv)
summary(wa_prop1)
car::vif(wa_prop1)
r.squaredGLMM(wa_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
wa_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
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

#########################################################
################# Species-Level Test - ab - Raw Basal Area - Excludes RE of site

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

#########################################################
################# Species-Level Test - nro - Raw Basal Area

##### Subset nro data
nro_surv <- subset(allSpDat, Sp == "nro")

## Null model - species-level
nro_15m_null <- glmer(Year_1_counts ~ 1 +
                       (1|Site) + (1|survYear),
                     family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_null)
r.squaredGLMM(nro_15m_null) 

## Testing Panel A conceptual figure - CNDD - species-level
nro_15m <- glmer(Year_1_counts ~ consp_15m +
                  logPrior + log0yr +
                  (1|Site) + (1|survYear),
                family = poisson(link = "log"), data = nro_surv)
summary(nro_15m)
car::vif(nro_15m)
r.squaredGLMM(nro_15m) 

## Testing Partial Panel C conceptual figure (with AM) - species-level
nro_15m1 <- glmer(Year_1_counts ~ consp_15m + 
                   AM_15m + ECM_15m + 
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
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
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_15m2)
car::vif(nro_15m2)
r.squaredGLMM(nro_15m2) 

## Model comparison - species-level
### Compare AICs - lowest is for nro_15m4 but nro_15m3 close 2nd
AIC(nro_15m_null, nro_15m, nro_15m1, nro_15m2)
# More complex model is NOT better
anova(nro_15m1, nro_15m2)

################# Species-Level Test - nro - Proportional Basal Area

## Null model - species-level (proportional BA)
nro_15m_prop_null <- glmer(Year_1_counts ~ 1 +
                            (1|Site) + (1|survYear),
                          family = poisson(link = "log"), data = nro_surv) 
summary(nro_15m_prop_null)
r.squaredGLMM(nro_15m_prop_null) 

## Testing Panel A conceptual figure - CNDD - species-level (proportional BA)
nro_prop <- glmer(Year_1_counts ~ propConsp +
                   logPrior + log0yr +
                   (1|Site) + (1|survYear),
                 family = poisson(link = "log"), data = nro_surv)
summary(nro_prop)
car::vif(nro_prop)
r.squaredGLMM(nro_prop) 

## Testing Partial Panel C conceptual figure (with AM) - species-level (proportional BA)
nro_prop1 <- glmer(Year_1_counts ~ propConsp + propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
                  family = poisson(link = "log"), data = nro_surv)
summary(nro_prop1)
car::vif(nro_prop1)
r.squaredGLMM(nro_prop1) 

## Testing Partial Panel D conceptual figure (with AM) - species-level (proportional BA)
nro_prop2 <- glmer(Year_1_counts ~ propConsp*propAM + 
                    logPrior + log0yr +
                    (1|Site) + (1|survYear),
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


################################################ 
################################################ 
################################################ 
##########  PLOTTING CODE FOR LATER

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


