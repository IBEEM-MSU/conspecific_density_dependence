
#################################################
####### Implementing an SEM to integrate biotic and abiotic!
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 12 August 2025
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
library(lavaan)
library(piecewiseSEM)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1
allSpDat$survYear = as.factor(allSpDat$survYear)

## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

# Add X column (quadrat by the 0.5)
allSpDat$X = allSpDat$quadrat - 0.5

# Add true site number
allSpDat$SiteNo = ifelse(allSpDat$SiteID == 1, "49.01",
                  ifelse(allSpDat$SiteID == 2, "59.12",
                  ifelse(allSpDat$SiteID == 3, "8.21",
                  ifelse(allSpDat$SiteID == 4, "39.21",
                  ifelse(allSpDat$SiteID == 5, "16.35",
                  ifelse(allSpDat$SiteID == 6, "60.35",
                  ifelse(allSpDat$SiteID == 7, "7.40",
                  ifelse(allSpDat$SiteID == 8, "40.40",
                  ifelse(allSpDat$SiteID == 9, "21.43",
                  ifelse(allSpDat$SiteID == 10, "25.43",
                  ifelse(allSpDat$SiteID == 11, "6.45", "24.45")))))))))))

# Create TransectID
allSpDat$TransectID = paste0(allSpDat$SiteNo, allSpDat$X)

### Left join to add climate variables to data
allSpDat_surv_clim <- allSpDat %>%
  left_join(clim_env_allSites, by = c("TransectID", "survYear"))

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

############################################
##### Most complex abiotic model
allSp_abiotic <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      percOpen + Base_CEC + log0yr +
                      (1|SiteID) + (1|survYear) + (1|Sp),
                    family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_abiotic)
car::vif(allSp_abiotic)
r.squaredGLMM(allSp_abiotic)

##### Most complex biotic model
allSp_biotic <- glmer(Year_1_counts ~ consp_15m + Seedling_myc + 
                       AM_15m + ECM_15m + 
                       consp_15m:AM_15m + consp_15m:ECM_15m + 
                       Seedling_myc:AM_15m + Seedling_myc:ECM_15m + 
                       consp_15m:Seedling_myc:AM_15m + 
                       consp_15m:Seedling_myc:ECM_15m + 
                       logPrior + log0yr + (1|SiteID) + (1|survYear) + (1|Sp),
                     family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_biotic)
car::vif(allSp_biotic)
r.squaredGLMM(allSp_biotic)




dat <- allSpDat_surv_clim
dat <- filter(allSpDat_surv_clim,Sp == "sm")

## 1) Biotic mediator models (Abiotic → Biotic)
m_AM <- lmer(AM_15m ~ mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC + (1|survYear)
              # (1|SiteID) #+ (1|survYear) + (1|Sp)
             , data = dat)

m_ECM <- lmer(ECM_15m ~ mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC + (1|survYear)
               # (1|SiteID) #+ (1|survYear) + (1|Sp)
              , data = dat)

m_CON <- lmer(consp_15m ~ mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC + (1|survYear)
                #(1|SiteID) #+ (1|survYear) + (1|Sp)
              , data = dat)

## 2) Seedling outcome model (Biotic + Abiotic + Trait + Interactions + Controls)
m_seedlings <- glmer(
  Year_1_counts ~
    # main effects
    consp_15m + #Seedling_myc + 
    AM_15m + ECM_15m + 
    #consp_15m:AM_15m + consp_15m:ECM_15m + 
    #Seedling_myc:AM_15m + Seedling_myc:ECM_15m + 
    #consp_15m:Seedling_myc:AM_15m + 
    #consp_15m:Seedling_myc:ECM_15m +
    mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC +
    logPrior + log0yr + (1|survYear),
    # random effects
   # (1|SiteID) +  + (1|Sp),
  family = poisson(link = "log"), data = dat)

## 3) Assemble the SEM
sem_mod <- psem(
  m_AM, m_ECM, m_CON,  # mediators
  m_seedlings          # outcome
)

## 4) Inspect fit and paths
summary(sem_mod, .progress = FALSE, .verbose = FALSE)

summary(sem_mod)       # path coefficients, directed separation (d-sep) tests
sem.coefs(sem_mod)     # standardized coefficients (where applicable)
plot(sem_mod)          # quick path plot (base plotting)




#########################################################
#########################################################

## Practice example - confirmatory factor analysis

HS.model <- ' visual  =~ x1 + x2 + x3 
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)

summary(fit, fit.measures = TRUE)


## Practice example - structural equation model
model <- '
  # measurement model
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit <- sem(model, data = PoliticalDemocracy)
summary(fit, standardized = TRUE)

#  all free parameters are automatically named according to a simple set of rules. This is convenient, for example, if equality constraints are needed
model <- '
  # latent variable definitions
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # residual (co)variances
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit <- sem(model, 
           data = PoliticalDemocracy)


fit <- cfa(HS.model, 
           data = HolzingerSwineford1939, 
           meanstructure = TRUE)
summary(fit, remove.unused = FALSE)

# coef() extracts the estimated values of the free parameters in the model, together with their names. Each name consists of three parts and reflects the part of the formula where the parameter was involved. The first part is the variable name that appears on the left-hand side (lhs) of the formula. The middle part is the operator type (op) of the formula, and the third part is the variable in the right-hand side (rhs) of the formula that corresponds with the parameter.
coef(fit)

