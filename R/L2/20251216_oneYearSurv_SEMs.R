
#################################################
####### Implementing an SEM to integrate biotic and abiotic!
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 16 December 2025
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
library(DiagrammeR)
library(performance)

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
  left_join(clim_env_allSites, by = c("TransectID", "survYear", "X"))

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

# ### Create binary (0/1) variable for Seedling_myc
# allSpDat_surv_clim$Sdlg_myc_bin <- ifelse(allSpDat_surv_clim$Seedling_myc 
#                                           == "AM", 1, 0)

### Recode Seedling_myc as +0.5 for AM and -0.5 for ECM
# makes the main slopes interpretable as average effects across mycorrhizal types, and the interactions interpretable as deviations from that average.
allSpDat_surv_clim$Sdlg_myc_bin <- ifelse(allSpDat_surv_clim$Seedling_myc 
                                          == "AM", 0.5, -0.5)





### Most complex model to fit (15-m radius)
###  Abiotic model of one-year survivorship 
# (#1 yr) ~ meanGS_PDSI + minDailyTemp + maxDailyTemp + percOpenness + baseCEC + log(#0 yr) + (1 | Site) + (1 | Year)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")
allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

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
allSpDat$survYear = as.factor(allSpDat$survYear)
allSpDat_surv_clim <- allSpDat %>%
  left_join(clim_env_allSites, by = c("X", "TransectID", "survYear"))

### Center min and max temps by the minimum of each
allSpDat_surv_clim$mintemp = allSpDat_surv_clim$tmin_c - 
  min(allSpDat_surv_clim$tmin_c) 
allSpDat_surv_clim$maxtemp = allSpDat_surv_clim$tmax_c - 
  min(allSpDat_surv_clim$tmax_c)
### Make Site a factor
allSpDat_surv_clim$SiteID = as.factor(allSpDat_surv_clim$SiteID)

##### Subset out only the relevant columns
dat <- allSpDat_surv_clim[,c(2,6,8,35:38,41:43,49,52,54,55,56)]
# Rename columns
colnames(dat) <- c("SiteID","Sp","Year_1_counts","Seedling_myc", 
                   "ConspBA","Log_Recruits","Log_Prior","survYear",
                   "AM_BA","ECM_BA", "Mean_PDSI","perc_Openness",
                   "Base_CEC", "Min_Temp", "Max_Temp") 

## Subset for AM
AM_dat = filter(dat, Seedling_myc == "AM")
## Subset for ECM
ECM_dat = filter(dat, Seedling_myc == "ECM")

############################################
############################################
############ AM SEM

## ================================
## 1) Biotic mediator models (Abiotic → Biotic)
## ================================

m_AM_AMBA <- lmer(AM_BA ~ Mean_PDSI + Min_Temp + Max_Temp +
             perc_Openness + Base_CEC + (1|Sp) + (1|survYear) + (1|SiteID),
             data = AM_dat)

m_AM_ECMBA <- lmer(
  hetECM_15m ~ mean_pdsi + mintemp + maxtemp +
    percOpen + Base_CEC +
    (1 | Sp) + (1 | survYear) + (1 | SiteID),
  data = allSpDat_AM
)

m_AM_CONS <- lmer(
  consp_15m ~ mean_pdsi + mintemp + maxtemp +
    percOpen + Base_CEC +
    (1 | Sp) + (1 | survYear) + (1 | SiteID),
  data = allSpDat_AM
)

## ================================
## 2) AM Seedling outcome model
## ================================

m_AM_seedlings <- glmmTMB(
  Year_1_counts ~
    # Biotic main effects
    consp_15m + hetAM_15m + hetECM_15m +
    # Biotic interactions (from best model)
    consp_15m:hetAM_15m +
    consp_15m:hetECM_15m +
    # Abiotic direct effects
    mean_pdsi + mintemp + maxtemp +
    percOpen + Base_CEC +
    # Controls
    logPrior + log0yr +
    # Random effects
    (1 | Sp) + (1 | survYear) + (1 | SiteID),
  family = nbinom2,
  data = allSpDat_AM
)

## ================================
## 3) Assemble AM SEM
## ================================

sem_AM <- psem(
  m_AM_AMBA,
  m_AM_ECMBA,
  m_AM_CONS,
  m_AM_seedlings
)

summary(sem_AM, .progress = FALSE)

# Extract coefficients from piecewise SEM
coef_table <- coefs(sem_mod) %>%
  as.data.frame() %>%
  select(Response, Predictor, Estimate, Std.Error, Std.Estimate, DF, P.Value) %>%
  arrange(Response, Predictor)

# # Optional: round numeric columns for readability
# coef_table <- coef_table %>%
#   mutate(across(c(Estimate, Std.Error, DF, P.Value), ~round(., 3))) %>%
#   mutate(across(c(Estimate, Std.Error, DF, P.Value), 
#                 ~ifelse(. == 0, "<0.001", as.character(.))))


# Identify which predictors exist in the data
predictors_real <- coef_table$Predictor[coef_table$Predictor %in% names(dat)]

# Compute SDs for real predictors and responses
SDs <- dat %>%
  select(all_of(predictors_real), all_of(unique(coef_table$Response))) %>%
  summarise(across(everything(), \(x) sd(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "SD")

# Convert to named vector for easy lookup
SD_lookup <- setNames(SDs$SD, SDs$Variable)

# Compute standardized SEs without touching other columns
coef_table <- coef_table %>%
  rowwise() %>%
  mutate(Std.Error.Std = ifelse(
    Predictor %in% names(SD_lookup) & Response %in% names(SD_lookup),
    Std.Error * (SD_lookup[Predictor] / SD_lookup[Response]),
    NA_real_  # NA for interactions or missing columns
  )) %>%
  ungroup()

# View
coef_table

write.csv(coef_table, "./R/L1/SEM_Full_Coefficients.csv", row.names = FALSE)

#################### Plotting the SEM

# Super basic plot
plot(sem_mod)

### More customized plot
plot(sem_mod, layout = "nicely", 
     node_attrs = list(shape = "oval",
  width = c(0.4, 0.4, 0.4, 0.7, 
            0.6, 0.6, 0.6, 0.6, 0.6, 
            0.6, 0.6, 0.6, 0.9, 0.9, 0.9),
  height = 0.18, # shrink node height
  fillcolor = c(rep("#9dc183", 3),   # Biotic - Neighborhood BA
    rep("white", 1),     # Response variable
    rep("orange", 5),    # Abiotic variables
    rep("lightblue", 1),      # Trait - Myc_type 
    rep("gray", 2),    # Controls - Log_Prior + Log_Recruits
    rep("pink", 3)       # Interactions
  ), outline = "black", color = "black",
  fontsize = 5.3, fontcolor = "black"),
  # edge_attrs = list(
  #   penwidth = 1,
  #   label = "Std.Estimate",
  #   fontsize = 6,fontcolor = "black",
  #   color = "black", labeldistance = 3 ))
  edge_attrs = list(penwidth = 1,
                    # scale factor to adjust thickness
    label = NA,fontsize = 0,fontcolor = "white",
    color = "black"))



# # Get coefficients including Std.Estimate
# sem_coefs <- coefs(sem_mod)
# 
# # Inspect
# head(sem_coefs)
# 
# # Select relevant columns for plotting
# edge_labels <- sem_coefs[, c("Response", "Predictor", "Std.Estimate")]
# 
# # Optionally, round for cleaner plotting
# edge_labels$Std.Estimate <- round(edge_labels$Std.Estimate, 3)
# 
# edge_labels




############################################
############################################

### calculations to get group-specific slopes (AM vs ECM) and put the multiplicative effects on the response scale.
# reported estimates
beta_AM  <- -0.0753
beta_M_AM <-  0.0289  # Myc_type:AM_BA

beta_ECM <-  0.0780
beta_M_ECM <- 0.0319  # Myc_type:ECM_BA

# slopes
slope_AM_for_AMseedling  <- beta_AM  + 0.5 * beta_M_AM
slope_AM_for_ECMseedling <- beta_AM  - 0.5 * beta_M_AM

slope_ECM_for_AMseedling  <- beta_ECM + 0.5 * beta_M_ECM
slope_ECM_for_ECMseedling <- beta_ECM - 0.5 * beta_M_ECM

slope_AM_for_AMseedling
slope_AM_for_ECMseedling
slope_ECM_for_AMseedling
slope_ECM_for_ECMseedling

# multiplicative effects:
exp(slope_AM_for_AMseedling)   # multiplicative factor per unit AM_BA for AM seedlings
exp(slope_AM_for_ECMseedling)  # for ECM seedlings

### Using emmeans / emtrends on your fitted glmer m_seedlings to estimate slopes by Myc_type

library(emmeans)

# estimated slope of AM_BA for each Myc_type
emtrends(m_seedlings, ~ Myc_type, var = "AM_BA")

# same for ECM_BA
emtrends(m_seedlings, ~ Myc_type, var = "ECM_BA")

# predicted line for AM_BA across a range for each Myc_type:
newdat <- with(dat, 
               expand.grid(AM_BA = seq(min(AM_BA), max(AM_BA), length=30),
                           ECM_BA = mean(ECM_BA, na.rm=TRUE),
                           Myc_type = c(-0.5,0.5),
                           Log_Prior = mean(Log_Prior, na.rm=TRUE),
                           Log_Recruits = mean(Log_Recruits, na.rm=TRUE),
                           Mean_PDSI = mean(Mean_PDSI, na.rm=TRUE),
                           Min_Temp = mean(Min_Temp, na.rm=TRUE),
                           Max_Temp = mean(Max_Temp, na.rm=TRUE),
                           perc_Openness = mean(perc_Openness, na.rm=TRUE),
                           Base_CEC = mean(Base_CEC, na.rm=TRUE),
                           Sp = NA, survYear = NA, SiteID = NA))

# get predictions (on response scale)
pred <- predict(m_seedlings, newdata = newdat, re.form = NA, type = "response")
newdat$pred <- pred

library(ggplot2)
ggplot(newdat, aes(x = AM_BA, y = pred, color = factor(Myc_type))) +
  geom_line() +
  scale_color_manual(values = c("-0.5"="blue","0.5"="red"),
                     labels = c("ECM","AM"),
                     name = "Myc type") +
  labs(y = "Predicted Year_1_counts", x = "AM_BA (units)") +
  theme_minimal()




############################################
############################################
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

## Testing for multicollinearity
m_seedlings_AMEM <- lm(Year_1_counts ~
                         # Biotic main effects
                         consp_15m + Sdlg_myc_bin + 
                         AM_15m + ECM_15m +
                         # Biotic interactions
                         consp_15m:AM_15m + consp_15m:ECM_15m +
                         # Seedling_myc:AM_15m + Seedling_myc:ECM_15m +
                         # Abiotic direct effects
                         mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC +
                         # Controls
                         logPrior + log0yr , data = dat)

m_seedlings_AM <- lm(Year_1_counts ~
                       # Biotic main effects
                       consp_15m + Sdlg_myc_bin + 
                       AM_15m + ECM_15m +
                       # Biotic interactions
                       consp_15m:AM_15m + 
                       # Seedling_myc:AM_15m + Seedling_myc:ECM_15m +
                       # Abiotic direct effects
                       mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC +
                       # Controls
                       logPrior + log0yr , data = dat)

m_seedlings_EM <- lm(Year_1_counts ~
                       # Biotic main effects
                       consp_15m + Sdlg_myc_bin + 
                       AM_15m + ECM_15m +
                       # Biotic interactions
                       consp_15m:ECM_15m +
                       # Seedling_myc:AM_15m + Seedling_myc:ECM_15m +
                       # Abiotic direct effects
                       mean_pdsi + mintemp + maxtemp + percOpen + Base_CEC +
                       # Controls
                       logPrior + log0yr , data = dat)

car::vif(m_seedlings_AMEM)
car::vif(m_seedlings_AM)
car::vif(m_seedlings_EM)

AIC(m_seedlings_AMEM, m_seedlings_AM, m_seedlings_EM)


