
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

## Full model to fit: 
# (#1 yr / #0 yr) ~ Conspecific BA * % EMF hosts + Con-myc BA + Hetero-myc BA + (total # RM present - #0 yr seedlings) + #0 yr seedlings + (1 | Site) + (1 | survYear)

#################################################
#################################################
##### Read in file that has every TransectID for every site
TI <- read.csv("./data/L1/all_TransectIDs_Sites.csv")

TI$Site = as.character(TI$Site)

##### Expand out to have each row repeated for each Obs.Year (2001-2024)    
# Expand with Obs.Year from 2001 to 2024
TI_ObsYear <- crossing(TI, Obs.Year = 2001:2024)
               
##### Read in previous three years of seedling counts for each species for each quadrat/obs year combination, and rename "Prior_3_year_counts" to include sp name
rm_dens <- read.csv( "./data/L1/rm_prior_cohort_densities.csv")
rm_dens <- rename(rm_dens, rm_Prior3_counts = Prior_3_year_counts)
sm_dens <- read.csv( "./data/L1/sm_prior_cohort_densities.csv")
sm_dens <- rename(sm_dens, sm_Prior3_counts = Prior_3_year_counts)
ab_dens <- read.csv( "./data/L1/ab_prior_cohort_densities.csv")
ab_dens <- rename(ab_dens, ab_Prior3_counts = Prior_3_year_counts)
wa_dens <- read.csv( "./data/L1/wa_prior_cohort_densities.csv")
wa_dens <- rename(wa_dens, wa_Prior3_counts = Prior_3_year_counts)
bc_dens <- read.csv( "./data/L1/bc_prior_cohort_densities.csv")
bc_dens <- rename(bc_dens, bc_Prior3_counts = Prior_3_year_counts)
nro_dens <- read.csv( "./data/L1/nro_prior_cohort_densities.csv")
nro_dens <- rename(nro_dens, nro_Prior3_counts = Prior_3_year_counts)
bonro_dens <- read.csv( "./data/L1/bonro_prior_cohort_densities.csv")
bonro_dens <- rename(bonro_dens, bonro_Prior3_counts = Prior_3_year_counts)
wo_dens <- read.csv( "./data/L1/wo_prior_cohort_densities.csv")
wo_dens <- rename(wo_dens, wo_Prior3_counts = Prior_3_year_counts)
iw_dens <- read.csv( "./data/L1/iw_prior_cohort_densities.csv")
iw_dens <- rename(iw_dens, iw_Prior3_counts = Prior_3_year_counts)
bs_dens <- read.csv( "./data/L1/bs_prior_cohort_densities.csv")
bs_dens <- rename(bs_dens, bs_Prior3_counts = Prior_3_year_counts)

#### Add the 3-year count columns for each species to the same data frame

# Step 1: Create a named list of all your species dataframes
species_dfs <- list(rm = rm_dens, sm = sm_dens,
  ab = ab_dens, wa = wa_dens, bc = bc_dens,
  nro = nro_dens, bonro = bonro_dens, wo = wo_dens,
  iw = iw_dens, bs = bs_dens)

# Step 2: Subset each dataframe to keep only the relevant merge columns and the species count
species_dfs_cleaned <- lapply(names(species_dfs), function(sp) {
  df <- species_dfs[[sp]]
  df_subset <- df[, c("TransectID", "Site", 
          "SiteID", "Obs.Year", paste0(sp, "_Prior3_counts"))]
  return(df_subset)
})

# Step 3: Merge all species-specific dataframes into TI_ObsYear
allSp_TI_ObsYear <- Reduce(function(x, y) {
  merge(x, y, by = c("TransectID", "Site", 
                     "SiteID", "Obs.Year"), all.x = TRUE)
}, c(list(TI_ObsYear), species_dfs_cleaned))

##### Write CSV with the counts by species for every combination
# write.csv(allSp_TI_ObsYear, "./data/L1/allSp_prior_cohort_densities.csv",
#           row.names = FALSE)

#################################################
#################################################
##### Red Maple

### Read in data file from the R/L1 folder 
rm_surv_temp <- read.csv( "./data/L1/rm_oneYrSurv_matureBA.csv", row.names = 1)

### Filter to include 2001 as first First.Year
rm_surv_temp <- filter(rm_surv_temp, First.Year > 2000)

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
# Heterospecific
rm_surv_temp$hetSp_Prior3_counts = 
  allSp_TI_ObsYear$hetSp_Prior3_counts[match(rm_surv_temp$TI_FY,
                                             allSp_TI_ObsYear$TI_FY)]

rm_surv <- rm_surv_temp

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
rm_surv$percEMF_5m = (rm_surv$hetMyc_5m / (rm_surv$hetMyc_5m + rm_surv$conMyc_5m))
# %EMF - 10m
rm_surv$percEMF_10m = (rm_surv$hetMyc_10m / (rm_surv$hetMyc_10m + rm_surv$conMyc_10m))
# %EMF - 15m
rm_surv$percEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))

###############################
###############################

## Full model to fit: 
# (total_count_yr1 / total_count_yr0) ~ rm_10m * % EMF hosts + conMycNS_10m + hetMyc_10m + Prior_3_year_counts + total_count_yr0 + (1 | Site) + (1 | survYear)

#################### 5-m radius
# cor(rm_surv[, c("rm_5m", "percEMF_5m", "conMyc_5m", "hetMyc_5m",
#                 "Prior_3_year_counts", "total_count_yr0")], use = "complete.obs")

### Most complex model
# Scaling the count covariates using log10() to account for difference in magnitude 
rm_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_5m * percEMF_5m + 
                      conMycNS_5m + hetMyc_5m + 
                      log10(rm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                      family = binomial(link = "logit"), 
                      data = rm_surv)
summary(rm_mod_5m)

car::vif(rm_mod_5m)

#################### 10-m radius
### Most complex model
# Scaling the count covariates using log10() to account for difference in magnitude 
rm_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                     rm_10m * percEMF_10m + 
                     conMycNS_10m + hetMyc_10m + 
                     log10(rm_Prior3_counts + 1) + 
                     log10(hetSp_Prior3_counts + 1) + 
                     log10(total_count_yr0 + 1) + 
                     (1 | Site) + (1 | survYear),
                   family = binomial(link = "logit"), 
                   data = rm_surv)
summary(rm_mod_10m)

car::vif(rm_mod_10m)


#################### 15-m radius
### Most complex model
# Scaling the count covariates using log10() to account for difference in magnitude 
rm_mod_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_15m * percEMF_15m + 
                      conMycNS_15m + hetMyc_15m + 
                      log10(rm_Prior3_counts + 1) + 
                      log10(hetSp_Prior3_counts + 1) + 
                      log10(total_count_yr0 + 1) + 
                      (1 | Site) + (1 | survYear),
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(rm_mod_15m)

car::vif(rm_mod_15m)







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



