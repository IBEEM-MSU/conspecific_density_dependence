
#################################################
####### Analyses of one-year survivorship of Red Maple WRT mature tree neighborhood (BA)
####### Author: Bailey McNichol
####### Date: 17 March 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)

# For these models of one-year survivorship, we elected to include both random effects (REs) of site and survivorship year
# (1 | Site) + (1 | survYear)

## Full model to fit: 
# (#1 yr / #0 yr) ~ Conspecific BA * % EMF hosts + Con-myc BA + Hetero-myc BA + (total # RM present - #0 yr seedlings) + #0 yr seedlings + (1 | Site) + (1 | survYear)


#################################################
#################################################
##### Red Maple

### Read in data file from the R/L1 folder 
rm_surv_temp <- read.csv( "./data/L1/rm_oneYrSurv_matureBA.csv", row.names = 1)
rm_dens <- read.csv( "./data/L1/rm_prior_cohort_densities.csv")

# Merging the dataframes by matching TransectID and Obs.Year (in rm_dens) with TransectID and survYear (in rm_surv)
rm_surv_temp$TI_FY = paste(rm_surv_temp$TransectID,rm_surv_temp$First.Year)
rm_dens$TI_FY = paste(rm_dens$TransectID,rm_dens$First.Year)

## Add new column with match
rm_surv_temp$Prior_3_year_counts = rm_dens$Prior_3_year_counts[match(rm_surv_temp$TI_FY,
                                                                     rm_dens$TI_FY)]

rm_surv <- rm_surv_temp

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Calculate the proportion of EMF hosts within each radius
# %EMF - 5m
rm_surv$percEMF_5m = (rm_surv$hetMyc_5m / (rm_surv$hetMyc_5m + rm_surv$conMyc_5m))
# rm_surv$percEMF_5m = (rm_surv$hetMyc_5m / (rm_surv$hetMyc_5m + rm_surv$conMycNS_5m))
# %EMF - 10m
rm_surv$percEMF_10m = (rm_surv$hetMyc_10m / (rm_surv$hetMyc_10m + rm_surv$conMyc_10m))
# %EMF - 15m
rm_surv$percEMF_15m = (rm_surv$hetMyc_15m / (rm_surv$hetMyc_15m + rm_surv$conMyc_15m))

#### Calculate the proportion conspecific BA, proportion con-mycorrhizAL species BA including the species, and proportion con-mycorrhizAL species BA excluding the species for each radius
# ## Conspecific - 5m
# rm_surv$propCS_5m = (rm_surv$rm_5m / (rm_surv$conMyc_5m + rm_surv$hetMyc_5m))
# ## Con-Mycorrhizal including species - 5m
# rm_surv$propCM_5m = (rm_surv$conMyc_5m / (rm_surv$conMyc_5m + rm_surv$hetMyc_5m))
# ## Con-Mycorrhizal excluding species - 5m
# rm_surv$propCM_NS_5m = (rm_surv$conMycNS_5m / (rm_surv$conMycNS_5m + rm_surv$hetMyc_5m))
# ## Conspecific - 10m
# rm_surv$propCS_10m = (rm_surv$rm_10m / (rm_surv$conMyc_10m + rm_surv$hetMyc_10m))
# ## Con-Mycorrhizal including species - 10m
# rm_surv$propCM_10m = (rm_surv$conMyc_10m / (rm_surv$conMyc_10m + rm_surv$hetMyc_10m))
# ## Con-Mycorrhizal excluding species - 10m
# rm_surv$propCM_NS_10m = (rm_surv$conMycNS_10m / (rm_surv$conMycNS_10m + rm_surv$hetMyc_10m))
# ## Conspecific - 15m
# rm_surv$propCS_15m = (rm_surv$rm_15m / (rm_surv$conMyc_15m + rm_surv$hetMyc_15m))
# ## Con-Mycorrhizal including species - 15m
# rm_surv$propCM_15m = (rm_surv$conMyc_15m / (rm_surv$conMyc_15m + rm_surv$hetMyc_15m))
# ## Con-Mycorrhizal excluding species - 15m
# rm_surv$propCM_NS_15m = (rm_surv$conMycNS_15m / (rm_surv$conMycNS_15m + rm_surv$hetMyc_15m))

###############################
###############################

## Full model to fit: 
# (total_count_yr1 / total_count_yr0) ~ rm_10m * % EMF hosts + conMyc_10m + hetMyc_10m + Prior_3_year_counts + total_count_yr0 + (1 | Site) + (1 | survYear)

## 5-m radius
cor(rm_surv[, c("rm_5m", "percEMF_5m", "conMyc_5m", "hetMyc_5m",
                "Prior_3_year_counts", "total_count_yr0")], use = "complete.obs")

rm_mod_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_5m * percEMF_5m + 
                      conMyc_5m + hetMyc_5m + 
                      # Prior_3_year_counts + 
                      # total_count_yr0 + 
                      (1 | Site) + (1 | survYear),
                    #control = glmerControl(optimizer = "nloptwrap"),
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(rm_mod_5m)


## 10-m radius
cor(rm_surv[, c("rm_10m", "percEMF_10m", "conMyc_10m", "hetMyc_10m",
                "Prior_3_year_counts", "total_count_yr0")], use = "complete.obs")

rm_mod_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_10m * percEMF_10m + 
                      conMyc_10m + hetMyc_10m + 
                     # Prior_3_year_counts + 
                     # total_count_yr0 + 
                      (1 | Site) + (1 | survYear),
               #control = glmerControl(optimizer = "nloptwrap"),
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(rm_mod_10m)









##### Candidate Models - 5-m radius - total BA
## Effect of conspecific DD in 5-m radius
cs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                 rm_5m + (1 | Site) + (1 | survYear), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cs_5m)
## Effect of con-mycorrhizal DD in 5-m radius - including the species
cm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                 conMyc_5m + (1 | Site) + (1 | survYear), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cm_5m)
## Effect of con-mycorrhizal DD in 5-m radius - excluding the species
cmNs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                   conMycNS_5m + (1 | Site) + (1 | survYear), 
                 family = binomial(link = "logit"), 
                 data = rm_surv)
summary(cmNs_5m)
## Effect of hetero-mycorrhizal DD in 5-m radius 
hm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                 hetMyc_5m + (1 | Site) + (1 | survYear), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hm_5m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 5-m radius - including species 
hmcm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                   hetMyc_5m*conMyc_5m  + (1 | Site) + (1 | survYear), 
                 family = binomial(link = "logit"), 
                 data = rm_surv)
summary(hmcm_5m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 5-m radius - excluding species
hmcmNS_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                     hetMyc_5m * conMycNS_5m  + (1 | Site) + (1 | survYear), 
                   family = binomial(link = "logit"), 
                   data = rm_surv)
summary(hmcmNS_5m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_5m <- tibble(Model = c("cs_5m", "cm_5m", "cmNs_5m", "hm_5m", 
                              "hmcm_5m", "hmcmNS_5m"),
                    AIC = sapply(list(cs_5m, cm_5m, cmNs_5m, hm_5m, 
                                      hmcm_5m, hmcmNS_5m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_5m)
# Pseudo R2 for best models
r.squaredGLMM(hmcm_5m)
r.squaredGLMM(cm_5m)

###############################
##### Candidate Models - 5-m radius - proportional BA
## Effect of proportional conspecific DD in 5-m radius
prop_cs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
               propCS_5m + (1 | Site) + (1 | survYear), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(prop_cs_5m)
## Effect of proportional Con-Mycorrhizal DD in 5-m radius - including the species
prop_cm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      propCM_5m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(prop_cm_5m)
## Effect of proportional Con-Mycorrhizal DD in 5-m radius - excluding the species
prop_cmNS_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                    propCM_NS_5m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(prop_cmNS_5m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_5m_prop <- tibble(Model = c("prop_cs_5m", "prop_cm_5m", "prop_cmNS_5m"),
                    AIC = sapply(list(prop_cs_5m, prop_cm_5m, prop_cmNS_5m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_5m_prop)
# Pseudo R2 for best models
r.squaredGLMM(prop_cmNS_5m)

###############################
##################################
##### Candidate Models - 10-m radius - total BA
## Effect of conspecific DD in 10-m radius
cs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  rm_10m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cs_10m)
## Effect of con-mycorrhizal DD in 10-m radius - including the species
cm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  conMyc_10m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cm_10m)
## Effect of con-mycorrhizal DD in 10-m radius - excluding the species
cmNs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                    conMycNS_10m + (1 | Site) + (1 | survYear), 
                  family = binomial(link = "logit"), 
                  data = rm_surv)
summary(cmNs_10m)
## Effect of hetero-mycorrhizal DD in 10-m radius 
hm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  hetMyc_10m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(hm_10m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 10-m radius - including species 
hmcm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                    hetMyc_10m*conMyc_10m  + (1 | Site) + (1 | survYear), 
                  family = binomial(link = "logit"), 
                  data = rm_surv)
summary(hmcm_10m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 10-m radius - excluding species
hmcmNS_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      hetMyc_10m * conMycNS_10m  + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(hmcmNS_10m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_10m <- tibble(Model = c("cs_10m", "cm_10m", "cmNs_10m", "hm_10m", 
                               "hmcm_10m", "hmcmNS_10m"),
                     AIC = sapply(list(cs_10m, cm_10m, cmNs_10m, hm_10m, 
                                       hmcm_10m, hmcmNS_10m), AIC)) %>%
  arrange(AIC)  # Sort by AIC

# Print the sorted AIC values
print(aic_rm_10m)

###############################
##### Candidate Models - 10-m radius - proportional BA
## Effect of proportional conspecific DD in 10-m radius
prop_cs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      propCS_10m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(prop_cs_10m)
## Effect of proportional Con-Mycorrhizal DD in 10-m radius - including the species
prop_cm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      propCM_10m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(prop_cm_10m)
## Effect of proportional Con-Mycorrhizal DD in 10-m radius - excluding the species
prop_cmNS_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                        propCM_NS_10m + (1 | Site) + (1 | survYear), 
                      family = binomial(link = "logit"), 
                      data = rm_surv)
summary(prop_cmNS_10m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_10m_prop <- tibble(Model = c("prop_cs_10m", "prop_cm_10m", "prop_cmNS_10m"),
                         AIC = sapply(list(prop_cs_10m, prop_cm_10m, prop_cmNS_10m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_10m_prop)

###############################
##################################
##### Candidate Models - 15-m radius - total BA
## Effect of conspecific DD in 15-m radius
cs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  rm_15m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cs_15m)
## Effect of con-mycorrhizal DD in 15-m radius - including the species
cm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  conMyc_15m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cm_15m)
## Effect of con-mycorrhizal DD in 15-m radius - excluding the species
cmNs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                    conMycNS_15m + (1 | Site) + (1 | survYear), 
                  family = binomial(link = "logit"), 
                  data = rm_surv)
summary(cmNs_15m)
## Effect of hetero-mycorrhizal DD in 15-m radius 
hm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                  hetMyc_15m + (1 | Site) + (1 | survYear), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(hm_15m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 15-m radius - including species 
hmcm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                    hetMyc_15m*conMyc_15m  + (1 | Site) + (1 | survYear), 
                  family = binomial(link = "logit"), 
                  data = rm_surv)
summary(hmcm_15m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 15-m radius - excluding species
hmcmNS_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                      hetMyc_15m * conMycNS_15m  + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(hmcmNS_15m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_15m <- tibble(Model = c("cs_15m", "cm_15m", "cmNs_15m", "hm_15m", 
                               "hmcm_15m", "hmcmNS_15m"),
                     AIC = sapply(list(cs_15m, cm_15m, cmNs_15m, hm_15m, 
                                       hmcm_15m, hmcmNS_15m), AIC)) %>%
  arrange(AIC)  # Sort by AIC

# Print the sorted AIC values
print(aic_rm_15m)


###############################
##### Candidate Models - 15-m radius - proportional BA
## Effect of proportional conspecific DD in 15-m radius
prop_cs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       propCS_15m + (1 | Site) + (1 | survYear), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(prop_cs_15m)
## Effect of proportional Con-Mycorrhizal DD in 15-m radius - including the species
prop_cm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       propCM_15m + (1 | Site) + (1 | survYear), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(prop_cm_15m)
## Effect of proportional Con-Mycorrhizal DD in 15-m radius - excluding the species
prop_cmNS_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                         propCM_NS_15m + (1 | Site) + (1 | survYear), 
                       family = binomial(link = "logit"), 
                       data = rm_surv)
summary(prop_cmNS_15m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_15m_prop <- tibble(Model = c("prop_cs_15m", "prop_cm_15m", "prop_cmNS_15m"),
                          AIC = sapply(list(prop_cs_15m, prop_cm_15m, prop_cmNS_15m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_15m_prop)


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


###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### 
###### Plot survivorship versus PROPORTIONAL conmycorrhizal BA with fit lines from:
# 1) prop_cm_5m ; 2) prop_cmNS_5m; 3) prop_cs_5m
# Extract model coefficients (logit scale)
models <- data.frame(
  Model = c("Con-Mycorrhizal - RM", "Con-Mycorrhizal - no RM", "RM"),
  Intercept = c(fixef(prop_cm_5m)[1], fixef(prop_cmNS_5m)[1], fixef(prop_cs_5m)[1]),
  Slope = c(fixef(prop_cm_5m)[2], fixef(prop_cmNS_5m)[2], fixef(prop_cs_5m)[2])
)

# Generate predictor values (common range for x-axis)
x_values_prop <- seq(min(rm_surv$propCM_5m, na.rm = TRUE), 
                     max(rm_surv$propCM_5m, na.rm = TRUE), length.out = 100)

# Compute back-transformed survival probabilities
fit_data_prop <- do.call(rbind, lapply(1:nrow(models), function(i) {
  logit_pred <- models$Intercept[i] + models$Slope[i] * x_values_prop
  prob_pred <- plogis(logit_pred)  # Transform logit to probability
  data.frame(propCM_5m = x_values_prop, survivorship = prob_pred, Model = models$Model[i])
}))

## geom_smooth() - add line and CIs
## r.squaredGLMM() - pseudo R2 - gives theoretical and delta values

## Re-order levels
fit_data_prop$Model <- factor(fit_data_prop$Model, 
                         levels = c("Con-Mycorrhizal - RM", "Con-Mycorrhizal - no RM", "RM"))

# Scatterplot with fitted lines from all models
ggplot(rm_surv, aes(x = propCM_5m, y = total_count_yr1 / total_count_yr0)) +
  geom_point(alpha = 0.5) +  # Observed survivorship
  geom_line(data = fit_data_prop, aes(x = propCM_5m, y = survivorship, 
                                      color = Model), linewidth = 1) +
  scale_color_manual(values = c("blue", "red", "gold")) +  # Custom colors
  labs(x = "Proportion of Con-Mycorrhizal BA in 5-m Neighborhood", 
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
