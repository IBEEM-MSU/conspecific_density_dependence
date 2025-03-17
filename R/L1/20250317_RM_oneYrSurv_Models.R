
#################################################
####### Analyses of one-year survivorship WRT mature tree neighborhood (BA)
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

# For these models of one-year survivorship, we elected to include both random effects (REs) of site and survivorship year, as, for every neighborhood basal area, likelihood ratio tests (run as: anova(cs_5m, cs_5m_Yr, test = "LRT") where the two models being compared were one with just RE of site and another with REs of site and year) indicated that the RE of year significantly improved the model fit (AIC values also much lower)

##### Potential random effect nesting to test: 
# (1 | survYear/Site)
# (1 | Site/survYear)
# (1 | Site) + (1 | survYear)

## Create color palette for sites
fertility = c("49.01" = "#FCECAEFF","59.12" = "#FECD90FF",
              "8.21" = "#FEAC76FF","39.21" = "#FC8B62FF",
              "16.35" = "#F56B5CFF","60.35" = "#D8456CFF",
              "7.4" = "#BD3977FF","40.4" = "#A02F7FFF",
              "21.43" = "#762181FF","25.43" = "#5B167EFF",
              "6.45" = "#3F0F72FF","24.45" = "#21114EFF")

#################################################
#################################################
##### Red Maple

### Read in data file from the R/L1 folder 
rm_surv <- read.csv( "./data/L1/rm_oneYrSurv_matureBA.csv", row.names = 1)

#### Make site and year into factors
rm_surv$Site = as.factor(rm_surv$Site)
rm_surv$survYear = as.factor(rm_surv$survYear)

#### Calculate the percent conspecific BA, percent con-mycorrhizAL species BA including the species, and percent con-mycorrhizAL species BA excluding the species for each radius
## Conspecific - 5m
rm_surv$percCS_5m = (rm_surv$rm_5m / (rm_surv$conMyc_5m + rm_surv$hetMyc_5m))
## Con-Mycorrhizal including species - 5m
rm_surv$percCM_5m = (rm_surv$conMyc_5m / (rm_surv$conMyc_5m + rm_surv$hetMyc_5m))
## Con-Mycorrhizal excluding species - 5m
rm_surv$percCM_NS_5m = (rm_surv$conMycNS_5m / (rm_surv$conMycNS_5m + rm_surv$hetMyc_5m))
## Conspecific - 10m
rm_surv$percCS_10m = (rm_surv$rm_10m / (rm_surv$conMyc_10m + rm_surv$hetMyc_10m))
## Con-Mycorrhizal including species - 10m
rm_surv$percCM_10m = (rm_surv$conMyc_10m / (rm_surv$conMyc_10m + rm_surv$hetMyc_10m))
## Con-Mycorrhizal excluding species - 10m
rm_surv$percCM_NS_10m = (rm_surv$conMycNS_10m / (rm_surv$conMycNS_10m + rm_surv$hetMyc_10m))
## Conspecific - 15m
rm_surv$percCS_15m = (rm_surv$rm_15m / (rm_surv$conMyc_15m + rm_surv$hetMyc_15m))
## Con-Mycorrhizal including species - 15m
rm_surv$percCM_15m = (rm_surv$conMyc_15m / (rm_surv$conMyc_15m + rm_surv$hetMyc_15m))
## Con-Mycorrhizal excluding species - 15m
rm_surv$percCM_NS_15m = (rm_surv$conMycNS_15m / (rm_surv$conMycNS_15m + rm_surv$hetMyc_15m))

###############################
###############################
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

###############################
##### Candidate Models - 5-m radius - proportional BA
## Effect of proportional conspecific DD in 5-m radius
perc_cs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
               percCS_5m + (1 | Site) + (1 | survYear), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(perc_cs_5m)
## Effect of proportional Con-Mycorrhizal DD in 5-m radius - including the species
perc_cm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      percCM_5m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(perc_cm_5m)
## Effect of proportional Con-Mycorrhizal DD in 5-m radius - excluding the species
perc_cmNS_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                    percCM_NS_5m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(perc_cmNS_5m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_5m_prop <- tibble(Model = c("perc_cs_5m", "perc_cm_5m", "perc_cmNS_5m"),
                    AIC = sapply(list(perc_cs_5m, perc_cm_5m, perc_cmNS_5m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_5m_prop)

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
perc_cs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      percCS_10m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(perc_cs_10m)
## Effect of proportional Con-Mycorrhizal DD in 10-m radius - including the species
perc_cm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      percCM_10m + (1 | Site) + (1 | survYear), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(perc_cm_10m)
## Effect of proportional Con-Mycorrhizal DD in 10-m radius - excluding the species
perc_cmNS_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                        percCM_NS_10m + (1 | Site) + (1 | survYear), 
                      family = binomial(link = "logit"), 
                      data = rm_surv)
summary(perc_cmNS_10m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_10m_prop <- tibble(Model = c("perc_cs_10m", "perc_cm_10m", "perc_cmNS_10m"),
                         AIC = sapply(list(perc_cs_10m, perc_cm_10m, perc_cmNS_10m), AIC)) %>%
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
perc_cs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       percCS_15m + (1 | Site) + (1 | survYear), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(perc_cs_15m)
## Effect of proportional Con-Mycorrhizal DD in 15-m radius - including the species
perc_cm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       percCM_15m + (1 | Site) + (1 | survYear), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(perc_cm_15m)
## Effect of proportional Con-Mycorrhizal DD in 15-m radius - excluding the species
perc_cmNS_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                         percCM_NS_15m + (1 | Site) + (1 | survYear), 
                       family = binomial(link = "logit"), 
                       data = rm_surv)
summary(perc_cmNS_15m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_rm_15m_prop <- tibble(Model = c("perc_cs_15m", "perc_cm_15m", "perc_cmNS_15m"),
                          AIC = sapply(list(perc_cs_15m, perc_cm_15m, perc_cmNS_15m), AIC)) %>%
  arrange(AIC)  # Sort by AIC
# Print the sorted AIC values
print(aic_rm_15m_prop)


###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### 
###### Plot survivorship versus conmycorrhizal BA with fit lines from:
# 1) cm_5m ; 2) cmNs_5m; 3) cs_5m
# Extract model coefficients
intercept_cm_5m <- fixef(cm_5m)[1]  # Intercept on logit scale
slope_cm_5m <- fixef(cm_5m)[2]      # Slope on logit scale

# Generate predictor values (range of conMycNS_5m)
x_values <- seq(min(rm_surv$conMyc_5m), max(rm_surv$conMyc_5m), length.out = 100)

# Convert logit-scale predictions to probability scale
logit_predictions <- intercept_cm_5m + slope_cm_5m * x_values
prob_predictions <- plogis(logit_predictions)  # Inverse logit function

# Create dataframe for plotting the fitted line
fit_data <- data.frame(conMycNS_5m = x_values, survivorship = prob_predictions)

# Create scatterplot with fitted line
ggplot(rm_surv, aes(x = conMycNS_5m, y = total_count_yr1 / total_count_yr0)) +
  geom_point(alpha = 0.5) +  # Observed survivorship
  geom_line(data = fit_data, aes(x = conMycNS_5m, y = survivorship), color = "red", linewidth = 1) +
  labs(x = "Basal Area of Trees (conMycNS_5m)", y = "Seedling Survivorship", title = "Seedling Survivorship vs. Basal Area") +
  theme_minimal()


ggplot(rm_surv, aes(x = conMyc_5m, y = Survivorship, color = Site)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = -0.78743, slope = -0.84482, color = "black") + # conMyc_5m
  scale_color_manual(values = fertility) +
  labs(x = "Con-mycorrhizal Basal Area (m²) in 5-m Neighborhood",
       y = "Survivorship") +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 10), 
        axis.title = element_text(size = 10), legend.key = element_blank(), 
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "none")
