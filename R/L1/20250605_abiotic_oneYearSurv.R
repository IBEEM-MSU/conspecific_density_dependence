
#################################################
####### Analyses of one-year survivorship WRT abiotic environment
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

### Most complex model to fit (15-m radius)
###  Abiotic model of one-year survivorship 
# (#1 yr) ~ meanGS_PDSI + minDailyTemp + maxDailyTemp + percOpenness + baseCEC + log(#0 yr) + (1 | Site) + (1 | Year)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")
## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1

# Add transformed variables
allSpDat <- allSpDat %>%
  mutate(
    Site = as.factor(SiteID),
    survYear = as.factor(survYear),
    log0yr = log(Year_0_counts),
    logPrior = log(Prior_3_year_cohorts_counts + 1))

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

#########################################################
#########################################################
################# All Species Model

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

#################### Effects of environment and climate
allSp_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                   percOpen + Base_CEC + #log0yr +
                   (1 | SiteID) + (1 | survYear) + (1 | Sp),
                 family = poisson(link = "log"), data = allSpDat_surv_clim)
summary(allSp_abio)
car::vif(allSp_abio)

#########################################################
#########################################################
################# Species-Level Test - RM

##### Subset RM data
rm_surv <- subset(allSpDat_surv_clim, Sp == "rm")

#################### Effects of environment and climate
rm_abio <- glmer(Year_1_counts ~ mean_pdsi + mintemp + maxtemp +
                      percOpen + Base_CEC + #log0yr +
                      (1 | SiteID) + (1 | survYear),
                    family = poisson(link = "log"), data = rm_surv)
summary(rm_abio)
car::vif(rm_abio)

# ## VIFs and pseudo r-squared values
# car::vif(rm_abio) # all < 2
# r.squaredGLMM(rm_abio) 
# 
# # Simulate residuals
# sim_res <- simulateResiduals(fittedModel = rm_abio, n = 1000)
# # QQ plot of residuals
# plotQQunif(sim_res)


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
models <- list(rm = rm_abio, sm = sm_abio, ab = ab_abio, wa = wa_abio,
               bc = bc_abio, nro = nro_abio, bonro = bonro_abio,
               wo = wo_abio, iw = iw_abio)

# Map terms by model
allSp_results <- purrr::imap_dfr(models, extract_model_summary)

# Add uniform names for covariates across species
allSp_results <- allSp_results %>%
  rowwise() %>%
  mutate(component = case_when(
    term == "mean_pdsi" ~ "Mean PDSI",
    term == "mintemp" ~ "Minimum Temp",
    term == "maxtemp" ~ "Maximum Temp",
    term == "percOpen" ~ "Canopy Openness (%)",
    term == "Base_CEC" ~ "Base CEC",
   # term == "log0yr" ~ "Log(Recruit Density)",
    TRUE ~ NA_character_
  )) %>%
  ungroup()
# Make into an ordered factor
allSp_results$component = factor(allSp_results$component,
                                 levels = c( #"Log(Recruit Density)", 
                                             "Base CEC",
                                             "Canopy Openness (%)",
                                             "Maximum Temp",
                                             "Minimum Temp",
                                             "Mean PDSI"))

# Add AM or ECM to each row by species
allSp_results$Myc_Type = ifelse(allSp_results$species == "rm", "AM",
                        ifelse(allSp_results$species == "sm", "AM",
                        ifelse(allSp_results$species == "wa", "AM",
                        ifelse(allSp_results$species == "bc", "AM","ECM"))))

# Add column for significance of parameter estimate
allSp_results$Significant = ifelse(allSp_results$p.value <= 0.05, "yes", "no")

######################
## Create a color vector with red shades for AM and blue shades for ECM
species_cols <- c(ab = "#B4D4DA", bonro = "#72BCDC", iw = "#3F8BBA",
                  nro = "#1C65A3", wo = "#26456E",
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
    component = factor(component, levels = c( #"Log(Recruit Density)", 
                                              "Base CEC",
                                              "Canopy Openness (%)",
                                              "Maximum Temp",
                                              "Minimum Temp",
                                              "Mean PDSI")),
    Significant = factor(Significant, levels = c("yes", "no")))

## Calculate mean parameter estimates for AM versus ECM
group_means <- allSp_results %>%
  group_by(Myc_Type, component) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    std.error = sd(estimate, na.rm = TRUE) / sqrt(n()),  
    .groups = "drop")

# Plot - All Species (ECM and AM)
ggplot(allSp_results, aes(x = estimate, y = component, 
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

### Plot - ECM only
# Make just an ECM group
ecm_group <- group_means %>% filter(Myc_Type == "ECM" )
                                    #& component != "Log(Recruit Density)")
# Plot
allSp_results %>% filter(., Myc_Type == "ECM" ) %>% 
                        # & component != "Log(Recruit Density)")  %>% 
  ggplot(., aes(x = estimate,  y = component, 
                color = species_label,  shape = Significant)) +
  geom_point(aes(size = log1p(percBA)), 
             position = position_dodge(width = 0.6)) +   # size only here
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error), 
                 height = 0, position = position_dodge(width = 0.6)) +
  ## Add ECM and AM group means
  geom_point(data = ecm_group, 
             aes(x = mean_estimate, y = component, shape = Myc_Type), 
             color = "black", size = 3, inherit.aes = FALSE) +
  geom_errorbarh(data = ecm_group,
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


### Plot - AM only
# Make just an AM group
am_group <- group_means %>% filter(Myc_Type == "AM" )
                                  # &   component != "Log(Recruit Density)")
# Plot
allSp_results %>% filter(., Myc_Type == "AM" ) %>% 
                        # &  component != "Log(Recruit Density)")  %>% 
  ggplot(., aes(x = estimate,  y = component, 
                color = species_label,  shape = Significant)) +
  geom_point(aes(size = log1p(percBA)), 
             position = position_dodge(width = 0.6)) +   # size only here
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
  geom_errorbarh(aes(xmin = estimate - std.error, 
                     xmax = estimate + std.error), 
                 height = 0, position = position_dodge(width = 0.6)) +
  ## Add ECM and AM group means
  geom_point(data = am_group, 
             aes(x = mean_estimate, y = component, shape = Myc_Type), 
             color = "black", size = 3, inherit.aes = FALSE) +
  geom_errorbarh(data = am_group,
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