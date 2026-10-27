
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
library(janitor)

#### Read in csv with all survivorship data for all species
allSpDat <- read.csv("./data/L3/final_seedling_demogr_modeling_data.csv")
allSpDat$Seedling_myc <- factor(allSpDat$Seedling_myc, levels = c("AM", "ECM"))

## Environmental/climate data
clim_env_allSites <- read.csv("./data/L1/clim_envVars_allSites.csv")
clim_env_allSites$Site = as.character(clim_env_allSites$Site)
clim_env_allSites$Site[clim_env_allSites$Site=="7.4"] <- "7.40"
clim_env_allSites$Site[clim_env_allSites$Site=="40.4"] <- "40.40"
clim_env_allSites$survYear = as.factor(clim_env_allSites$survYear)

# Add X column (quadrat by the 0.5)
allSpDat$survYear = as.factor(allSpDat$survYear)
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

m_AM_AMBA <- lmer(AM_BA ~ Mean_PDSI + Min_Temp + Max_Temp +
             perc_Openness + Base_CEC + 
             (1|Sp) + (1|survYear) + (1|SiteID),
             data = AM_dat)

m_AM_ECMBA <- lmer(ECM_BA ~ Mean_PDSI + Min_Temp + Max_Temp +
              perc_Openness + Base_CEC + 
              (1|Sp) + (1|survYear) + (1|SiteID),
              data = AM_dat)

m_AM_CONS <- lmer(ConspBA ~ Mean_PDSI + Min_Temp + Max_Temp +
                  perc_Openness + Base_CEC + 
                  (1|Sp) + (1|survYear) + (1|SiteID),
                  data = AM_dat)

## ================================
## 2) AM Seedling outcome model

m_AM_seedlings <- glmmTMB(Year_1_counts ~
    # Biotic main effects
      ConspBA + AM_BA + ECM_BA +
    # Biotic interactions (from best model)
      ConspBA:AM_BA +
      ConspBA:ECM_BA +
    # Abiotic direct effects
      Mean_PDSI + Min_Temp + Max_Temp +
      perc_Openness + Base_CEC +
    # Controls
      Log_Prior + Log_Recruits +
    # Random effects
      (1|Sp) + (1|survYear) + (1|SiteID),
  family = nbinom2, data = AM_dat)

## ================================
## 3) Assemble AM SEM

sem_AM <- psem(m_AM_AMBA, m_AM_ECMBA,
  m_AM_CONS, m_AM_seedlings,
  AM_BA %~~% ECM_BA,
  AM_BA %~~% ConspBA,
  ECM_BA %~~% ConspBA)

summary(sem_AM, .progress = FALSE)

## ================================
## 4) Extract coefficients from piecewise SEM
coef_table_AM <- coefs(sem_AM) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  filter( !grepl("^~~", response),
    !grepl("^~~", predictor)) %>%
  select(response, predictor, estimate,
    std_error, std_estimate, df,p_value) %>%
   arrange(response, predictor)

## ================================
## 5) Identify numeric variables that actually exist in data
numeric_vars <- AM_dat %>%
  select(where(is.numeric)) %>%
  names()

vars_for_sd <- union(coef_table_AM$predictor,
  coef_table_AM$response) %>%
  intersect(numeric_vars)

## ================================
## 6) Compute SDs for valid predictors and responses
SDs <- AM_dat %>%
  select(all_of(vars_for_sd)) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(),
    names_to = "variable", values_to = "sd")

SD_lookup <- setNames(SDs$sd, SDs$variable)

coef_table_AM <- coef_table_AM %>%
  mutate(
    estimate     = as.numeric(estimate),
    std_error    = as.numeric(std_error),
    std_estimate = as.numeric(std_estimate),
    df           = as.numeric(df),
    p_value      = as.numeric(p_value))

## ================================
## 7) Compute standardized standard errors
coef_table_AM <- coef_table_AM %>%
  rowwise() %>%
  mutate(std_error_std = if_else(
      predictor %in% names(SD_lookup) &
        response  %in% names(SD_lookup),
      std_error * (SD_lookup[predictor] / SD_lookup[response]),
      NA_real_ )) %>%
  ungroup()

## ================================
## 8) View final table
View(coef_table_AM)

############################################
############################################
############ ECM SEM

## ================================
## 1) Biotic mediator models (Abiotic → Biotic)

m_ECM_AMBA <- lmer(AM_BA ~ perc_Openness + Base_CEC +
    (1 | Sp) + (1 | survYear) + (1 | SiteID),
    data = ECM_dat)

m_ECM_ECMBA <- lmer(ECM_BA ~ perc_Openness + Base_CEC +
                      (1 | Sp) + (1 | survYear) + (1 | SiteID),
                    data = ECM_dat)

m_ECM_CONS <- lmer(ConspBA ~ perc_Openness + Base_CEC +
                     (1 | Sp) + (1 | survYear) + (1 | SiteID),
                   data = ECM_dat)

## ================================
## 2) ECM Seedling outcome model

m_ECM_seedlings <- glmer(Year_1_counts ~
    # Biotic main effects
      ConspBA + AM_BA + ECM_BA +
    # Biotic interactions (supported)
      ConspBA:AM_BA +
      ConspBA:ECM_BA +
    # Abiotic direct effects (best model)
      perc_Openness + Base_CEC +
      # Controls
      Log_Prior + Log_Recruits +
      # Random effects
      (1|Sp) + (1|survYear) + (1|SiteID),
  family = poisson(link = "log"),
  data = ECM_dat)

## ================================
## 3) Assemble AM SEM

sem_ECM <- psem(m_ECM_AMBA, m_ECM_ECMBA,
                m_ECM_CONS, m_ECM_seedlings,
               AM_BA %~~% ECM_BA,
               AM_BA %~~% ConspBA,
               ECM_BA %~~% ConspBA)

summary(sem_ECM, .progress = FALSE)

## ================================
## 4) Extract coefficients from piecewise SEM
coef_table_ECM <- coefs(sem_ECM) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  filter( !grepl("^~~", response),
          !grepl("^~~", predictor)) %>%
  select(response, predictor, estimate,
         std_error, std_estimate, df,p_value) %>%
  arrange(response, predictor)

## ================================
## 5) Identify numeric variables that actually exist in data
numeric_vars <- ECM_dat %>%
  select(where(is.numeric)) %>%
  names()

vars_for_sd <- union(coef_table_ECM$predictor,
                     coef_table_ECM$response) %>%
  intersect(numeric_vars)

## ================================
## 6) Compute SDs for valid predictors and responses
SDs <- ECM_dat %>%
  select(all_of(vars_for_sd)) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable", values_to = "sd")

SD_lookup <- setNames(SDs$sd, SDs$variable)

coef_table_ECM <- coef_table_ECM %>%
  mutate(
    estimate     = as.numeric(estimate),
    std_error    = as.numeric(std_error),
    std_estimate = as.numeric(std_estimate),
    df           = as.numeric(df),
    p_value      = as.numeric(p_value))

## ================================
## 7) Compute standardized standard errors
coef_table_ECM <- coef_table_ECM %>%
  rowwise() %>%
  mutate(std_error_std = if_else(
    predictor %in% names(SD_lookup) &
      response  %in% names(SD_lookup),
    std_error * (SD_lookup[predictor] / SD_lookup[response]),
    NA_real_ )) %>%
  ungroup()

## ================================
## 8) View final table
coef_table_ECM

#################### #################### #################### 
#################### Plotting the SEMs

#################### AM
## Define nodes and spacing
nodes <- data.frame(
  name = c("Mean_PDSI", "Min_Temp", "Max_Temp", "perc_Openness", "Base_CEC",
    "AM_BA", "ECM_BA", "ConspBA",
    "Year_1_counts",
    "Log_Prior", "Log_Recruits",
    "ConspBA:AM_BA", "ConspBA:ECM_BA"),
  x = c(rep(0.10, 5), rep(0.45, 3), 0.90, rep(0.70, 2), rep(0.70, 2)),
  y = c(seq(1.05, 0.15, length.out = 5),
    c(0.80, 0.50, 0.20),
    0.50,
    c(0.15, 0.35),
    c(0.6, 0.8)),
  fill = c(
    rep("orange", 5),
    rep("#9dc183", 3),
    "white",
    rep("gray", 2),
    rep("pink", 2)))
## Define edges
edges <- data.frame(
  from = c("Mean_PDSI","Min_Temp","Max_Temp","perc_Openness","Base_CEC",
    "Mean_PDSI","Min_Temp","Max_Temp","perc_Openness","Base_CEC",
    "Mean_PDSI","Min_Temp","Max_Temp","perc_Openness","Base_CEC",
    "Mean_PDSI","Min_Temp","Max_Temp","perc_Openness","Base_CEC",
    "AM_BA","ECM_BA","ConspBA",
    "Log_Prior","Log_Recruits",
    "ConspBA:AM_BA","ConspBA:ECM_BA"),
  to = c(rep("AM_BA",5), rep("ECM_BA",5),
    rep("ConspBA",5), rep("Year_1_counts",5), rep("Year_1_counts",3),
    rep("Year_1_counts",2),"Year_1_counts","Year_1_counts"))
## Define node radius
r_node <- 0.065
## 
clip_to_circle <- function(x0, y0, x1, y1, r) {
  dx <- x1 - x0
  dy <- y1 - y0
  d  <- sqrt(dx^2 + dy^2)
  c(x = x0 + dx / d * r,
    y = y0 + dy / d * r)
}
## Define plot space
plot(NA, xlim = c(0, 1), ylim = c(0, 1.1),
  axes = FALSE, xlab = "",ylab = "")
## Specify arrow plotting
for (i in seq_len(nrow(edges))) {
  
  from <- nodes[nodes$name == edges$from[i], ]
  to   <- nodes[nodes$name == edges$to[i], ]
  
  start <- clip_to_circle(
    x0 = from$x, y0 = from$y,
    x1 = to$x,   y1 = to$y,
    r  = r_node)
  
  end <- clip_to_circle(
    x0 = to$x, y0 = to$y,
    x1 = from$x, y1 = from$y,
    r  = r_node)
  
  arrows(x0 = start["x"], y0 = start["y"],
    x1 = end["x"],   y1 = end["y"],
    length = 0.08,
    lwd = 1)
}
## Add circle nodes
symbols(x = nodes$x,
  y = nodes$y,circles = rep(r_node, nrow(nodes)),
  inches = FALSE,bg = nodes$fill,
  fg = "black",add = TRUE)
## Add text labels
text(nodes$x, nodes$y, labels = nodes$name, cex = 0.65)


#################### ECM
## Define nodes and spacing
nodes <- data.frame(
  name = c(# Abiotic drivers (ECM model only)
    "perc_Openness", "Base_CEC",
    # Mediators
    "AM_BA", "ECM_BA", "ConspBA",
    # Response
    "Year_1_counts",
    # Controls
    "Log_Prior", "Log_Recruits",
    # Interactions
    "ConspBA:AM_BA", "ConspBA:ECM_BA"),
  x = c(rep(0.10, 2),   # abiotic
    rep(0.45, 3),   # mediators
    0.90,           # response
    rep(0.70, 2),   # controls
    rep(0.70, 2) ),    # interactions
  y = c(c(0.70, 0.30),      # abiotic (spread vertically)
    c(0.80, 0.50, 0.20),# mediators
    0.50,               # response
    c(0.15, 0.35),      # controls
    c(0.65, 0.85)),      # interactions
  fill = c(
    rep("orange", 2),   # abiotic
    rep("#9dc183", 3),  # mediators
    "white",            # response
    rep("gray", 2),     # controls
    rep("pink", 2)      # interactions
  )
)
## Define edges
edges <- data.frame(
  from = c(
    # Abiotic → mediators
    "perc_Openness","Base_CEC",
    "perc_Openness","Base_CEC",
    "perc_Openness","Base_CEC",
    # Mediators → response
    "AM_BA","ECM_BA","ConspBA",
    # Abiotic → response
    "perc_Openness","Base_CEC",
    # Controls → response
    "Log_Prior","Log_Recruits",
    # Interactions → response
    "ConspBA:AM_BA","ConspBA:ECM_BA"),
  to = c( "AM_BA","AM_BA",
    "ECM_BA","ECM_BA",
    "ConspBA","ConspBA",
    rep("Year_1_counts", 9)))

## Define node radius
r_node <- 0.065
## 
clip_to_circle <- function(x0, y0, x1, y1, r) {
  dx <- x1 - x0
  dy <- y1 - y0
  d  <- sqrt(dx^2 + dy^2)
  c(x = x0 + dx / d * r,
    y = y0 + dy / d * r)
}
## Define plot space
plot(NA, xlim = c(0, 1), ylim = c(0, 1.1),
     axes = FALSE, xlab = "",ylab = "")
## Specify arrow plotting
for (i in seq_len(nrow(edges))) {
  
  from <- nodes[nodes$name == edges$from[i], ]
  to   <- nodes[nodes$name == edges$to[i], ]
  
  start <- clip_to_circle(
    x0 = from$x, y0 = from$y,
    x1 = to$x,   y1 = to$y,
    r  = r_node)
  
  end <- clip_to_circle(
    x0 = to$x, y0 = to$y,
    x1 = from$x, y1 = from$y,
    r  = r_node)
  
  arrows(x0 = start["x"], y0 = start["y"],
    x1 = end["x"],   y1 = end["y"],
    length = 0.08,
    lwd = 1)
}

## Add circle nodes
symbols(x = nodes$x, y = nodes$y,
  circles = rep(r_node, nrow(nodes)),
  inches = FALSE,bg = nodes$fill,
  fg = "black",add = TRUE)
# Text labels
text(nodes$x, nodes$y, labels = nodes$name, cex = 0.65)











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





