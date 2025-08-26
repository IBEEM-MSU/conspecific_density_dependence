
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
library(DiagrammeR)

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

### Create binary (0/1) variable for Seedling_myc
allSpDat_surv_clim$Sdlg_myc_bin <- ifelse(allSpDat_surv_clim$Seedling_myc 
                                          == "AM", 1, 0)

############################################
############################################
############ Fitting a structural equation model integrating biotic and abiotic components!

##### FULL VERSION WE WANT TO FIT
dat <- allSpDat_surv_clim[,c(1,5:7,14,15,34,38:40,
                             44,50,53,55:58)]
# Rename columns
colnames(dat) <- c("SiteID","Sp","Year_0_counts",
            "Year_1_counts","AM_BA","ECM_BA",
            "Seedling_myc", "Consp_BA", "Log_Recruits",
            "Log_Prior","survYear", "Mean_PDSI","perc_Openness",
             "Base_CEC", "Min_Temp", "Max_Temp","Myc_type") 

## 1) Biotic mediator models (Abiotic → Biotic)
m_AM <- lmer(AM_BA ~ Mean_PDSI + Min_Temp + Max_Temp + 
               perc_Openness + Base_CEC +
              (1|Sp) + (1|survYear) + (1|SiteID), data = dat)

m_ECM <- lmer(ECM_BA ~ Mean_PDSI + Min_Temp + Max_Temp + 
                perc_Openness + Base_CEC +
               (1|Sp) + (1|survYear) + (1|SiteID), data = dat)

m_CON <- lmer(Consp_BA ~ Mean_PDSI + Min_Temp + Max_Temp + 
                perc_Openness + Base_CEC  +
                (1|Sp) + (1|survYear) + (1|SiteID), data = dat)

## 2) Seedling outcome model (Biotic + Abiotic + Trait + Interactions + Controls)
#### Outcome model?
m_seedlings <- glmer(Year_1_counts ~
    # Biotic main effects
    Consp_BA + Myc_type +  AM_BA + ECM_BA +
    # Biotic interactions
    #consp_15m:AM_15m + 
      Consp_BA:ECM_BA + Myc_type:AM_BA + Myc_type:ECM_BA +
    # Abiotic direct effects
    Mean_PDSI + Min_Temp + Max_Temp + perc_Openness + Base_CEC +
    # Controls
      Log_Prior + Log_Recruits + 
    # Random effects
      (1|Sp) + (1|survYear) + (1|SiteID), 
  family = poisson(link = "log"), data = dat)

## 3) Assemble the SEM
sem_mod <- psem(
  m_AM, m_ECM, m_CON,  # mediators
  m_seedlings)        # outcome

## 4) Inspect fit and paths
summary(sem_mod, .progress = FALSE, .verbose = FALSE)

# Super basic plot
plot(sem_mod)

### More customized plot

plot(sem_mod, node_attrs = list(
  shape = "oval",
  width = 0.7,  # shrink node width
  height = 0.28, # shrink node height
  fillcolor = c(
    rep("#9dc183", 3),   # Biotic - Neighborhood BA
    rep("white", 1),      # Response variable
    rep("orange", 5),      # Abiotic variables
    rep("gray", 3),      # Controls - logPrior and log0yr
    rep("pink", 3)        # Interactions
  ), outline = "black", color = "black",
  fontsize = 5, fontcolor = "black"))



data(keeley)

mod <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  lm(firesev ~ age, data=keeley),
  data = keeley
)

plot(mod)

### More customized plot

plot(mod, node_attrs = list(
  shape = "rectangle", color = "black",
  fillcolor = "orange", #x = 3, y=1:4
  ))




## Fancier plot!

# 1. Extract all significant paths
paths <- coefs(sem_mod)  # this gives you a data frame of all paths with estimates and p-values

# 2. Create edge list (only significant paths if desired)
edges <- paths[, c("Predictor", "Response", "Estimate")]
#edges <- edges[edges$P.Value < 0.05, ]  # optional: filter by significance

# 3. Make igraph object
g <- graph_from_data_frame(edges, directed = TRUE)

# 4. Set node colors by type (you can define a vector)
V(g)$color <- "lightblue"
V(g)$size <- 30
V(g)$label.cex <- 0.5

# 5. Edge width scaled by coefficient magnitude
E(g)$width <- abs(E(g)$Estimate) * 5

# 6. Plot
plot(g, vertex.label.color = "black", layout = layout_as_tree(g))

######## Maybe??

# --- 1) Define nodes ---
nodes <- create_node_df(
  n = 12,
  label = c("mean_pdsi", "mintemp", "maxtemp", "percOpen", "Base_CEC", 
            "AM_15m", "ECM_15m", "consp_15m", 
            "Sdlg_myc_bin", "logPrior", "log0yr",
            "Year_1_counts"),
  type = c(rep("Abiotic", 5), rep("Mediator", 3), 
           "Trait", rep("Control", 2), "Outcome")
)

# --- 2) Define edges ---
edges <- create_edge_df(
  from = c(1,2,3,4,5, 1,2,3,4,5, 1,2,3,4,5, 6,7,8,9, 8,7,7, 1,2,3,4,5, 10,11),
  to   = c(6,6,6,6,6, 7,7,7,7,7, 8,8,8,8,8, 12,12,12,12, 12,12,12, 12,12,12,12,12, 12,12)
)

# --- 3) Build graph ---
graph <- create_graph(
  nodes_df = nodes,
  edges_df = edges,
  directed = TRUE
)

# --- 4) Add colors to nodes ---
graph <- graph %>%
  set_node_attrs(node_attr = "color",
                 values = c(rep("lightblue",5), rep("orange",3), "green", rep("purple",2), "red")) %>%
  set_node_attrs(node_attr = "style", values = "filled")

# --- 5) Render ---
render_graph(graph, layout = "tree")



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


