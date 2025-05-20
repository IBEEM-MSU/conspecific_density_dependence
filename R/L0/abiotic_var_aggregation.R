
#################################################
####### Derivation of environmental covariates needed for abiotic models of survivorship
####### Author: Bailey McNichol 
####### Date: 20 May 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(dplyr)
library(ggplot2)

##### Read in csv files with currently derived environmental covariates
envVars <- read.csv("./data/L0/quad_envVars_allSites.csv", row.names = 1)
## Not variables we need because they are already pre-summarized
climVars <- read.csv("./data/L0/quad_climateData_allSites.csv", row.names = 1)

# Subset to only include the relevant columns
envVars_mod = envVars[,c(1:3,5:6,13)]
# Multiple percOpen by 100 so it's actually a percentage
envVars_mod$percOpen = envVars_mod$percOpen*100

head(climVars)

