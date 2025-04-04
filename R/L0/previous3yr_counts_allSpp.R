
#################################################
####### Extracting previous 3 yrs of counts per quadrat/per species
####### Author: Bailey McNichol and Andrew Eagar
####### Date: 04 April 2025
#################################################

# Load Packages
library(tidyverse) # For data manipulation & graphing. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823


############################################################################
#####################  Red Maple
### Read in data file from the R/L0 folder 
rm_demography <- read.csv( "./data/L0/rm_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

rm_demography_prev3yr_dens <- rm_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(rm_demography_prev3yr_dens, "./data/L1/rm_prior_cohort_densities.csv",
#           row.names = FALSE)


############################################################################
#####################  Sugar Maple
### Read in data file from the R/L0 folder 
sm_demography <- read.csv( "./data/L0/sm_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

sm_demography_prev3yr_dens <- sm_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(sm_demography_prev3yr_dens, "./data/L1/sm_prior_cohort_densities.csv",
#           row.names = FALSE)


############################################################################
#####################  American Beech
### Read in data file from the R/L0 folder 
ab_demography <- read.csv( "./data/L0/ab_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

ab_demography_prev3yr_dens <- ab_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(ab_demography_prev3yr_dens, "./data/L1/ab_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  White Ash
### Read in data file from the R/L0 folder 
wa_demography <- read.csv( "./data/L0/wa_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

wa_demography_prev3yr_dens <- wa_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

 # write.csv(wa_demography_prev3yr_dens, "./data/L1/wa_prior_cohort_densities.csv",
 #          row.names = FALSE)

############################################################################
#####################  Black Cherry
### Read in data file from the R/L0 folder 
bc_demography <- read.csv( "./data/L0/bc_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

bc_demography_prev3yr_dens <- bc_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(bc_demography_prev3yr_dens, "./data/L1/bc_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  Northern Red Oak
### Read in data file from the R/L0 folder 
nro_demography <- read.csv( "./data/L0/nro_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

nro_demography_prev3yr_dens <- nro_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(nro_demography_prev3yr_dens, "./data/L1/nro_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  Black Oak
### Read in data file from the R/L0 folder 
bonro_demography <- read.csv( "./data/L0/bonro_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

bonro_demography_prev3yr_dens <- bonro_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(bonro_demography_prev3yr_dens, "./data/L1/bonro_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  White Oak
### Read in data file from the R/L0 folder 
wo_demography <- read.csv( "./data/L0/wo_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

wo_demography_prev3yr_dens <- wo_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(wo_demography_prev3yr_dens, "./data/L1/wo_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  Ironwood
### Read in data file from the R/L0 folder 
iw_demography <- read.csv( "./data/L0/iw_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

iw_demography_prev3yr_dens <- iw_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(iw_demography_prev3yr_dens, "./data/L1/iw_prior_cohort_densities.csv",
#           row.names = FALSE)

############################################################################
#####################  Basswood
### Read in data file from the R/L0 folder 
bs_demography <- read.csv( "./data/L0/bs_allCohortsSites.csv", row.names = 1) %>%
  replace(is.na(.), 0) # replace the NAs with 0s.

bs_demography_prev3yr_dens <- bs_demography %>% # Calculate existing seedling densities over the previous 3 years 
  filter(sumAll > 0) %>% # Remove quadrats with no observed seedlings
  gather(c(11:37), key = "Obs.Year", value = "Count") %>% # just pull columns of counts
  separate(Obs.Year, into = c("X", "Obs.Year"), sep = "X") %>% # drop the leading X
  filter(First.Year == as.numeric(Obs.Year) - 1 | 
           First.Year == as.numeric(Obs.Year) - 2 | 
           First.Year == as.numeric(Obs.Year) - 3) %>% # pull previous 3 years
  group_by(Site, quadrat, Obs.Year) %>% # group into sites, quadrats, and observed year
  mutate(Prior_3_year_counts = sum(Count)) %>% # sum the previous 3 years counts
  filter(Obs.Year > 2000) %>% # exclude years before 2001 (cannot sum 3 prev years)
  select(c("TransectID", "Site", "SiteID",
           "quadrat", "Obs.Year", "Prior_3_year_counts")) %>% # just columns we want
  unique()

# write.csv(bs_demography_prev3yr_dens, "./data/L1/bs_prior_cohort_densities.csv",
#           row.names = FALSE)

