
data_25.43 <- read.csv("./data/L0/all_SpYears_25.43.csv", row.names = 1)

head(data_25.43)

rm <- read.csv("./data/L0/rm_allCohortsSites.csv", row.names = 1)

head(rm)

envVars <- read.csv("./data/L0/quad_envVars_allSites.csv", row.names = 1)

head(envVars)

climVars <- read.csv("./data/L0/quad_climateData_allSites.csv", row.names = 1)

head(climVars)

BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

head(BA_data)
