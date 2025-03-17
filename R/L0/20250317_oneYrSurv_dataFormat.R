
#################################################
####### Creation of data files for analyses of one-year survivorship WRT mature neighborhood
####### Author: Bailey McNichol
####### Date: 17 March 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(dplyr)
library(fuzzyjoin) # for merging data

#################################################
#################################################
##### Red Maple
# Read in data 
rm <- read.csv("./data/L0/rm_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
rm <- rm %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
rm$total_count_yr0 <- rm[cbind(1:nrow(rm), 
                  match(paste0("X", rm$First.Year), names(rm)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
rm$total_count_yr1 <- rm[cbind(1:nrow(rm),  
                               match(paste0("X", rm$First.Year + 1), names(rm)))]

## Convert to numeric
rm$total_count_yr0 = as.numeric(rm$total_count_yr0)
rm$total_count_yr1 = as.numeric(rm$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
rm$total_count_yr1 <- ifelse(rm$total_count_yr1 > rm$total_count_yr0,
                             rm$total_count_yr0, rm$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
rm[,c(11:37)] <- rm[,c(11:37)]/apply(rm[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
rm$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(rm)) {
  fy_value <- rm$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    rm$Survivorship[i] <- rm[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
rm_data <- rm[rm$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
rm_data <- rm_data[rm_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
rm_data <- subset(rm_data, FY != 27)
# Check dimensions of cleaned data
dim(rm_data)

## Add a column that denotes the year for proportion of survivorship
rm_data$survYear = as.character(rm_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
rm_data$baYear = ifelse(rm_data$survYear == "1999", "1999",
                        ifelse(rm_data$survYear == "2000", "1999",
                        ifelse(rm_data$survYear == "2001", "1999",
                        ifelse(rm_data$survYear == "2002", "1999",
                        ifelse(rm_data$survYear == "2003", "2005",
                        ifelse(rm_data$survYear == "2004", "2005",
                        ifelse(rm_data$survYear == "2005", "2005",
                        ifelse(rm_data$survYear == "2006", "2005",
                        ifelse(rm_data$survYear == "2007", "2007",
                        ifelse(rm_data$survYear == "2008", "2007",
                        ifelse(rm_data$survYear == "2009", "2009",
                        ifelse(rm_data$survYear == "2010", "2009",
                        ifelse(rm_data$survYear == "2011", "2011",
                        ifelse(rm_data$survYear == "2012", "2011",
                        ifelse(rm_data$survYear == "2013", "2013",
                        ifelse(rm_data$survYear == "2014", "2013",
                        ifelse(rm_data$survYear == "2015", "2016",
                        ifelse(rm_data$survYear == "2016", "2016",
                        ifelse(rm_data$survYear == "2017", "2016",
                        ifelse(rm_data$survYear == "2018", "2019",
                        ifelse(rm_data$survYear == "2019", "2019",
                        ifelse(rm_data$survYear == "2020", "2019",
                        ifelse(rm_data$survYear == "2021", "2022",
                        ifelse(rm_data$survYear == "2022", "2022",
                        ifelse(rm_data$survYear == "2023", "2022",
                               "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                        "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                         "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                          "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                          "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                          "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "sm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "sm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                       "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "sm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into rm seedling data
columns_to_merge <- c("rm_5m","rm_10m","rm_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
rm_merged <- rm_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(rm_merged)
rm_surv <- rm_merged[,c(1:7,9,39:58)]

# Make Site a character 
rm_surv$Site = as.character(rm_surv$Site)

# Convert NAs in the basal area data to 0's
rm_surv[,14:28][is.na(rm_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(rm_surv, "./data/L1/rm_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### Sugar Maple
# Read in data 
sm <- read.csv("./data/L0/sm_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
sm <- sm %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
sm$total_count_yr0 <- sm[cbind(1:nrow(sm), 
                               match(paste0("X", sm$First.Year), names(sm)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
sm$total_count_yr1 <- sm[cbind(1:nrow(sm),  
                               match(paste0("X", sm$First.Year + 1), names(sm)))]

## Convert to numeric
sm$total_count_yr0 = as.numeric(sm$total_count_yr0)
sm$total_count_yr1 = as.numeric(sm$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
sm$total_count_yr1 <- ifelse(sm$total_count_yr1 > sm$total_count_yr0,
                             sm$total_count_yr0, sm$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
sm[,c(11:37)] <- sm[,c(11:37)]/apply(sm[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
sm$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(sm)) {
  fy_value <- sm$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    sm$Survivorship[i] <- sm[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
sm_data <- sm[sm$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
sm_data <- sm_data[sm_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
sm_data <- subset(sm_data, FY != 27)
# Check dimensions of cleaned data
dim(sm_data)

## Add a column that denotes the year for proportion of survivorship
sm_data$survYear = as.character(sm_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
sm_data$baYear = ifelse(sm_data$survYear == "1999", "1999",
                 ifelse(sm_data$survYear == "2000", "1999",
                  ifelse(sm_data$survYear == "2001", "1999",
                 ifelse(sm_data$survYear == "2002", "1999",
               ifelse(sm_data$survYear == "2003", "2005",
              ifelse(sm_data$survYear == "2004", "2005",
               ifelse(sm_data$survYear == "2005", "2005",
               ifelse(sm_data$survYear == "2006", "2005",
               ifelse(sm_data$survYear == "2007", "2007",
            ifelse(sm_data$survYear == "2008", "2007",
           ifelse(sm_data$survYear == "2009", "2009",
             ifelse(sm_data$survYear == "2010", "2009",
          ifelse(sm_data$survYear == "2011", "2011",
   ifelse(sm_data$survYear == "2012", "2011",
       ifelse(sm_data$survYear == "2013", "2013",
               ifelse(sm_data$survYear == "2014", "2013",
                   ifelse(sm_data$survYear == "2015", "2016",
                             ifelse(sm_data$survYear == "2016", "2016",
                                  ifelse(sm_data$survYear == "2017", "2016",
                                        ifelse(sm_data$survYear == "2018", "2019",
                             ifelse(sm_data$survYear == "2019", "2019",
                                    ifelse(sm_data$survYear == "2020", "2019",
                               ifelse(sm_data$survYear == "2021", "2022",
                                     ifelse(sm_data$survYear == "2022", "2022",
               ifelse(sm_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("rm_5m", "ae_5m", "wa_5m", 
                                          "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("rm_10m", "ae_10m", "wa_10m", 
                                           "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("rm_15m", "ae_15m", "wa_15m", 
                                           "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into sm seedling data
columns_to_merge <- c("sm_5m","sm_10m","sm_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
sm_merged <- sm_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(sm_merged)
sm_surv <- sm_merged[,c(1:7,9,39:58)]

# Make Site a character 
sm_surv$Site = as.character(sm_surv$Site)

# Convert NAs in the basal area data to 0's
sm_surv[,14:28][is.na(sm_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(sm_surv, "./data/L1/sm_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### American Beech
# Read in data 
ab <- read.csv("./data/L0/ab_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
ab <- ab %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
ab$total_count_yr0 <- ab[cbind(1:nrow(ab), 
                               match(paste0("X", ab$First.Year), names(ab)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
ab$total_count_yr1 <- ab[cbind(1:nrow(ab),  
                               match(paste0("X", ab$First.Year + 1), names(ab)))]

## Convert to numeric
ab$total_count_yr0 = as.numeric(ab$total_count_yr0)
ab$total_count_yr1 = as.numeric(ab$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
ab$total_count_yr1 <- ifelse(ab$total_count_yr1 > ab$total_count_yr0,
                             ab$total_count_yr0, ab$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
ab[,c(11:37)] <- ab[,c(11:37)]/apply(ab[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
ab$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(ab)) {
  fy_value <- ab$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    ab$Survivorship[i] <- ab[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
ab_data <- ab[ab$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
ab_data <- ab_data[ab_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
ab_data <- subset(ab_data, FY != 27)
# Check dimensions of cleaned data
dim(ab_data)

## Add a column that denotes the year for proportion of survivorship
ab_data$survYear = as.character(ab_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
ab_data$baYear = ifelse(ab_data$survYear == "1999", "1999",
           ifelse(ab_data$survYear == "2000", "1999",
       ifelse(ab_data$survYear == "2001", "1999",
           ifelse(ab_data$survYear == "2002", "1999",
        ifelse(ab_data$survYear == "2003", "2005",
           ifelse(ab_data$survYear == "2004", "2005",
       ifelse(ab_data$survYear == "2005", "2005",
       ifelse(ab_data$survYear == "2006", "2005",
     ifelse(ab_data$survYear == "2007", "2007",
       ifelse(ab_data$survYear == "2008", "2007",
         ifelse(ab_data$survYear == "2009", "2009",
       ifelse(ab_data$survYear == "2010", "2009",
     ifelse(ab_data$survYear == "2011", "2011",
      ifelse(ab_data$survYear == "2012", "2011",
         ifelse(ab_data$survYear == "2013", "2013",
           ifelse(ab_data$survYear == "2014", "2013",
           ifelse(ab_data$survYear == "2015", "2016",
            ifelse(ab_data$survYear == "2016", "2016",
              ifelse(ab_data$survYear == "2017", "2016",
               ifelse(ab_data$survYear == "2018", "2019",
                 ifelse(ab_data$survYear == "2019", "2019",
                 ifelse(ab_data$survYear == "2020", "2019",
                  ifelse(ab_data$survYear == "2021", "2022",
            ifelse(ab_data$survYear == "2022", "2022",
         ifelse(ab_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into ab seedling data
columns_to_merge <- c("ab_5m","ab_10m","ab_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
ab_merged <- ab_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(ab_merged)
ab_surv <- ab_merged[,c(1:7,9,39:58)]

# Make Site a character 
ab_surv$Site = as.character(ab_surv$Site)

# Convert NAs in the basal area data to 0's
ab_surv[,14:28][is.na(ab_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(ab_surv, "./data/L1/ab_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### White Ash
# Read in data 
wa <- read.csv("./data/L0/wa_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
wa <- wa %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
wa$total_count_yr0 <- wa[cbind(1:nrow(wa), 
                               match(paste0("X", wa$First.Year), names(wa)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
wa$total_count_yr1 <- wa[cbind(1:nrow(wa),  
                               match(paste0("X", wa$First.Year + 1), names(wa)))]

## Convert to numeric
wa$total_count_yr0 = as.numeric(wa$total_count_yr0)
wa$total_count_yr1 = as.numeric(wa$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
wa$total_count_yr1 <- ifelse(wa$total_count_yr1 > wa$total_count_yr0,
                             wa$total_count_yr0, wa$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
wa[,c(11:37)] <- wa[,c(11:37)]/apply(wa[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
wa$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(wa)) {
  fy_value <- wa$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    wa$Survivorship[i] <- wa[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
wa_data <- wa[wa$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
wa_data <- wa_data[wa_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
wa_data <- subset(wa_data, FY != 27)
# Check dimensions of cleaned data
dim(wa_data)

## Add a column that denotes the year for proportion of survivorship
wa_data$survYear = as.character(wa_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
wa_data$baYear = ifelse(wa_data$survYear == "1999", "1999",
       ifelse(wa_data$survYear == "2000", "1999",
        ifelse(wa_data$survYear == "2001", "1999",
       ifelse(wa_data$survYear == "2002", "1999",
     ifelse(wa_data$survYear == "2003", "2005",
      ifelse(wa_data$survYear == "2004", "2005",
   ifelse(wa_data$survYear == "2005", "2005",
   ifelse(wa_data$survYear == "2006", "2005",
    ifelse(wa_data$survYear == "2007", "2007",
    ifelse(wa_data$survYear == "2008", "2007",
   ifelse(wa_data$survYear == "2009", "2009",
     ifelse(wa_data$survYear == "2010", "2009",
      ifelse(wa_data$survYear == "2011", "2011",
      ifelse(wa_data$survYear == "2012", "2011",
        ifelse(wa_data$survYear == "2013", "2013",
           ifelse(wa_data$survYear == "2014", "2013",
           ifelse(wa_data$survYear == "2015", "2016",
         ifelse(wa_data$survYear == "2016", "2016",
          ifelse(wa_data$survYear == "2017", "2016",
          ifelse(wa_data$survYear == "2018", "2019",
          ifelse(wa_data$survYear == "2019", "2019",
         ifelse(wa_data$survYear == "2020", "2019",
            ifelse(wa_data$survYear == "2021", "2022",
              ifelse(wa_data$survYear == "2022", "2022",
             ifelse(wa_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("wa_5m", "ae_5m", "sm_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("wa_10m", "ae_10m", "sm_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("wa_15m", "ae_15m", "sm_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("rm_5m", "ae_5m", "sm_5m", 
                                          "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("rm_10m", "ae_10m", "sm_10m", 
                                           "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("rm_15m", "ae_15m", "sm_15m", 
                                           "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "sm_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "sm_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "sm_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into wa seedling data
columns_to_merge <- c("wa_5m","wa_10m","wa_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
wa_merged <- wa_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(wa_merged)
wa_surv <- wa_merged[,c(1:7,9,39:58)]

# Make Site a character 
wa_surv$Site = as.character(wa_surv$Site)

# Convert NAs in the basal area data to 0's
wa_surv[,14:28][is.na(wa_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(wa_surv, "./data/L1/wa_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### Black Cherry
# Read in data 
bc <- read.csv("./data/L0/bc_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
bc <- bc %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
bc$total_count_yr0 <- bc[cbind(1:nrow(bc), 
                               match(paste0("X", bc$First.Year), names(bc)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
bc$total_count_yr1 <- bc[cbind(1:nrow(bc),  
                               match(paste0("X", bc$First.Year + 1), names(bc)))]

## Convert to numeric
bc$total_count_yr0 = as.numeric(bc$total_count_yr0)
bc$total_count_yr1 = as.numeric(bc$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
bc$total_count_yr1 <- ifelse(bc$total_count_yr1 > bc$total_count_yr0,
                             bc$total_count_yr0, bc$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
bc[,c(11:37)] <- bc[,c(11:37)]/apply(bc[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
bc$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(bc)) {
  fy_value <- bc$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    bc$Survivorship[i] <- bc[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
bc_data <- bc[bc$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
bc_data <- bc_data[bc_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
bc_data <- subset(bc_data, FY != 27)
# Check dimensions of cleaned data
dim(bc_data)

## Add a column that denotes the year for proportion of survivorship
bc_data$survYear = as.character(bc_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
bc_data$baYear = ifelse(bc_data$survYear == "1999", "1999",
             ifelse(bc_data$survYear == "2000", "1999",
              ifelse(bc_data$survYear == "2001", "1999",
               ifelse(bc_data$survYear == "2002", "1999",
            ifelse(bc_data$survYear == "2003", "2005",
              ifelse(bc_data$survYear == "2004", "2005",
             ifelse(bc_data$survYear == "2005", "2005",
              ifelse(bc_data$survYear == "2006", "2005",
            ifelse(bc_data$survYear == "2007", "2007",
               ifelse(bc_data$survYear == "2008", "2007",
               ifelse(bc_data$survYear == "2009", "2009",
               ifelse(bc_data$survYear == "2010", "2009",
               ifelse(bc_data$survYear == "2011", "2011",
              ifelse(bc_data$survYear == "2012", "2011",
              ifelse(bc_data$survYear == "2013", "2013",
               ifelse(bc_data$survYear == "2014", "2013",
                ifelse(bc_data$survYear == "2015", "2016",
                 ifelse(bc_data$survYear == "2016", "2016",
                 ifelse(bc_data$survYear == "2017", "2016",
                     ifelse(bc_data$survYear == "2018", "2019",
                   ifelse(bc_data$survYear == "2019", "2019",
                    ifelse(bc_data$survYear == "2020", "2019",
                      ifelse(bc_data$survYear == "2021", "2022",
                      ifelse(bc_data$survYear == "2022", "2022",
                ifelse(bc_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("bc_5m", "ae_5m", "sm_5m", 
                                        "rm_5m", "wa_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("bc_10m", "ae_10m", "sm_10m", 
                                         "rm_10m", "wa_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("bc_15m", "ae_15m", "sm_15m", 
                                         "rm_15m", "wa_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("rm_5m", "ae_5m", "sm_5m", 
                                          "wa_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("rm_10m", "ae_10m", "sm_10m", 
                                           "wa_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("rm_15m", "ae_15m", "sm_15m", 
                                           "wa_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "sm_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","wa_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "sm_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "wa_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "sm_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "wa_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into bc seedling data
columns_to_merge <- c("bc_5m","bc_10m","bc_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
bc_merged <- bc_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(bc_merged)
bc_surv <- bc_merged[,c(1:7,9,39:58)]

# Make Site a character 
bc_surv$Site = as.character(bc_surv$Site)

# Convert NAs in the basal area data to 0's
bc_surv[,14:28][is.na(bc_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(bc_surv, "./data/L1/bc_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### White Oak
# Read in data 
wo <- read.csv("./data/L0/wo_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
wo <- wo %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
wo$total_count_yr0 <- wo[cbind(1:nrow(wo), 
                                 match(paste0("X", wo$First.Year), names(wo)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
wo$total_count_yr1 <- wo[cbind(1:nrow(wo),  
                                 match(paste0("X", wo$First.Year + 1), names(wo)))]

## Convert to numeric
wo$total_count_yr0 = as.numeric(wo$total_count_yr0)
wo$total_count_yr1 = as.numeric(wo$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
wo$total_count_yr1 <- ifelse(wo$total_count_yr1 > wo$total_count_yr0,
                              wo$total_count_yr0, wo$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
wo[,c(11:37)] <- wo[,c(11:37)]/apply(wo[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
wo$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(wo)) {
  fy_value <- wo$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    wo$Survivorship[i] <- wo[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
wo_data <- wo[wo$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
wo_data <- wo_data[wo_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
wo_data <- subset(wo_data, FY != 27)
# Check dimensions of cleaned data
dim(wo_data)

## Add a column that denotes the year for proportion of survivorship
wo_data$survYear = as.character(wo_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
wo_data$baYear = ifelse(wo_data$survYear == "1999", "1999",
                ifelse(wo_data$survYear == "2000", "1999",
                  ifelse(wo_data$survYear == "2001", "1999",
                  ifelse(wo_data$survYear == "2002", "1999",
                  ifelse(wo_data$survYear == "2003", "2005",
                   ifelse(wo_data$survYear == "2004", "2005",
                    ifelse(wo_data$survYear == "2005", "2005",
                    ifelse(wo_data$survYear == "2006", "2005",
                     ifelse(wo_data$survYear == "2007", "2007",
                    ifelse(wo_data$survYear == "2008", "2007",
                    ifelse(wo_data$survYear == "2009", "2009",
                      ifelse(wo_data$survYear == "2010", "2009",
                    ifelse(wo_data$survYear == "2011", "2011",
                      ifelse(wo_data$survYear == "2012", "2011",
                  ifelse(wo_data$survYear == "2013", "2013",
                   ifelse(wo_data$survYear == "2014", "2013",
                       ifelse(wo_data$survYear == "2015", "2016",
                      ifelse(wo_data$survYear == "2016", "2016",
                       ifelse(wo_data$survYear == "2017", "2016",
                    ifelse(wo_data$survYear == "2018", "2019",
                     ifelse(wo_data$survYear == "2019", "2019",
                      ifelse(wo_data$survYear == "2020", "2019",
                      ifelse(wo_data$survYear == "2021", "2022",
                      ifelse(wo_data$survYear == "2022", "2022",
                   ifelse(wo_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("bs_5m", "qa_5m","bta_5m", 
                                          "ab_5m","iw_5m", "bo_5m","bonro_5m",
                                          "nro_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("bs_10m", "qa_10m","bta_10m", 
                                           "ab_10m","iw_10m", "bo_10m","bonro_10m",
                                           "nro_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("bs_15m", "qa_15m","bta_15m", 
                                           "ab_15m","iw_15m", "bo_15m","bonro_15m",
                                           "nro_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "ab_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","nro_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "ab_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","nro_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "ab_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","nro_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into wo seedling data
columns_to_merge <- c("wo_5m","wo_10m","wo_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
wo_merged <- wo_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(wo_merged)
wo_surv <- wo_merged[,c(1:7,9,39:58)]

# Make Site a character 
wo_surv$Site = as.character(wo_surv$Site)

# Convert NAs in the basal area data to 0's
wo_surv[,14:28][is.na(wo_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(wo_surv, "./data/L1/wo_oneYrSurv_matureBA.csv")

#################################################
#################################################
##### Northern Red Oak
# Read in data 
nro <- read.csv("./data/L0/nro_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
nro <- nro %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
nro$total_count_yr0 <- nro[cbind(1:nrow(nro), 
                               match(paste0("X", nro$First.Year), names(nro)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
nro$total_count_yr1 <- nro[cbind(1:nrow(nro),  
                               match(paste0("X", nro$First.Year + 1), names(nro)))]

## Convert to numeric
nro$total_count_yr0 = as.numeric(nro$total_count_yr0)
nro$total_count_yr1 = as.numeric(nro$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
nro$total_count_yr1 <- ifelse(nro$total_count_yr1 > nro$total_count_yr0,
                             nro$total_count_yr0, nro$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
nro[,c(11:37)] <- nro[,c(11:37)]/apply(nro[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
nro$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(nro)) {
  fy_value <- nro$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    nro$Survivorship[i] <- nro[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
nro_data <- nro[nro$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
nro_data <- nro_data[nro_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
nro_data <- subset(nro_data, FY != 27)
# Check dimensions of cleaned data
dim(nro_data)

## Add a column that denotes the year for proportion of survivorship
nro_data$survYear = as.character(nro_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
nro_data$baYear = ifelse(nro_data$survYear == "1999", "1999",
                        ifelse(nro_data$survYear == "2000", "1999",
                               ifelse(nro_data$survYear == "2001", "1999",
                                      ifelse(nro_data$survYear == "2002", "1999",
                         ifelse(nro_data$survYear == "2003", "2005",
                     ifelse(nro_data$survYear == "2004", "2005",
                  ifelse(nro_data$survYear == "2005", "2005",
             ifelse(nro_data$survYear == "2006", "2005",
     ifelse(nro_data$survYear == "2007", "2007",
      ifelse(nro_data$survYear == "2008", "2007",
        ifelse(nro_data$survYear == "2009", "2009",
         ifelse(nro_data$survYear == "2010", "2009",
       ifelse(nro_data$survYear == "2011", "2011",
       ifelse(nro_data$survYear == "2012", "2011",
          ifelse(nro_data$survYear == "2013", "2013",
           ifelse(nro_data$survYear == "2014", "2013",
             ifelse(nro_data$survYear == "2015", "2016",
             ifelse(nro_data$survYear == "2016", "2016",
              ifelse(nro_data$survYear == "2017", "2016",
             ifelse(nro_data$survYear == "2018", "2019",
             ifelse(nro_data$survYear == "2019", "2019",
           ifelse(nro_data$survYear == "2020", "2019",
           ifelse(nro_data$survYear == "2021", "2022",
              ifelse(nro_data$survYear == "2022", "2022",
           ifelse(nro_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("bs_5m", "qa_5m","bta_5m", 
                                          "ab_5m","iw_5m", "bo_5m","bonro_5m",
                                          "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("bs_10m", "qa_10m","bta_10m", 
                                           "ab_10m","iw_10m", "bo_10m","bonro_10m",
                                           "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("bs_15m", "qa_15m","bta_15m", 
                                           "ab_15m","iw_15m", "bo_15m","bonro_15m",
                                           "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "ab_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "ab_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "ab_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into nro seedling data
columns_to_merge <- c("nro_5m","nro_10m","nro_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
nro_merged <- nro_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(nro_merged)
nro_surv <- nro_merged[,c(1:7,9,39:58)]

# Make Site a character 
nro_surv$Site = as.character(nro_surv$Site)

# Convert NAs in the basal area data to 0's
nro_surv[,14:28][is.na(nro_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(nro_surv, "./data/L1/nro_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### Black Oak/Northern Red Oak
# Read in data 
bonro <- read.csv("./data/L0/bonro_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
bonro <- bonro %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
bonro$total_count_yr0 <- bonro[cbind(1:nrow(bonro), 
                                 match(paste0("X", bonro$First.Year), names(bonro)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
bonro$total_count_yr1 <- bonro[cbind(1:nrow(bonro),  
                                 match(paste0("X", bonro$First.Year + 1), names(bonro)))]

## Convert to numeric
bonro$total_count_yr0 = as.numeric(bonro$total_count_yr0)
bonro$total_count_yr1 = as.numeric(bonro$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
bonro$total_count_yr1 <- ifelse(bonro$total_count_yr1 > bonro$total_count_yr0,
                              bonro$total_count_yr0, bonro$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
bonro[,c(11:37)] <- bonro[,c(11:37)]/apply(bonro[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
bonro$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(bonro)) {
  fy_value <- bonro$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    bonro$Survivorship[i] <- bonro[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
bonro_data <- bonro[bonro$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
bonro_data <- bonro_data[bonro_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
bonro_data <- subset(bonro_data, FY != 27)
# Check dimensions of cleaned data
dim(bonro_data)

## Add a column that denotes the year for proportion of survivorship
bonro_data$survYear = as.character(bonro_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
bonro_data$baYear = ifelse(bonro_data$survYear == "1999", "1999",
                         ifelse(bonro_data$survYear == "2000", "1999",
                                ifelse(bonro_data$survYear == "2001", "1999",
                                       ifelse(bonro_data$survYear == "2002", "1999",
                   ifelse(bonro_data$survYear == "2003", "2005",
                    ifelse(bonro_data$survYear == "2004", "2005",
                    ifelse(bonro_data$survYear == "2005", "2005",
                   ifelse(bonro_data$survYear == "2006", "2005",
                    ifelse(bonro_data$survYear == "2007", "2007",
                   ifelse(bonro_data$survYear == "2008", "2007",
                   ifelse(bonro_data$survYear == "2009", "2009",
                    ifelse(bonro_data$survYear == "2010", "2009",
               ifelse(bonro_data$survYear == "2011", "2011",
               ifelse(bonro_data$survYear == "2012", "2011",
               ifelse(bonro_data$survYear == "2013", "2013",
               ifelse(bonro_data$survYear == "2014", "2013",
                ifelse(bonro_data$survYear == "2015", "2016",
                  ifelse(bonro_data$survYear == "2016", "2016",
                  ifelse(bonro_data$survYear == "2017", "2016",
                  ifelse(bonro_data$survYear == "2018", "2019",
                    ifelse(bonro_data$survYear == "2019", "2019",
                    ifelse(bonro_data$survYear == "2020", "2019",
                     ifelse(bonro_data$survYear == "2021", "2022",
                      ifelse(bonro_data$survYear == "2022", "2022",
          ifelse(bonro_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "bonro_5m","iw_5m", "bo_5m","nro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                         "bonro_10m","iw_10m", "bo_10m","nro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "bonro_15m","iw_15m", "bo_15m","nro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("bs_5m", "qa_5m","bta_5m", 
                                          "ab_5m","iw_5m", "bo_5m","nro_5m",
                                          "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("bs_10m", "qa_10m","bta_10m", 
                                           "ab_10m","iw_10m", "bo_10m","nro_10m",
                                           "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("bs_15m", "qa_15m","bta_15m", 
                                           "ab_15m","iw_15m", "bo_15m","nro_15m",
                                           "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "bs_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "ab_5m","iw_5m", "bo_5m","bc_5m",
                                       "nro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "bs_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "ab_10m","iw_10m", "bo_10m",
                                        "bc_10m","nro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "bs_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "ab_15m","iw_15m", "bo_15m",
                                        "bc_15m","nro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into bonro seedling data
columns_to_merge <- c("bonro_5m","bonro_10m","bonro_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
bonro_merged <- bonro_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(bonro_merged)
bonro_surv <- bonro_merged[,c(1:7,9,39:58)]

# Make Site a character 
bonro_surv$Site = as.character(bonro_surv$Site)

# Convert NAs in the basal area data to 0's
bonro_surv[,14:28][is.na(bonro_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(bonro_surv, "./data/L1/bonro_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### American Basswood
# Read in data 
bs <- read.csv("./data/L0/bs_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
bs <- bs %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
bs$total_count_yr0 <- bs[cbind(1:nrow(bs), 
                               match(paste0("X", bs$First.Year), names(bs)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
bs$total_count_yr1 <- bs[cbind(1:nrow(bs),  
                               match(paste0("X", bs$First.Year + 1), names(bs)))]

## Convert to numeric
bs$total_count_yr0 = as.numeric(bs$total_count_yr0)
bs$total_count_yr1 = as.numeric(bs$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
bs$total_count_yr1 <- ifelse(bs$total_count_yr1 > bs$total_count_yr0,
                             bs$total_count_yr0, bs$total_count_yr1)

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
bs[,c(11:37)] <- bs[,c(11:37)]/apply(bs[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
bs$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(bs)) {
  fy_value <- bs$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    bs$Survivorship[i] <- bs[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
bs_data <- bs[bs$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
bs_data <- bs_data[bs_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
bs_data <- subset(bs_data, FY != 27)
# Check dimensions of cleaned data
dim(bs_data)

## Add a column that denotes the year for proportion of survivorship
bs_data$survYear = as.character(bs_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
bs_data$baYear = ifelse(bs_data$survYear == "1999", "1999",
                        ifelse(bs_data$survYear == "2000", "1999",
                               ifelse(bs_data$survYear == "2001", "1999",
                                      ifelse(bs_data$survYear == "2002", "1999",
                      ifelse(bs_data$survYear == "2003", "2005",
                    ifelse(bs_data$survYear == "2004", "2005",
           ifelse(bs_data$survYear == "2005", "2005",
              ifelse(bs_data$survYear == "2006", "2005",
               ifelse(bs_data$survYear == "2007", "2007",
                 ifelse(bs_data$survYear == "2008", "2007",
                   ifelse(bs_data$survYear == "2009", "2009",
                    ifelse(bs_data$survYear == "2010", "2009",
               ifelse(bs_data$survYear == "2011", "2011",
            ifelse(bs_data$survYear == "2012", "2011",
              ifelse(bs_data$survYear == "2013", "2013",
              ifelse(bs_data$survYear == "2014", "2013",
                ifelse(bs_data$survYear == "2015", "2016",
                  ifelse(bs_data$survYear == "2016", "2016",
                 ifelse(bs_data$survYear == "2017", "2016",
                   ifelse(bs_data$survYear == "2018", "2019",
                     ifelse(bs_data$survYear == "2019", "2019",
                      ifelse(bs_data$survYear == "2020", "2019",
                       ifelse(bs_data$survYear == "2021", "2022",
                      ifelse(bs_data$survYear == "2022", "2022",
                     ifelse(bs_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("bs_5m", "ab_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("bs_10m", "ab_10m", "qa_10m","bta_10m", 
                                         "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("bs_15m", "ab_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("ab_5m", "qa_5m","bta_5m", 
                                          "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                          "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("ab_10m", "qa_10m","bta_10m", 
                                           "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                           "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("ab_15m", "qa_15m","bta_15m", 
                                           "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                           "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "ab_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "nro_5m","iw_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "ab_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "nro_10m","iw_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "ab_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "nro_15m","iw_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into bs seedling data
columns_to_merge <- c("bs_5m","bs_10m","bs_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
bs_merged <- bs_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(bs_merged)
bs_surv <- bs_merged[,c(1:7,9,39:58)]

# Make Site a character 
bs_surv$Site = as.character(bs_surv$Site)

# Convert NAs in the basal area data to 0's
bs_surv[,14:28][is.na(bs_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(bs_surv, "./data/L1/bs_oneYrSurv_matureBA.csv")


#################################################
#################################################
##### Ironwood
# Read in data 
iw <- read.csv("./data/L0/iw_allCohortsSites.csv", row.names = 1)

# Replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
iw <- iw %>%
  replace(is.na(.), 0) 

## Add in the total count in the starting year of the cohort to use as weights in glm
iw$total_count_yr0 <- iw[cbind(1:nrow(iw), 
                               match(paste0("X", iw$First.Year), names(iw)))]

## Add in the total count after 1 year for each cohort to calculate survivorship
iw$total_count_yr1 <- iw[cbind(1:nrow(iw),  
                               match(paste0("X", iw$First.Year + 1), names(iw)))]

## Convert to numeric
iw$total_count_yr0 = as.numeric(iw$total_count_yr0)
iw$total_count_yr1 = as.numeric(iw$total_count_yr1)

## If the count in yr1 is greater than the count in yr0, set to yr0
iw$total_count_yr1 <- ifelse(iw$total_count_yr1 > iw$total_count_yr0,
                             iw$total_count_yr0, iw$total_count_yr1)

# Convert time by cohort oiwervation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row
iw[,c(11:37)] <- iw[,c(11:37)]/apply(iw[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
iw$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(iw)) {
  fy_value <- iw$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    iw$Survivorship[i] <- iw[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
iw_data <- iw[iw$sumAll != 0, ]
# Remove all rows that had zeroes in the first year per cohort
iw_data <- iw_data[iw_data$total_count_yr0 != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
iw_data <- subset(iw_data, FY != 27)
# Check dimensions of cleaned data
dim(iw_data)

## Add a column that denotes the year for proportion of survivorship
iw_data$survYear = as.character(iw_data$First.Year + 1)

## Use an ifelse statement to add a column with the closest BA census year
iw_data$baYear = ifelse(iw_data$survYear == "1999", "1999",
                        ifelse(iw_data$survYear == "2000", "1999",
                               ifelse(iw_data$survYear == "2001", "1999",
                  ifelse(iw_data$survYear == "2002", "1999",
                  ifelse(iw_data$survYear == "2003", "2005",
                  ifelse(iw_data$survYear == "2004", "2005",
                  ifelse(iw_data$survYear == "2005", "2005",
                 ifelse(iw_data$survYear == "2006", "2005",
                  ifelse(iw_data$survYear == "2007", "2007",
                 ifelse(iw_data$survYear == "2008", "2007",
                    ifelse(iw_data$survYear == "2009", "2009",
                     ifelse(iw_data$survYear == "2010", "2009",
         ifelse(iw_data$survYear == "2011", "2011",
          ifelse(iw_data$survYear == "2012", "2011",
             ifelse(iw_data$survYear == "2013", "2013",
               ifelse(iw_data$survYear == "2014", "2013",
              ifelse(iw_data$survYear == "2015", "2016",
                 ifelse(iw_data$survYear == "2016", "2016",
                ifelse(iw_data$survYear == "2017", "2016",
                  ifelse(iw_data$survYear == "2018", "2019",
                ifelse(iw_data$survYear == "2019", "2019",
               ifelse(iw_data$survYear == "2020", "2019",
                ifelse(iw_data$survYear == "2021", "2022",
                 ifelse(iw_data$survYear == "2022", "2022",
         ifelse(iw_data$survYear == "2023", "2022",  "2022")))))))))))))))))))))))))

## Read in basal area data
BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

## Calculate species-specific conmycorrhizal, heterospecific, and heteromycorrhizal BA neighborhoods
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH species included
BA_data$conMyc_5m = rowSums(BA_data[, c("iw_5m", "ab_5m", "qa_5m","bta_5m", 
                                        "nro_5m","bs_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("iw_10m", "ab_10m", "qa_10m","bta_10m", 
                                         "nro_10m","bs_10m", "bo_10m","bonro_10m",
                                         "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("iw_15m", "ab_15m", "qa_15m","bta_15m", 
                                         "nro_15m","bs_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITHOUT species included
BA_data$conMycNS_5m = rowSums(BA_data[, c("ab_5m", "qa_5m","bta_5m", 
                                          "nro_5m","bs_5m", "bo_5m","bonro_5m",
                                          "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
BA_data$conMycNS_10m = rowSums(BA_data[, c("ab_10m", "qa_10m","bta_10m", 
                                           "nro_10m","bs_10m", "bo_10m","bonro_10m",
                                           "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
BA_data$conMycNS_15m = rowSums(BA_data[, c("ab_15m", "qa_15m","bta_15m", 
                                           "nro_15m","bs_15m", "bo_15m","bonro_15m",
                                           "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
BA_data$hetMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$hetMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                         "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$hetMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heterospecific BA in a 5-m, 10-m, and 15-m radius
BA_data$hetSp_5m = rowSums(BA_data[, c("sm_5m", "ab_5m", "rm_5m", "qa_5m","ae_5m", "wa_5m",
                                       "bta_5m", "nro_5m","bs_5m", "bo_5m","bc_5m",
                                       "bonro_5m","wo_5m", "hem_5m","wp_5m",
                                       "jp_5m","svb_5m" )],  na.rm = TRUE)
BA_data$hetSp_10m = rowSums(BA_data[, c("sm_10m", "ab_10m", "rm_10m", "qa_10m","ae_10m",
                                        "wa_10m","bta_10m", "nro_10m","bs_10m", "bo_10m",
                                        "bc_10m","bonro_10m","wo_10m", "hem_10m","wp_10m",
                                        "jp_10m","svb_10m" )],  na.rm = TRUE)
BA_data$hetSp_15m = rowSums(BA_data[, c("sm_15m", "ab_15m", "rm_15m", "qa_15m","ae_15m",
                                        "wa_15m","bta_15m", "nro_15m","bs_15m", "bo_15m",
                                        "bc_15m","bonro_15m","wo_15m", "hem_15m","wp_15m",
                                        "jp_15m","svb_15m" )],  na.rm = TRUE)

# Make year a character
BA_data$Year = as.character(BA_data$Year)

## Merge the BA in the neighborhood with the seedling data using TransectID and baYear
# Specify columns to merge from BA data into iw seedling data
columns_to_merge <- c("iw_5m","iw_10m","iw_15m",
                      "conMyc_5m", "conMyc_10m", "conMyc_15m", 
                      "conMycNS_5m", "conMycNS_10m", "conMycNS_15m",
                      "hetMyc_5m", "hetMyc_10m", "hetMyc_15m",
                      "hetSp_5m", "hetSp_10m", "hetSp_15m")

# Perform the merge
iw_merged <- iw_data %>%
  left_join(BA_data %>% select(TransectID, Year, all_of(columns_to_merge)), 
            by = c("TransectID" = "TransectID", "baYear" = "Year"))

## Filter to only include columns we need for analyses
names(iw_merged)
iw_surv <- iw_merged[,c(1:7,9,39:58)]

# Make Site a character 
iw_surv$Site = as.character(iw_surv$Site)

# Convert NAs in the basal area data to 0's
iw_surv[,14:28][is.na(iw_surv[,14:28])] <- 0

### Write data file for species to the L1 folder
write.csv(iw_surv, "./data/L1/iw_oneYrSurv_matureBA.csv")

