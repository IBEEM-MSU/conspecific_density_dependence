
## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)

# Read in data for rm
rm <- read.csv("./data/L0/rm_allCohortsSites.csv", row.names = 1)
head(rm)

# This replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
rm <- rm %>%
  replace(is.na(.), 0) 

# Convert time by cohort observation matrix into % of seedlings remaining in each quadrat based on the highest number of seedlings recorded in that row. This helps account for missed observations in early cohort years (yr 1, 10 seedlings -> yr 2, 11 seedlings) and frees the data from its more linear structure.
rm[,c(11:37)] <- rm[,c(11:37)]/apply(rm[,c(11:37)], 1, max)

# Create column called Survivorship and fill with NAs
rm$Survivorship <- NA  # Initialize the new column

# Loop over rows and fill the Survivorship column based on FY
# For each row
for(i in 1:nrow(rm)) {
  fy_value <- rm$FY[i]  # Get the FY value for the current row
  if (fy_value >= 1 & fy_value <= 26) { # Only do this for survivorship up to 2024 (can't calculate for 2024 data)
    # Map FY to the correct year column (X1999 to X2024)
    year_column <- paste0("X", 1998 + fy_value)  # Construct the column name based on FY
    rm$Survivorship[i] <- rm[i, year_column]  # Assign the proportion from the correct year column
  }
}

# Remove all rows that never had seedlings across a cohort
rm_data <- rm[rm$sumAll != 0, ]
# Remove all rows where FY = 27 (seedlings born in 2024 - can't look at survivorship)
rm_data <- subset(rm_data, FY != 27)
# Check dimensions of cleaned data
dim(rm_data)


BA_data <- read.csv("./data/L0/MNF_ba_radii_allCensuses.csv", row.names = 1)

head(BA_data)
