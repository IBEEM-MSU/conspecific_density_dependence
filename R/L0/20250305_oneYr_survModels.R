
## Read in packages
library(tidyverse) # For data manipulation & graphing
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823.
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.
library(dplyr)
library(fuzzyjoin) # for merging data

# Read in data for rm
rm <- read.csv("./data/L0/rm_allCohortsSites.csv", row.names = 1)

# This replaces the NA parings in the matrix with 0s (to account for years a cohort did not exist)
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

## Calculate conmycorrhizal, heterospecific, and heteromycorrhizal BA for RM
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH rm
BA_data$conMyc_5m = rowSums(BA_data[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
BA_data$conMyc_10m = rowSums(BA_data[, c("sm_10m", "ae_10m", "wa_10m", 
                                        "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
BA_data$conMyc_15m = rowSums(BA_data[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - NO rm
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
rm_surv <- rm_merged[,c(1:7,9,39:58)]

# Make Site a character 
rm_surv$Site = as.character(rm_surv$Site)

# Convert NAs in the basal area data to 0's
rm_surv[,14:28][is.na(rm_surv[,14:28])] <- 0

########################################
##### Candidate Models - 5-m radius
## Effect of conspecific DD in 5-m radius
cs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ rm_5m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cs_5m)
## Effect of con-mycorrhizal DD in 5-m radius
cm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMyc_5m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cm_5m)
## Effect of con-mycorrhizal DD in 5-m radius - excluding RM
cmNs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMycNS_5m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cmNs_5m)
## Effect of hetero-mycorrhizal DD in 5-m radius 
hm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetMyc_5m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hm_5m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 5-m radius - excluding RM
hmcmNS_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                     hetMyc_5m * conMycNS_5m  + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hmcmNS_5m)
## Effect of hetero-mycorrhizal * con-mycorrhizal DD in 5-m radius 
hmcm_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~
                     hetMyc_5m*conMyc_5m  + (1 | Site), 
                   family = binomial(link = "logit"), 
                   data = rm_surv)
summary(hmcm_5m)
## Effect of hetero-specific DD in 5-m radius 
hs_5m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetSp_5m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hs_5m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_df_5m <- tibble(Model = c("cs_5m", "cm_5m", "cmNs_5m", "hm_5m", 
                           "hmcmNS_5m", "hmcm_5m", "hs_5m"),
  AIC = sapply(list(cs_5m, cm_5m, cmNs_5m, hm_5m, 
                    hmcmNS_5m, hmcm_5m, hs_5m), AIC)) %>%
  arrange(AIC)  # Sort by AIC

# Print the sorted AIC values
print(aic_df_5m)

########################################
##### Candidate Models - 10-m radius
## Effect of conspecific DD in 10-m radius
cs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ rm_10m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cs_10m)
## Effect of con-mycorrhizal DD in 10-m radius
cm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMyc_10m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(cm_10m)
## Effect of con-mycorrhizal DD in 10-m radius - excluding RM
cmNs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMycNS_10m + (1 | Site), 
                 family = binomial(link = "logit"), 
                 data = rm_surv)
summary(cmNs_10m)
## Interactive effect between conspecific and mycorrhizal dominance
csHmInt_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_10m*hetMyc_10m + (1 | Site), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(csHmInt_10m)
## Additive effect between conspecific and mycorrhizal dominance
csHmAdd_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                      rm_10m + hetMyc_10m + (1 | Site), 
                    family = binomial(link = "logit"), 
                    data = rm_surv)
summary(csHmAdd_10m)
## Effect of hetero-mycorrhizal DD in 10-m radius 
hm_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetMyc_10m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hm_10m)
## Effect of hetero-specific DD in 10-m radius 
hs_10m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetSp_10m + (1 | Site), 
               family = binomial(link = "logit"), 
               data = rm_surv)
summary(hs_10m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_df_10m <- tibble(Model = c("cs_10m", "cm_10m", "cmNs_10m", "csHmInt_10m", 
                           "csHmAdd_10m", "hm_10m", "hs_10m"),
                 AIC = sapply(list(cs_10m, cm_10m, cmNs_10m, csHmInt_10m, 
                                   csHmAdd_10m, hm_10m, hs_10m), AIC)) %>%
  arrange(AIC)  # Sort by AIC

# Print the sorted AIC values
print(aic_df_10m)

########################################
##### Candidate Models - 15-m radius
## Effect of conspecific DD in 15-m radius
cs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ rm_15m + (1 | Site), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cs_15m)
## Effect of con-mycorrhizal DD in 15-m radius
cm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMyc_15m + (1 | Site), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(cm_15m)
## Effect of con-mycorrhizal DD in 15-m radius - excluding RM
cmNs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ conMycNS_15m + (1 | Site), 
                  family = binomial(link = "logit"), 
                  data = rm_surv)
summary(cmNs_15m)
## Interactive effect between conspecific and mycorrhizal dominance
csHmInt_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       rm_15m*hetMyc_15m + (1 | Site), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(csHmInt_15m)
## Additive effect between conspecific and mycorrhizal dominance
csHmAdd_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ 
                       rm_15m + hetMyc_15m + (1 | Site), 
                     family = binomial(link = "logit"), 
                     data = rm_surv)
summary(csHmAdd_15m)
## Effect of hetero-mycorrhizal DD in 15-m radius 
hm_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetMyc_15m + (1 | Site), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(hm_15m)
## Effect of hetero-specific DD in 15-m radius 
hs_15m <- glmer(cbind(total_count_yr1, total_count_yr0 - total_count_yr1) ~ hetSp_15m + (1 | Site), 
                family = binomial(link = "logit"), 
                data = rm_surv)
summary(hs_15m)

## Model selection/evaluation using AIC
# Create a data frame with model names and AIC values
aic_df_15m <- tibble(Model = c("cs_15m", "cm_15m", "cmNs_15m", "csHmInt_15m", 
                               "csHmAdd_15m", "hm_15m", "hs_15m"),
                     AIC = sapply(list(cs_15m, cm_15m, cmNs_15m, csHmInt_15m, 
                                       csHmAdd_15m, hm_15m, hs_15m), AIC)) %>%
  arrange(AIC)  # Sort by AIC

# Print the sorted AIC values
print(aic_df_15m)

