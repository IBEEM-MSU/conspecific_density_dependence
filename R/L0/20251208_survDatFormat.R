
#################################################
####### Adding in corrected proportional and heterospecific BA variables
####### Author: Bailey McNichol
####### Date: 8 December 2025
#################################################

## Read in packages
library(tidyverse) # For data manipulation & graphing
library(dplyr)
library(fuzzyjoin) # for merging data

#################################################
#################################################
# Read in current data structure
allSpDat <- read.csv("./data/L3/Combined_demography_modeling_data.csv")

## Add survYear
allSpDat$survYear = allSpDat$Cohort.Year + 1

## Create two new columns, hetAM_15m and hetECM_15m, that exclude conspecific BA from calculations of AM and ECM associating basal area

# Define species as AM or ECM associating
AM_species  <- c("bc", "rm", "sm", "wa")
ECM_species <- c("ab", "bonro", "bs", "iw", "nro", "wo")

# Create placeholder columns
allSpDat$hetAM_15m  <- NA
allSpDat$hetECM_15m <- NA

# Match and substract
allSpDat <- allSpDat %>%
  rowwise() %>%
  mutate(
    # Identify the column containing the focal species' basal area
    focal_col = paste0(Sp, "_15m"),
    focal_value = get(focal_col),
    # Calculate heterospecific AM and ECM basal area
    hetAM_15m = ifelse(Sp %in% AM_species,
                       AM_15m - focal_value,   # subtract conspecific
                       AM_15m),  # all AM are heterospecific for ECM species
    hetECM_15m = ifelse(Sp %in% ECM_species,
                        ECM_15m - focal_value,  # subtract conspecific
                        ECM_15m)) %>%    # all ECM are heterospecific for AM species
  ungroup() %>%
  select(-focal_col, -focal_value)   # cleanup temporary columns

# make sure data.frame
allSpDat <- as.data.frame(allSpDat)

# define species BA columns (exact names in your DF)
species_cols <- c("ab_15m","ae_15m","bc_15m","bo_15m","bonro_15m","bs_15m",
                  "bta_15m","hem_15m","iw_15m","jp_15m","nro_15m","qa_15m",
                  "rm_15m","sm_15m","svb_15m","wa_15m","wo_15m","wp_15m")

# get numeric matrix of BA and row sums (total BA)
ba_mat <- as.matrix(allSpDat[, species_cols])
# If these columns might be factors/characters, coerce to numeric safely:
ba_mat <- apply(ba_mat, 2, function(x) as.numeric(as.character(x)))
total_ba <- rowSums(ba_mat, na.rm = TRUE)       # use na.rm just in case

# handle rows with total_ba == 0 by setting totals to NA (avoids dividing by zero)
total_ba_zero <- total_ba == 0
total_ba[total_ba_zero] <- NA_real_

# compute proportions (each row sums to 1 unless total_ba was 0)
ba_props_mat <- sweep(ba_mat, 1, total_ba, "/")  # NA where total_ba was 0

# conspecific BA proportion:
# construct expected conspecific column names from Sp (e.g., 'ab' -> 'ab_15m')
cons_col_names <- paste0(allSpDat$Sp, "_15m")
cons_idx <- match(cons_col_names, species_cols)   # NA if species code not in list
# extract conspecific BA (numeric), 0 when species column missing
n <- nrow(ba_mat)
cons_ba <- rep(0, n)
valid_rows <- !is.na(cons_idx)
cons_ba[valid_rows] <- ba_mat[cbind(which(valid_rows), cons_idx[valid_rows])]
cons_ba_prop <- cons_ba / rowSums(ba_mat, na.rm = TRUE)   # will be NA when total_ba was 0

# define AM species columns (based on your mapping)
AM_species_cols <- c("sm_15m","ae_15m","wa_15m","rm_15m","bc_15m","svb_15m")
AM_idx <- match(AM_species_cols, species_cols)

# AM total proportion (AM BA / total BA)
AM_total_ba <- rowSums(ba_mat[, AM_idx, drop = FALSE], na.rm = TRUE)
AM_total_prop <- AM_total_ba / rowSums(ba_mat, na.rm = TRUE)   # NA when total_ba == 0

# conspecific BA that is in AM pool (if conspecific is AM; else 0)
cons_in_AM <- ifelse(cons_col_names %in% AM_species_cols, cons_ba, 0)
cons_in_AM_prop <- cons_in_AM / rowSums(ba_mat, na.rm = TRUE)

# numerator and denominator for heterospecific AM proportion
denom_het <- 1 - cons_ba_prop  # heterospecific BA proportion
AM_het_numer <- AM_total_prop - cons_in_AM_prop  # AM BA excluding conspecific (if any)

# final: proportion of heterospecific BA that is AM
propHetAM <- ifelse(is.na(denom_het) | denom_het == 0, NA_real_, AM_het_numer / denom_het)

# add proportiron of heterospecific AM BA to allSpDat
allSpDat$propHetAM   <- propHetAM




## Trim the unnecessary columns
allSpDat_Final <- allSpDat[,c(1:34,38:42,44:47)]

## Write the final data file for analyses
write.csv(allSpDat_Final, "./data/L3/final_seedling_demogr_modeling_data.csv")



