
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

# species definitions 
AM_species  <- c("bc", "rm", "sm", "wa")
ECM_species <- c("ab", "bonro", "bs", "iw", "nro", "wo")

# BA column names in DF (exact)
species_cols <- c("ab_15m","ae_15m","bc_15m","bo_15m","bonro_15m","bs_15m",
                  "bta_15m","hem_15m","iw_15m","jp_15m","nro_15m","qa_15m",
                  "rm_15m","sm_15m","svb_15m","wa_15m","wo_15m","wp_15m")

# Coerce BA columns to numeric (safe)
allSpDat[species_cols] <- lapply(allSpDat[species_cols],
                                 function(x) as.numeric(as.character(x)))

# --- 0/0 SAFE DIVISION FUNCTION ---
safe_div <- function(num, den) {
  out <- num / den
  out[den == 0] <- 0     # 0/0 → 0
  out
}

# Compute raw BA matrix and totals
ba_mat <- as.matrix(allSpDat[, species_cols])
ba_mat_zeroNA <- ba_mat
ba_mat_zeroNA[is.na(ba_mat_zeroNA)] <- 0

total_ba <- rowSums(ba_mat_zeroNA)
rows_all_NA <- apply(ba_mat, 1, function(r) all(is.na(r)))

# Build conspecific column names and find indices
cons_col_names <- paste0(allSpDat$Sp, "_15m")
cons_idx <- match(cons_col_names, species_cols)

# Extract conspecific BA safely
n <- nrow(allSpDat)
cons_ba <- numeric(n)
cons_missing_flag <- logical(n)

valid_rows <- which(!is.na(cons_idx))
if(length(valid_rows)>0) {
  cons_ba_vals <- ba_mat[cbind(valid_rows, cons_idx[valid_rows])]
  cons_ba[valid_rows] <- cons_ba_vals
  cons_missing_flag[valid_rows] <- is.na(cons_ba_vals)
}
cons_ba[is.na(cons_ba)] <- 0

# Compute guild totals (AM and ECM)
AM_species_cols  <- paste0(AM_species, "_15m")
ECM_species_cols <- paste0(ECM_species, "_15m")

AM_idx  <- match(AM_species_cols, species_cols)
ECM_idx <- match(ECM_species_cols, species_cols)

AM_ba <- if(all(is.na(AM_idx))) rep(0, n) else rowSums(ba_mat_zeroNA[, AM_idx, drop = FALSE], na.rm = TRUE)
ECM_ba <- if(all(is.na(ECM_idx))) rep(0, n) else rowSums(ba_mat_zeroNA[, ECM_idx, drop = FALSE], na.rm = TRUE)

# Raw heterospecific BA totals
focal_value <- cons_ba
focal_value_na <- cons_missing_flag

focal_value_for_calc <- focal_value
focal_value_for_calc[focal_value_na] <- 0

hetAM_raw  <- ifelse(allSpDat$Sp %in% AM_species, AM_ba - focal_value_for_calc, AM_ba)
hetECM_raw <- ifelse(allSpDat$Sp %in% ECM_species, ECM_ba - focal_value_for_calc, ECM_ba)

# Prevent small negative values
hetAM_raw  <- pmax(hetAM_raw, 0)
hetECM_raw <- pmax(hetECM_raw, 0)

allSpDat$hetAM_15m  <- hetAM_raw
allSpDat$hetECM_15m <- hetECM_raw

# Total BA for proportions
total_ba_for_prop <- total_ba
total_ba_for_prop[total_ba_for_prop == 0] <- NA_real_

# --- PROPORTIONS WITH 0/0 PROTECTION ---

# Conspecific BA proportion
cons_ba_prop <- safe_div(cons_ba, total_ba_for_prop)

# AM total proportion
AM_total_prop <- safe_div(AM_ba, total_ba_for_prop)

# Conspecific BA in AM guild
cons_in_AM <- ifelse(allSpDat$Sp %in% AM_species, cons_ba, 0)
cons_in_AM_prop <- safe_div(cons_in_AM, total_ba_for_prop)

# Heterospecific BA proportion denominator
denom_het <- 1 - cons_ba_prop   # already safe because cons_ba_prop handled

# AM heterospecific numerator
AM_het_numer <- AM_total_prop - cons_in_AM_prop

# Final proportion of heterospecific BA that is AM
propHetAM <- safe_div(AM_het_numer, denom_het)

allSpDat$propHetAM <- propHetAM

## Trim the unnecessary columns
allSpDat_Final <- allSpDat[,c(1:34,38:42,44:47)]

## Write the final data file for analyses
write.csv(allSpDat_Final, "./data/L3/final_seedling_demogr_modeling_data.csv")


