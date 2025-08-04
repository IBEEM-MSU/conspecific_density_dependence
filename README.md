# Repository for IBEEM Conspecific Density Dependence Project

[*For more information on how to use this repository, see: [AboutThisTemplate.md](https://github.com/NEON-biodiversity/template_repository/blob/main/AboutThisTemplate.md)*]

## Table of Contents
- [Introduction](#Introduction)
- [Workflow](#Workflow)
- [Location of data](#Location-of-data)
- [Spatiotemporal extent and resolution](#Spatiotemporal-extent-and-resolution)
- [Usage](#Usage)
- [File naming conventions](#File-naming-conventions)
- [License](#License)
- [Contributors](#Contributors)
- [Contact Information](#Contact-information)

## Introduction

The overarching goal of our project is to investigate how mycorrhizal associations mediate density-dependent interactions between adult trees and seedlings along resource gradients under climate change. 

## Workflow

Here are the steps that we will take to analyze the tree data!

[*Provide a high-level overview of the steps in the workflow (e.g., Download raw data files, merge to monthly time scale, and calculate geodiversity statistics).*] 

## Location of data 

The seedling data files for analyses, for the 12 sites (49.01, 59.12, 8.21, 39.21, 16.35, 60.35, 7.40, 40.40, 21.43, 25.43, 6.45, and 24.45) in the Manistee National Forest (6 excel files per site), are archived in the data/L0 folder.

Some of the relevant covariates are also now in this folder - specifically,
1) "quad_envVars_allSites.csv", which contains information on the light availability, soil nutrients and moisture for each quadrat
2) "quad_climateData_allSites.csv", which contains climatic variables for each quadrat-year combination (from 1998-2023)
3) "MNF_ba_radii_allCensuses.csv", which contains the species-level BA in a 5-m, 10-m, and 15-m radius around each quadrat for each of the nine years that adult trees were censused (1999, 2005, 2007, 2009, 2011, 2013, 2016, 2019, and 2022).

## Spatiotemporal extent and resolution 

[*Describe the spatial and temporal extent and resolution of the data sets resulting from the workflow. For example:*]  
- Spatial extent: 12 sites across the Manistee National Forest in the northernwestern Lower Peninsula, Michigan, USA.
- Spatial resolution: Seedling data are censused in each 1 m2 quadrat along a 200 m x 1 m transect situated in a 240 m x 41 m forest inventory plot.
- Temporal extent: 1998-2024
- Temporal resolution: Annual seedling recruitment and survival data (1 census per year, between June-July).

## Usage

[*Specify which software and versions were used to create the scripts (e.g., R version 4.3.0).*]

### File Naming Conventions

**Seedling Data Files**: Naming for the seedling data files follows the convention "SiteNumber_FirstYear_LastYear.xlsx", with each site having 6 excel files covering subsets of the 26 years of data, with every age of seedling (all suriviving cohorts) censused in a given year included in the excel file with that year. Here's some additional information on the column meanings in these files:
- **uniqueID**: concatenation of Transect ID and Julian (5-digit) date
- **TransectID**: concatenation of Site Name and Quadrat (X) Number
- **Site**: Site Name
- **AgeofSeedling**: Age of seedlings censused, with one row of data per age per quadrat per year
- **Date**: exact date of seedling census
- **Year**: Year of census
- **X**: Quadrat number (quadrats 1-200 identified as 0.5 to 199.5, with labeling by the 0.5-m mark to avoid confusion)
-  **rm through ae**: Abbreviations for seedling species' names; RM = red maple, SM = sugar maple, AB = American beech, WA = white ash,  BC = black cherry, WO = white oak, NR = northern red oak, BO = black oak, BONRO/BNR = black oak/northern red oak, BS = American basswood, wp = white pine, IW = ironwood, ae = American elm (note we do not have sufficient replication to model HEM. WP, or AE)
-  **Notes**: any pertinent notes about the quadrat/censusing that year

**Environmental (non-climatic) Covariates**: "quad_envVars_allSites.csv"
Here's some additional information on the column meanings in this file:
- **TransectID**: concatenation of Site Name and Quadrat (X) Number
- **Site**: Site Name
- **X**: Quadrat number (quadrats 1-200 identified as 0.5 to 199.5, with labeling by the 0.5-m mark to avoid confusion)
- **dateOpen**: Exact date that percent canopy openness (percOpen) was measured
- **percOpen**: percent canopy openness (%/100) - derived from diffuse non-interceptance measured with a hemispherical lens; proportion representing understory light availability (percentage divided by 100) where 1 = 100% open (no shade); measured every 2 m and interpolated to the quadrat-level
- **mean_VWC**: Mean volumetric water content (%), based on an average of 2 measurements taken in early April and late May 1999; measurements were taken every 10 m and interpolated to the quadrat-level; represents a relative site indicator of dryness (only 1 year of measurements)
- **Ca_ppm**: Concentrations of total extractable Calcium, in ppm; measured in 2011 every 10 m and interpolated to the quadrat-level
- **K_ppm**: Concentrations of total extractable Potassium, in ppm; measured in 2011 every 10 m and interpolated to the quadrat-level
- **Mg_ppm**: Concentrations of total extractable Magnesium, in ppm; measured in 2011 every 10 m and interpolated to the quadrat-level
- **NH4_ppm**: Concentrations of ammonium, in ppm, determined with acolorimetric assay using an absorbance microplate reader; measured in 2011 every 10 m and interpolated to the quadrat-level
- **NO3_ppm**: Concentrations of nitrate, in ppm, determined with acolorimetric assay using an absorbance microplate reader; measured in 2011 every 10 m and interpolated to the quadrat-level
- **Base_CEC**: Base cation exchange capacity (milliequivalent, meq, per 100 g soil); integrated index of soil fertility calculates as Bas_CEC = (ppm Ca/200) + (ppm K/390) + (ppm Mg/120), where 200, 390, and 120 represent the meq/100 g soil for base cations of Ca, K, and Mg, respectively. Calculated every 10 m and interpolated to the quadrat-level.

**Climatic Covariates**: "quad_climateData_allSites.csv"
Here's some additional information on the column meanings in this file:
- **TransectID**: concatenation of Site Name and Quadrat Number
- **Year**:Year of climate variables (1998-2023; will add in 2024 if we need it)
- **Site**: Site Name
- **precipTotGS_mm**: Total accumulated precipitation (mm) during the growing season (May to September)
- **sradTotGS_wm2**: Total solar radiation (W/m2) during the growing season 
- **max_tempGS_c**: Maximum daily temperature (C) reached during the growing season 
- **min_tempGS_c**: Minimum daily temperature (C) reached during the growing season 
- **mean_pdsiGS**: Mean monthly growing season Palmer Drought Severity Index (PDSI) (average across May-September monthly averages)
- **min_pdsiGS**: Minimum monthly growing season Palmer Drought Severity Index (PDSI) (minimum of the May-September monthly averages)
	
**Species' Basal Area in 5-m, 10-m, and 15-m radii**: "MNF_ba_radii_allCensuses.csv"
Here's some additional information on the column meanings in this file:
- **TransectID**: concatenation of Site Name and Quadrat (X) Number
- **Site**: Site Name
- **X**: X coordinate denoting the position of quadrat along the transect
- **Y**: Y coordinate, denoting that the quadrats are along transect (always = 0)
- **Midpoint_X**: X coordinate denoting the center of the quadrat
- **Midpoint_Y**: Y coordinate denoting the vertical center point of each quadrat are along transect (always = 0.5)
- **Year**: Year of mature tree census (1999, 2005, 2007, 2009, 2011, 2013, 2016, 2019, and 2022).
- **All other columns**: take the form xx_ym and hold the basal area values (m2 per radius), where xx is a 2-letter species abbreviation, and ym is one of 5m, 10m, or 15m representing the radius for which BA was calculated in meters
ab = American beech, bs = basswood, sm = sugar maple, qa = quaking aspen, ae = American elm, wa = white ash, bta = bigtooth aspen, nro = Northern red oak, iw = ironwood, rm = red maple, bo = black oak, bc = black cherry, bonro = black oak/nro hybrid, wo = white oak, hem = eastern hemlock, wp = white pine, jp = jack pine, svb = serviceberry

**Combined_demography_modeling_data**: This master file was created after cleaning up and manipulating existing raw data files. It includes information for each species' cohort of seedlings in a given year at a given site in a given quadrat and information about the surrounding adult community. Unless otherwise noted, all adult community data was sampled in a 15m radius plot area surrounding each quadrat. Column names refer to the following:
- **SiteID**: ID of each of the 12 sites.
- **quadrat**: Each of the 200 quadrats per site
- **Cohort.Year**: Year each cohort in each quadrat at each site was censused.
- **GLM.Slope.Estimate**: Slope parameter estimate from a simple linear model estimating seedling count reductions over 3 years. More negative values indicate greater rates of attrition through time for a given cohort in a given quadrat at a given site. These values are used as the response variable in the permutation test.
- **Sp**: Species ID of recorded seedlings.
- **Year_0_counts**: Starting number of seedlings in a cohort.
- **Year_1_counts**: Number of seedlings remaining for a given cohort after 1 year.
- **Year_2_counts**: Number of seedlings remaining for a given cohort after 2 years.
- **Year_3_counts**: Number of seedlings remaining for a given cohort after 3 years.
- **Avg_3_year_prior_counts**: Average number of other seedlings a cohort contends with over a 3 year span.
- **Prior_1_year_cohorts_and_1st_year_and_new_2nd_year_counts**: Number of seedlings from the previous year and following two years of new seedlings that emerged that challenge a given cohort of seedlings.
- **Prior_2_year_cohorts_and_new_1st_year_counts**: Number of seedlings from the previous 2 years and new seedlings that emerged that challenge a given cohort of seedlings.
- **Prior_3_year_cohorts_counts**: Number of seedlings from the previous 3 years that challenge an emerging cohort of seedlings.
- **AM_15m**: Total AM adult tree basal area.
- **ECM_15m**: Total ECM adult tree basal area.
- **ab_15m**: Total American Beech adult tree basal area.
- **ae_15m**: Total American Elm adult tree basal area.
- **bc_15m**: Total Black Cherry adult tree basal area.
- **bo_15m**: Total Black Oak adult tree basal area.
- **bonro_15m**: Total Black Oak-Northern Red Oak Hybrid adult tree basal area.
- **bs_15m**: Total Basswood adult tree basal area.
- **bta_15m**: Total Bigtooth Aspen adult tree basal area.
- **hem_15m**: Total Eastern Hemlock adult tree basal area.
- **iw_15m**: Total Ironwood adult tree basal area.
- **jp_15m**: Total Jack Pine adult tree basal area.
- **nro_15m**: Total Northern Red Oak adult tree basal area.
- **qa_15m**: Total Quaking Aspen adult tree basal area.
- **rm_15m**: Total Red Maple adult tree basal area.
- **sm_15m**: Total Sugar Mapel adult tree basal area.
- **svb_15m**: Total Serviceberry adult shrub basal area.
- **wa_15m**: Total White Ash adult tree basal area.
- **wo_15m**: Total White Oak adult tree basal area.
- **wp_15m**: Total White Pine adult tree basal area.
- **Seedling_myc**: Mycorrhizal type of seedling species.
- **conMyc_15m**: Total conmycorrhizal adult tree basal area, based on Sp column species id with conspecifics included.
- **conMycNS_15m**: Total conmycorrhizal adult tree basal area, based on Sp column species id with conspecifics removed.
- **hetMyc_15m**: Total heterospecific adult tree basal area, based on Sp column species id.
- **consp_15m**: Total conspecific adult tree basal area, based on Sp column species id.
- **log0yr**: Log10-transformed counts of a specific cohort of seedlings in a quadrat.
- **logPrior**: Log10-transformed counts of all prior seedlings in a quadrat that are not the cohort of interest.
- **propConsp**: Proportion of total adult tree basal area that are conspecific trees based on Sp column species id.
- **propConmyc**: Proportion of total adult tree basal area that are conmycorrhizal trees based on Sp column species id.
- **propAM**: Proportion of total adult tree basal area that are AM-associated trees.

**Scripts**: [*Use descriptive names for scripts and follow a consistent naming convention, such as `process_data.py` or `cleaning_script.R`. If scripts are designed to be run in sequence, it may be helpful to number them in order (e.g., `1_cleaning_script.R`, `2_process_data.R`)*]

## Scripts

[*Provide descriptions of the scripts in the `/L0` folder, including their purpose, inputs, and outputs. Specify whether scripts are standalone or designed to be run in sequence.*] 

### [*`preprocess_data.py`*]

- **Purpose**: [*Cleans and preprocesses raw Level 0 data.*]
- **Inputs**: [*Raw data files in the `/data` folder (`YYYY-MM-DD_sensorID_raw.ext`).*]
- **Outputs**: [*Processed data files in a new folder (`/processed_data`).*]

### [*`merge_datasets.R`*]

- **Purpose**: [*Merges multiple preprocessed datasets.*]
- **Inputs**: [*Processed data files in the `/processed_data` folder.*]
- **Outputs**: [*Merged dataset saved as `merged_data.csv` in the `/output` folder.*]

## Contributors

[*List individuals or organizations involved in collecting, processing, and maintaining the Level 0 data.*]

## Contact Information

[*Provide contact information for inquiries related to the data and scripts.*]

## Edits 1/23/2025

This is a test for testing things that need tested. Test.
