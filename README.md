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

Trees hate trees. A lot.

## Workflow

Here are the steps that we will take to analyze the tree data!

[*Provide a high-level overview of the steps in the workflow (e.g., Download raw data files, merge to monthly time scale, and calculate geodiversity statistics).*] 

## Location of data 

The preliminary seedling data files for analyses, for the 12 sites (49.01, 59.12, 8.21, 39.21, 16.35, 60.35, 7.40, 40.40, 21.43, 25.43, 6.45, and 24.45) in the Manistee National Forest (6 excel files per site) are archived in the L0 folder.

## Spatiotemporal extent and resolution 

[*Describe the spatial and temporal extent and resolution of the data sets resulting from the workflow. For example:*]  
- Spatial extent: 12 sites across the Manistee National Forest in the northernwestern Lower Peninsula, Michigan, USA.
- Spatial resolution: Seedling data are censused in each 1 m2 quadrat along a 200 m x 1 m transect situated in a 240 m x 41 m forest inventory plot.
- Temporal extent: 1998-2023
- Temporal resolution: Annual seedling recruitment and survival data (1 census per year, between June-July).

## Usage

[*Specify which software and versions were used to create the scripts (e.g., R version 4.3.0).*]

### File Naming Conventions

**Data Files**: Naming for the seedling data files follows the convention "SiteNumber_FirstYear_LastYear.xlsx", with each site having 6 excel files covering subsets of the 26 years of data, with every age of seedling (all suriviving cohorts) censused in a given year included in the excel file with that year. Here's some additional information on the column meanings in these files:
- **uniqueID**: concatenation of Transect ID and Julian (5-digit) date
- **TransectID**: concatenation of Site Name and Quadrat (X) Number
- **Site**: Site Name
- **AgeofSeedling**: Age of seedlings censused, with one row of data per age per quadrat per year
- **Date**: exact date of seedling census
- **Year**: Year of census
- **X**: Quadrat number (quadrats 1-200 identified as 0.5 to 199.5, with labeling by the 0.5-m mark to avoid confusion)
-  **rm through ae**: Abbreviations for seedling species' names; RM = red maple, SM = sugar maple, AB = American beech, WA = white ash,  BC = black cherry, WO = white oak, NR = northern red oak, BO = black oak, BONRO/BNR = black oak/northern red oak, BS = American basswood, HEM = hemlock, wp = white pine, IW = ironwood, ae = American elm (note we do not have sufficient replication to model HEM. WP, or AE)
-  **Notes**: any pertinent notes about the quadrat/censusing that year
- **Scripts**: [*Use descriptive names for scripts and follow a consistent naming convention, such as `process_data.py` or `cleaning_script.R`. If scripts are designed to be run in sequence, it may be helpful to number them in order (e.g., `1_cleaning_script.R`, `2_process_data.R`)*]

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
