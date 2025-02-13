

library(dplyr)
library(readxl)

# Find all the excel files in the L0 data folder 
temp <- list.files("./data/L0", pattern = ".xlsx", full.names=T) 
# Read them into a list (Use indexing if you don't want all of the files!)
t <- lapply(temp, readxl::read_xlsx)
# Bind them all together 
u <- do.call(dplyr::bind_rows, t)

# Take a look 
head(u)
tail(u)