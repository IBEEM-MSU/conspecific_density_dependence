

library(dplyr)
library(readxl)


# Find all the excel files in the L0 data folder 
temp <- list.files("./data/L0/Site_Year_Data", pattern = ".xlsx", full.names=T) 
# Read them into a list (Use indexing if you don't want all of the files!)
t <- lapply(temp, readxl::read_excel, col_types = c("text", "text", "text",
                                                    "text",  "text","text",
                                                    "text","numeric","numeric",
                                                    "numeric","numeric","numeric",
                                                    "numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","text"))
# Bind them all together 
u <- do.call(dplyr::bind_rows, t)

# Take a look 
head(u)
tail(u)

u = as.data.frame(u)

# Write csv
write.csv(u, "./data/L0/all_Seedlings_SpSiteYear.csv")


# Find all the excel files in the L0 data folder 
temp <- list.files("./data/L0/Species_Site_Data", pattern = "^wp.*\\.xlsx$", full.names = TRUE)

# Read them into a list (Use indexing if you don't want all of the files!)
t <- lapply(temp, readxl::read_excel, col_types = c("text", "text", "numeric","numeric",
                                                    "numeric","numeric","text",
                                                    "text","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric",
                                                    "numeric","numeric","numeric","numeric"))
# Bind them all together 
u <- do.call(dplyr::bind_rows, t)

str(u)
u$Site[u$Site == "8.210000000000001"] <- "8.21"
u$Site = factor(u$Site, levels = c("49.01","59.12","8.21","39.21",
                                   "16.35","60.35","7.40","40.40",
                                   "21.43","25.43","6.45","24.45"))


u = as.data.frame(u)

# Write csv
write.csv(u, "./data/L0/wp_allCohortsSites.csv")

