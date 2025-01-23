
library(readxl)

dat.test <- read_excel("./data/L0/rm_16.35.xlsx")

plot(dat.test$'1998', dat.test$'1999')
