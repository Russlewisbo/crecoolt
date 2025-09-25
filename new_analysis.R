###############################
#### CRECOLT_2025 Data analysis 
###############################

library (dplyr)
library (janitor)
library(summarytools)

# Read the CSV file
data <- read.csv("CRECOOLT_overall.csv", stringsAsFactors = FALSE)

# Filter rows where retro_or_pros is not missing or empty
data_new <- subset(data, !is.na(retro_or_pros) & retro_or_pros != "")

# Check results
dim(data_new)
head(data_new)

