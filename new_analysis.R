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

# Define the base variable stems you are interested in
vars <- c(
  "multisite_colonization",
  "post_olt_compli",
  "cre_colon",
  "date_positive_post_olt",
  "date_of_transplant",
  "type_of_sample",
  "mec_of_carbapenem_resi",
  "cre_colonization",
  "preemptive",
  "cre_infection",
  "cre_infection_date",
  "infection_source",
  "cre_isolate",
  "bsi_specify",
  "source_control",
  "event_type",
  "score",
  "score_impact",
  "drug",
  "combination_treat",
  "start_date",
  "end_date",
  "indication",
  "death",
  "death_date",
  "discharge_date",
  "cre_clearence",   # spelling from your file
  "clearence_date",  # spelling from your file
  "cre_recurrence",
  "relapse_date"
)

# Build regex pattern to catch both base variables and checkbox expansions
pattern <- paste0("^(", paste(vars, collapse = "|"), ")")

# Get all matching variables
matched_vars <- names(data_new)[grepl(pattern, names(data_new))]

# Subset into a new dataframe
data_extract <- data_new[, matched_vars, drop = FALSE]

# Check what you got
names(data_extract)
dim(data_extract)
saveRDS(data_extract, "data_extract.rds")
