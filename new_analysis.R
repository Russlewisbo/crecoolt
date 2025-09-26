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
  "retro_or_pros",
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

# Make sure dates are in Date format (assuming they are in "YYYY-MM-DD")
data_extract$date_of_transplant <- as.Date(data_extract$date_of_transplant)
data_extract$discharge_date     <- as.Date(data_extract$discharge_date)

# Create new variable: days from OLT to hospital discharge
data_extract$discharge_d <- as.numeric(
  data_extract$discharge_date - data_extract$date_of_transplant
)

# Quick check
summary(data_extract$discharge_d)
head(data_extract[, c("date_of_transplant", "discharge_date", "discharge_d")])

# Ensure both are in Date format
data_extract$date_of_transplant <- as.Date(data_extract$date_of_transplant)
data_extract$death_date         <- as.Date(data_extract$death_date)

# Calculate days from OLT to death
data_extract$death_d <- as.numeric(
  data_extract$death_date - data_extract$date_of_transplant
)

# Quick check
summary(data_extract$death_d)
head(data_extract[, c("date_of_transplant", "death_date", "death_d")])

# Make sure both are Date format
data_extract$date_of_transplant  <- as.Date(data_extract$date_of_transplant)
data_extract$cre_infection_date  <- as.Date(data_extract$cre_infection_date)

# Calculate days from OLT to CRE infection
data_extract$cre_d <- as.numeric(
  data_extract$cre_infection_date - data_extract$date_of_transplant
)

# Quick check
summary(data_extract$cre_d)
head(data_extract[, c("date_of_transplant", "cre_infection_date", "cre_d")])

# Ensure date variables are in Date format
data_extract$date_of_transplant   <- as.Date(data_extract$date_of_transplant)
data_extract$cre_infection_date   <- as.Date(data_extract$cre_infection_date)
data_extract$discharge_date       <- as.Date(data_extract$discharge_date)
data_extract$death_date           <- as.Date(data_extract$death_date)

# Define 180-day censoring date
data_extract$day180 <- data_extract$date_of_transplant + 180

# Find earliest event among CRE, discharge, death
event_date <- pmin(data_extract$cre_infection_date,
                   data_extract$discharge_date,
                   data_extract$death_date,
                   data_extract$day180,
                   na.rm = TRUE)

# Time from OLT to earliest event (capped at 180 days)
data_extract$time_c <- as.numeric(event_date - data_extract$date_of_transplant)

## preparation of competing risk variable ##

# Outcome coding:
# 0 = censored (no event or after 180d)
# 1 = CRE infection
# 2 = Discharge
# 3 = Death

# 1. Make sure all dates are Date format
data_extract$date_of_transplant   <- as.Date(data_extract$date_of_transplant)
data_extract$cre_infection_date   <- as.Date(data_extract$cre_infection_date)
data_extract$discharge_date       <- as.Date(data_extract$discharge_date)
data_extract$death_date           <- as.Date(data_extract$death_date)

# 2. Define censoring cutoff (180 days after transplant)
data_extract$day180 <- data_extract$date_of_transplant + 180

# 3. Compute earliest observed date among CRE, discharge, death, and 180-day cutoff
event_date <- pmin(data_extract$cre_infection_date,
                   data_extract$discharge_date,
                   data_extract$death_date,
                   data_extract$day180,
                   na.rm = TRUE)

# 4. Time from OLT to earliest event (capped at 180 days)
data_extract$time_c <- as.numeric(event_date - data_extract$date_of_transplant)

# 5. Initialize outcome as censored (0)
data_extract$outcome_comp <- 0

# 6. Assign events based on earliest event type
data_extract$outcome_comp[!is.na(data_extract$cre_infection_date) &
                          data_extract$cre_infection_date <= data_extract$day180 &
                          data_extract$cre_infection_date == event_date] <- 1

data_extract$outcome_comp[!is.na(data_extract$discharge_date) &
                          data_extract$discharge_date <= data_extract$day180 &
                          data_extract$discharge_date == event_date] <- 2

data_extract$outcome_comp[!is.na(data_extract$death_date) &
                          data_extract$death_date <= data_extract$day180 &
                          data_extract$death_date == event_date] <- 3

# 7. Quick sanity check
table(data_extract$outcome_comp, useNA = "ifany")
summary(data_extract$time_c)
head(data_extract[, c("date_of_transplant", "cre_infection_date",
                      "discharge_date", "death_date",
                      "time_c", "outcome_comp")])

write.csv(data_extract, file = here::here("data", "data_extract.csv"), row.names = FALSE)


## select variables for risk model

data_extract$ARF <- ifelse(data_extract$post_olt_compli == 1, 1, 0)


## Baseline Fine-Gray Model



# Outcome coding in your data:
# 0 = censored
# 1 = CRE infection (event of interest)
# 2 = discharge (competing risk)
# 3 = death (competing risk)

# -------------------------------------------------
# 1. Unadjusted Fine–Gray (CRE infection as event)
# -------------------------------------------------
library(prodlim)
library(riskRegression)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)   # for reproducibility

# ----------------------------
# 1. Prepare data
# ----------------------------
df_model <- data.frame(
  time_c = as.numeric(data_extract$time_c),
  outcome_comp = as.integer(data_extract$outcome_comp)
) %>%
  filter(!is.na(time_c) & time_c >= 0)

# ----------------------------
# 2. Base CIF fit
# ----------------------------
ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_model)

# ----------------------------
# 3. Times to evaluate
# ----------------------------
times_seq <- seq(0, 180, by = 5)

# ----------------------------
# 4. Bootstrap CIF estimates
# ----------------------------
B <- 200   # number of bootstrap resamples (increase to 1000 for publication)

boot_cif <- function(data, times, cause) {
  fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = data)
  predict(fit, cause = cause, times = times)
}

boot_results <- lapply(1:3, function(cause) {
  mat <- replicate(B, {
    idx <- sample(nrow(df_model), replace = TRUE)
    boot_data <- df_model[idx, ]
    boot_cif(boot_data, times_seq, cause)
  })
  # mean and CI across bootstrap runs
  est <- predict(ci_fit, cause = cause, times = times_seq)
  lower <- apply(mat, 1, quantile, 0.025, na.rm = TRUE)
  upper <- apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  data.frame(times = times_seq, CIF = est, lower = lower, upper = upper,
             cause = factor(cause))
})

ci_tidy <- bind_rows(boot_results) %>%
  mutate(cause = recode(cause,
                        "1" = "CRE infection",
                        "2" = "Discharge",
                        "3" = "Death"))

# ----------------------------
# 5. NEJM-style plot
# ----------------------------
nejm_colors <- c("CRE infection" = "#D81B60",
                 "Discharge"     = "#1E88E5",
                 "Death"         = "black")

p_nejm <- ggplot(ci_tidy, aes(x = times, y = CIF, color = cause, fill = cause)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(size = 1.2) +
  labs(x = "Days since OLT",
       y = "Cumulative incidence",
       color = NULL, fill = NULL) +
  scale_color_manual(values = nejm_colors) +
  scale_fill_manual(values = nejm_colors) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p_nejm)


#############################################
## Multivariable Fine–Gray + Calibration Plot
## with Jittered Dots + Legend for Outcomes
#############################################

library(dplyr)
library(cmprsk)
library(riskRegression)
library(prodlim)
library(ggplot2)

# ============================================================
# 1. Prepare dataset with covariates
# ============================================================
df_model <- data_extract %>%
  mutate(
    post_olt_compli1 = ifelse(post_olt_compli___1 == 1, 1, 0),
    post_olt_compli3 = ifelse(post_olt_compli___3 == 1, 1, 0),
    post_olt_compli5 = ifelse(post_olt_compli___5 == 1, 1, 0),
    multisite_col    = ifelse(multisite_colonization == 1, 1, 0),
    cre_col1         = ifelse(cre_colonization___1 == 1, 1, 0),
    cre_col2         = ifelse(cre_colonization___2 == 1, 1, 0),
    mec1             = ifelse(mec_of_carbapenem_resi___1 == 1, 1, 0)
  ) %>%
  select(time_c, outcome_comp,
         post_olt_compli1, post_olt_compli3, post_olt_compli5,
         multisite_col, cre_col1, cre_col2, mec1) %>%
  filter(!is.na(time_c) & time_c >= 0)

# ============================================================
# 2. Fit Fine–Gray model (cause = CRE infection)
# ============================================================
fg_risk <- FGR(
  Hist(time_c, outcome_comp) ~ 
    post_olt_compli1 + post_olt_compli3 + post_olt_compli5 +
    multisite_col + cre_col1 + cre_col2 + mec1,
  data  = df_model,
  cause = 1
)

# ============================================================
# 3. Predictions vs Observed CIF
# ============================================================
times_seq <- seq(0, 180, by = 5)

pred_cif <- predict(
  fg_risk,
  newdata = df_model,
  cause = 1,
  times = times_seq,
  se = FALSE,
  keep.newdata = FALSE
)

pred_df <- data.frame(
  times = times_seq,
  CIF_pred = colMeans(pred_cif, na.rm = TRUE)
)

ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_model)
obs_df <- data.frame(
  times = times_seq,
  CIF_obs = predict(ci_fit, cause = 1, times = times_seq)
)

plot_df <- left_join(obs_df, pred_df, by = "times")

# ============================================================
# 4. Dot dataset with event type
# ============================================================
dot_df <- df_model %>%
  mutate(
    dot_type = case_when(
      outcome_comp == 1 ~ "CRE infection",
      outcome_comp == 2 ~ "Discharge",
      outcome_comp == 3 ~ "Death",
      outcome_comp == 0 ~ "Censored"
    ),
    # assign strip positions
    y_pos = case_when(
      outcome_comp == 1 ~ -0.01,   # infections (bottom)
      outcome_comp == 2 ~ 0.26,    # discharge
      outcome_comp == 3 ~ 0.27,    # death
      outcome_comp == 0 ~ 0.28     # censored
    )
  ) %>%
  filter(!is.na(dot_type))

# ============================================================
# 5. NEJM-style plot with jittered dots + legend
# ============================================================
nejm_colors <- c("Observed" = "black",
                 "Predicted" = "#D81B60")

dot_colors <- c("CRE infection" = "red",
                "Discharge"    = "green3",
                "Death"        = "black",
                "Censored"     = "grey50")

p_calib <- ggplot(plot_df, aes(x = times)) +
  # Calibration curves
  geom_line(aes(y = CIF_obs, color = "Observed"), size = 1.2) +
  geom_line(aes(y = CIF_pred, color = "Predicted"),
            size = 1.2, linetype = "dashed") +
  # Jittered dots with legend
  geom_jitter(
    data = dot_df,
    aes(x = time_c, y = y_pos, color = dot_type),
    inherit.aes = FALSE,
    width = 0, height = 0.002,
    shape = 16, alpha = 0.6, size = 1.5
  ) +
  labs(x = "Days since OLT",
       y = "Cumulative incidence of CRE infection",
       title = "Observed vs Predicted CIF with Patient Distributions") +
  scale_color_manual(values = c(nejm_colors, dot_colors)) +
  coord_cartesian(ylim = c(0, max(plot_df$CIF_obs, plot_df$CIF_pred, na.rm = TRUE) * 1.1)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p_calib)



#############################################
## Stratified Calibration Plot (Retrospective vs Prospective)
## with Separated Outcome Dots
#############################################

#############################################
## Stratified calibration with separated outcome dots
#############################################

library(dplyr)
library(riskRegression)
library(prodlim)
library(ggplot2)

## 1) Prep with study design
df_model <- data_extract %>%
  mutate(
    post_olt_compli1 = as.integer(post_olt_compli___1 == 1),
    post_olt_compli3 = as.integer(post_olt_compli___3 == 1),
    post_olt_compli5 = as.integer(post_olt_compli___5 == 1),
    multisite_col    = as.integer(multisite_colonization == 1),
    cre_col1         = as.integer(cre_colonization___1 == 1),
    cre_col2         = as.integer(cre_colonization___2 == 1),
    mec1             = as.integer(mec_of_carbapenem_resi___1 == 1),
    retro_group      = factor(ifelse(retro_or_pros == 1, "Retrospective", "Prospective"))
  ) %>%
  select(time_c, outcome_comp, retro_group,
         post_olt_compli1, post_olt_compli3, post_olt_compli5,
         multisite_col, cre_col1, cre_col2, mec1) %>%
  filter(!is.na(time_c) & time_c >= 0 & !is.na(retro_group))

## 2) Fit FG once on full data
fg_risk <- FGR(
  Hist(time_c, outcome_comp) ~ 
    post_olt_compli1 + post_olt_compli3 + post_olt_compli5 +
    multisite_col + cre_col1 + cre_col2 + mec1,
  data  = df_model,
  cause = 1
)

## 3) Time grid
times_seq <- seq(0, 180, by = 5)

## 4) Helper to get obs + pred per subgroup
get_obs_pred <- function(dat, group_label){
  if (nrow(dat) < 2) return(NULL)
  ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = dat)
  obs_df <- data.frame(
    times   = times_seq,
    CIF_obs = predict(ci_fit, cause = 1, times = times_seq)
  )
  pred_cif <- predict(
    fg_risk, newdata = dat, cause = 1, times = times_seq,
    se = FALSE, keep.newdata = FALSE
  )
  pred_df <- data.frame(
    times    = times_seq,
    CIF_pred = colMeans(pred_cif, na.rm = TRUE)
  )
  dplyr::left_join(obs_df, pred_df, by = "times") %>%
    mutate(group = group_label)
}

## 5) Build calibration data for facets
plot_df <- bind_rows(
  get_obs_pred(filter(df_model, retro_group == "Retrospective"), "Retrospective"),
  get_obs_pred(filter(df_model, retro_group == "Prospective"),  "Prospective")
)

## --- Dynamic strip positions so dots never get clipped ----
y_max <- max(plot_df$CIF_obs, plot_df$CIF_pred, na.rm = TRUE)
y_bottom  <- -0.03 * y_max                  # CRE dots (bottom)
y_top_base <-  1.02 * y_max                 # base just above curves
y_discharge <- y_top_base + 0.01 * y_max
y_death     <- y_top_base + 0.03 * y_max
y_censored  <- y_top_base + 0.05 * y_max
y_min_plot  <- y_bottom - 0.01 * y_max
y_max_plot  <- y_top_base + 0.07 * y_max

## 6) Dot data (with facet variable)
dot_df <- df_model %>%
  mutate(
    group   = retro_group,  # match facet variable
    dot_type = case_when(
      outcome_comp == 1 ~ "CRE infection",
      outcome_comp == 2 ~ "Discharge",
      outcome_comp == 3 ~ "Death",
      outcome_comp == 0 ~ "Censored",
      TRUE ~ NA_character_
    ),
    y_pos = case_when(
      outcome_comp == 1 ~ y_bottom,
      outcome_comp == 2 ~ y_discharge,
      outcome_comp == 3 ~ y_death,
      outcome_comp == 0 ~ y_censored,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(dot_type))

## 7) Plot
nejm_colors <- c("Observed" = "black",
                 "Predicted" = "#D81B60")

dot_colors <- c("CRE infection" = "red",
                "Discharge"    = "green3",
                "Death"        = "black",
                "Censored"     = "grey50")

p_calib_strat <- ggplot(plot_df, aes(x = times)) +
  # curves
  geom_line(aes(y = CIF_obs, color = "Observed"), size = 1.2) +
  geom_line(aes(y = CIF_pred, color = "Predicted"),
            size = 1.2, linetype = "dashed") +
  # dots (with facet-aware data)
  geom_jitter(
    data = dot_df,
    aes(x = time_c, y = y_pos, color = dot_type),
    inherit.aes = FALSE,
    width = 0, height = 0.002,
    shape = 16,
    alpha = 0.7,
    size = ifelse(dot_df$dot_type == "Death", 1.9, 1.4)  # make deaths a touch larger
  ) +
  facet_wrap(~group) +
  labs(x = "Days since OLT",
       y = "Cumulative incidence of CRE infection",
       title = "Observed vs Predicted CIF by Study Design") +
  scale_color_manual(values = c(nejm_colors, dot_colors)) +
  coord_cartesian(ylim = c(y_min_plot, y_max_plot)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text  = element_text(size = 11),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p_calib_strat)
