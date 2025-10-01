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

#############################################
## Bar plot of covariate distribution
## Retrospective vs Prospective cohorts
## with custom labels
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Select the covariates of interest
# ============================================================
covariates <- c("post_olt_compli1", "post_olt_compli3", "post_olt_compli5",
                "multisite_col", "cre_col1", "cre_col2", "mec1")

df_cov <- df_model %>%
  select(retro_group, all_of(covariates))

# ============================================================
# 2. Reshape to long format
# ============================================================
df_long <- df_cov %>%
  pivot_longer(
    cols = all_of(covariates),
    names_to = "covariate",
    values_to = "value"
  )

# ============================================================
# 3. Summarize % positive per group
# ============================================================
df_summary <- df_long %>%
  group_by(retro_group, covariate) %>%
  summarise(
    n = n(),
    n_pos = sum(value == 1, na.rm = TRUE),
    pct_pos = 100 * n_pos / n,
    .groups = "drop"
  ) %>%
  mutate(
    covariate_label = recode(covariate,
                             "cre_col1"          = "Prior",
                             "cre_col2"          = "Post",
                             "mec1"              = "KPC",
                             "multisite_col"     = "Multisite",
                             "post_olt_compli1"  = "ARF",
                             "post_olt_compli3"  = "Vent",
                             "post_olt_compli5"  = "Reint"
    )
  )

# ============================================================
# 4. Bar plot (NEJM-style)
# ============================================================
p_cov <- ggplot(df_summary, aes(x = covariate_label, y = pct_pos, fill = retro_group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Covariate", y = "Percent positive (%)",
       title = "Distribution of Covariates by Study Design") +
  scale_fill_manual(values = c("Retrospective" = "#1E88E5", "Prospective" = "#D81B60")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

print(p_cov)

#############################################
## Bar plot of covariate distribution
## Previous study (CRE3 dataset)
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Load and clean CRE3 dataset
# ============================================================
CRE3 <- read.csv("CRE3.csv", stringsAsFactors = FALSE)

CRE3_clean <- CRE3 %>%
  mutate(
    time        = as.numeric(time),
    status      = as.numeric(status),
    reint       = as.numeric(reint),
    mv          = as.numeric(mv),
    arf         = as.numeric(arf),
    crepre_60   = as.numeric(crepre_60),
    crepost_60  = as.numeric(crepost_60),
    multipost   = as.numeric(multipost)
  ) %>%
  filter(!is.na(time), !is.na(status), time > 0)

# ============================================================
# 2. Select predictors and reshape to long format
# ============================================================
predictors <- c("reint", "mv", "arf", "crepre_60", "crepost_60", "multipost")

df_long <- CRE3_clean %>%
  select(all_of(predictors)) %>%
  pivot_longer(
    cols = all_of(predictors),
    names_to = "covariate",
    values_to = "value"
  )

# ============================================================
# 3. Summarize % positive for each covariate
# ============================================================
df_summary <- df_long %>%
  group_by(covariate) %>%
  summarise(
    n = n(),
    n_pos = sum(value == 1, na.rm = TRUE),
    pct_pos = 100 * n_pos / n,
    .groups = "drop"
  ) %>%
  mutate(
    covariate_label = recode(covariate,
                             "reint"       = "Reint",
                             "mv"          = "Vent",
                             "arf"         = "ARF",
                             "crepre_60"   = "Coln prior",
                             "crepost_60"  = "Coln post",
                             "multipost"   = "Multisite"
    )
  )

# ============================================================
# 4. Bar plot (NEJM-style)
# ============================================================
p_cov_CRE3 <- ggplot(df_summary, aes(x = covariate_label, y = pct_pos)) +
  geom_col(fill = "#1E88E5", width = 0.7) +
  labs(x = "Covariate", y = "Percent positive (%)",
       title = "Covariate Distribution (Previous Study, CRE3)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p_cov_CRE3)

#############################################
## Risk Stratification and Survival Analysis
## Based on Predicted CRE Infection Risk
#############################################

library(dplyr)
library(cmprsk)
library(riskRegression)
library(prodlim)
library(survival)
library(ggplot2)
library(tidyr)

# ============================================================
# 1. Prepare dataset with covariates (matching your analysis)
# ============================================================
df_model_full <- data_extract %>%
  mutate(
    post_olt_compli1 = ifelse(post_olt_compli___1 == 1, 1, 0),
    post_olt_compli3 = ifelse(post_olt_compli___3 == 1, 1, 0),
    post_olt_compli5 = ifelse(post_olt_compli___5 == 1, 1, 0),
    multisite_col    = ifelse(multisite_colonization == 1, 1, 0),
    cre_col1         = ifelse(cre_colonization___1 == 1, 1, 0),
    cre_col2         = ifelse(cre_colonization___2 == 1, 1, 0),
    mec1             = ifelse(mec_of_carbapenem_resi___1 == 1, 1, 0),
    retro_group      = factor(ifelse(retro_or_pros == 1, "Retrospective", "Prospective"))
  ) %>%
  select(time_c, outcome_comp, retro_group,
         post_olt_compli1, post_olt_compli3, post_olt_compli5,
         multisite_col, cre_col1, cre_col2, mec1)

# Create complete cases dataset for modeling
df_model <- df_model_full %>%
  filter(!is.na(time_c) & time_c >= 0) %>%
  filter(complete.cases(.))  # Only keep rows with no missing values

cat("Total rows in original data:", nrow(df_model_full), "\n")
cat("Rows with complete data for modeling:", nrow(df_model), "\n")

# ============================================================
# 2. Refit Fine-Gray model on complete cases
# ============================================================
fg_risk <- FGR(
  Hist(time_c, outcome_comp) ~ 
    post_olt_compli1 + post_olt_compli3 + post_olt_compli5 +
    multisite_col + cre_col1 + cre_col2 + mec1,
  data  = df_model,
  cause = 1
)

# ============================================================
# 3. Get predicted risks at specific time
# ============================================================
# Choose evaluation time for risk prediction
eval_time <- 30  # Can adjust (e.g., 14, 30, 60 days)

# Get predicted CRE infection risk for patients with complete data
pred_risk_30d <- predict(
  fg_risk,
  newdata = df_model,
  cause = 1,
  times = eval_time,
  se = FALSE,
  keep.newdata = FALSE
)

# Check dimensions
cat("\nPrediction dimensions:", dim(pred_risk_30d), "\n")
cat("Data dimensions:", nrow(df_model), "\n")

# Add predicted risk to dataframe
if(is.matrix(pred_risk_30d)) {
  df_model$pred_cre_risk <- as.numeric(pred_risk_30d[,1])
} else {
  df_model$pred_cre_risk <- as.numeric(pred_risk_30d)
}

# ============================================================
# 4. Define risk groups based on tertiles
# ============================================================
# Method 1: Tertiles (3 equal-sized groups)
risk_tertiles <- quantile(df_model$pred_cre_risk, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
df_model$risk_group_tertile <- cut(
  df_model$pred_cre_risk,
  breaks = risk_tertiles,
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

# Method 2: Clinical cutpoints (you can adjust these)
# First check the range of predicted risks
cat("\n=== Predicted Risk Range ===\n")
summary(df_model$pred_cre_risk)

# Adjust clinical cutpoints based on actual risk distribution
risk_quantiles <- quantile(df_model$pred_cre_risk, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
cat("\nRisk quartiles:", round(risk_quantiles, 3), "\n")

# Use data-driven cutpoints if clinical ones don't work
if(max(df_model$pred_cre_risk, na.rm = TRUE) < 0.15) {
  clinical_cuts <- c(0, 
                     quantile(df_model$pred_cre_risk, 0.33, na.rm = TRUE),
                     quantile(df_model$pred_cre_risk, 0.67, na.rm = TRUE),
                     1)
  df_model$risk_group_clinical <- cut(
    df_model$pred_cre_risk,
    breaks = clinical_cuts,
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )
} else {
  clinical_cuts <- c(0, 0.05, 0.15, 1)
  df_model$risk_group_clinical <- cut(
    df_model$pred_cre_risk,
    breaks = clinical_cuts,
    labels = c("Low (<5%)", "Medium (5-15%)", "High (>15%)"),
    include.lowest = TRUE
  )
}

# Check distribution
cat("\n=== Risk Group Distribution (Tertiles) ===\n")
table(df_model$risk_group_tertile)
cat("\nRisk ranges by tertile:\n")
df_model %>%
  group_by(risk_group_tertile) %>%
  summarise(
    n = n(),
    min_risk = min(pred_cre_risk, na.rm = TRUE),
    median_risk = median(pred_cre_risk, na.rm = TRUE),
    max_risk = max(pred_cre_risk, na.rm = TRUE)
  ) %>%
  print()

# ============================================================
# 5. Calculate 180-day survival by risk group
# ============================================================
# Create survival outcome (death = 1, alive = 0)
df_model$death_180d <- ifelse(
  df_model$outcome_comp == 3 & df_model$time_c <= 180, 
  1, 
  0
)

# Summary statistics by risk group
survival_summary <- df_model %>%
  group_by(risk_group_tertile) %>%
  summarise(
    n = n(),
    n_deaths = sum(death_180d),
    death_rate = n_deaths / n,
    survival_rate = 1 - death_rate,
    n_cre = sum(outcome_comp == 1),
    cre_rate = n_cre / n,
    n_discharge = sum(outcome_comp == 2),
    discharge_rate = n_discharge / n
  ) %>%
  mutate(across(c(death_rate, survival_rate, cre_rate, discharge_rate), 
                ~ round(. * 100, 1)))

cat("\n=== 180-Day Outcomes by Risk Group (Tertiles) ===\n")
print(survival_summary)

# ============================================================
# 6. Kaplan-Meier Survival Curves by Risk Group
# ============================================================
# Prepare data for survival analysis
df_surv <- df_model %>%
  filter(!is.na(risk_group_tertile)) %>%
  mutate(
    # For survival analysis: event = death, censor others
    surv_event = ifelse(outcome_comp == 3, 1, 0),
    surv_time = time_c
  )

# Fit KM curves
km_fit <- survfit(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_surv)

# Print survival at key timepoints
cat("\n=== Survival Estimates at Key Timepoints ===\n")
summary(km_fit, times = c(30, 60, 90, 180))

# Create tidy data for plotting
km_tidy <- broom::tidy(km_fit) %>%
  mutate(risk_group = gsub("risk_group_tertile=", "", strata))

# NEJM-style survival plot
p_survival <- ggplot(km_tidy, aes(x = time, y = estimate, color = risk_group)) +
  geom_step(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = risk_group), 
              alpha = 0.15, linetype = 0) +
  labs(
    x = "Days since OLT",
    y = "Survival probability",
    title = "180-Day Survival by Predicted CRE Infection Risk",
    color = "Risk Group",
    fill = "Risk Group"
  ) +
  scale_color_manual(values = c("Low" = "green3", "Medium" = "#FFA500", "High" = "red")) +
  scale_fill_manual(values = c("Low" = "green3", "Medium" = "#FFA500", "High" = "red")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 180)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

print(p_survival)

# ============================================================
# 7. Log-rank test for survival differences
# ============================================================
logrank_test <- survdiff(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_surv)
cat("\n=== Log-rank Test for Survival Differences ===\n")
print(logrank_test)

# ============================================================
# 8. Cox Proportional Hazards Model
# ============================================================
# Compare hazard ratios between risk groups
cox_model <- coxph(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_surv)
cat("\n=== Cox Proportional Hazards Model ===\n")
print(summary(cox_model))

# Extract hazard ratios with confidence intervals
hr_table <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(., 2)),
         p.value = round(p.value, 4))

cat("\n=== Hazard Ratios (Reference: Low Risk) ===\n")
print(hr_table)

# ============================================================
# 9. Competing Risks Analysis: Death as Outcome
# ============================================================
# Calculate cumulative incidence of death by risk group
times_eval <- seq(0, 180, by = 10)

ci_death_by_risk <- lapply(levels(df_model$risk_group_tertile), function(rg) {
  df_subset <- filter(df_model, risk_group_tertile == rg)
  if(nrow(df_subset) > 0) {
    ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_subset)
    data.frame(
      times = times_eval,
      CI_death = predict(ci_fit, cause = 3, times = times_eval),
      risk_group = rg
    )
  }
})

ci_death_df <- bind_rows(ci_death_by_risk)

# Plot cumulative incidence of death
p_ci_death <- ggplot(ci_death_df, aes(x = times, y = CI_death, color = risk_group)) +
  geom_line(size = 1.2) +
  labs(
    x = "Days since OLT",
    y = "Cumulative incidence of death",
    title = "Cumulative Incidence of Death by CRE Risk Group",
    color = "Risk Group"
  ) +
  scale_color_manual(values = c("Low" = "green3", "Medium" = "#FFA500", "High" = "red")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

print(p_ci_death)

# ============================================================
# 10. Detailed Table of Outcomes at Specific Timepoints
# ============================================================
timepoints <- c(30, 60, 90, 180)

detailed_outcomes <- lapply(levels(df_model$risk_group_tertile), function(rg) {
  df_subset <- filter(df_model, risk_group_tertile == rg)
  if(nrow(df_subset) > 0) {
    ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_subset)
    
    # Get CI for each outcome at each timepoint
    results <- lapply(timepoints, function(t) {
      data.frame(
        risk_group = rg,
        time = t,
        n = nrow(df_subset),
        CI_CRE = predict(ci_fit, cause = 1, times = t),
        CI_discharge = predict(ci_fit, cause = 2, times = t),
        CI_death = predict(ci_fit, cause = 3, times = t),
        survival = 1 - predict(ci_fit, cause = 3, times = t)
      )
    })
    bind_rows(results)
  }
})

detailed_table <- bind_rows(detailed_outcomes) %>%
  mutate(across(c(starts_with("CI_"), survival), ~ round(. * 100, 1)))

cat("\n=== Cumulative Incidence and Survival (%) by Risk Group and Time ===\n")
print(detailed_table)

# ============================================================
# 11. Statistical Comparison Between Groups
# ============================================================
# Pairwise comparisons
if(length(unique(df_model$risk_group_tertile)) > 2) {
  cat("\n=== Pairwise Log-rank Tests ===\n")
  
  # Low vs Medium
  df_low_med <- filter(df_surv, risk_group_tertile %in% c("Low", "Medium"))
  lr_low_med <- survdiff(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_low_med)
  cat("Low vs Medium: p =", round(1 - pchisq(lr_low_med$chisq, 1), 4), "\n")
  
  # Low vs High
  df_low_high <- filter(df_surv, risk_group_tertile %in% c("Low", "High"))
  lr_low_high <- survdiff(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_low_high)
  cat("Low vs High: p =", round(1 - pchisq(lr_low_high$chisq, 1), 4), "\n")
  
  # Medium vs High
  df_med_high <- filter(df_surv, risk_group_tertile %in% c("Medium", "High"))
  lr_med_high <- survdiff(Surv(surv_time, surv_event) ~ risk_group_tertile, data = df_med_high)
  cat("Medium vs High: p =", round(1 - pchisq(lr_med_high$chisq, 1), 4), "\n")
}

# ============================================================
# 12. Create Combined Summary Plot
# ============================================================
# Prepare data for all outcomes
all_outcomes <- lapply(levels(df_model$risk_group_tertile), function(rg) {
  df_subset <- filter(df_model, risk_group_tertile == rg)
  if(nrow(df_subset) > 0) {
    ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_subset)
    bind_rows(
      data.frame(
        times = times_eval,
        CI = predict(ci_fit, cause = 1, times = times_eval),
        risk_group = rg,
        outcome = "CRE infection"
      ),
      data.frame(
        times = times_eval,
        CI = predict(ci_fit, cause = 3, times = times_eval),
        risk_group = rg,
        outcome = "Death"
      ),
      data.frame(
        times = times_eval,
        CI = 1 - predict(ci_fit, cause = 3, times = times_eval),
        risk_group = rg,
        outcome = "Survival"
      )
    )
  }
})

combined_df <- bind_rows(all_outcomes)

p_combined <- ggplot(combined_df, aes(x = times, y = CI, color = risk_group)) +
  geom_line(size = 1.2) +
  facet_wrap(~outcome, scales = "free_y") +
  labs(
    x = "Days since OLT",
    y = "Probability / Cumulative incidence",
    title = "Outcomes by Predicted CRE Risk Group",
    color = "Risk Group"
  ) +
  scale_color_manual(values = c("Low" = "green3", "Medium" = "#FFA500", "High" = "red")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

print(p_combined)

# ============================================================
# 13. Risk Score Performance Metrics
# ============================================================
# C-statistic for discrimination
library(survAUC)

# For survival outcome
c_stat_surv <- concordance(Surv(surv_time, surv_event) ~ pred_cre_risk, data = df_surv)
cat("\n=== Model Performance ===\n")
cat("C-statistic for mortality prediction:", round(c_stat_surv$concordance, 3), "\n")
cat("Standard error:", round(sqrt(c_stat_surv$var), 3), "\n")

# ============================================================
# 14. Export Key Results
# ============================================================
# Save summary statistics
write.csv(survival_summary, "survival_by_risk_group.csv", row.names = FALSE)
write.csv(detailed_table, "detailed_outcomes_by_timepoint.csv", row.names = FALSE)
write.csv(hr_table, "hazard_ratios.csv", row.names = FALSE)

# Save patient-level data with risk scores
df_export <- df_model %>%
  select(time_c, outcome_comp, pred_cre_risk, risk_group_tertile, 
         death_180d, retro_group) %>%
  mutate(patient_id = row_number())

write.csv(df_export, "patient_risk_scores.csv", row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Key findings have been saved to CSV files.\n")
cat("Total patients analyzed:", nrow(df_model), "\n")
cat("Plots show survival and competing outcomes stratified by CRE risk.\n")



