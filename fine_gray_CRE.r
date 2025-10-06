#############################################
## FINE-GRAY COMPETING RISKS ANALYSIS
## Outcome: CRE Infection 
## Exposure: CRE Management Strategy
## Competing risks: Death without CRE, Discharge without CRE
## Based on original CRECOLT analysis framework
#############################################

library(dplyr)
library(cmprsk)
library(riskRegression)
library(prodlim)
library(ggplot2)
library(tidyr)

cat("============================================\n")
cat("COMPETING RISKS ANALYSIS: CRE INFECTION\n")
cat("Impact of Management Strategies\n")
cat("============================================\n")

# Load data
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}

# Use data_extract if baseline_data doesn't have the variables we need
if(!"cre_colonization___1" %in% names(baseline_data)) {
  if(file.exists("data_extract.rds")) {
    data_extract <- readRDS("data_extract.rds")
    data_analysis <- data_extract
    cat("Using data_extract.rds\n")
  } else if(exists("data_extract")) {
    data_analysis <- data_extract
    cat("Using existing data_extract object\n")
  } else {
    stop("Cannot find required variables. Please ensure data_extract is loaded.")
  }
} else {
  data_analysis <- baseline_data
  cat("Using baseline_data\n")
}

# ============================================================
# DIAGNOSTIC SECTION: IDENTIFY AVAILABLE VARIABLES
# ============================================================

cat("\n=== Variable Availability Check ===\n")

# Check for key variables
key_vars <- c(
  "cre_colonization___1", "cre_colonization___2", "cre_colonization___3",  # Management strategies
  "post_olt_compli___1", "post_olt_compli___3", "post_olt_compli___5",     # Complications
  "multisite_colonization",                                                 # Multisite
  "cre_colon", "cre_colon___1", "cre_colon___2",                          # Colonization timing
  "mec_of_carbapenem_resi___1",                                            # Resistance mechanism
  "time_c", "outcome_comp",                                                # Outcome variables
  "cre_infection_date", "death_date", "discharge_date",                    # Date variables
  "retro_or_pros",                                                         # Study design
  "age", "age_clean", "meld", "meld_clean"                                # Patient characteristics
)

cat("\nChecking for key variables:\n")
for(v in key_vars) {
  if(v %in% names(data_analysis)) {
    n_valid <- sum(!is.na(data_analysis[[v]]))
    cat("  ✓ ", v, " (n=", n_valid, ")\n", sep="")
  }
}

cat("\n--- Variables NOT found (but checked for): ---\n")
for(v in key_vars) {
  if(!v %in% names(data_analysis)) {
    cat("  ✗ ", v, "\n", sep="")
  }
}

cat("\n--- First 20 actual variable names in dataset: ---\n")
cat(paste(head(names(data_analysis), 20), collapse=", "), "\n")

# ============================================================
# STEP 1: CREATE VARIABLES AND PREPARE DATA
# ============================================================

cat("\n=== Preparing Data ===\n")

# Create management strategy variable
data_analysis$cre_management_strategy <- NA
data_analysis$cre_management_strategy[data_analysis$cre_colonization___1 == 1] <- 1
data_analysis$cre_management_strategy[data_analysis$cre_colonization___2 == 1] <- 2
data_analysis$cre_management_strategy[data_analysis$cre_colonization___3 == 1] <- 3

# Check for multiple strategies
multiple_strategies <- rowSums(data_analysis[, c("cre_colonization___1", 
                                                 "cre_colonization___2", 
                                                 "cre_colonization___3")], na.rm = TRUE)
n_multiple <- sum(multiple_strategies > 1, na.rm = TRUE)
if(n_multiple > 0) {
  cat("WARNING: ", n_multiple, " patients have multiple strategies recorded\n", sep="")
  cat("Using highest strategy number for these patients\n")
}

# Create factor with meaningful labels
data_analysis$cre_strategy <- factor(
  data_analysis$cre_management_strategy,
  levels = c(1, 2, 3),
  labels = c("No Strategy", "Targeted Prophylaxis", "Pre-emptive Strategy")
)

cat("\n--- Management Strategy Distribution ---\n")
print(table(data_analysis$cre_strategy, useNA = "ifany"))

# Create covariates (matching original analysis)
# First check what variables are available
cat("\nChecking available variables...\n")
available_vars <- names(data_analysis)

# Create safe versions of variables (check if they exist first)
df_model <- data_analysis

# Post-OLT complications
if("post_olt_compli___1" %in% available_vars) {
  df_model$post_olt_compli1 <- ifelse(data_analysis$post_olt_compli___1 == 1, 1, 0)
} else {
  df_model$post_olt_compli1 <- NA
  cat("Warning: post_olt_compli___1 not found\n")
}

if("post_olt_compli___3" %in% available_vars) {
  df_model$post_olt_compli3 <- ifelse(data_analysis$post_olt_compli___3 == 1, 1, 0)
} else {
  df_model$post_olt_compli3 <- NA
  cat("Warning: post_olt_compli___3 not found\n")
}

if("post_olt_compli___5" %in% available_vars) {
  df_model$post_olt_compli5 <- ifelse(data_analysis$post_olt_compli___5 == 1, 1, 0)
} else {
  df_model$post_olt_compli5 <- NA
  cat("Warning: post_olt_compli___5 not found\n")
}

# Multisite colonization
if("multisite_colonization" %in% available_vars) {
  df_model$multisite_col <- ifelse(data_analysis$multisite_colonization == 1, 1, 0)
} else {
  df_model$multisite_col <- NA
  cat("Warning: multisite_colonization not found\n")
}

# CRE colonization timing - try different possible names
if("cre_colon" %in% available_vars) {
  # If cre_colon is a single variable with values 1,2,3
  df_model$cre_col1 <- ifelse(data_analysis$cre_colon == 1, 1, 0)  # Prior
  df_model$cre_col2 <- ifelse(data_analysis$cre_colon == 2, 1, 0)  # Post
} else if("cre_colonization___1" %in% available_vars) {
  # These might be the management strategies, not timing
  df_model$cre_col1 <- NA
  df_model$cre_col2 <- NA
  cat("Note: cre_colon variables not found, skipping\n")
} else {
  df_model$cre_col1 <- NA
  df_model$cre_col2 <- NA
  cat("Warning: CRE colonization timing variables not found\n")
}

# Mechanism of resistance
if("mec_of_carbapenem_resi___1" %in% available_vars) {
  df_model$mec1 <- ifelse(data_analysis$mec_of_carbapenem_resi___1 == 1, 1, 0)
} else {
  df_model$mec1 <- NA
  cat("Warning: mec_of_carbapenem_resi___1 not found\n")
}

# Study design
if("retro_or_pros" %in% available_vars) {
  df_model$retro_group <- factor(ifelse(data_analysis$retro_or_pros == "retro", 
                                        "Retrospective", "Prospective"))
} else {
  df_model$retro_group <- NA
}

# Patient characteristics
if("age_clean" %in% available_vars) {
  df_model$age_years <- data_analysis$age_clean
} else if("age" %in% available_vars) {
  df_model$age_years <- data_analysis$age
} else {
  df_model$age_years <- NA
}

if("meld_clean" %in% available_vars) {
  df_model$meld_score <- data_analysis$meld_clean
} else if("meld" %in% available_vars) {
  df_model$meld_score <- data_analysis$meld
} else {
  df_model$meld_score <- NA
}

cat("\nVariables created for model:\n")
model_vars <- c("post_olt_compli1", "post_olt_compli3", "post_olt_compli5",
                "multisite_col", "cre_col1", "cre_col2", "mec1")
for(v in model_vars) {
  if(v %in% names(df_model)) {
    n_valid <- sum(!is.na(df_model[[v]]))
    if(n_valid > 0) {
      cat("  ", v, ": ", n_valid, " valid values\n", sep="")
    }
  }
}

# Use existing outcome_comp and time_c if available, otherwise create them
if(!"outcome_comp" %in% names(df_model) || !"time_c" %in% names(df_model)) {
  cat("\nCreating competing risk outcomes...\n")
  
  # Ensure dates are in Date format
  df_model$date_of_transplant <- as.Date(df_model$date_of_transplant)
  df_model$cre_infection_date <- as.Date(df_model$cre_infection_date)
  df_model$discharge_date <- as.Date(df_model$discharge_date)
  df_model$death_date <- as.Date(df_model$death_date)
  
  # Define 180-day cutoff
  df_model$day180 <- df_model$date_of_transplant + 180
  
  # Find earliest event
  event_date <- pmin(df_model$cre_infection_date,
                     df_model$discharge_date,
                     df_model$death_date,
                     df_model$day180,
                     na.rm = TRUE)
  
  # Time to event
  df_model$time_c <- as.numeric(event_date - df_model$date_of_transplant)
  
  # Initialize outcome
  df_model$outcome_comp <- 0
  
  # Assign events based on earliest event type
  df_model$outcome_comp[!is.na(df_model$cre_infection_date) &
                          df_model$cre_infection_date <= df_model$day180 &
                          df_model$cre_infection_date == event_date] <- 1
  
  df_model$outcome_comp[!is.na(df_model$discharge_date) &
                          df_model$discharge_date <= df_model$day180 &
                          df_model$discharge_date == event_date] <- 2
  
  df_model$outcome_comp[!is.na(df_model$death_date) &
                          df_model$death_date <= df_model$day180 &
                          df_model$death_date == event_date] <- 3
}

# Filter to complete cases for key variables
df_model_complete <- df_model %>%
  filter(!is.na(time_c) & time_c >= 0 & !is.na(cre_strategy))

cat("\n--- Outcome Distribution ---\n")
outcome_table <- table(df_model_complete$outcome_comp, useNA = "ifany")
cat("0 = Censored: ", outcome_table[1], "\n")
cat("1 = CRE infection: ", outcome_table[2], "\n")
cat("2 = Discharge: ", outcome_table[3], "\n")
cat("3 = Death: ", outcome_table[4], "\n")

cat("\nTotal patients with complete strategy data: ", nrow(df_model_complete), "\n")

# ============================================================
# STEP 2: CUMULATIVE INCIDENCE BY STRATEGY (UNADJUSTED)
# ============================================================

cat("\n=== Cumulative Incidence of CRE Infection by Strategy ===\n")

# Calculate cumulative incidence for each strategy
times_seq <- seq(0, 180, by = 5)

ci_by_strategy <- lapply(levels(df_model_complete$cre_strategy), function(strat) {
  df_subset <- filter(df_model_complete, cre_strategy == strat)
  if(nrow(df_subset) > 10) {
    ci_fit <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_subset)
    data.frame(
      times = times_seq,
      CIF_cre = predict(ci_fit, cause = 1, times = times_seq),
      CIF_death = predict(ci_fit, cause = 3, times = times_seq),
      CIF_discharge = predict(ci_fit, cause = 2, times = times_seq),
      strategy = strat,
      n = nrow(df_subset)
    )
  }
})

ci_strategy_df <- bind_rows(ci_by_strategy)

# Summary at key timepoints
key_times <- c(14, 30, 60, 90, 180)
summary_table <- ci_strategy_df %>%
  filter(times %in% key_times) %>%
  select(strategy, times, CIF_cre, n) %>%
  mutate(CIF_cre_pct = round(CIF_cre * 100, 1)) %>%
  pivot_wider(names_from = times, values_from = CIF_cre_pct, names_prefix = "Day_")

cat("\n--- CRE Infection Cumulative Incidence (%) by Strategy ---\n")
print(as.data.frame(summary_table))

# Plot cumulative incidence curves
strategy_colors <- c("No Strategy" = "#E41A1C", 
                     "Targeted Prophylaxis" = "#377EB8", 
                     "Pre-emptive Strategy" = "#4DAF4A")

p_ci_strategy <- ggplot(ci_strategy_df, aes(x = times, y = CIF_cre, color = strategy)) +
  geom_line(size = 1.2) +
  labs(
    x = "Days since OLT",
    y = "Cumulative incidence of CRE infection",
    title = "CRE Infection by Management Strategy (Unadjusted)",
    color = "Strategy"
  ) +
  scale_color_manual(values = strategy_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

print(p_ci_strategy)
ggsave("cre_infection_by_strategy_unadjusted.png", p_ci_strategy, 
       width = 10, height = 7, dpi = 300)

# ============================================================
# STEP 3: GRAY'S TEST FOR EQUALITY OF CIF
# ============================================================

cat("\n=== Gray's Test for Equality of CIF ===\n")

# Prepare data for cmprsk
time <- df_model_complete$time_c
status <- df_model_complete$outcome_comp
group <- as.numeric(df_model_complete$cre_strategy)

# Gray's test
gray_test <- cuminc(time, status, group, cencode = 0)

cat("\nTesting equality of CRE infection incidence across strategies:\n")
print(gray_test$Tests)

# Extract p-value
gray_p <- gray_test$Tests[1, "pv"]
cat("\nGray's test p-value for CRE infection: ", format.pval(gray_p), "\n")

if(gray_p < 0.05) {
  cat("SIGNIFICANT: Strategies differ in CRE infection incidence\n")
} else {
  cat("NOT SIGNIFICANT: No evidence of difference between strategies\n")
}

# ============================================================
# STEP 4: FINE-GRAY MODEL - UNADJUSTED
# ============================================================

cat("\n=== Fine-Gray Model: Unadjusted ===\n")

# Fit Fine-Gray model with strategy as only predictor
fg_unadjusted <- FGR(
  Hist(time_c, outcome_comp) ~ cre_strategy,
  data = df_model_complete,
  cause = 1
)

# Extract coefficients
fg_coef_unadj <- summary(fg_unadjusted)
print(fg_coef_unadj)

# Calculate subdistribution hazard ratios
cat("\n--- Subdistribution Hazard Ratios (Unadjusted) ---\n")
cat("Reference: No Strategy\n")

coef_vals <- fg_unadjusted$crrFit$coef
se_vals <- sqrt(diag(fg_unadjusted$crrFit$var))

for(i in 1:length(coef_vals)) {
  hr <- exp(coef_vals[i])
  ci_low <- exp(coef_vals[i] - 1.96 * se_vals[i])
  ci_high <- exp(coef_vals[i] + 1.96 * se_vals[i])
  z_val <- coef_vals[i] / se_vals[i]
  p_val <- 2 * (1 - pnorm(abs(z_val)))
  
  cat(names(coef_vals)[i], ":\n")
  cat("  SHR = ", round(hr, 2), " (95% CI: ", round(ci_low, 2), 
      "-", round(ci_high, 2), "), p = ", format.pval(p_val), "\n", sep="")
}

# ============================================================
# STEP 5: FINE-GRAY MODEL - ADJUSTED
# ============================================================

cat("\n=== Fine-Gray Model: Adjusted for Risk Factors ===\n")

# Select covariates for adjustment (only those that exist and have data)
potential_covariates <- c("post_olt_compli1", "post_olt_compli3", "post_olt_compli5",
                          "multisite_col", "cre_col1", "cre_col2", "mec1", 
                          "age_years", "meld_score")

# Check which covariates are available with sufficient data
available_covars <- c()
for(covar in potential_covariates) {
  if(covar %in% names(df_model_complete)) {
    n_valid <- sum(!is.na(df_model_complete[[covar]]))
    if(n_valid > 10) {  # Only use if we have at least 10 non-missing values
      available_covars <- c(available_covars, covar)
    }
  }
}

cat("Available covariates for adjustment:\n")
if(length(available_covars) > 0) {
  print(available_covars)
} else {
  cat("  None - will perform unadjusted analysis only\n")
}

# Create formula and fit model if we have covariates
if(length(available_covars) > 0) {
  formula_adjusted <- as.formula(paste("Hist(time_c, outcome_comp) ~ cre_strategy +", 
                                       paste(available_covars, collapse = " + ")))
  
  # Filter to complete cases for all variables
  vars_needed <- c("time_c", "outcome_comp", "cre_strategy", available_covars)
  df_model_adjusted <- df_model_complete %>%
    select(all_of(vars_needed)) %>%
    filter(complete.cases(.))
  
  cat("\nPatients with complete data for adjusted model: ", nrow(df_model_adjusted), "\n")
  
  if(nrow(df_model_adjusted) > 50) {  # Need reasonable sample size
    # Fit adjusted Fine-Gray model
    fg_adjusted <- tryCatch({
      FGR(formula_adjusted, data = df_model_adjusted, cause = 1)
    }, error = function(e) {
      cat("Error fitting adjusted model: ", e$message, "\n")
      cat("Falling back to strategy-only model\n")
      FGR(Hist(time_c, outcome_comp) ~ cre_strategy, data = df_model_adjusted, cause = 1)
    })
    
    # Extract coefficients
    fg_coef_adj <- summary(fg_adjusted)
    print(fg_coef_adj)
    
    # Extract strategy effects
    cat("\n--- Subdistribution Hazard Ratios for Strategies (Adjusted) ---\n")
    cat("Reference: No Strategy\n")
    
    coef_vals_adj <- fg_adjusted$crrFit$coef
    se_vals_adj <- sqrt(diag(fg_adjusted$crrFit$var))
    
    # Find strategy coefficients
    strategy_indices <- grep("cre_strategy", names(coef_vals_adj))
    
    if(length(strategy_indices) > 0) {
      for(i in strategy_indices) {
        hr <- exp(coef_vals_adj[i])
        ci_low <- exp(coef_vals_adj[i] - 1.96 * se_vals_adj[i])
        ci_high <- exp(coef_vals_adj[i] + 1.96 * se_vals_adj[i])
        z_val <- coef_vals_adj[i] / se_vals_adj[i]
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        
        strategy_name <- gsub("cre_strategy", "", names(coef_vals_adj)[i])
        cat(strategy_name, ":\n")
        cat("  SHR = ", round(hr, 2), " (95% CI: ", round(ci_low, 2), 
            "-", round(ci_high, 2), "), p = ", format.pval(p_val), "\n", sep="")
      }
    }
  } else {
    cat("Insufficient data for adjusted model\n")
    df_model_adjusted <- df_model_complete
    fg_adjusted <- fg_unadjusted
  }
} else {
  cat("\nNo covariates available for adjustment - using unadjusted model\n")
  df_model_adjusted <- df_model_complete
  fg_adjusted <- fg_unadjusted
}

# ============================================================
# STEP 6: PREDICTED VS OBSERVED CALIBRATION
# ============================================================

cat("\n=== Model Calibration ===\n")

# Predict CIF for all patients
pred_cif <- predict(
  fg_adjusted,
  newdata = df_model_adjusted,
  cause = 1,
  times = times_seq,
  se = FALSE,
  keep.newdata = FALSE
)

# Average predicted CIF
pred_df <- data.frame(
  times = times_seq,
  CIF_pred = colMeans(pred_cif, na.rm = TRUE)
)

# Observed CIF
ci_fit_all <- prodlim(Hist(time_c, outcome_comp) ~ 1, data = df_model_adjusted)
obs_df <- data.frame(
  times = times_seq,
  CIF_obs = predict(ci_fit_all, cause = 1, times = times_seq)
)

calib_df <- left_join(obs_df, pred_df, by = "times")

# Calibration plot
p_calibration <- ggplot(calib_df, aes(x = times)) +
  geom_line(aes(y = CIF_obs, color = "Observed"), size = 1.2) +
  geom_line(aes(y = CIF_pred, color = "Predicted"), 
            size = 1.2, linetype = "dashed") +
  labs(
    x = "Days since OLT",
    y = "Cumulative incidence of CRE infection",
    title = "Model Calibration: Observed vs Predicted"
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "#D81B60")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(size = 14, face = "bold")
  )

print(p_calibration)
ggsave("cre_model_calibration.png", p_calibration, 
       width = 10, height = 7, dpi = 300)

# ============================================================
# STEP 7: STRATIFIED ANALYSIS BY COLLECTION METHOD
# ============================================================

if("retro_group" %in% names(df_model_adjusted) && !all(is.na(df_model_adjusted$retro_group))) {
  
  cat("\n=== Stratified Analysis by Study Design ===\n")
  
  # Retrospective cohort
  df_retro <- filter(df_model_adjusted, retro_group == "Retrospective")
  if(nrow(df_retro) > 30) {
    cat("\n--- RETROSPECTIVE COHORT (n=", nrow(df_retro), ") ---\n", sep="")
    
    fg_retro <- FGR(
      formula_adjusted,
      data = df_retro,
      cause = 1
    )
    
    # Extract strategy effects
    coef_retro <- fg_retro$crrFit$coef
    se_retro <- sqrt(diag(fg_retro$crrFit$var))
    strategy_idx_retro <- grep("cre_strategy", names(coef_retro))
    
    for(i in strategy_idx_retro) {
      hr <- exp(coef_retro[i])
      ci_low <- exp(coef_retro[i] - 1.96 * se_retro[i])
      ci_high <- exp(coef_retro[i] + 1.96 * se_retro[i])
      p_val <- 2 * (1 - pnorm(abs(coef_retro[i] / se_retro[i])))
      
      strategy_name <- gsub("cre_strategy", "", names(coef_retro)[i])
      cat(strategy_name, ": SHR = ", round(hr, 2), 
          " (", round(ci_low, 2), "-", round(ci_high, 2), 
          "), p = ", format.pval(p_val), "\n", sep="")
    }
  }
  
  # Prospective cohort
  df_prosp <- filter(df_model_adjusted, retro_group == "Prospective")
  if(nrow(df_prosp) > 30) {
    cat("\n--- PROSPECTIVE COHORT (n=", nrow(df_prosp), ") ---\n", sep="")
    
    fg_prosp <- FGR(
      formula_adjusted,
      data = df_prosp,
      cause = 1
    )
    
    # Extract strategy effects
    coef_prosp <- fg_prosp$crrFit$coef
    se_prosp <- sqrt(diag(fg_prosp$crrFit$var))
    strategy_idx_prosp <- grep("cre_strategy", names(coef_prosp))
    
    for(i in strategy_idx_prosp) {
      hr <- exp(coef_prosp[i])
      ci_low <- exp(coef_prosp[i] - 1.96 * se_prosp[i])
      ci_high <- exp(coef_prosp[i] + 1.96 * se_prosp[i])
      p_val <- 2 * (1 - pnorm(abs(coef_prosp[i] / se_prosp[i])))
      
      strategy_name <- gsub("cre_strategy", "", names(coef_prosp)[i])
      cat(strategy_name, ": SHR = ", round(hr, 2), 
          " (", round(ci_low, 2), "-", round(ci_high, 2), 
          "), p = ", format.pval(p_val), "\n", sep="")
    }
  }
}

# ============================================================
# STEP 8: NUMBER NEEDED TO TREAT CALCULATION
# ============================================================

cat("\n=== Number Needed to Treat (NNT) Analysis ===\n")

# Calculate 30, 60, and 180-day CRE incidence by strategy
nnt_times <- c(30, 60, 180)

nnt_data <- ci_strategy_df %>%
  filter(times %in% nnt_times) %>%
  select(strategy, times, CIF_cre) %>%
  pivot_wider(names_from = strategy, values_from = CIF_cre)

cat("\n--- Cumulative Incidence by Strategy ---\n")
print(nnt_data)

# Calculate ARR and NNT compared to no strategy
for(t in nnt_times) {
  row_data <- nnt_data[nnt_data$times == t, ]
  
  if(nrow(row_data) > 0 && "No Strategy" %in% names(row_data)) {
    no_strategy_risk <- row_data[["No Strategy"]]
    
    # Check if no_strategy_risk is valid
    if(!is.na(no_strategy_risk)) {
      
      # Targeted Prophylaxis
      if("Targeted Prophylaxis" %in% names(row_data)) {
        targeted_risk <- row_data[["Targeted Prophylaxis"]]
        
        if(!is.na(targeted_risk)) {
          arr_targeted <- no_strategy_risk - targeted_risk
          
          cat("\nDay ", t, " - Targeted Prophylaxis:\n", sep="")
          cat("  No Strategy risk: ", round(no_strategy_risk * 100, 1), "%\n", sep="")
          cat("  Targeted risk: ", round(targeted_risk * 100, 1), "%\n", sep="")
          cat("  ARR = ", round(arr_targeted * 100, 1), "%\n", sep="")
          
          if(!is.na(arr_targeted) && arr_targeted > 0.001) {  # Small positive threshold
            nnt_targeted <- 1 / arr_targeted
            cat("  NNT = ", round(nnt_targeted), "\n", sep="")
          } else if(!is.na(arr_targeted) && arr_targeted < -0.001) {  # Negative threshold
            nnh_targeted <- -1 / arr_targeted
            cat("  NNH = ", round(nnh_targeted), " (increases risk)\n", sep="")
          } else {
            cat("  No meaningful difference\n")
          }
        } else {
          cat("\nDay ", t, " - Targeted Prophylaxis: No data available\n", sep="")
        }
      }
      
      # Pre-emptive Strategy
      if("Pre-emptive Strategy" %in% names(row_data)) {
        preemptive_risk <- row_data[["Pre-emptive Strategy"]]
        
        if(!is.na(preemptive_risk)) {
          arr_preemptive <- no_strategy_risk - preemptive_risk
          
          cat("\nDay ", t, " - Pre-emptive Strategy:\n", sep="")
          cat("  No Strategy risk: ", round(no_strategy_risk * 100, 1), "%\n", sep="")
          cat("  Pre-emptive risk: ", round(preemptive_risk * 100, 1), "%\n", sep="")
          cat("  ARR = ", round(arr_preemptive * 100, 1), "%\n", sep="")
          
          if(!is.na(arr_preemptive) && arr_preemptive > 0.001) {  # Small positive threshold
            nnt_preemptive <- 1 / arr_preemptive
            cat("  NNT = ", round(nnt_preemptive), "\n", sep="")
          } else if(!is.na(arr_preemptive) && arr_preemptive < -0.001) {  # Negative threshold
            nnh_preemptive <- -1 / arr_preemptive
            cat("  NNH = ", round(nnh_preemptive), " (increases risk)\n", sep="")
          } else {
            cat("  No meaningful difference\n")
          }
        } else {
          cat("\nDay ", t, " - Pre-emptive Strategy: No data available\n", sep="")
        }
      }
    } else {
      cat("\nDay ", t, ": No Strategy baseline risk not available\n", sep="")
    }
  } else {
    cat("\nDay ", t, ": Insufficient data for NNT calculation\n", sep="")
  }
}

# ============================================================
# STEP 9: SUMMARY AND CLINICAL IMPLICATIONS
# ============================================================

cat("\n============================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("============================================\n")

cat("\n1. OVERALL EFFECT ON CRE INFECTION:\n")
cat("   Gray's test p-value: ", format.pval(gray_p), "\n")

cat("\n2. UNADJUSTED SUBDISTRIBUTION HAZARD RATIOS:\n")
cat("   (See results above)\n")

cat("\n3. ADJUSTED SUBDISTRIBUTION HAZARD RATIOS:\n")
cat("   (See results above)\n")

cat("\n4. CLINICAL INTERPRETATION:\n")

# Determine effectiveness based on results
targeted_effective <- FALSE
preemptive_effective <- FALSE

if(length(strategy_indices) >= 1) {
  # Check if Targeted Prophylaxis reduces risk (HR < 1)
  targeted_idx <- grep("Targeted", names(coef_vals_adj))
  if(length(targeted_idx) > 0) {
    targeted_hr <- exp(coef_vals_adj[targeted_idx[1]])
    if(targeted_hr < 1) {
      targeted_effective <- TRUE
      cat("   - Targeted Prophylaxis REDUCES CRE infection risk\n")
    } else {
      cat("   - Targeted Prophylaxis does NOT reduce CRE infection risk\n")
    }
  }
  
  # Check if Pre-emptive reduces risk (HR < 1)
  preemptive_idx <- grep("Pre-emptive", names(coef_vals_adj))
  if(length(preemptive_idx) > 0) {
    preemptive_hr <- exp(coef_vals_adj[preemptive_idx[1]])
    if(preemptive_hr < 1) {
      preemptive_effective <- TRUE
      cat("   - Pre-emptive Strategy REDUCES CRE infection risk\n")
    } else {
      cat("   - Pre-emptive Strategy does NOT reduce CRE infection risk\n")
    }
  }
}

if(!targeted_effective && !preemptive_effective) {
  cat("\n   CONCLUSION: Neither management strategy significantly reduces CRE infection\n")
  cat("   despite being designed for this purpose.\n")
} else if(targeted_effective || preemptive_effective) {
  cat("\n   CONCLUSION: Some strategies show benefit for preventing CRE infection.\n")
}

cat("\n5. PARADOX:\n")
cat("   If strategies don't prevent infection but worsen mortality,\n")
cat("   this suggests potential harm from the interventions themselves.\n")

# ============================================================
# STEP 10: SAVE RESULTS
# ============================================================

# Save key results
results_list <- list(
  gray_test = gray_test,
  fg_unadjusted = fg_unadjusted,
  fg_adjusted = fg_adjusted,
  cumulative_incidence = ci_strategy_df,
  nnt_analysis = nnt_data,
  n_patients = nrow(df_model_complete),
  n_cre_infections = sum(df_model_complete$outcome_comp == 1),
  n_deaths = sum(df_model_complete$outcome_comp == 3)
)

save(results_list, file = "cre_infection_competing_risks_results.RData")

cat("\n============================================\n")
cat("Analysis complete. Results saved.\n")
cat("Plots saved: cre_infection_by_strategy_unadjusted.png\n")
cat("            cre_model_calibration.png\n")

