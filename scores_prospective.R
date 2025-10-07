#############################################
## SCORE TRAJECTORY ANALYSIS - REDCAP STRUCTURE
## Longitudinal Analysis Using Filter Variable
## Main records: filter_$ == 1
## Score records: filter_$ < 1
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(survival)
library(viridis)
library(patchwork)
library(haven)

cat("============================================\n")
cat("SCORE TRAJECTORY ANALYSIS\n")
cat("REDCap Longitudinal Structure\n")
cat("============================================\n")

# Load data
data <- read_sav("CRECOOLT_overall.sav")

# Check the filter variable structure
cat("\n=== Understanding Filter Structure ===\n")
cat("Distribution of filter_$ variable:\n")
print(table(data$`filter_$`, useNA = "ifany"))

# ============================================================
# STEP 1: SEPARATE BASELINE AND LONGITUDINAL DATA
# ============================================================

cat("\n=== Separating Baseline and Longitudinal Data ===\n")

# Extract baseline patient records
baseline_data <- data %>%
  filter(`filter_$` == 1)

cat("Baseline records (filter_$ == 1):", nrow(baseline_data), "\n")
cat("Unique patients in baseline:", n_distinct(baseline_data$record_id), "\n")

# Extract longitudinal score records
score_data <- data %>%
  filter(`filter_$` < 1)

cat("\nScore records (filter_$ < 1):", nrow(score_data), "\n")
cat("Unique patients with scores:", n_distinct(score_data$record_id), "\n")

# Check score data structure
score_summary <- score_data %>%
  group_by(record_id) %>%
  summarise(
    n_scores = sum(!is.na(score)),
    n_score_dates = sum(!is.na(score_date)),
    n_records = n(),
    .groups = "drop"
  )

cat("\nScore data structure:\n")
cat("  Mean score records per patient:", round(mean(score_summary$n_records), 2), "\n")
cat("  Median score records per patient:", median(score_summary$n_records), "\n")
cat("  Range:", min(score_summary$n_records), "-", max(score_summary$n_records), "\n")
cat("  Patients with actual score values:", sum(score_summary$n_scores > 0), "\n")

# ============================================================
# STEP 2: PREPARE BASELINE DATA
# ============================================================

cat("\n=== Preparing Baseline Data ===\n")

# Select key baseline variables
baseline_processed <- baseline_data %>%
  select(
    record_id,
    retro_or_pros,
    date_of_transplant,
    cre_infection,
    cre_infection_date,
    death,
    death_date,
    discharge_date,
    cre_colonization___1,
    cre_colonization___2,
    cre_colonization___3,
    multisite_colonization,
    post_olt_compli___1,
    post_olt_compli___3,
    post_olt_compli___5
  ) %>%
  mutate(
    # Convert dates
    date_of_transplant = as.Date(date_of_transplant, format = "%Y-%m-%d"),
    cre_infection_date = as.Date(cre_infection_date, format = "%Y-%m-%d"),
    death_date = as.Date(death_date, format = "%Y-%m-%d"),
    discharge_date = as.Date(discharge_date, format = "%Y-%m-%d"),
    
    # Calculate days to events
    days_to_cre = as.numeric(cre_infection_date - date_of_transplant),
    days_to_death = as.numeric(death_date - date_of_transplant),
    days_to_discharge = as.numeric(discharge_date - date_of_transplant),
    
    # Define primary outcome
    primary_outcome = case_when(
      cre_infection == 1 & !is.na(days_to_cre) & days_to_cre <= 180 ~ "CRE Infection",
      death == 1 & !is.na(days_to_death) & days_to_death <= 180 ~ "Death",
      !is.na(days_to_discharge) & days_to_discharge <= 180 ~ "Discharge",
      TRUE ~ "Censored"
    ),
    
    # Time to primary outcome
    time_to_outcome = case_when(
      primary_outcome == "CRE Infection" ~ days_to_cre,
      primary_outcome == "Death" ~ days_to_death,
      primary_outcome == "Discharge" ~ days_to_discharge,
      TRUE ~ 180  # Censoring at 180 days
    ),
    
    # Treatment strategy
    treatment_strategy = case_when(
      cre_colonization___1 == 1 ~ "No Prophylaxis",
      cre_colonization___2 == 1 ~ "Prophylaxis",
      cre_colonization___3 == 1 ~ "Pre-emptive",
      TRUE ~ "None"
    )
  )

cat("Baseline cohort distribution:\n")
table_cohort <- table(baseline_processed$retro_or_pros)
cat("  Retrospective (1):", table_cohort["1"], "\n")
cat("  Prospective (2):", table_cohort["2"], "\n")

# ============================================================
# STEP 3: PREPARE LONGITUDINAL SCORE DATA
# ============================================================

cat("\n=== Preparing Longitudinal Score Data ===\n")

score_processed <- score_data %>%
  select(record_id, score, score_date, events, score_impact___1, score_impact___2, score_impact___3) %>%
  filter(!is.na(score) & !is.na(score_date)) %>%
  mutate(
    score_date = as.Date(score_date, format = "%Y-%m-%d")
  )

cat("Score records with valid data:", nrow(score_processed), "\n")

# ============================================================
# STEP 4: JOIN BASELINE AND SCORE DATA
# ============================================================

cat("\n=== Joining Baseline and Score Data ===\n")

# Join datasets
longitudinal_data <- score_processed %>%
  inner_join(baseline_processed, by = "record_id") %>%
  mutate(
    # Calculate days from transplant for each score
    days_from_olt = as.numeric(score_date - date_of_transplant),
    
    # Flag if score led to treatment change
    treatment_impact = case_when(
      score_impact___1 == 1 ~ "Changed therapy",
      score_impact___2 == 1 ~ "Changed diagnostics",
      score_impact___3 == 1 ~ "No change",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(days_from_olt >= 0 & days_from_olt <= 180)  # Focus on first 180 days

cat("Final longitudinal dataset:\n")
cat("  Total score measurements:", nrow(longitudinal_data), "\n")
cat("  Unique patients:", n_distinct(longitudinal_data$record_id), "\n")

# Check cohort distribution in longitudinal data
cohort_dist <- longitudinal_data %>%
  group_by(record_id, retro_or_pros) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(retro_or_pros)

cat("\nPatients with longitudinal scores by cohort:\n")
print(cohort_dist)

# ============================================================
# STEP 5: ANALYZE SCORE TRAJECTORIES
# ============================================================

cat("\n=== Analyzing Score Trajectories ===\n")

# Summary by patient
patient_trajectories <- longitudinal_data %>%
  group_by(record_id, primary_outcome, treatment_strategy, retro_or_pros) %>%
  summarise(
    n_scores = n(),
    baseline_score = first(score),
    final_score = last(score),
    max_score = max(score, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    score_trend = final_score - baseline_score,
    cv_score = sd(score, na.rm = TRUE) / mean_score,
    
    # Early warning signal
    early_max = max(score[days_from_olt <= 30], na.rm = TRUE),
    
    # Treatment changes
    n_treatment_changes = sum(treatment_impact == "Changed therapy", na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    early_max = ifelse(is.infinite(early_max), NA, early_max)
  )

cat("\nTrajectory summary statistics:\n")
cat("  Patients analyzed:", nrow(patient_trajectories), "\n")
cat("  Mean scores per patient:", round(mean(patient_trajectories$n_scores), 2), "\n")

# ============================================================
# STEP 6: VISUALIZATIONS
# ============================================================

cat("\n=== Creating Visualizations ===\n")

# 1. Spaghetti plot by outcome
p_trajectories <- ggplot(longitudinal_data, 
                         aes(x = days_from_olt, y = score, group = record_id)) +
  geom_line(alpha = 0.3, size = 0.5) +
  geom_smooth(aes(group = NULL), method = "loess", se = TRUE, 
              color = "red", size = 1.2) +
  facet_wrap(~ primary_outcome, scales = "free_y") +
  labs(
    title = "Score Trajectories by Patient Outcome",
    subtitle = "Individual trajectories with smoothed average",
    x = "Days from Transplant",
    y = "Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 11, face = "bold")
  )

print(p_trajectories)
ggsave("score_trajectories_by_outcome.png", p_trajectories, 
       width = 12, height = 8, dpi = 300)

# 2. Treatment impact over time
treatment_impact_summary <- longitudinal_data %>%
  filter(treatment_impact != "Unknown") %>%
  mutate(
    time_bin = cut(days_from_olt, 
                   breaks = c(0, 30, 60, 90, 120, 150, 180),
                   include.lowest = TRUE,
                   labels = c("0-30", "31-60", "61-90", "91-120", "121-150", "151-180"))
  ) %>%
  group_by(time_bin, treatment_impact) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    .groups = "drop"
  )

p_treatment_impact <- ggplot(treatment_impact_summary, 
                             aes(x = time_bin, y = n, fill = treatment_impact)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Treatment Decisions Based on Score Over Time",
    x = "Days from Transplant",
    y = "Number of Decisions",
    fill = "Decision Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_treatment_impact)
ggsave("treatment_impact_timeline.png", p_treatment_impact, 
       width = 10, height = 6, dpi = 300)

# 3. Score trajectories by treatment strategy
p_by_strategy <- ggplot(longitudinal_data %>% filter(treatment_strategy != "None"), 
                        aes(x = days_from_olt, y = score, color = treatment_strategy)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2) +
  facet_wrap(~ primary_outcome) +
  labs(
    title = "Mean Score Trajectories by Treatment Strategy and Outcome",
    x = "Days from Transplant",
    y = "Score",
    color = "Treatment Strategy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

print(p_by_strategy)
ggsave("trajectories_by_strategy.png", p_by_strategy, 
       width = 12, height = 8, dpi = 300)

# ============================================================
# STEP 7: STATISTICAL ANALYSIS
# ============================================================

cat("\n=== Statistical Analysis ===\n")

# Test if early high scores predict CRE infection
if(sum(patient_trajectories$primary_outcome == "CRE Infection") > 5) {
  
  # Logistic regression: early scores → CRE infection
  model_data <- patient_trajectories %>%
    filter(!is.na(early_max)) %>%
    mutate(
      cre_outcome = as.integer(primary_outcome == "CRE Infection"),
      early_max_scaled = scale(early_max),
      mean_score_scaled = scale(mean_score),
      cv_score_scaled = scale(cv_score)
    )
  
  logit_model <- glm(
    cre_outcome ~ early_max_scaled + mean_score_scaled + cv_score_scaled + 
      treatment_strategy + n_treatment_changes,
    data = model_data,
    family = binomial()
  )
  
  cat("\n--- Logistic Regression: Score Features → CRE Infection ---\n")
  print(summary(logit_model))
  
  # Calculate odds ratios
  or_table <- exp(cbind(OR = coef(logit_model), confint(logit_model)))
  cat("\nOdds Ratios (95% CI):\n")
  print(round(or_table, 3))
}

# Compare scores between those who had treatment changes vs not
cat("\n--- Score Comparison: Treatment Change vs No Change ---\n")

change_comparison <- longitudinal_data %>%
  group_by(treatment_impact) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

print(change_comparison)

# Test difference
if(sum(longitudinal_data$treatment_impact == "Changed therapy", na.rm = TRUE) > 10) {
  therapy_scores <- longitudinal_data$score[longitudinal_data$treatment_impact == "Changed therapy"]
  no_change_scores <- longitudinal_data$score[longitudinal_data$treatment_impact == "No change"]
  
  wilcox_test <- wilcox.test(therapy_scores, no_change_scores)
  cat("\nWilcoxon test (Changed therapy vs No change):\n")
  cat("  p-value:", format.pval(wilcox_test$p.value), "\n")
  
  if(wilcox_test$p.value < 0.05) {
    cat("  SIGNIFICANT: Scores differ between treatment change decisions\n")
  }
}

# ============================================================
# STEP 8: COHORT COMPARISON
# ============================================================

cat("\n=== Retrospective vs Prospective Comparison ===\n")

cohort_comparison <- patient_trajectories %>%
  group_by(retro_or_pros) %>%
  summarise(
    n_patients = n(),
    mean_n_scores = mean(n_scores),
    mean_baseline = mean(baseline_score, na.rm = TRUE),
    mean_max = mean(max_score, na.rm = TRUE),
    prop_treatment_change = mean(n_treatment_changes > 0, na.rm = TRUE),
    
    # Outcomes
    prop_cre = mean(primary_outcome == "CRE Infection"),
    prop_death = mean(primary_outcome == "Death"),
    
    .groups = "drop"
  ) %>%
  mutate(
    cohort = ifelse(retro_or_pros == 1, "Retrospective", "Prospective")
  )

cat("\nCohort characteristics:\n")
print(cohort_comparison)

# ============================================================
# STEP 9: TIME-VARYING ANALYSIS
# ============================================================

cat("\n=== Time-Varying Score Analysis ===\n")

# Landmark analysis at key time points
landmark_times <- c(14, 30, 60)

landmark_analysis <- function(landmark_day) {
  # Get most recent score before landmark
  landmark_data <- longitudinal_data %>%
    filter(days_from_olt <= landmark_day) %>%
    group_by(record_id) %>%
    arrange(desc(days_from_olt)) %>%
    slice(1) %>%
    ungroup() %>%
    filter(time_to_outcome > landmark_day) %>%  # Still at risk
    mutate(
      cre_within_30 = as.integer(
        primary_outcome == "CRE Infection" & 
          time_to_outcome <= (landmark_day + 30)
      )
    )
  
  if(nrow(landmark_data) > 20 & sum(landmark_data$cre_within_30) > 2) {
    return(data.frame(
      landmark = landmark_day,
      n_at_risk = nrow(landmark_data),
      n_cre_30d = sum(landmark_data$cre_within_30),
      mean_score_cre = mean(landmark_data$score[landmark_data$cre_within_30 == 1], na.rm = TRUE),
      mean_score_no_cre = mean(landmark_data$score[landmark_data$cre_within_30 == 0], na.rm = TRUE)
    ))
  }
  return(NULL)
}

landmark_results <- bind_rows(lapply(landmark_times, landmark_analysis))

if(nrow(landmark_results) > 0) {
  cat("\nLandmark Analysis Results:\n")
  cat("Score at landmark → CRE infection within 30 days\n")
  print(landmark_results)
}

# ============================================================
# STEP 10: SUMMARY
# ============================================================

cat("\n============================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("============================================\n")

cat("\n1. DATA STRUCTURE:\n")
cat("   - ", nrow(baseline_processed), " patients in baseline\n")
cat("   - ", n_distinct(longitudinal_data$record_id), " patients with longitudinal scores\n")
cat("   - ", nrow(longitudinal_data), " total score measurements\n")

cat("\n2. TREATMENT IMPACT:\n")
impact_table <- table(longitudinal_data$treatment_impact)
for(impact in names(impact_table)) {
  cat("   - ", impact, ": ", impact_table[impact], " (", 
      round(100 * impact_table[impact] / sum(impact_table), 1), "%)\n", sep="")
}

cat("\n3. KEY INSIGHTS:\n")
cat("   - Score trajectories can be linked to treatment decisions\n")
cat("   - Early high scores may identify high-risk patients\n")
cat("   - Treatment changes occur in response to score patterns\n")

# Save results
trajectory_results <- list(
  baseline_data = baseline_processed,
  longitudinal_data = longitudinal_data,
  patient_trajectories = patient_trajectories,
  cohort_comparison = cohort_comparison,
  landmark_results = landmark_results
)

save(trajectory_results, file = "trajectory_analysis_results.RData")

cat("\n============================================\n")
cat("Analysis complete. Results and plots saved.\n")
cat("============================================\n")


############# Limit to prospectively studied patients ################

#############################################
## SIMPLIFIED SCORE TRAJECTORY ANALYSIS
## PROSPECTIVE COHORT ONLY
## Avoiding filter column naming issues
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)

cat("============================================\n")
cat("SCORE TRAJECTORY ANALYSIS - PROSPECTIVE ONLY\n")
cat("Simplified Approach\n")
cat("============================================\n")

# Load data
data <- read_sav("CRECOOLT_overall.sav")
cat("Total records loaded:", nrow(data), "\n")

# ============================================================
# IDENTIFY FILTER COLUMN
# ============================================================

cat("\n=== Identifying Data Structure ===\n")

# Method 1: Check common filter column names
filter_col <- NULL

# Try different possible names
if("filter_$" %in% names(data)) {
  filter_col <- "filter_$"
  cat("Found filter column: filter_$\n")
} else if("filter" %in% names(data)) {
  filter_col <- "filter"
  cat("Found filter column: filter\n")
} else if("filter_" %in% names(data)) {
  filter_col <- "filter_"
  cat("Found filter column: filter_\n")
} else {
  # Search for any column with filter in the name
  possible <- grep("filter", names(data), ignore.case = TRUE, value = TRUE)
  if(length(possible) > 0) {
    filter_col <- possible[1]
    cat("Found filter column:", filter_col, "\n")
  }
}

if(is.null(filter_col)) {
  cat("\nNo filter column found. Looking for alternative data structure indicators...\n")
  
  # Check for REDCap repeat structure
  if("redcap_repeat_instance" %in% names(data)) {
    cat("Found redcap_repeat_instance - using this for longitudinal structure\n")
    
    # Baseline: no repeat instance or instance = 1
    baseline_data <- subset(data, 
                            (is.na(redcap_repeat_instance) | redcap_repeat_instance == 1) & 
                              retro_or_pros == 2)
    
    # Longitudinal: repeat instance > 1
    score_data <- subset(data, 
                         redcap_repeat_instance > 1 & 
                           record_id %in% baseline_data$record_id)
  } else {
    stop("Cannot identify longitudinal data structure. Please check your data.")
  }
} else {
  # Use the filter column we found
  cat("\nUsing filter column:", filter_col, "\n")
  
  # Show distribution using base R to avoid issues
  cat("Filter column distribution:\n")
  print(table(data[[filter_col]], useNA = "ifany"))
  
  # ============================================================
  # SEPARATE DATA - PROSPECTIVE ONLY
  # ============================================================
  
  cat("\n=== Extracting Prospective Data ===\n")
  
  # Baseline: filter == 1 and prospective
  baseline_data <- subset(data, 
                          data[[filter_col]] == 1 & 
                            retro_or_pros == 2)
  
  cat("Prospective baseline records:", nrow(baseline_data), "\n")
  cat("Unique prospective patients:", length(unique(baseline_data$record_id)), "\n")
  
  # Get prospective patient IDs
  prosp_ids <- unique(baseline_data$record_id)
  
  # Longitudinal: filter < 1 for prospective patients
  score_data <- subset(data, 
                       data[[filter_col]] < 1 & 
                         record_id %in% prosp_ids)
  
  cat("\nProspective score records:", nrow(score_data), "\n")
  cat("Unique patients with scores:", length(unique(score_data$record_id)), "\n")
}

# ============================================================
# PROCESS BASELINE DATA
# ============================================================

cat("\n=== Processing Baseline Data ===\n")

baseline_processed <- baseline_data %>%
  select(
    record_id,
    date_of_transplant,
    cre_infection,
    cre_infection_date,
    death,
    death_date,
    discharge_date,
    any_of(c("cre_colonization___1", "cre_colonization___2", "cre_colonization___3"))
  ) %>%
  mutate(
    date_of_transplant = as.Date(date_of_transplant),
    cre_infection_date = as.Date(cre_infection_date),
    death_date = as.Date(death_date),
    discharge_date = as.Date(discharge_date),
    
    # Calculate outcomes
    days_to_cre = as.numeric(cre_infection_date - date_of_transplant),
    days_to_death = as.numeric(death_date - date_of_transplant),
    days_to_discharge = as.numeric(discharge_date - date_of_transplant),
    
    # Primary outcome
    primary_outcome = case_when(
      cre_infection == 1 & !is.na(days_to_cre) ~ "CRE Infection",
      death == 1 & !is.na(days_to_death) ~ "Death",
      !is.na(days_to_discharge) ~ "Discharge",
      TRUE ~ "Censored"
    )
  )

cat("Baseline processing complete\n")
cat("Outcome distribution:\n")
print(table(baseline_processed$primary_outcome, useNA = "ifany"))

# ============================================================
# PROCESS SCORE DATA
# ============================================================

cat("\n=== Processing Score Data ===\n")

if(nrow(score_data) > 0) {
  score_processed <- score_data %>%
    select(record_id, score, score_date, any_of(c("events", "score_impact___1", "score_impact___2", "score_impact___3"))) %>%
    filter(!is.na(score) & !is.na(score_date)) %>%
    mutate(
      score_date = as.Date(score_date)
    )
  
  cat("Valid score records:", nrow(score_processed), "\n")
  
  if(nrow(score_processed) > 0) {
    # ============================================================
    # JOIN DATA
    # ============================================================
    
    cat("\n=== Joining Data ===\n")
    
    longitudinal_data <- score_processed %>%
      inner_join(baseline_processed, by = "record_id") %>%
      mutate(
        days_from_olt = as.numeric(score_date - date_of_transplant)
      ) %>%
      filter(days_from_olt >= 0 & days_from_olt <= 180)
    
    cat("Final longitudinal records:", nrow(longitudinal_data), "\n")
    cat("Unique patients:", n_distinct(longitudinal_data$record_id), "\n")
    
    # ============================================================
    # TRAJECTORY ANALYSIS
    # ============================================================
    
    cat("\n=== Trajectory Analysis ===\n")
    
    patient_summary <- longitudinal_data %>%
      group_by(record_id, primary_outcome) %>%
      summarise(
        n_scores = n(),
        baseline_score = first(score),
        final_score = last(score),
        max_score = max(score),
        mean_score = mean(score),
        trend = final_score - baseline_score,
        .groups = "drop"
      )
    
    cat("\nPatient trajectory summary:\n")
    cat("  Patients analyzed:", nrow(patient_summary), "\n")
    cat("  Mean scores per patient:", round(mean(patient_summary$n_scores), 2), "\n")
    
    # ============================================================
    # VISUALIZATIONS
    # ============================================================
    
    cat("\n=== Creating Visualizations ===\n")
    
    # 1. Trajectories by outcome
    p1 <- ggplot(longitudinal_data, aes(x = days_from_olt, y = score, group = record_id)) +
      geom_line(alpha = 0.3) +
      geom_smooth(aes(group = NULL), method = "loess", se = TRUE, color = "red") +
      facet_wrap(~ primary_outcome) +
      labs(
        title = "Score Trajectories by Outcome - Prospective Cohort",
        x = "Days from Transplant",
        y = "Score"
      ) +
      theme_minimal()
    
    print(p1)
    ggsave("prospective_trajectories.png", p1, width = 10, height = 6, dpi = 300)
    
    # 2. Mean scores by outcome
    outcome_summary <- patient_summary %>%
      group_by(primary_outcome) %>%
      summarise(
        n = n(),
        mean_baseline = mean(baseline_score),
        mean_max = mean(max_score),
        mean_final = mean(final_score),
        .groups = "drop"
      )
    
    cat("\nOutcome comparison:\n")
    print(outcome_summary)
    
    # ============================================================
    # STATISTICAL TESTS
    # ============================================================
    
    cat("\n=== Statistical Analysis ===\n")
    
    # Test if baseline scores differ by outcome
    if(length(unique(patient_summary$primary_outcome)) > 1) {
      kw_test <- kruskal.test(baseline_score ~ primary_outcome, data = patient_summary)
      cat("\nKruskal-Wallis test (baseline score by outcome):\n")
      cat("  p-value:", format.pval(kw_test$p.value), "\n")
      
      if(kw_test$p.value < 0.05) {
        cat("  SIGNIFICANT: Baseline scores differ by outcome\n")
      }
    }
    
    # Save results
    results <- list(
      baseline = baseline_processed,
      longitudinal = longitudinal_data,
      patient_summary = patient_summary,
      outcome_summary = outcome_summary
    )
    
    save(results, file = "prospective_trajectory_results.RData")
    cat("\nResults saved to prospective_trajectory_results.RData\n")
    
  } else {
    cat("\nNo valid score data found after filtering\n")
  }
} else {
  cat("\nNo score data found for prospective patients\n")
  cat("This may indicate:\n")
  cat("  1. Scores haven't been entered yet\n")
  cat("  2. Scores are in a different format\n")
  cat("  3. The filter structure is different than expected\n")
}

cat("\n============================================\n")
cat("Analysis complete\n")
cat("============================================\n")

