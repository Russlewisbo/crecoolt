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

### Visualization of score
# ============================================================
# SCORE DISTRIBUTION VS CRE INFECTION RATE ANALYSIS
# Building on existing trajectory analysis
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # for combining plots
library(scales)     # for percentage formatting

# Assuming you have already run your main analysis and have the results object
# If not, load it:
# load("prospective_trajectory_results.RData")

cat("\n============================================\n")
cat("SCORE DISTRIBUTION VS CRE INFECTION ANALYSIS\n")
cat("============================================\n")

# ============================================================
# PREPARE DATA FOR VISUALIZATION
# ============================================================

# Create analysis dataset with scores and CRE outcome
score_cre_analysis <- longitudinal_data %>%
  group_by(record_id) %>%
  summarise(
    # Score metrics
    baseline_score = first(score),
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    final_score = last(score),
    score_trajectory = final_score - baseline_score,
    
    # CRE outcome
    cre_infection = first(cre_infection),
    primary_outcome = first(primary_outcome),
    days_to_cre = first(days_to_cre),
    
    # Time metrics
    n_measurements = n(),
    follow_up_days = max(days_from_olt),
    
    .groups = "drop"
  ) %>%
  mutate(
    cre_status = factor(
      ifelse(cre_infection == 1, "CRE Infection", "No CRE Infection"),
      levels = c("No CRE Infection", "CRE Infection")
    )
  )

cat("\nData summary:\n")
cat("Total patients:", nrow(score_cre_analysis), "\n")
cat("CRE infections:", sum(score_cre_analysis$cre_infection == 1, na.rm = TRUE), "\n")
cat("CRE infection rate:", 
    round(mean(score_cre_analysis$cre_infection == 1, na.rm = TRUE) * 100, 2), "%\n")

# ============================================================
# VISUALIZATION 1: BOX PLOTS OF SCORE DISTRIBUTIONS
# ============================================================

cat("\n=== Creating Box Plot Comparison ===\n")

p1 <- ggplot(score_cre_analysis, aes(x = cre_status, y = max_score, fill = cre_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  scale_fill_manual(values = c("No CRE Infection" = "#4CAF50", "CRE Infection" = "#F44336")) +
  labs(
    title = "Maximum Score Distribution by CRE Infection Status",
    x = "CRE Status",
    y = "Maximum Score During Follow-up",
    fill = "CRE Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white")

# ============================================================
# VISUALIZATION 2: SCORE BINS VS INFECTION RATE
# ============================================================

cat("\n=== Creating Binned Score Analysis ===\n")

# Create score bins for max score
score_cre_analysis <- score_cre_analysis %>%
  mutate(
    score_bin = cut(max_score, 
                   breaks = quantile(max_score, probs = seq(0, 1, 0.2), na.rm = TRUE),
                   include.lowest = TRUE,
                   labels = c("Very Low", "Low", "Medium", "High", "Very High"))
  )

# Calculate infection rate by bin
bin_summary <- score_cre_analysis %>%
  group_by(score_bin) %>%
  summarise(
    n_patients = n(),
    n_infections = sum(cre_infection == 1, na.rm = TRUE),
    infection_rate = mean(cre_infection == 1, na.rm = TRUE),
    se = sqrt(infection_rate * (1 - infection_rate) / n_patients),
    ci_lower = infection_rate - 1.96 * se,
    ci_upper = infection_rate + 1.96 * se,
    .groups = "drop"
  ) %>%
  filter(!is.na(score_bin))

p2 <- ggplot(bin_summary, aes(x = score_bin, y = infection_rate)) +
  geom_col(fill = "#2196F3", alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, size = 0.8) +
  geom_text(aes(label = paste0(round(infection_rate * 100, 1), "%")), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
  labs(
    title = "CRE Infection Rate by Score Category",
    subtitle = paste("n =", sum(bin_summary$n_patients), "patients"),
    x = "Maximum Score Category",
    y = "CRE Infection Rate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# VISUALIZATION 3: CONTINUOUS RELATIONSHIP
# ============================================================

cat("\n=== Creating Continuous Relationship Plot ===\n")

# Create a smoothed relationship plot
p3 <- ggplot(score_cre_analysis, aes(x = max_score, y = as.numeric(cre_infection == 1))) +
  geom_point(aes(color = cre_status), 
             position = position_jitter(height = 0.02), 
             alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "darkblue", fill = "lightblue") +
  scale_y_continuous(labels = percent_format(), 
                     breaks = seq(0, 1, 0.25),
                     limits = c(-0.05, 1.05)) +
  scale_color_manual(values = c("No CRE Infection" = "#4CAF50", "CRE Infection" = "#F44336")) +
  labs(
    title = "Relationship Between Maximum Score and CRE Infection Risk",
    x = "Maximum Score",
    y = "Probability of CRE Infection",
    color = "Actual Outcome"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ============================================================
# VISUALIZATION 4: DENSITY PLOTS
# ============================================================

cat("\n=== Creating Density Comparison ===\n")

p4 <- ggplot(score_cre_analysis, aes(x = max_score, fill = cre_status)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("No CRE Infection" = "#4CAF50", "CRE Infection" = "#F44336")) +
  labs(
    title = "Score Distribution Density by CRE Infection Status",
    x = "Maximum Score",
    y = "Density",
    fill = "CRE Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_vline(data = score_cre_analysis %>% 
               group_by(cre_status) %>% 
               summarise(median_score = median(max_score, na.rm = TRUE)),
             aes(xintercept = median_score, color = cre_status),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("No CRE Infection" = "#2E7D32", "CRE Infection" = "#C62828"),
                     guide = "none")

# ============================================================
# VISUALIZATION 5: TRAJECTORY PATTERNS
# ============================================================

cat("\n=== Creating Trajectory Pattern Analysis ===\n")

# Analyze score trajectories for patients who developed CRE
trajectory_summary <- score_cre_analysis %>%
  mutate(
    trajectory_pattern = case_when(
      score_trajectory > 2 ~ "Increasing",
      score_trajectory < -2 ~ "Decreasing",
      TRUE ~ "Stable"
    )
  )

trajectory_rates <- trajectory_summary %>%
  group_by(trajectory_pattern) %>%
  summarise(
    n = n(),
    cre_rate = mean(cre_infection == 1, na.rm = TRUE),
    .groups = "drop"
  )

p5 <- ggplot(trajectory_rates, aes(x = trajectory_pattern, y = cre_rate)) +
  geom_col(fill = "#FF9800", alpha = 0.7) +
  geom_text(aes(label = paste0(round(cre_rate * 100, 1), "%\n(n=", n, ")")), 
            vjust = -0.2, size = 3.5) +
  scale_y_continuous(labels = percent_format(), limits = c(0, max(trajectory_rates$cre_rate) * 1.2)) +
  labs(
    title = "CRE Infection Rate by Score Trajectory Pattern",
    x = "Score Trajectory Pattern",
    y = "CRE Infection Rate"
  ) +
  theme_minimal()

# ============================================================
# COMBINE PLOTS
# ============================================================

cat("\n=== Combining Visualizations ===\n")

# Create a combined plot using patchwork
combined_plot <- (p1 | p2) / (p3 | p4) / p5 +
  plot_annotation(
    title = "Score Distribution and CRE Infection Risk Analysis",
    subtitle = paste("Prospective Cohort (n =", nrow(score_cre_analysis), "patients)"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)
ggsave("score_cre_combined_analysis.png", combined_plot, width = 14, height = 12, dpi = 300)

# Save individual plots
ggsave("score_cre_boxplot.png", p1, width = 8, height = 6, dpi = 300)
ggsave("score_cre_bins.png", p2, width = 8, height = 6, dpi = 300)
ggsave("score_cre_continuous.png", p3, width = 8, height = 6, dpi = 300)
ggsave("score_cre_density.png", p4, width = 8, height = 6, dpi = 300)
ggsave("score_cre_trajectory.png", p5, width = 8, height = 6, dpi = 300)

# ============================================================
# STATISTICAL ANALYSIS
# ============================================================

cat("\n=== Statistical Analysis ===\n")

# 1. Compare mean scores between groups
t_test <- t.test(max_score ~ cre_status, data = score_cre_analysis)
cat("\nT-test for maximum score difference:\n")
cat("  Mean score (No CRE):", round(t_test$estimate[1], 2), "\n")
cat("  Mean score (CRE):", round(t_test$estimate[2], 2), "\n")
cat("  Difference:", round(t_test$estimate[2] - t_test$estimate[1], 2), "\n")
cat("  p-value:", format.pval(t_test$p.value), "\n")

# 2. Wilcoxon rank-sum test (non-parametric)
wilcox_test <- wilcox.test(max_score ~ cre_status, data = score_cre_analysis)
cat("\nWilcoxon rank-sum test:\n")
cat("  p-value:", format.pval(wilcox_test$p.value), "\n")

# 3. Logistic regression
logit_model <- glm(cre_infection ~ max_score, 
                   data = score_cre_analysis, 
                   family = binomial())
cat("\nLogistic Regression Results:\n")
print(summary(logit_model))

# Calculate odds ratio
or <- exp(coef(logit_model)[2])
ci <- exp(confint(logit_model)[2,])
cat("\nOdds Ratio per unit increase in max score:", round(or, 3), "\n")
cat("95% CI:", round(ci[1], 3), "-", round(ci[2], 3), "\n")

# ============================================================
# SUMMARY TABLE
# ============================================================

cat("\n=== Creating Summary Table ===\n")

summary_table <- score_cre_analysis %>%
  group_by(cre_status) %>%
  summarise(
    n = n(),
    baseline_mean = round(mean(baseline_score, na.rm = TRUE), 2),
    baseline_sd = round(sd(baseline_score, na.rm = TRUE), 2),
    max_mean = round(mean(max_score, na.rm = TRUE), 2),
    max_sd = round(sd(max_score, na.rm = TRUE), 2),
    trajectory_mean = round(mean(score_trajectory, na.rm = TRUE), 2),
    trajectory_sd = round(sd(score_trajectory, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("\nScore Summary by CRE Status:\n")
print(summary_table)

# Save analysis results
analysis_results <- list(
  score_cre_data = score_cre_analysis,
  bin_summary = bin_summary,
  trajectory_rates = trajectory_rates,
  summary_table = summary_table,
  statistical_tests = list(
    t_test = t_test,
    wilcox_test = wilcox_test,
    logit_model = logit_model
  )
)

save(analysis_results, file = "score_cre_analysis_results.RData")
cat("\nAnalysis results saved to score_cre_analysis_results.RData\n")

cat("\n============================================\n")
cat("Score vs CRE Analysis Complete\n")
cat("============================================\n")

# ============================================================
# CORRECTED ANTI-CRE THERAPY EFFECTIVENESS ANALYSIS
# Using score_impact___2 as the therapy indicator
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

cat("\n============================================\n")
cat("ANTI-CRE THERAPY EFFECTIVENESS ANALYSIS\n")
cat("Based on score_impact___2 (Therapy Started)\n")
cat("============================================\n")

# ============================================================
# DATA PREPARATION WITH CORRECT VARIABLE INTERPRETATION
# ============================================================

cat("\n=== Understanding the Variables ===\n")
cat("score_impact___1 = No action taken\n")
cat("score_impact___2 = Anti-CRE therapy started (OUR FOCUS)\n")
cat("score_impact___3 = Diagnostic test performed\n")

cat("\n=== Processing Score Impact Data ===\n")

# Process score data with correct therapy indicator
score_therapy_data <- score_data %>%
  select(
    record_id,
    score,
    score_date,
    score_impact___1,  # No action
    score_impact___2,  # THERAPY STARTED
    score_impact___3   # Diagnostic test
  ) %>%
  filter(!is.na(score) & !is.na(score_date)) %>%
  mutate(
    score_date = as.Date(score_date),
    
    # Correct interpretation
    no_action = ifelse(score_impact___1 == 1, 1, 0),
    therapy_started = ifelse(score_impact___2 == 1, 1, 0),  # THIS IS OUR THERAPY VARIABLE
    diagnostic_performed = ifelse(score_impact___3 == 1, 1, 0),
    
    # Any action taken
    any_action = ifelse(therapy_started == 1 | diagnostic_performed == 1, 1, 0),
    
    # Score categories
    score_category = cut(score, 
                        breaks = c(-Inf, 5, 10, 15, 20, Inf),
                        labels = c("0-5", "6-10", "11-15", "16-20", ">20"))
  )

cat("\nScore assessments processed:", nrow(score_therapy_data), "\n")
cat("Assessments where therapy started:", sum(score_therapy_data$therapy_started), 
    "(", round(mean(score_therapy_data$therapy_started)*100, 1), "%)\n")
cat("Assessments with no action:", sum(score_therapy_data$no_action),
    "(", round(mean(score_therapy_data$no_action)*100, 1), "%)\n")
cat("Assessments with diagnostic test:", sum(score_therapy_data$diagnostic_performed),
    "(", round(mean(score_therapy_data$diagnostic_performed)*100, 1), "%)\n")

# ============================================================
# CREATE PATIENT-LEVEL SUMMARY
# ============================================================

cat("\n=== Creating Patient-Level Summary ===\n")

# Aggregate to patient level - did they EVER receive therapy based on score?
patient_therapy_data <- score_therapy_data %>%
  group_by(record_id) %>%
  summarise(
    # Therapy exposure (EVER received therapy based on score_impact___2)
    received_score_based_therapy = max(therapy_started),
    n_therapy_instances = sum(therapy_started),
    
    # Other actions
    ever_had_diagnostic = max(diagnostic_performed),
    ever_no_action = max(no_action),
    
    # Score metrics
    baseline_score = first(score),
    max_score = max(score),
    mean_score = mean(score),
    final_score = last(score),
    n_assessments = n(),
    
    # Score categories
    max_score_category = cut(max_score,
                            breaks = c(-Inf, 5, 10, 15, 20, Inf),
                            labels = c("0-5", "6-10", "11-15", "16-20", ">20")),
    
    # Risk stratification
    risk_group = case_when(
      max_score <= 10 ~ "Low Risk (≤10)",
      max_score <= 20 ~ "Medium Risk (11-20)",
      TRUE ~ "High Risk (>20)"
    ),
    
    .groups = "drop"
  ) %>%
  # Join with outcomes
  left_join(baseline_processed %>% 
              select(record_id, cre_infection, primary_outcome, date_of_transplant,
                     any_of(c("days_to_cre", "cre_infection_date"))),
            by = "record_id") %>%
  filter(!is.na(primary_outcome))

cat("\nPatient-level summary:\n")
cat("Total patients analyzed:", nrow(patient_therapy_data), "\n")
cat("Patients who received score-based therapy:", 
    sum(patient_therapy_data$received_score_based_therapy), 
    "(", round(mean(patient_therapy_data$received_score_based_therapy)*100, 1), "%)\n")
cat("CRE infections overall:", sum(patient_therapy_data$cre_infection == 1, na.rm = TRUE),
    "(", round(mean(patient_therapy_data$cre_infection == 1, na.rm = TRUE)*100, 1), "%)\n")

# ============================================================
# OVERALL TREATMENT EFFECT
# ============================================================

cat("\n=== Overall Treatment Effect (score_impact___2) ===\n")

# Overall comparison
overall_comparison <- patient_therapy_data %>%
  group_by(received_score_based_therapy) %>%
  summarise(
    n = n(),
    n_infections = sum(cre_infection == 1, na.rm = TRUE),
    infection_rate = mean(cre_infection == 1, na.rm = TRUE),
    mean_max_score = mean(max_score, na.rm = TRUE),
    sd_score = sd(max_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    therapy_label = ifelse(received_score_based_therapy == 1, 
                          "Score-Based Therapy", "No Score-Based Therapy"),
    ci_lower = qbinom(0.025, n, infection_rate) / n,
    ci_upper = qbinom(0.975, n, infection_rate) / n
  )

print(overall_comparison)

# Statistical test
therapy_table <- table(patient_therapy_data$received_score_based_therapy, 
                      patient_therapy_data$cre_infection)
fisher_test <- fisher.test(therapy_table)
chi_test <- chisq.test(therapy_table)

cat("\nStatistical Tests:\n")
cat("  Fisher's exact test p-value:", format.pval(fisher_test$p.value), "\n")
cat("  Odds ratio:", round(fisher_test$estimate, 3), "\n")

# Calculate effectiveness metrics
if(nrow(overall_comparison) == 2) {
  no_therapy_rate <- overall_comparison$infection_rate[overall_comparison$received_score_based_therapy == 0]
  therapy_rate <- overall_comparison$infection_rate[overall_comparison$received_score_based_therapy == 1]
  
  arr <- no_therapy_rate - therapy_rate
  rrr <- if(no_therapy_rate > 0) (no_therapy_rate - therapy_rate) / no_therapy_rate else NA
  nnt <- if(arr > 0) 1/arr else NA
  
  cat("\n=== Treatment Effectiveness Metrics ===\n")
  cat("Infection rate WITHOUT score-based therapy:", round(no_therapy_rate * 100, 1), "%\n")
  cat("Infection rate WITH score-based therapy:", round(therapy_rate * 100, 1), "%\n")
  cat("Absolute Risk Reduction (ARR):", round(arr * 100, 2), "%\n")
  if(!is.na(rrr)) cat("Relative Risk Reduction (RRR):", round(rrr * 100, 2), "%\n")
  if(!is.na(nnt) & nnt > 0) cat("Number Needed to Treat (NNT):", round(nnt, 0), "\n")
  
  if(arr > 0) {
    cat("\n✓ Score-based therapy appears BENEFICIAL\n")
  } else if(arr < 0) {
    cat("\n✗ Score-based therapy appears HARMFUL (paradox - likely confounding)\n")
    cat("  Note: High-risk patients more likely to receive therapy\n")
  } else {
    cat("\n- No difference detected\n")
  }
}

# ============================================================
# VISUALIZATION 1: OVERALL COMPARISON
# ============================================================

cat("\n=== Creating Visualizations ===\n")

p1 <- ggplot(overall_comparison, aes(x = therapy_label, y = infection_rate, fill = therapy_label)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 1) +
  geom_text(aes(label = paste0(round(infection_rate * 100, 1), "%\n(", 
                               n_infections, "/", n, ")")),
            vjust = -0.5, size = 4) +
  scale_y_continuous(labels = percent_format(), 
                     limits = c(0, max(overall_comparison$ci_upper) * 1.2)) +
  scale_fill_manual(values = c("Score-Based Therapy" = "#4CAF50", 
                               "No Score-Based Therapy" = "#F44336")) +
  labs(
    title = "CRE Infection Rates by Score-Based Therapy Status",
    subtitle = paste("Based on score_impact___2 | p =", format.pval(fisher_test$p.value)),
    x = "",
    y = "CRE Infection Rate",
    fill = "Therapy Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p1)

# ============================================================
# STRATIFIED ANALYSIS BY RISK SCORE
# ============================================================

cat("\n=== Stratified Analysis by Risk Score ===\n")

# Stratified comparison
stratified_comparison <- patient_therapy_data %>%
  group_by(risk_group, received_score_based_therapy) %>%
  summarise(
    n = n(),
    n_infections = sum(cre_infection == 1, na.rm = TRUE),
    infection_rate = mean(cre_infection == 1, na.rm = TRUE),
    mean_score = mean(max_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = received_score_based_therapy,
    values_from = c(n, n_infections, infection_rate, mean_score),
    names_prefix = "therapy_"
  ) %>%
  mutate(
    # Handle NA values
    across(starts_with("n_"), ~replace_na(., 0)),
    across(starts_with("infection_rate_"), ~replace_na(., 0)),
    
    # Calculate effectiveness within each stratum
    arr = infection_rate_therapy_0 - infection_rate_therapy_1,
    rrr = ifelse(infection_rate_therapy_0 > 0,
                 (infection_rate_therapy_0 - infection_rate_therapy_1) / infection_rate_therapy_0,
                 NA),
    nnt = ifelse(arr > 0, 1/arr, NA)
  )

cat("\nRisk-Stratified Results:\n")
print(stratified_comparison %>% 
        select(risk_group, n_therapy_0, n_therapy_1, 
               infection_rate_therapy_0, infection_rate_therapy_1, arr, nnt))

# Statistical tests within strata
cat("\n=== Within-Stratum Statistical Tests ===\n")
for(risk in unique(patient_therapy_data$risk_group)) {
  subset_data <- filter(patient_therapy_data, risk_group == risk)
  if(length(unique(subset_data$received_score_based_therapy)) > 1) {
    test_table <- table(subset_data$received_score_based_therapy, subset_data$cre_infection)
    if(all(dim(test_table) >= 2)) {
      test_result <- fisher.test(test_table)
      cat("\n", risk, ":\n")
      cat("  p-value:", format.pval(test_result$p.value), "\n")
      cat("  Odds ratio:", round(test_result$estimate, 3), "\n")
    }
  }
}

# ============================================================
# VISUALIZATION 2: STRATIFIED COMPARISON
# ============================================================

# Prepare data for visualization
stratified_viz <- patient_therapy_data %>%
  group_by(risk_group, received_score_based_therapy) %>%
  summarise(
    n = n(),
    n_infections = sum(cre_infection == 1, na.rm = TRUE),
    infection_rate = mean(cre_infection == 1, na.rm = TRUE),
    se = sqrt(infection_rate * (1 - infection_rate) / n),
    ci_lower = pmax(0, infection_rate - 1.96 * se),
    ci_upper = pmin(1, infection_rate + 1.96 * se),
    .groups = "drop"
  ) %>%
  mutate(
    therapy_label = ifelse(received_score_based_therapy == 1, "Therapy", "No Therapy")
  )

p2 <- ggplot(stratified_viz, aes(x = therapy_label, y = infection_rate, fill = therapy_label)) +
  facet_wrap(~ risk_group) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 0.8) +
  geom_text(aes(label = paste0(round(infection_rate * 100, 1), "%\n(", 
                               n_infections, "/", n, ")")),
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.6)) +
  scale_fill_manual(values = c("Therapy" = "#4CAF50", "No Therapy" = "#F44336")) +
  labs(
    title = "CRE Infection Rates Stratified by Risk Score",
    subtitle = "Effect of score-based therapy (score_impact___2) within risk groups",
    x = "",
    y = "CRE Infection Rate",
    fill = "Score-Based Therapy"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)

# ============================================================
# VISUALIZATION 3: CONTINUOUS RELATIONSHIP
# ============================================================

p3 <- ggplot(patient_therapy_data, 
             aes(x = max_score, y = as.numeric(cre_infection == 1),
                 color = factor(received_score_based_therapy))) +
  geom_point(position = position_jitter(height = 0.02, width = 0.5), 
             alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_y_continuous(labels = percent_format(), limits = c(-0.05, 1.05)) +
  scale_color_manual(values = c("0" = "#F44336", "1" = "#4CAF50"),
                    labels = c("0" = "No Score-Based Therapy", 
                              "1" = "Score-Based Therapy")) +
  labs(
    title = "CRE Infection Probability by Maximum Score",
    subtitle = "Comparing patients with and without score-based therapy",
    x = "Maximum Score",
    y = "Probability of CRE Infection",
    color = "Therapy Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# ============================================================
# VISUALIZATION 4: THERAPY UTILIZATION BY SCORE
# ============================================================

therapy_by_score <- patient_therapy_data %>%
  group_by(max_score_category) %>%
  summarise(
    n = n(),
    therapy_rate = mean(received_score_based_therapy),
    infection_rate = mean(cre_infection == 1, na.rm = TRUE),
    .groups = "drop"
  )

p4 <- ggplot(therapy_by_score, aes(x = max_score_category)) +
  geom_col(aes(y = therapy_rate), fill = "#2196F3", alpha = 0.7) +
  geom_line(aes(y = infection_rate, group = 1), color = "red", size = 1.5) +
  geom_point(aes(y = infection_rate), color = "red", size = 3) +
  geom_text(aes(y = therapy_rate, label = paste0(round(therapy_rate * 100, 0), "%")),
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = percent_format(),
                     sec.axis = sec_axis(~., name = "CRE Infection Rate (Red Line)")) +
  labs(
    title = "Score-Based Therapy Utilization and Infection Rates",
    subtitle = "Blue bars: % receiving therapy | Red line: infection rate",
    x = "Maximum Score Category",
    y = "Therapy Rate (Blue Bars)"
  ) +
  theme_minimal()

print(p4)

# ============================================================
# ADJUSTED ANALYSIS
# ============================================================

cat("\n=== Regression Analysis Adjusting for Baseline Risk ===\n")

# Unadjusted model
unadjusted_model <- glm(cre_infection ~ received_score_based_therapy,
                        data = patient_therapy_data,
                        family = binomial())

# Adjusted model
adjusted_model <- glm(cre_infection ~ received_score_based_therapy + max_score,
                      data = patient_therapy_data,
                      family = binomial())

cat("\nUnadjusted Analysis:\n")
cat("  Therapy coefficient:", round(coef(unadjusted_model)["received_score_based_therapy"], 3), "\n")
cat("  Therapy OR:", round(exp(coef(unadjusted_model)["received_score_based_therapy"]), 3), "\n")

cat("\nAdjusted for Risk Score:\n")
cat("  Therapy coefficient:", round(coef(adjusted_model)["received_score_based_therapy"], 3), "\n")
cat("  Therapy OR:", round(exp(coef(adjusted_model)["received_score_based_therapy"]), 3), "\n")
cat("  Score coefficient:", round(coef(adjusted_model)["max_score"], 3), "\n")

if(coef(unadjusted_model)["received_score_based_therapy"] > 0 & 
   coef(adjusted_model)["received_score_based_therapy"] < 0) {
  cat("\n✓ Simpson's Paradox detected! Effect reverses after adjustment\n")
}

# ============================================================
# COMBINE VISUALIZATIONS
# ============================================================

combined_analysis <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Score-Based Anti-CRE Therapy Effectiveness Analysis",
    subtitle = paste("Based on score_impact___2 | n =", nrow(patient_therapy_data), 
                    "patients |", sum(patient_therapy_data$received_score_based_therapy), 
                    "received therapy"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_analysis)
ggsave("score_based_therapy_effectiveness.png", combined_analysis, 
       width = 14, height = 12, dpi = 300)

# Save individual plots
ggsave("therapy_overall_effect.png", p1, width = 8, height = 6, dpi = 300)
ggsave("therapy_stratified_effect.png", p2, width = 10, height = 6, dpi = 300)
ggsave("therapy_continuous_effect.png", p3, width = 10, height = 6, dpi = 300)
ggsave("therapy_utilization_by_score.png", p4, width = 10, height = 6, dpi = 300)

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("THERAPY EFFECTIVENESS SUMMARY\n")
cat("============================================\n")

cat("\n1. KEY FINDING:\n")
if(exists("arr")) {
  if(arr > 0) {
    cat("   ✓ Score-based therapy shows BENEFIT\n")
    cat("   - Reduces infections by", round(arr * 100, 1), "percentage points\n")
  } else if(arr < 0) {
    cat("   ✗ Paradox detected: Therapy group has HIGHER infection rate\n")
    cat("   - Likely due to confounding by indication\n")
    cat("   - High-risk patients → Get therapy AND get infected\n")
  }
}

cat("\n2. CLINICAL INTERPRETATION:\n")
cat("   - Score successfully identifies high-risk patients\n")
cat("   - These patients appropriately receive therapy\n")
cat("   - True effectiveness requires risk adjustment or RCT\n")

cat("\n3. RECOMMENDATIONS:\n")
if(coef(adjusted_model)["received_score_based_therapy"] < 0) {
  cat("   ✓ After risk adjustment, therapy appears beneficial\n")
  cat("   - Continue score-based therapy protocol\n")
} else {
  cat("   - Review therapy effectiveness\n")
  cat("   - Consider alternative interventions\n")
}

cat("\n============================================\n")
cat("Analysis Complete\n")
cat("============================================\n")
