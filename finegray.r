#############################################
## CORRECTED ADJUSTED ANALYSIS FOR CRE INFECTION
## Accounting for Confounding by Indication
#############################################

library(dplyr)
library(cmprsk)
library(riskRegression)
library(prodlim)
library(ggplot2)
library(survival)

# ============================================================
# STEP 1: Examine baseline differences between groups
# ============================================================

print("=========== BASELINE CHARACTERISTICS BY MANAGEMENT ===========")

# Create full dataset with all covariates
data_full <- data_analysis %>%
  mutate(
    # Create binary indicators for all risk factors
    post_olt_compli1 = ifelse(post_olt_compli___1 == 1, 1, 0),
    post_olt_compli3 = ifelse(post_olt_compli___3 == 1, 1, 0),
    post_olt_compli5 = ifelse(post_olt_compli___5 == 1, 1, 0),
    multisite_col = ifelse(multisite_colonization == 1, 1, 0),
    cre_col1 = ifelse(cre_colonization___1 == 1, 1, 0),
    cre_col2 = ifelse(cre_colonization___2 == 1, 1, 0),
    mec1 = ifelse(mec_of_carbapenem_resi___1 == 1, 1, 0),
    Management = ifelse(management == 1, "Pre-emptive", "Standard care")
  )

# Compare baseline characteristics
baseline_comparison <- data_full %>%
  group_by(Management) %>%
  summarise(
    n = n(),
    n_cre_infection = sum(outcome_comp == 1),
    infection_rate = mean(outcome_comp == 1) * 100,
    
    # Risk factors
    pct_ARF = mean(post_olt_compli1, na.rm = TRUE) * 100,
    pct_ventilation = mean(post_olt_compli3, na.rm = TRUE) * 100,
    pct_reintervention = mean(post_olt_compli5, na.rm = TRUE) * 100,
    pct_multisite = mean(multisite_col, na.rm = TRUE) * 100,
    pct_cre_pre = mean(cre_col1, na.rm = TRUE) * 100,
    pct_cre_post = mean(cre_col2, na.rm = TRUE) * 100,
    pct_KPC = mean(mec1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

print("Baseline characteristics showing selection bias:")
print(baseline_comparison)

# ============================================================
# STEP 2: Fit ADJUSTED Fine-Gray Model
# ============================================================

print("\n=========== ADJUSTED FINE-GRAY MODEL ===========")

# Prepare data for modeling (remove NAs)
model_data <- data_full %>%
  select(time_c, outcome_comp, management,
         post_olt_compli1, post_olt_compli3, post_olt_compli5,
         multisite_col, cre_col1, cre_col2, mec1) %>%
  na.omit()

print(paste("Patients in adjusted analysis:", nrow(model_data)))

# Try fitting the adjusted model with error handling
tryCatch({
  # Fit adjusted Fine-Gray model
  fg_adjusted <- FGR(
    Hist(time_c, outcome_comp) ~ 
      management +  # Treatment effect
      post_olt_compli1 + post_olt_compli3 + post_olt_compli5 +  # Complications
      multisite_col + cre_col1 + cre_col2 + mec1,  # Colonization factors
    data = model_data,
    cause = 1  # CRE infection
  )
  
  # Extract results
  coef_adj <- fg_adjusted$crrFit$coef
  se_adj <- sqrt(diag(fg_adjusted$crrFit$var))
  hr_adj <- exp(coef_adj)
  ci_lower <- exp(coef_adj - 1.96*se_adj)
  ci_upper <- exp(coef_adj + 1.96*se_adj)
  p_values <- 2 * pnorm(-abs(coef_adj/se_adj))
  
  # Create results table
  results_adjusted <- data.frame(
    Variable = names(coef_adj),
    HR = round(hr_adj, 2),
    CI_95 = paste0("[", round(ci_lower, 2), "-", round(ci_upper, 2), "]"),
    p_value = round(p_values, 4)
  )
  
  print("\nADJUSTED Hazard Ratios:")
  print(results_adjusted)
  
  # Highlight treatment effect
  mgmt_idx <- which(names(coef_adj) == "management")
  if(length(mgmt_idx) > 0) {
    print("\n========== ADJUSTED TREATMENT EFFECT ==========")
    print(results_adjusted[mgmt_idx, ])
    
    if(hr_adj[mgmt_idx] < 1) {
      print(paste0("After adjustment: Pre-emptive therapy reduces CRE infection risk by ",
                   round((1 - hr_adj[mgmt_idx]) * 100, 1), "%"))
    } else {
      print(paste0("After adjustment: Pre-emptive therapy increases CRE infection risk by ",
                   round((hr_adj[mgmt_idx] - 1) * 100, 1), "%"))
    }
  }
  
}, error = function(e) {
  print(paste("Full model failed:", e$message))
  print("Trying reduced model...")
})

# ============================================================
# STEP 3: Create ADJUSTED Cumulative Incidence Curves
# ============================================================

print("\n=========== CREATING ADJUSTED CIF CURVES ===========")

# Method 1: Use predictions for all patients and average by group
times_seq <- seq(0, 180, by = 5)

# Get predictions for all patients
pred_matrix <- predict(fg_adjusted, 
                       newdata = model_data,
                       times = times_seq, 
                       cause = 1, 
                       se = FALSE)

# Calculate average CIF by management group
adjusted_cif_data <- model_data %>%
  mutate(patient_id = row_number()) %>%
  bind_cols(as.data.frame(pred_matrix)) %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "time_idx", 
               values_to = "cif") %>%
  mutate(
    time_idx = as.numeric(gsub("V", "", time_idx)),
    time = times_seq[time_idx],
    Management = ifelse(management == 1, "Pre-emptive therapy", "Standard care")
  ) %>%
  group_by(time, Management) %>%
  summarise(
    mean_cif = mean(cif, na.rm = TRUE),
    se_cif = sd(cif, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower = pmax(0, mean_cif - 1.96 * se_cif),
    ci_upper = pmin(1, mean_cif + 1.96 * se_cif)
  )

# ADJUSTED cumulative incidence plot
p_adjusted <- ggplot(adjusted_cif_data, 
                     aes(x = time, y = mean_cif, color = Management, fill = Management)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Standard care" = "#D81B60", 
                                "Pre-emptive therapy" = "#1E88E5")) +
  scale_fill_manual(values = c("Standard care" = "#D81B60", 
                               "Pre-emptive therapy" = "#1E88E5")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
  labs(
    title = "ADJUSTED Cumulative Incidence of CRE Infection",
    subtitle = "Adjusted for baseline risk factors",
    x = "Days since liver transplantation",
    y = "Adjusted cumulative incidence",
    color = "Management",
    fill = "Management"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )

print(p_adjusted)

# ============================================================
# STEP 4: Side-by-side comparison - Unadjusted vs Adjusted
# ============================================================

print("\n=========== UNADJUSTED VS ADJUSTED COMPARISON ===========")

# Get unadjusted CIF using prodlim
unadj_model <- prodlim(Hist(time_c, outcome_comp) ~ management, data = model_data)

unadj_cif_0 <- predict(unadj_model, 
                       newdata = data.frame(management = 0),
                       times = times_seq, cause = 1)

unadj_cif_1 <- predict(unadj_model, 
                       newdata = data.frame(management = 1),
                       times = times_seq, cause = 1)

# Combine unadjusted and adjusted data
comparison_data <- rbind(
  # Unadjusted
  data.frame(
    time = rep(times_seq, 2),
    cif = c(unadj_cif_0, unadj_cif_1),
    Management = rep(c("Standard care", "Pre-emptive therapy"), each = length(times_seq)),
    Analysis = "Unadjusted"
  ),
  # Adjusted
  adjusted_cif_data %>%
    select(time, cif = mean_cif, Management) %>%
    mutate(Analysis = "Adjusted")
)

# Faceted plot
p_comparison <- ggplot(comparison_data, 
                       aes(x = time, y = cif, color = Management)) +
  geom_line(size = 1.2) +
  facet_wrap(~Analysis, ncol = 2) +
  scale_color_manual(values = c("Standard care" = "#D81B60", 
                                "Pre-emptive therapy" = "#1E88E5")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
  labs(
    title = "Unadjusted vs Adjusted Analysis",
    subtitle = "Impact of adjusting for baseline risk factors",
    x = "Days since liver transplantation",
    y = "Cumulative incidence of CRE infection"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )

print(p_comparison)

# ============================================================
# STEP 5: Calculate risk at key time points
# ============================================================

key_times <- c(30, 60, 90, 180)

risk_summary <- adjusted_cif_data %>%
  filter(time %in% key_times) %>%
  select(time, Management, mean_cif) %>%
  pivot_wider(names_from = Management, values_from = mean_cif) %>%
  mutate(
    risk_difference = `Standard care` - `Pre-emptive therapy`,
    relative_reduction = risk_difference / `Standard care` * 100
  )

print("\n=========== ADJUSTED RISK AT KEY TIME POINTS ===========")
print(risk_summary)

# ============================================================
# STEP 6: Simple alternative - stratified analysis
# ============================================================

print("\n=========== STRATIFIED ANALYSIS BY RISK LEVEL ===========")

# Create risk score based on baseline factors
model_data <- model_data %>%
  mutate(
    risk_score = post_olt_compli1 + post_olt_compli3 + post_olt_compli5 +
      multisite_col + cre_col1 + cre_col2 + mec1,
    risk_category = case_when(
      risk_score <= 1 ~ "Low risk",
      risk_score <= 3 ~ "Medium risk",
      TRUE ~ "High risk"
    )
  )

# Stratified results
stratified_results <- model_data %>%
  group_by(risk_category, management) %>%
  summarise(
    n = n(),
    n_cre = sum(outcome_comp == 1),
    cre_rate = n_cre / n * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Management = ifelse(management == 1, "Pre-emptive", "Standard")
  ) %>%
  pivot_wider(names_from = Management, 
              values_from = c(n, cre_rate),
              values_fill = list(n = 0, cre_rate = 0))

print("CRE infection rates stratified by baseline risk:")
print(stratified_results)

# Visualization of stratified results
stratified_plot_data <- model_data %>%
  mutate(
    Management = ifelse(management == 1, "Pre-emptive therapy", "Standard care")
  ) %>%
  group_by(risk_category, Management) %>%
  summarise(
    n = n(),
    n_cre = sum(outcome_comp == 1),
    cre_rate = n_cre / n,
    se = sqrt(cre_rate * (1 - cre_rate) / n),
    .groups = "drop"
  )

p_stratified <- ggplot(stratified_plot_data, 
                       aes(x = risk_category, y = cre_rate, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = cre_rate - 1.96*se, ymax = cre_rate + 1.96*se),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("Standard care" = "#D81B60", 
                               "Pre-emptive therapy" = "#1E88E5")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "CRE Infection Rates Stratified by Baseline Risk",
    x = "Baseline risk category",
    y = "CRE infection rate",
    fill = "Management"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )

print(p_stratified)

# Save plots
ggsave("adjusted_cif_plot.png", p_adjusted, width = 10, height = 6, dpi = 300)
ggsave("comparison_plot.png", p_comparison, width = 12, height = 6, dpi = 300)
ggsave("stratified_analysis.png", p_stratified, width = 10, height = 6, dpi = 300)

print("\n=========== ADJUSTED ANALYSIS COMPLETE ===========")

