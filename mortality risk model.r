#############################################
## FIXED RISK ANALYSIS WITH NA HANDLING
## Properly handles missing values throughout
#############################################

library(dplyr)
library(ggplot2)

cat("============================================\n")
cat("SIMPLE RISK ANALYSIS WITH NA HANDLING\n")
cat("============================================\n")

# ============================================================
# STEP 1: LOAD DATA
# ============================================================

cat("\n=== Loading Data ===\n")

data_full <- read.csv("CRECOOLT_overall.csv", stringsAsFactors = FALSE)
cat("Total records:", nrow(data_full), "\n")

# Basic subset
data_work <- data_full %>%
  filter(!is.na(retro_or_pros))

cat("After filtering:", nrow(data_work), "\n")

# ============================================================
# STEP 2: CREATE OUTCOMES
# ============================================================

cat("\n=== Creating Outcomes ===\n")

# Convert dates
data_work$date_of_transplant <- as.Date(data_work$date_of_transplant)
data_work$death_date <- as.Date(data_work$death_date)
data_work$discharge_date <- as.Date(data_work$discharge_date)

# Create death outcome
data_work$days_to_death <- as.numeric(data_work$death_date - data_work$date_of_transplant)
data_work$death_180 <- ifelse(!is.na(data_work$days_to_death) & 
                                data_work$days_to_death <= 180 & 
                                data_work$days_to_death > 0, 1, 0)
data_work$death_180[is.na(data_work$death_180)] <- 0

# Create time variable
data_work$time_simple <- 180
data_work$time_simple[data_work$death_180 == 1] <- data_work$days_to_death[data_work$death_180 == 1]

# Filter valid cases
data_work <- data_work %>%
  filter(!is.na(time_simple) & time_simple > 0 & time_simple <= 180)

cat("Valid patients:", nrow(data_work), "\n")
cat("Deaths within 180 days:", sum(data_work$death_180), "\n")
cat("Crude mortality:", round(mean(data_work$death_180) * 100, 1), "%\n")

# ============================================================
# STEP 3: PREPARE PREDICTORS
# ============================================================

cat("\n=== Preparing Predictors ===\n")

# Initialize all predictor columns
data_work$age_clean <- NA
data_work$meld_clean <- NA
data_work$cci_clean <- NA

# Age
if("age" %in% names(data_work)) {
  data_work$age_clean <- suppressWarnings(as.numeric(as.character(data_work$age)))
  cat("Age: ", sum(!is.na(data_work$age_clean)), " valid values, range ", 
      round(min(data_work$age_clean, na.rm=T)), "-", 
      round(max(data_work$age_clean, na.rm=T)), "\n", sep="")
}

# MELD
if("meld_score" %in% names(data_work)) {
  data_work$meld_clean <- suppressWarnings(as.numeric(as.character(data_work$meld_score)))
  cat("MELD: ", sum(!is.na(data_work$meld_clean)), " valid values, range ", 
      round(min(data_work$meld_clean, na.rm=T)), "-", 
      round(max(data_work$meld_clean, na.rm=T)), "\n", sep="")
}

# CCI
if("cci" %in% names(data_work)) {
  data_work$cci_clean <- suppressWarnings(as.numeric(as.character(data_work$cci)))
  cat("CCI: ", sum(!is.na(data_work$cci_clean)), " valid values\n", sep="")
}

# Binary complications
cat("\nComplications:\n")
for(i in 1:11) {
  var_name <- paste0("post_olt_compli___", i)
  if(var_name %in% names(data_work)) {
    data_work[[paste0("compli_", i)]] <- ifelse(data_work[[var_name]] == 1, 1, 0)
    data_work[[paste0("compli_", i)]][is.na(data_work[[paste0("compli_", i)]])] <- 0
    n_cases <- sum(data_work[[paste0("compli_", i)]])
    if(n_cases > 0) {
      cat("  Complication ", i, ": ", n_cases, " cases\n", sep="")
    }
  }
}

# ============================================================
# STEP 4: CREATE SIMPLE RISK SCORE
# ============================================================

cat("\n=== Creating Risk Score ===\n")

# Initialize risk score for all patients
data_work$risk_points <- 0

# Add age points (if available)
if(sum(!is.na(data_work$age_clean)) > 30) {
  age_median <- median(data_work$age_clean, na.rm = TRUE)
  data_work$age_risk <- ifelse(!is.na(data_work$age_clean) & data_work$age_clean > age_median, 1, 0)
  data_work$risk_points <- data_work$risk_points + data_work$age_risk
  cat("Added age component (median = ", round(age_median), ")\n", sep="")
}

# Add MELD points (if available)
if(sum(!is.na(data_work$meld_clean)) > 30) {
  meld_median <- median(data_work$meld_clean, na.rm = TRUE)
  data_work$meld_risk <- ifelse(!is.na(data_work$meld_clean) & data_work$meld_clean > meld_median, 1, 0)
  data_work$risk_points <- data_work$risk_points + data_work$meld_risk
  cat("Added MELD component (median = ", round(meld_median), ")\n", sep="")
}

# Add major complications
major_complications <- c("compli_1", "compli_3", "compli_5")  # ARF, Vent, Reintervention
for(comp in major_complications) {
  if(comp %in% names(data_work)) {
    data_work$risk_points <- data_work$risk_points + data_work[[comp]]
    cat("Added ", comp, " (n=", sum(data_work[[comp]]), ")\n", sep="")
  }
}

# Check risk score distribution
cat("\nRisk score distribution:\n")
print(table(data_work$risk_points))

# ============================================================
# STEP 5: CREATE RISK GROUPS
# ============================================================

cat("\n=== Creating Risk Groups ===\n")

# Method 1: Try tertiles if enough variation
if(length(unique(data_work$risk_points)) >= 3) {
  # Use tertiles with NA handling
  risk_tertiles <- quantile(data_work$risk_points, 
                            probs = c(0, 1/3, 2/3, 1), 
                            na.rm = TRUE)  # This is the fix!
  
  # Handle ties by adjusting breaks
  if(length(unique(risk_tertiles)) < 4) {
    # If ties, use simple equal groups
    data_work$risk_group <- cut(rank(data_work$risk_points, ties.method = "random"),
                                breaks = 3,
                                labels = c("Low", "Medium", "High"))
  } else {
    data_work$risk_group <- cut(data_work$risk_points,
                                breaks = risk_tertiles,
                                labels = c("Low", "Medium", "High"),
                                include.lowest = TRUE)
  }
} else {
  # Method 2: Simple grouping based on score
  cat("Using simple grouping due to limited variation\n")
  data_work$risk_group <- factor(
    ifelse(data_work$risk_points == 0, "Low",
           ifelse(data_work$risk_points <= 1, "Medium", "High")),
    levels = c("Low", "Medium", "High")
  )
}

# Check distribution
cat("\nRisk group distribution:\n")
print(table(data_work$risk_group, useNA = "ifany"))

# ============================================================
# STEP 6: ADD CRE INFORMATION
# ============================================================

cat("\n=== Adding CRE Status ===\n")

# Initialize CRE variable
data_work$cre_180 <- 0

# Check different CRE variable names
if("cre_infection_date" %in% names(data_work)) {
  data_work$cre_infection_date <- as.Date(data_work$cre_infection_date)
  cre_days <- as.numeric(data_work$cre_infection_date - data_work$date_of_transplant)
  data_work$cre_180 <- ifelse(!is.na(cre_days) & cre_days <= 180 & cre_days > 0, 1, 0)
  cat("CRE from dates: ", sum(data_work$cre_180), " infections\n", sep="")
} else if("cre_infection" %in% names(data_work)) {
  data_work$cre_180 <- ifelse(data_work$cre_infection == 1, 1, 0)
  data_work$cre_180[is.na(data_work$cre_180)] <- 0
  cat("CRE from binary variable: ", sum(data_work$cre_180), " infections\n", sep="")
} else if("cre_180d" %in% names(data_work)) {
  data_work$cre_180 <- ifelse(data_work$cre_180d == 1, 1, 0)
  data_work$cre_180[is.na(data_work$cre_180)] <- 0
  cat("CRE from existing variable: ", sum(data_work$cre_180), " infections\n", sep="")
}

cat("Total CRE infections within 180 days: ", sum(data_work$cre_180), "\n", sep="")

# ============================================================
# STEP 7: ANALYSIS
# ============================================================

cat("\n=== ANALYSIS RESULTS ===\n")

# Overall mortality
cat("\nOverall 180-day mortality: ", 
    round(mean(data_work$death_180) * 100, 1), "%\n", sep="")

# By risk group
risk_table <- data_work %>%
  group_by(risk_group) %>%
  summarise(
    n = n(),
    deaths = sum(death_180),
    mortality_pct = round(deaths/n * 100, 1),
    cre_n = sum(cre_180),
    cre_pct = round(cre_n/n * 100, 1)
  )

cat("\n--- Mortality by Risk Group ---\n")
print(as.data.frame(risk_table))

# By CRE status
cre_table <- data_work %>%
  mutate(CRE_status = ifelse(cre_180 == 1, "CRE+", "CRE-")) %>%
  group_by(CRE_status) %>%
  summarise(
    n = n(),
    deaths = sum(death_180),
    mortality_pct = round(deaths/n * 100, 1)
  )

cat("\n--- Mortality by CRE Status ---\n")
print(as.data.frame(cre_table))

# Combined
combined_table <- data_work %>%
  mutate(CRE_status = ifelse(cre_180 == 1, "CRE+", "CRE-")) %>%
  group_by(risk_group, CRE_status) %>%
  summarise(
    n = n(),
    deaths = sum(death_180),
    mortality_pct = round(deaths/n * 100, 1),
    .groups = "drop"
  )

cat("\n--- Mortality by Risk Group and CRE Status ---\n")
print(as.data.frame(combined_table))

# ============================================================
# STEP 8: VISUALIZATION
# ============================================================

# Create bar plot
p <- ggplot(combined_table, aes(x = risk_group, y = mortality_pct, fill = CRE_status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(mortality_pct, "%\n(n=", n, ")")), 
            position = position_dodge(width = 0.8),
            vjust = -0.2, size = 3) +
  scale_fill_manual(values = c("CRE-" = "steelblue", "CRE+" = "red")) +
  labs(title = "180-Day Mortality by Risk Group and CRE Status",
       x = "Risk Group",
       y = "Mortality (%)",
       fill = "CRE Status") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

print(p)

# Save results
baseline_data <- data_work
save(baseline_data, file = "baseline_data_final.RData")
write.csv(combined_table, "mortality_by_risk_and_cre.csv", row.names = FALSE)

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("\nData saved as 'baseline_data' and 'baseline_data_final.RData'\n")
cat("Results saved to 'mortality_by_risk_and_cre.csv'\n")


#############################################
## COX REGRESSION ANALYSIS BY CRE STATUS
## Extension of baseline risk analysis
#############################################

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(broom)

cat("============================================\n")
cat("COX REGRESSION ANALYSIS - CRE STATUS\n")
cat("============================================\n")

# ============================================================
# LOAD THE PROCESSED DATA
# ============================================================

# Assuming baseline_data is already in workspace from previous script
# If not, load it:
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}

data_cox <- baseline_data

cat("\nData loaded: ", nrow(data_cox), " patients\n", sep="")
cat("Deaths: ", sum(data_cox$death_180), "\n", sep="")
cat("CRE infections: ", sum(data_cox$cre_180), "\n", sep="")

# ============================================================
# STEP 1: CREATE SURVIVAL OBJECT
# ============================================================

cat("\n=== Creating Survival Object ===\n")

# Create survival object
surv_obj <- Surv(time = data_cox$time_simple, 
                 event = data_cox$death_180)

cat("Survival object created\n")
cat("Median follow-up: ", round(median(data_cox$time_simple)), " days\n", sep="")

# ============================================================
# STEP 2: KAPLAN-MEIER CURVES BY CRE STATUS
# ============================================================

cat("\n=== Kaplan-Meier Analysis ===\n")

# Fit KM curves
km_fit <- survfit(surv_obj ~ cre_180, data = data_cox)

# Print summary
print(summary(km_fit, times = c(30, 60, 90, 120, 150, 180)))

# Create KM plot with risk table
km_plot <- ggsurvplot(
  km_fit,
  data = data_cox,
  pval = TRUE,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  xlab = "Days from Transplant",
  ylab = "Survival Probability",
  title = "Survival by CRE Infection Status",
  legend.title = "CRE Status",
  legend.labs = c("CRE Negative", "CRE Positive"),
  palette = c("steelblue", "red"),
  risk.table = TRUE,
  risk.table.title = "Number at Risk",
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.25,
  ggtheme = theme_bw(),
  break.time.by = 30,
  xlim = c(0, 180)
)

print(km_plot)

# Save KM plot
ggsave("km_curve_cre.png", 
       plot = km_plot$plot, 
       width = 10, height = 8, dpi = 300)

# ============================================================
# STEP 3: UNIVARIATE COX REGRESSION
# ============================================================

cat("\n=== Univariate Cox Regression ===\n")

# Fit univariate Cox model
cox_uni <- coxph(surv_obj ~ cre_180, data = data_cox)

# Print results
cat("\n--- CRE Status (Univariate) ---\n")
print(summary(cox_uni))

# Extract HR and CI
hr_uni <- exp(coef(cox_uni))
ci_uni <- exp(confint(cox_uni))

cat("\nHazard Ratio (CRE+): ", round(hr_uni, 2), 
    " (95% CI: ", round(ci_uni[1], 2), "-", round(ci_uni[2], 2), ")\n", sep="")

# ============================================================
# STEP 4: MULTIVARIATE COX REGRESSION
# ============================================================

cat("\n=== Multivariate Cox Regression ===\n")

# Build multivariate model with available covariates
formula_parts <- c("cre_180")

# Add age if available
if(sum(!is.na(data_cox$age_clean)) > 30) {
  formula_parts <- c(formula_parts, "age_clean")
}

# Add MELD if available
if(sum(!is.na(data_cox$meld_clean)) > 30) {
  formula_parts <- c(formula_parts, "meld_clean")
}

# Add major complications if available
if("compli_1" %in% names(data_cox)) {
  formula_parts <- c(formula_parts, "compli_1")  # ARF
}
if("compli_3" %in% names(data_cox)) {
  formula_parts <- c(formula_parts, "compli_3")  # Ventilation
}

# Create formula
cox_formula <- as.formula(paste("surv_obj ~", paste(formula_parts, collapse = " + ")))

# Fit multivariate model
cox_multi <- coxph(cox_formula, data = data_cox)

cat("\n--- Multivariate Model ---\n")
print(summary(cox_multi))

# ============================================================
# STEP 5: FOREST PLOT
# ============================================================

cat("\n=== Creating Forest Plot ===\n")

# Extract coefficients for forest plot
cox_tidy <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE)

# Clean variable names for display
cox_tidy$term <- gsub("_clean", "", cox_tidy$term)
cox_tidy$term <- gsub("cre_180", "CRE Infection", cox_tidy$term)
cox_tidy$term <- gsub("age", "Age", cox_tidy$term)
cox_tidy$term <- gsub("meld", "MELD Score", cox_tidy$term)
cox_tidy$term <- gsub("compli_1", "Acute Renal Failure", cox_tidy$term)
cox_tidy$term <- gsub("compli_3", "Mechanical Ventilation", cox_tidy$term)

# Create forest plot
forest_plot <- ggplot(cox_tidy, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.5, 1, 2, 4, 8),
                     labels = c("0.5", "1", "2", "4", "8")) +
  labs(x = "Hazard Ratio (95% CI)",
       y = "",
       title = "Multivariate Cox Regression - Risk Factors for 180-Day Mortality") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(forest_plot)

# Save forest plot
ggsave("forest_plot_cox.png", plot = forest_plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================
# STEP 6: STRATIFIED ANALYSIS BY RISK GROUP
# ============================================================

cat("\n=== Stratified Analysis by Risk Group ===\n")

# Fit stratified Cox model
cox_strat <- coxph(surv_obj ~ cre_180 + strata(risk_group), 
                   data = data_cox)

cat("\n--- CRE Effect Stratified by Risk Group ---\n")
print(summary(cox_strat))

# Alternative: Separate KM curves for each risk group
for(risk_lvl in levels(data_cox$risk_group)) {
  cat("\n--- Risk Group: ", risk_lvl, " ---\n", sep="")
  
  # Subset data
  data_subset <- data_cox[data_cox$risk_group == risk_lvl, ]
  
  # Fit KM for this group
  km_subset <- survfit(Surv(time_simple, death_180) ~ cre_180, 
                       data = data_subset)
  
  # Create plot
  km_plot_subset <- ggsurvplot(
    km_subset,
    data = data_subset,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Days from Transplant",
    ylab = "Survival Probability",
    title = paste("Survival by CRE Status -", risk_lvl, "Risk Group"),
    legend.title = "CRE Status",
    legend.labs = c("CRE-", "CRE+"),
    palette = c("steelblue", "red"),
    ggtheme = theme_bw()
  )
  
  print(km_plot_subset)
  
  # Save individual plots
  filename <- paste0("km_curve_", tolower(risk_lvl), "_risk.png")
  ggsave(filename, plot = km_plot_subset$plot, 
         width = 8, height = 6, dpi = 300)
}

# Try faceted plot (may not work with all versions)
tryCatch({
  km_strat <- survfit(surv_obj ~ cre_180 + risk_group, data = data_cox)
  
  # Plot stratified curves - corrected legend labels
  strat_plot <- ggsurvplot(
    km_strat,
    data = data_cox,
    facet.by = "risk_group",
    conf.int = TRUE,
    pval = TRUE,
    xlab = "Days from Transplant",
    ylab = "Survival Probability",
    title = "Survival by CRE Status Stratified by Risk Group",
    legend.title = "CRE Status",
    legend.labs = c("CRE-", "CRE+"),  # Fixed: only 2 labels needed
    palette = c("steelblue", "red"),   # Fixed: only 2 colors needed
    ggtheme = theme_bw()
  )
  
  print(strat_plot)
}, error = function(e) {
  cat("Note: Faceted plot not available - individual plots created instead\n")
})

# ============================================================
# STEP 7: ADJUSTED SURVIVAL CURVES
# ============================================================

cat("\n=== Adjusted Survival Curves ===\n")

# Create adjusted survival curves using survminer
# Wrap in tryCatch in case of convergence issues
tryCatch({
  adj_curves <- ggadjustedcurves(
    cox_multi,
    data = data_cox,
    variable = "cre_180",
    method = "average",
    ylab = "Survival Probability",
    xlab = "Days from Transplant",
    title = "Adjusted Survival Curves by CRE Status",
    palette = c("steelblue", "red"),
    legend = "top",
    legend.title = "CRE Status",
    legend.labs = c("CRE-", "CRE+"),
    ggtheme = theme_bw()
  )
  
  print(adj_curves)
  
  # Save adjusted curves
  ggsave("adjusted_survival_curves.png", plot = adj_curves, 
         width = 8, height = 6, dpi = 300)
}, error = function(e) {
  cat("Note: Could not create adjusted curves - ", e$message, "\n")
})

# ============================================================
# STEP 8: CUMULATIVE HAZARD PLOT
# ============================================================

cat("\n=== Cumulative Hazard Analysis ===\n")

# Create cumulative hazard plot
cumhaz_plot <- ggsurvplot(
  km_fit,
  data = data_cox,
  fun = "cumhaz",
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Days from Transplant",
  ylab = "Cumulative Hazard",
  title = "Cumulative Hazard by CRE Status",
  legend.title = "CRE Status",
  legend.labs = c("CRE Negative", "CRE Positive"),
  palette = c("steelblue", "red"),
  risk.table = TRUE,
  ggtheme = theme_bw()
)

print(cumhaz_plot)

# ============================================================
# STEP 9: PROPORTIONAL HAZARDS ASSUMPTION
# ============================================================

cat("\n=== Testing Proportional Hazards Assumption ===\n")

# Test PH assumption
ph_test <- cox.zph(cox_multi)
print(ph_test)

# Plot Schoenfeld residuals
par(mfrow = c(2, 2))
plot(ph_test)
par(mfrow = c(1, 1))

# ============================================================
# STEP 10: TIME-VARYING EFFECT (if PH violated)
# ============================================================

if(ph_test$table["cre_180", "p"] < 0.05) {
  cat("\n=== Time-Varying Effect Analysis ===\n")
  cat("Proportional hazards assumption violated for CRE\n")
  
  # Split follow-up time
  data_cox$time_period <- cut(data_cox$time_simple, 
                              breaks = c(0, 30, 90, 180),
                              labels = c("0-30d", "31-90d", "91-180d"))
  
  # Interaction with time
  cox_time <- coxph(surv_obj ~ cre_180 * time_period, data = data_cox)
  print(summary(cox_time))
}

# ============================================================
# STEP 11: SUMMARY TABLE
# ============================================================

cat("\n=== Creating Summary Table ===\n")

# Create summary table
summary_table <- data.frame(
  Variable = c("CRE Negative", "CRE Positive"),
  N = c(sum(data_cox$cre_180 == 0), sum(data_cox$cre_180 == 1)),
  Deaths = c(sum(data_cox$death_180[data_cox$cre_180 == 0]),
             sum(data_cox$death_180[data_cox$cre_180 == 1])),
  Mortality_Rate = NA,
  Median_Survival = NA,
  HR_Univariate = c(1, round(hr_uni, 2)),
  HR_Multivariate = c(1, round(exp(coef(cox_multi)["cre_180"]), 2)),
  P_Value = c(NA, round(summary(cox_multi)$coefficients["cre_180", "Pr(>|z|)"], 3))
)

# Calculate mortality rates
summary_table$Mortality_Rate <- paste0(
  round(summary_table$Deaths / summary_table$N * 100, 1), "%"
)

# Get median survival
km_summary <- summary(km_fit)
summary_table$Median_Survival[1] <- ifelse(is.na(km_summary$table["cre_180=0", "median"]),
                                           ">180", 
                                           round(km_summary$table["cre_180=0", "median"]))
summary_table$Median_Survival[2] <- ifelse(is.na(km_summary$table["cre_180=1", "median"]),
                                           ">180", 
                                           round(km_summary$table["cre_180=1", "median"]))

cat("\n--- Summary Table ---\n")
print(summary_table)

# Save summary table
write.csv(summary_table, "cox_regression_summary.csv", row.names = FALSE)

# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n============================================\n")
cat("COX REGRESSION ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("\nFiles saved:\n")
cat("- km_curve_cre.png: Kaplan-Meier survival curves\n")
cat("- forest_plot_cox.png: Forest plot of risk factors\n")
cat("- adjusted_survival_curves.png: Adjusted survival curves\n")
cat("- cox_regression_summary.csv: Summary statistics table\n")



#############################################
## CALIBRATION PLOT FOR MORTALITY MODEL
## Assessing model calibration and discrimination
#############################################

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(pROC)
library(rms)
library(ResourceSelection)

cat("============================================\n")
cat("MODEL CALIBRATION ANALYSIS\n")
cat("============================================\n")

# ============================================================
# LOAD DATA AND REFIT MODEL
# ============================================================

# Load the processed data
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}

data_cal <- baseline_data

cat("\nData loaded: ", nrow(data_cal), " patients\n", sep="")

# Create survival object
surv_obj <- Surv(time = data_cal$time_simple, 
                 event = data_cal$death_180)

# Build the model formula (same as Cox analysis)
formula_parts <- c("cre_180")

if(sum(!is.na(data_cal$age_clean)) > 30) {
  formula_parts <- c(formula_parts, "age_clean")
}

if(sum(!is.na(data_cal$meld_clean)) > 30) {
  formula_parts <- c(formula_parts, "meld_clean")
}

if("compli_1" %in% names(data_cal)) {
  formula_parts <- c(formula_parts, "compli_1")
}

if("compli_3" %in% names(data_cal)) {
  formula_parts <- c(formula_parts, "compli_3")
}

# Fit Cox model
cox_formula <- as.formula(paste("surv_obj ~", paste(formula_parts, collapse = " + ")))
cox_model <- coxph(cox_formula, data = data_cal)

cat("\nModel fitted with", length(formula_parts), "predictors\n")

# ============================================================
# STEP 1: CALCULATE PREDICTED PROBABILITIES
# ============================================================

cat("\n=== Calculating Predicted Probabilities ===\n")

# First, identify complete cases for the model
model_vars <- all.vars(cox_formula)
model_vars <- model_vars[model_vars != "surv_obj"]
complete_rows <- complete.cases(data_cal[, model_vars])

cat("Complete cases for prediction: ", sum(complete_rows), " out of ", nrow(data_cal), "\n", sep="")

# Get baseline survival at 180 days
basehaz_180 <- survfit(cox_model)
time_points <- basehaz_180$time
surv_probs <- basehaz_180$surv

# Find survival at 180 days (or closest time point)
idx_180 <- which.min(abs(time_points - 180))
baseline_surv_180 <- surv_probs[idx_180]

cat("Baseline survival at 180 days: ", round(baseline_surv_180, 3), "\n", sep="")

# Initialize prediction column with NAs
data_cal$pred_mort_180 <- NA

# Calculate linear predictor only for complete cases
lp <- predict(cox_model, newdata = data_cal[complete_rows, ], type = "lp")

# Calculate predicted mortality probability at 180 days
# P(death by 180) = 1 - S(180)^exp(lp)
data_cal$pred_mort_180[complete_rows] <- 1 - baseline_surv_180^exp(lp)

# Alternative method using survfit for each individual
cat("\nCalculating individual survival curves...\n")
data_cal$pred_mort_alt <- NA

# Only calculate for complete cases
complete_idx <- which(complete_rows)
for(i in complete_idx) {
  tryCatch({
    surv_i <- survfit(cox_model, newdata = data_cal[i,])
    # Find survival at 180 days
    idx <- which.min(abs(surv_i$time - 180))
    if(length(idx) > 0) {
      data_cal$pred_mort_alt[i] <- 1 - surv_i$surv[idx]
    }
  }, error = function(e) {
    # Keep as NA if prediction fails
  })
}

# Use the alternative method if available, otherwise use the first method
if(sum(!is.na(data_cal$pred_mort_alt)) > sum(complete_rows) * 0.8) {
  data_cal$pred_mort <- data_cal$pred_mort_alt
} else {
  data_cal$pred_mort <- data_cal$pred_mort_180
}

# Filter to only complete cases for analysis
data_cal_complete <- data_cal[complete_rows & !is.na(data_cal$pred_mort), ]
cat("Final dataset for calibration: ", nrow(data_cal_complete), " patients\n", sep="")

cat("Predicted mortality range: ", 
    round(min(data_cal_complete$pred_mort, na.rm=T)*100, 1), "% - ",
    round(max(data_cal_complete$pred_mort, na.rm=T)*100, 1), "%\n", sep="")

# ============================================================
# STEP 2: CREATE RISK GROUPS FOR CALIBRATION
# ============================================================

cat("\n=== Creating Risk Groups ===\n")

# Use the complete dataset from here on
data_cal <- data_cal_complete

# Create deciles of predicted risk
data_cal$risk_decile <- cut(data_cal$pred_mort,
                            breaks = quantile(data_cal$pred_mort, 
                                              probs = seq(0, 1, 0.1),
                                              na.rm = TRUE),
                            include.lowest = TRUE,
                            labels = paste0("D", 1:10))

# For smaller datasets, use quintiles instead
if(nrow(data_cal) < 200) {
  data_cal$risk_group_cal <- cut(data_cal$pred_mort,
                                 breaks = quantile(data_cal$pred_mort, 
                                                   probs = seq(0, 1, 0.2),
                                                   na.rm = TRUE),
                                 include.lowest = TRUE,
                                 labels = paste0("Q", 1:5))
  cat("Using quintiles due to sample size (n=", nrow(data_cal), ")\n", sep="")
} else {
  data_cal$risk_group_cal <- data_cal$risk_decile
  cat("Using deciles (n=", nrow(data_cal), ")\n", sep="")
}

# ============================================================
# STEP 3: CALCULATE OBSERVED VS EXPECTED
# ============================================================

cat("\n=== Calculating Calibration Statistics ===\n")

# Calculate observed and expected by risk group
calib_data <- data_cal %>%
  group_by(risk_group_cal) %>%
  summarise(
    n = n(),
    observed_events = sum(death_180),
    expected_events = sum(pred_mort),
    observed_prop = observed_events / n,
    expected_prop = expected_events / n,
    pred_mean = mean(pred_mort),
    .groups = "drop"
  ) %>%
  filter(!is.na(risk_group_cal))

cat("\nCalibration by risk group:\n")
print(as.data.frame(calib_data))

# ============================================================
# STEP 4: BASIC CALIBRATION PLOT
# ============================================================

cat("\n=== Creating Calibration Plot ===\n")

# Calculate confidence intervals for observed proportions
calib_data$obs_lower <- NA
calib_data$obs_upper <- NA
for(i in 1:nrow(calib_data)) {
  if(calib_data$n[i] > 0) {
    ci <- binom.test(calib_data$observed_events[i], 
                     calib_data$n[i])$conf.int
    calib_data$obs_lower[i] <- ci[1]
    calib_data$obs_upper[i] <- ci[2]
  }
}

# Create calibration plot
calib_plot <- ggplot(calib_data, aes(x = expected_prop, y = observed_prop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = obs_lower, ymax = obs_upper), 
                width = 0.01, alpha = 0.5) +
  geom_point(aes(size = n), color = "darkblue", alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "red", size = 1) +
  scale_x_continuous(limits = c(0, max(c(calib_data$expected_prop, 
                                         calib_data$observed_prop)) * 1.1),
                     labels = scales::percent) +
  scale_y_continuous(limits = c(0, max(c(calib_data$expected_prop, 
                                         calib_data$observed_prop)) * 1.1),
                     labels = scales::percent) +
  labs(
    title = "Calibration Plot: 180-Day Mortality Model",
    x = "Predicted Mortality",
    y = "Observed Mortality",
    size = "N patients"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  coord_equal()

print(calib_plot)

# Save calibration plot
ggsave("calibration_plot_basic.png", plot = calib_plot, 
       width = 8, height = 8, dpi = 300)

# ============================================================
# STEP 5: CALIBRATION BELT PLOT
# ============================================================

cat("\n=== Creating Calibration Belt Plot ===\n")

# Create more granular calibration data
n_bins <- 20
data_cal$pred_bin <- cut(data_cal$pred_mort,
                         breaks = seq(0, max(data_cal$pred_mort, na.rm=T), 
                                      length.out = n_bins + 1),
                         include.lowest = TRUE)

# Calculate for each bin
belt_data <- data_cal %>%
  group_by(pred_bin) %>%
  summarise(
    n = n(),
    observed = sum(death_180),
    expected = sum(pred_mort),
    pred_mean = mean(pred_mort),
    obs_prop = observed/n,
    .groups = "drop"
  ) %>%
  filter(!is.na(pred_bin) & n >= 5)  # Remove bins with too few observations

# Create smooth calibration curve
cal_smooth <- loess(obs_prop ~ pred_mean, data = belt_data, 
                    weights = belt_data$n)

# Predict over range
pred_range <- seq(min(belt_data$pred_mean), 
                  max(belt_data$pred_mean), 
                  length.out = 100)
cal_pred <- predict(cal_smooth, pred_range, se = TRUE)

# Create data frame for plotting
smooth_df <- data.frame(
  predicted = pred_range,
  observed = cal_pred$fit,
  se = cal_pred$se.fit
)

smooth_df$lower <- smooth_df$observed - 1.96 * smooth_df$se
smooth_df$upper <- smooth_df$observed + 1.96 * smooth_df$se

# Create belt plot
belt_plot <- ggplot() +
  geom_ribbon(data = smooth_df,
              aes(x = predicted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "blue") +
  geom_line(data = smooth_df,
            aes(x = predicted, y = observed),
            color = "blue", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "red") +
  geom_point(data = belt_data,
             aes(x = pred_mean, y = obs_prop, size = n),
             alpha = 0.5) +
  scale_x_continuous(limits = c(0, max(smooth_df$predicted) * 1.1),
                     labels = scales::percent) +
  scale_y_continuous(limits = c(0, max(smooth_df$upper, na.rm=T) * 1.1),
                     labels = scales::percent) +
  labs(
    title = "Calibration Belt: 180-Day Mortality Model",
    subtitle = "With 95% confidence band",
    x = "Predicted Mortality",
    y = "Observed Mortality",
    size = "N patients"
  ) +
  theme_bw(base_size = 12)

print(belt_plot)

# Save belt plot
ggsave("calibration_belt_plot.png", plot = belt_plot, 
       width = 8, height = 8, dpi = 300)

# ============================================================
# STEP 6: HOSMER-LEMESHOW TEST
# ============================================================

cat("\n=== Hosmer-Lemeshow Test ===\n")

# Perform Hosmer-Lemeshow test
hl_test <- hoslem.test(data_cal$death_180, 
                       data_cal$pred_mort, 
                       g = 10)

cat("\nHosmer-Lemeshow Test Results:\n")
cat("Chi-squared:", round(hl_test$statistic, 2), "\n")
cat("df:", hl_test$parameter, "\n")
cat("p-value:", round(hl_test$p.value, 4), "\n")

if(hl_test$p.value > 0.05) {
  cat("Conclusion: Model appears well calibrated (p > 0.05)\n")
} else {
  cat("Conclusion: Evidence of miscalibration (p < 0.05)\n")
}

# ============================================================
# STEP 7: DISCRIMINATION MEASURES
# ============================================================

cat("\n=== Discrimination Measures ===\n")

# Calculate C-statistic (AUC)
roc_obj <- roc(data_cal$death_180, data_cal$pred_mort)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

cat("\nC-statistic (AUC): ", round(auc_val, 3), 
    " (95% CI: ", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n", sep="")

# Create ROC curve
roc_plot <- ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "gray50") +
  geom_path(aes(x = 1 - roc_obj$specificities, 
                y = roc_obj$sensitivities), 
            color = "blue", size = 1.5) +
  labs(
    title = "ROC Curve: 180-Day Mortality Model",
    subtitle = paste0("AUC = ", round(auc_val, 3), 
                      " (95% CI: ", round(auc_ci[1], 3), 
                      "-", round(auc_ci[3], 3), ")"),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_bw(base_size = 12) +
  coord_equal()

print(roc_plot)

# Save ROC curve
ggsave("roc_curve.png", plot = roc_plot, 
       width = 7, height = 7, dpi = 300)

# ============================================================
# STEP 8: CALIBRATION BY CRE STATUS
# ============================================================

cat("\n=== Calibration by CRE Status ===\n")

# Separate calibration for CRE+ and CRE-
calib_cre <- data_cal %>%
  mutate(CRE_status = ifelse(cre_180 == 1, "CRE+", "CRE-")) %>%
  group_by(CRE_status, risk_group_cal) %>%
  summarise(
    n = n(),
    observed_events = sum(death_180),
    expected_events = sum(pred_mort),
    observed_prop = observed_events / n,
    expected_prop = expected_events / n,
    .groups = "drop"
  ) %>%
  filter(!is.na(risk_group_cal))

# Plot calibration by CRE status
calib_cre_plot <- ggplot(calib_cre, 
                         aes(x = expected_prop, y = observed_prop, 
                             color = CRE_status)) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, size = 1) +
  scale_color_manual(values = c("CRE-" = "steelblue", "CRE+" = "red")) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Calibration by CRE Status",
    x = "Predicted Mortality",
    y = "Observed Mortality",
    color = "CRE Status",
    size = "N patients"
  ) +
  theme_bw(base_size = 12) +
  coord_equal()

print(calib_cre_plot)

# Save CRE calibration plot
ggsave("calibration_plot_cre.png", plot = calib_cre_plot, 
       width = 9, height = 7, dpi = 300)

# ============================================================
# STEP 9: CALIBRATION METRICS SUMMARY
# ============================================================

cat("\n=== Calibration Metrics Summary ===\n")

# Calculate calibration slope
calib_model <- glm(death_180 ~ pred_mort, 
                   data = data_cal, 
                   family = binomial)
calib_slope <- coef(calib_model)[2]

# Calculate calibration intercept
calib_intercept <- coef(calib_model)[1]

# Brier score
brier_score <- mean((data_cal$pred_mort - data_cal$death_180)^2, na.rm = TRUE)

# Expected/Observed ratio
eo_ratio <- sum(data_cal$pred_mort, na.rm = TRUE) / sum(data_cal$death_180)

# Print summary
cat("\n--- Overall Calibration Metrics ---\n")
cat("Calibration intercept: ", round(calib_intercept, 3), "\n", sep="")
cat("Calibration slope: ", round(calib_slope, 3), "\n", sep="")
cat("Expected/Observed ratio: ", round(eo_ratio, 3), "\n", sep="")
cat("Brier score: ", round(brier_score, 4), "\n", sep="")
cat("C-statistic: ", round(auc_val, 3), "\n", sep="")
cat("Hosmer-Lemeshow p-value: ", round(hl_test$p.value, 4), "\n", sep="")

# Create summary table
calib_summary <- data.frame(
  Metric = c("Calibration Intercept", "Calibration Slope", 
             "E/O Ratio", "Brier Score", "C-statistic", 
             "Hosmer-Lemeshow p-value"),
  Value = c(round(calib_intercept, 3), round(calib_slope, 3),
            round(eo_ratio, 3), round(brier_score, 4),
            round(auc_val, 3), round(hl_test$p.value, 4)),
  Interpretation = c(
    ifelse(abs(calib_intercept) < 0.1, "Good", "Poor"),
    ifelse(abs(calib_slope - 1) < 0.1, "Good", "Poor"),
    ifelse(abs(eo_ratio - 1) < 0.1, "Good", "Poor"),
    ifelse(brier_score < 0.25, "Good", "Poor"),
    ifelse(auc_val > 0.7, "Good", "Fair"),
    ifelse(hl_test$p.value > 0.05, "Good calibration", "Poor calibration")
  )
)

cat("\n--- Calibration Summary Table ---\n")
print(calib_summary)

# Save summary table
write.csv(calib_summary, "calibration_metrics.csv", row.names = FALSE)

# ============================================================
# STEP 10: DECISION CURVE ANALYSIS
# ============================================================

cat("\n=== Decision Curve Analysis ===\n")

# Calculate net benefit across threshold probabilities
thresholds <- seq(0, 0.5, by = 0.01)
net_benefit <- numeric(length(thresholds))
net_benefit_all <- numeric(length(thresholds))
net_benefit_none <- rep(0, length(thresholds))

for(i in 1:length(thresholds)) {
  thresh <- thresholds[i]
  
  # Treat if predicted risk > threshold
  treat <- data_cal$pred_mort > thresh
  
  # True positives and false positives
  tp <- sum(treat & data_cal$death_180 == 1, na.rm = TRUE)
  fp <- sum(treat & data_cal$death_180 == 0, na.rm = TRUE)
  n <- sum(!is.na(data_cal$pred_mort))
  
  # Net benefit
  net_benefit[i] <- (tp/n) - (fp/n) * (thresh/(1-thresh))
  
  # Net benefit treat all
  event_rate <- mean(data_cal$death_180, na.rm = TRUE)
  net_benefit_all[i] <- event_rate - (1-event_rate) * (thresh/(1-thresh))
}

# Create decision curve data
dc_data <- data.frame(
  threshold = rep(thresholds, 3),
  net_benefit = c(net_benefit, net_benefit_all, net_benefit_none),
  strategy = factor(rep(c("Model", "Treat All", "Treat None"), 
                        each = length(thresholds)))
)

# Plot decision curve
dc_plot <- ggplot(dc_data, aes(x = threshold, y = net_benefit, 
                               color = strategy, linetype = strategy)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Model" = "blue", 
                                "Treat All" = "darkgreen", 
                                "Treat None" = "red")) +
  scale_linetype_manual(values = c("Model" = "solid", 
                                   "Treat All" = "dashed", 
                                   "Treat None" = "dotted")) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title = "Decision Curve Analysis",
    subtitle = "Net benefit across threshold probabilities",
    x = "Threshold Probability",
    y = "Net Benefit",
    color = "Strategy",
    linetype = "Strategy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

print(dc_plot)

# Save decision curve
ggsave("decision_curve.png", plot = dc_plot, 
       width = 9, height = 6, dpi = 300)

# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n============================================\n")
cat("CALIBRATION ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("\nFiles saved:\n")
cat("- calibration_plot_basic.png: Standard calibration plot\n")
cat("- calibration_belt_plot.png: Calibration with confidence belt\n")
cat("- calibration_plot_cre.png: Calibration by CRE status\n")
cat("- roc_curve.png: ROC curve with AUC\n")
cat("- decision_curve.png: Decision curve analysis\n")
cat("- calibration_metrics.csv: Summary metrics table\n")


