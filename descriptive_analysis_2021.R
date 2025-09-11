library (readr)
data_2025 <- read_csv("data/CRECOOLT_2025.csv")
data_2021 <- read_csv("data/CRE3.csv")

# ============================================================================
# CRE3 Dataset - Simple Descriptive Analysis by Death Status (Base R Version)
# ============================================================================
# This version uses primarily base R functions with minimal dependencies

# Load the data
cre3 <- read.csv("CRE3.csv", stringsAsFactors = FALSE)

# ============================================================================
# 1. BASIC DATASET INFORMATION
# ============================================================================

cat("============================================\n")
cat("DATASET OVERVIEW\n")
cat("============================================\n")
cat("Number of observations:", nrow(cre3), "\n")
cat("Number of variables:", ncol(cre3), "\n")
cat("Variable names:", paste(names(cre3), collapse = ", "), "\n\n")

# Death status distribution
cat("DEATH STATUS DISTRIBUTION:\n")
death_table <- table(cre3$death)
death_prop <- prop.table(death_table)
cat("Death = 0 (Survived):", death_table[1], "(", round(death_prop[1]*100, 2), "%)\n")
cat("Death = 1 (Died):", death_table[2], "(", round(death_prop[2]*100, 2), "%)\n\n")

# ============================================================================
# 2. BINARY VARIABLES ANALYSIS
# ============================================================================

cat("============================================\n")
cat("BINARY VARIABLES BY DEATH STATUS\n")
cat("============================================\n\n")

# List of binary variables (excluding death itself)
binary_vars <- c("atg", "multisite", "kpc", "arf", "mv", "reint", 
                 "crepre_60", "multipre", "multipost", "crepost_60", 
                 "status", "crestat")

# Create a results dataframe
binary_results <- data.frame(
  Variable = character(),
  Death_0_No = integer(),
  Death_0_Yes = integer(),
  Death_0_Yes_Pct = numeric(),
  Death_1_No = integer(),
  Death_1_Yes = integer(),
  Death_1_Yes_Pct = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Analyze each binary variable
for(var in binary_vars) {
  if(var %in% names(cre3)) {
    cat("---", toupper(var), "---\n")
    
    # Create contingency table
    cont_table <- table(cre3[[var]], cre3$death)
    
    # Calculate percentages
    col_pct <- prop.table(cont_table, margin = 2) * 100
    
    # Display results
    cat("Death = 0 (Survived):\n")
    cat("  ", var, "= 0:", cont_table[1,1], "(", round(col_pct[1,1], 1), "%)\n")
    cat("  ", var, "= 1:", cont_table[2,1], "(", round(col_pct[2,1], 1), "%)\n")
    
    cat("Death = 1 (Died):\n")
    cat("  ", var, "= 0:", cont_table[1,2], "(", round(col_pct[1,2], 1), "%)\n")
    cat("  ", var, "= 1:", cont_table[2,2], "(", round(col_pct[2,2], 1), "%)\n")
    
    # Chi-square test
    chi_test <- tryCatch(
      chisq.test(cont_table),
      error = function(e) list(p.value = NA)
    )
    cat("P-value (Chi-square):", format(chi_test$p.value, digits = 4), "\n\n")
    
    # Store results
    new_row <- data.frame(
      Variable = var,
      Death_0_No = cont_table[1,1],
      Death_0_Yes = cont_table[2,1],
      Death_0_Yes_Pct = round(col_pct[2,1], 1),
      Death_1_No = cont_table[1,2],
      Death_1_Yes = cont_table[2,2],
      Death_1_Yes_Pct = round(col_pct[2,2], 1),
      P_Value = chi_test$p.value,
      stringsAsFactors = FALSE
    )
    binary_results <- rbind(binary_results, new_row)
  }
}

# ============================================================================
# 3. CONTINUOUS VARIABLES ANALYSIS
# ============================================================================

cat("============================================\n")
cat("CONTINUOUS VARIABLES BY DEATH STATUS\n")
cat("============================================\n\n")

# List of continuous variables
continuous_vars <- c("time", "cretime", "futime")

# Create a results dataframe
continuous_results <- data.frame(
  Variable = character(),
  Death_0_N = integer(),
  Death_0_Mean = numeric(),
  Death_0_SD = numeric(),
  Death_0_Median = numeric(),
  Death_1_N = integer(),
  Death_1_Mean = numeric(),
  Death_1_SD = numeric(),
  Death_1_Median = numeric(),
  P_Value_Ttest = numeric(),
  stringsAsFactors = FALSE
)

# Analyze each continuous variable
for(var in continuous_vars) {
  if(var %in% names(cre3)) {
    cat("---", toupper(var), "---\n")
    
    # Split by death status
    group0 <- cre3[[var]][cre3$death == 0]
    group1 <- cre3[[var]][cre3$death == 1]
    
    # Remove NAs
    group0 <- group0[!is.na(group0)]
    group1 <- group1[!is.na(group1)]
    
    # Calculate statistics
    cat("Death = 0 (Survived):\n")
    cat("  N:", length(group0), "\n")
    cat("  Mean (SD):", round(mean(group0), 2), "(", round(sd(group0), 2), ")\n")
    cat("  Median [Q1, Q3]:", round(median(group0), 2), 
        "[", round(quantile(group0, 0.25), 2), ",", 
        round(quantile(group0, 0.75), 2), "]\n")
    cat("  Range:", round(min(group0), 2), "-", round(max(group0), 2), "\n")
    
    cat("Death = 1 (Died):\n")
    cat("  N:", length(group1), "\n")
    cat("  Mean (SD):", round(mean(group1), 2), "(", round(sd(group1), 2), ")\n")
    cat("  Median [Q1, Q3]:", round(median(group1), 2), 
        "[", round(quantile(group1, 0.25), 2), ",", 
        round(quantile(group1, 0.75), 2), "]\n")
    cat("  Range:", round(min(group1), 2), "-", round(max(group1), 2), "\n")
    
    # T-test
    if(length(group0) > 1 & length(group1) > 1) {
      t_test <- t.test(group0, group1)
      cat("P-value (T-test):", format(t_test$p.value, digits = 4), "\n")
      p_val <- t_test$p.value
    } else {
      cat("P-value: Not calculated (insufficient data)\n")
      p_val <- NA
    }
    cat("\n")
    
    # Store results
    new_row <- data.frame(
      Variable = var,
      Death_0_N = length(group0),
      Death_0_Mean = round(mean(group0), 2),
      Death_0_SD = round(sd(group0), 2),
      Death_0_Median = round(median(group0), 2),
      Death_1_N = length(group1),
      Death_1_Mean = round(mean(group1), 2),
      Death_1_SD = round(sd(group1), 2),
      Death_1_Median = round(median(group1), 2),
      P_Value_Ttest = p_val,
      stringsAsFactors = FALSE
    )
    continuous_results <- rbind(continuous_results, new_row)
  }
}

# ============================================================================
# 4. EXPORT RESULTS TO CSV
# ============================================================================

cat("============================================\n")
cat("EXPORTING RESULTS\n")
cat("============================================\n\n")

# Export binary variables results
write.csv(binary_results, "CRE3_binary_analysis.csv", row.names = FALSE)
cat("Binary variables analysis saved to: CRE3_binary_analysis.csv\n")

# Export continuous variables results
write.csv(continuous_results, "CRE3_continuous_analysis.csv", row.names = FALSE)
cat("Continuous variables analysis saved to: CRE3_continuous_analysis.csv\n")

# ============================================================================
# 5. SUMMARY OF SIGNIFICANT FINDINGS
# ============================================================================

cat("\n============================================\n")
cat("SUMMARY OF SIGNIFICANT FINDINGS (p < 0.05)\n")
cat("============================================\n\n")

# Significant binary variables
sig_binary <- binary_results[!is.na(binary_results$P_Value) & 
                               binary_results$P_Value < 0.05, ]
if(nrow(sig_binary) > 0) {
  cat("Significant Binary Variables:\n")
  for(i in 1:nrow(sig_binary)) {
    cat("-", sig_binary$Variable[i], 
        "(p =", format(sig_binary$P_Value[i], digits = 4), ")\n")
  }
} else {
  cat("No significant binary variables found.\n")
}

cat("\n")

# Significant continuous variables
sig_continuous <- continuous_results[!is.na(continuous_results$P_Value_Ttest) & 
                                       continuous_results$P_Value_Ttest < 0.05, ]
if(nrow(sig_continuous) > 0) {
  cat("Significant Continuous Variables:\n")
  for(i in 1:nrow(sig_continuous)) {
    cat("-", sig_continuous$Variable[i], 
        "(p =", format(sig_continuous$P_Value_Ttest[i], digits = 4), ")\n")
  }
} else {
  cat("No significant continuous variables found.\n")
}

# ============================================================================
# 6. CREATE SIMPLE VISUALIZATIONS (if possible)
# ============================================================================

# Try to create plots if graphics are available
if(capabilities("png")) {
  cat("\n============================================\n")
  cat("CREATING VISUALIZATIONS\n")
  cat("============================================\n\n")
  
  # Set up plotting area
  pdf("CRE3_simple_plots.pdf", width = 10, height = 8)
  
  # Plot 1: Death distribution
  par(mfrow = c(2, 2))
  barplot(death_table, 
          main = "Death Status Distribution",
          names.arg = c("Survived", "Died"),
          col = c("lightblue", "coral"),
          ylab = "Frequency")
  
  # Plot 2-4: Top 3 binary variables with smallest p-values
  top_vars <- head(binary_results[order(binary_results$P_Value), ], 3)
  
  for(i in 1:min(3, nrow(top_vars))) {
    var <- top_vars$Variable[i]
    if(var %in% names(cre3)) {
      cont_table <- table(cre3[[var]], cre3$death)
      barplot(cont_table,
              main = paste(var, "by Death Status"),
              xlab = "Death Status",
              ylab = "Frequency",
              names.arg = c("Survived", "Died"),
              col = c("lightgray", "darkgray"),
              legend = c("0", "1"),
              beside = TRUE)
    }
  }
  
  # Plot continuous variables
  par(mfrow = c(2, 2))
  for(var in continuous_vars) {
    if(var %in% names(cre3)) {
      boxplot(cre3[[var]] ~ cre3$death,
              main = paste(var, "by Death Status"),
              xlab = "Death Status",
              ylab = var,
              names = c("Survived", "Died"),
              col = c("lightblue", "coral"))
    }
  }
  
  dev.off()
  cat("Plots saved to: CRE3_simple_plots.pdf\n")
}

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")

