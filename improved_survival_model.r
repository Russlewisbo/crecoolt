#############################################
## SIMPLIFIED COMPLEX MODEL BUILDING
## Focus on practical model improvement
#############################################

library(survival)
library(dplyr)
library(ggplot2)
library(splines)
library(pROC)

cat("============================================\n")
cat("BUILDING BETTER CALIBRATED MODELS\n")
cat("============================================\n")

# Load data
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}

data_model <- baseline_data

cat("\nTotal patients: ", nrow(data_model), "\n", sep="")
cat("Outcome events: ", sum(data_model$death_180), "\n", sep="")
cat("Event rate: ", round(mean(data_model$death_180) * 100, 1), "%\n", sep="")

# ============================================================
# STEP 1: CREATE ENHANCED FEATURES
# ============================================================

cat("\n=== Creating Enhanced Features ===\n")

# A. Complication scores
cat("\n--- Complication Scores ---\n")

# Major complications (more serious)
major_compli_vars <- c("compli_1", "compli_3", "compli_5")  # ARF, Vent, Reintervention
data_model$major_complications <- 0
for(v in major_compli_vars) {
  if(v %in% names(data_model)) {
    data_model$major_complications <- data_model$major_complications + data_model[[v]]
  }
}

# Total complications
all_compli_vars <- names(data_model)[grep("^compli_", names(data_model))]
data_model$total_complications <- 0
for(v in all_compli_vars) {
  data_model$total_complications <- data_model$total_complications + data_model[[v]]
}

cat("Major complications range: 0-", max(data_model$major_complications, na.rm=T), "\n", sep="")
cat("Total complications range: 0-", max(data_model$total_complications, na.rm=T), "\n", sep="")

# B. Create interaction terms (CRITICAL for calibration)
cat("\n--- Interaction Terms ---\n")

# CRE × Age
if(sum(!is.na(data_model$age_clean)) > 30) {
  # Center age for better interpretation
  age_centered <- data_model$age_clean - mean(data_model$age_clean, na.rm=T)
  data_model$cre_age_int <- data_model$cre_180 * age_centered
  cat("Created CRE × Age interaction\n")
}

# CRE × MELD
if(sum(!is.na(data_model$meld_clean)) > 30) {
  # Center MELD
  meld_centered <- data_model$meld_clean - mean(data_model$meld_clean, na.rm=T)
  data_model$cre_meld_int <- data_model$cre_180 * meld_centered
  cat("Created CRE × MELD interaction\n")
}

# CRE × Major complications
data_model$cre_compli_int <- data_model$cre_180 * data_model$major_complications
cat("Created CRE × Complications interaction\n")

# Age × MELD (synergistic effect)
if(sum(!is.na(data_model$age_clean)) > 30 & sum(!is.na(data_model$meld_clean)) > 30) {
  age_std <- scale(data_model$age_clean)[,1]
  meld_std <- scale(data_model$meld_clean)[,1]
  data_model$age_meld_int <- age_std * meld_std
  cat("Created Age × MELD interaction\n")
}

# C. Non-linear terms (for capturing threshold effects)
cat("\n--- Non-linear Terms ---\n")

# Age splines (3 knots for simplicity)
if(sum(!is.na(data_model$age_clean)) > 50) {
  age_spline <- ns(data_model$age_clean, df=2)  # 2 df = 3 knots
  data_model$age_spline1 <- age_spline[,1]
  data_model$age_spline2 <- age_spline[,2]
  cat("Created age splines (2 df)\n")
}

# MELD splines
if(sum(!is.na(data_model$meld_clean)) > 50) {
  meld_spline <- ns(data_model$meld_clean, df=2)
  data_model$meld_spline1 <- meld_spline[,1]
  data_model$meld_spline2 <- meld_spline[,2]
  cat("Created MELD splines (2 df)\n")
}

# ============================================================
# STEP 2: BUILD PROGRESSIVE MODELS
# ============================================================

cat("\n=== Building Progressive Models ===\n")

# Create survival object
surv_obj <- Surv(time = data_model$time_simple, event = data_model$death_180)

# Initialize all models to NULL
model1 <- model2 <- model3 <- model4 <- NULL
vars1 <- vars2 <- vars3 <- vars4 <- NULL

# Model 1: ORIGINAL (Simple)
cat("\n--- Model 1: Original Simple ---\n")
formula1 <- surv_obj ~ cre_180 + age_clean + meld_clean + compli_1 + compli_3

# Only keep variables that exist
vars1 <- all.vars(formula1)[-1]  # Remove surv_obj
vars1 <- vars1[vars1 %in% names(data_model)]

if(length(vars1) > 0) {
  formula1 <- as.formula(paste("surv_obj ~", paste(vars1, collapse = " + ")))
  
  tryCatch({
    model1 <- coxph(formula1, data = data_model)
    n_events1 <- model1$nevent
    cat("Variables: ", length(vars1), ", Events: ", n_events1, "\n", sep="")
  }, error = function(e) {
    cat("Error fitting Model 1: ", e$message, "\n")
    # Fallback to minimal model
    model1 <<- coxph(surv_obj ~ cre_180, data = data_model)
    cat("Using minimal model with CRE only\n")
  })
} else {
  model1 <- coxph(surv_obj ~ cre_180, data = data_model)
  cat("Using minimal model with CRE only\n")
}

# Initialize vars2 with vars1 as fallback
vars2 <- vars1

# Model 2: ADD COMPLICATIONS SUMMARY
cat("\n--- Model 2: With Complication Scores ---\n")
formula2 <- surv_obj ~ cre_180 + age_clean + meld_clean + 
  major_complications + total_complications

vars2_temp <- all.vars(formula2)[-1]
vars2_temp <- vars2_temp[vars2_temp %in% names(data_model)]

if(length(vars2_temp) > length(vars1)) {
  vars2 <- vars2_temp
  formula2 <- as.formula(paste("surv_obj ~", paste(vars2, collapse = " + ")))
  
  tryCatch({
    model2 <- coxph(formula2, data = data_model)
    cat("Variables: ", length(vars2), "\n", sep="")
    if(!is.null(model1) && AIC(model2) < AIC(model1)) {
      cat("AIC improved by: ", round(AIC(model1) - AIC(model2), 1), "\n", sep="")
    }
  }, error = function(e) {
    cat("Model 2 failed, using Model 1\n")
    model2 <<- model1
    vars2 <<- vars1
  })
} else {
  model2 <- model1
  cat("No additional complication variables, using Model 1\n")
}

# Model 3: ADD KEY INTERACTIONS
cat("\n--- Model 3: With Interactions ---\n")
int_vars <- c("cre_age_int", "cre_meld_int", "cre_compli_int", "age_meld_int")
int_vars <- int_vars[int_vars %in% names(data_model)]

# Initialize vars3 with vars2 as fallback
vars3 <- vars2

if(length(int_vars) > 0) {
  vars3 <- c(vars2, int_vars)
  formula3 <- as.formula(paste("surv_obj ~", paste(vars3, collapse = " + ")))
  
  tryCatch({
    model3 <- coxph(formula3, data = data_model)
    cat("Variables: ", length(vars3), ", Added ", length(int_vars), " interactions\n", sep="")
    if(!is.null(model1) && AIC(model3) < AIC(model1)) {
      cat("AIC improved by: ", round(AIC(model1) - AIC(model3), 1), "\n", sep="")
    }
  }, error = function(e) {
    cat("Model 3 failed, using Model 2\n")
    model3 <<- model2
    vars3 <<- vars2
  })
} else {
  model3 <- model2
  cat("No interaction variables available, using Model 2\n")
}

# Model 4: ADD NON-LINEAR TERMS
cat("\n--- Model 4: With Non-linear Terms ---\n")
nl_vars <- c("age_spline1", "age_spline2", "meld_spline1", "meld_spline2")
nl_vars <- nl_vars[nl_vars %in% names(data_model)]

# Initialize vars4 with vars3 as fallback
vars4 <- vars3

if(length(nl_vars) > 0) {
  # Use splines instead of linear terms
  vars4 <- vars3[!vars3 %in% c("age_clean", "meld_clean")]
  vars4 <- c(vars4, nl_vars)
  formula4 <- as.formula(paste("surv_obj ~", paste(vars4, collapse = " + ")))
  
  tryCatch({
    model4 <- coxph(formula4, data = data_model)
    cat("Variables: ", length(vars4), "\n", sep="")
    if(!is.null(model1) && AIC(model4) < AIC(model1)) {
      cat("AIC improved by: ", round(AIC(model1) - AIC(model4), 1), "\n", sep="")
    }
  }, error = function(e) {
    cat("Model 4 failed, using Model 3\n")
    model4 <<- model3
    vars4 <<- vars3
  })
} else {
  model4 <- model3
  cat("No spline variables available, using Model 3\n")
}

# ============================================================
# STEP 3: EVALUATE CALIBRATION
# ============================================================

cat("\n=== Evaluating Model Calibration ===\n")

evaluate_calibration <- function(model, data, model_name) {
  
  # Check if model is NULL
  if(is.null(model)) {
    return(data.frame(
      Model = model_name,
      N = NA,
      Events = NA,
      Cal_Intercept = NA,
      Cal_Slope = NA,
      Slope_Diff = NA,
      Brier = NA,
      AUC = NA,
      Pred_SD = NA
    ))
  }
  
  # Get complete cases for this model's variables
  tryCatch({
    model_vars <- names(model$coefficients)
    complete_rows <- complete.cases(data[, model_vars])
    data_eval <- data[complete_rows, ]
    
    # Get linear predictor
    lp <- predict(model, newdata = data_eval, type = "lp")
    
    # Get baseline survival
    basehaz <- survfit(model)
    idx_180 <- which.min(abs(basehaz$time - 180))
    if(length(idx_180) > 0) {
      baseline_surv <- basehaz$surv[idx_180]
    } else {
      baseline_surv <- 0.8  # Default if 180 not reached
    }
    
    # Calculate predicted probabilities
    pred_prob <- 1 - baseline_surv^exp(lp)
    
    # Ensure valid range
    pred_prob <- pmax(pmin(pred_prob, 0.999), 0.001)
    
    # Get outcomes
    outcome <- data_eval$death_180
    
    # Calculate calibration metrics
    result <- tryCatch({
      # Logit transform
      pred_logit <- qlogis(pred_prob)
      
      # Fit calibration model
      if(length(unique(outcome)) > 1 & sd(pred_logit) > 0) {
        cal_model <- glm(outcome ~ pred_logit, family = binomial)
        cal_int <- round(coef(cal_model)[1], 3)
        cal_slope <- round(coef(cal_model)[2], 3)
      } else {
        cal_int <- NA
        cal_slope <- NA
      }
      
      # Other metrics
      brier <- round(mean((pred_prob - outcome)^2), 4)
      auc_val <- round(as.numeric(auc(roc(outcome, pred_prob, quiet = TRUE))), 3)
      pred_sd <- round(sd(pred_prob), 4)
      
      data.frame(
        Model = model_name,
        N = length(outcome),
        Events = sum(outcome),
        Cal_Intercept = cal_int,
        Cal_Slope = cal_slope,
        Slope_Diff = ifelse(is.na(cal_slope), NA, round(abs(cal_slope - 1), 3)),
        Brier = brier,
        AUC = auc_val,
        Pred_SD = pred_sd
      )
      
    }, error = function(e) {
      data.frame(
        Model = model_name,
        N = length(outcome),
        Events = sum(outcome),
        Cal_Intercept = NA,
        Cal_Slope = NA,
        Slope_Diff = NA,
        Brier = round(mean((pred_prob - outcome)^2), 4),
        AUC = NA,
        Pred_SD = round(sd(pred_prob), 4)
      )
    })
    
    return(result)
    
  }, error = function(e) {
    return(data.frame(
      Model = model_name,
      N = NA,
      Events = NA,
      Cal_Intercept = NA,
      Cal_Slope = NA,
      Slope_Diff = NA,
      Brier = NA,
      AUC = NA,
      Pred_SD = NA
    ))
  })
}

# Evaluate all models
cal_results <- rbind(
  evaluate_calibration(model1, data_model, "1. Original"),
  evaluate_calibration(model2, data_model, "2. +Complications"),
  evaluate_calibration(model3, data_model, "3. +Interactions"),
  evaluate_calibration(model4, data_model, "4. +Non-linear")
)

cat("\n--- Calibration Results ---\n")
print(cal_results)

# Find best model
valid_results <- cal_results[!is.na(cal_results$Slope_Diff), ]
if(nrow(valid_results) > 0) {
  best_idx <- which.min(cal_results$Slope_Diff)
  best_model_name <- cal_results$Model[best_idx]
} else {
  # If no valid calibration metrics, use AIC
  cat("\nNo valid calibration slopes, selecting based on AIC\n")
  aics <- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4))
  best_idx <- which.min(aics)
  best_model_name <- cal_results$Model[best_idx]
}
cat("\nBest calibrated model: ", best_model_name, "\n", sep="")

# Select best model object
if(best_idx == 1) {
  best_model <- model1
} else if(best_idx == 2) {
  best_model <- model2
} else if(best_idx == 3) {
  best_model <- model3
} else if(best_idx == 4) {
  best_model <- model4
} else {
  # Fallback to model1 if something went wrong
  best_model <- model1
  cat("Using Model 1 as fallback\n")
}

# ============================================================
# STEP 4: VISUALIZE IMPROVEMENT
# ============================================================

cat("\n=== Creating Calibration Plots ===\n")

# Function to create calibration plot
plot_calibration <- function(model, data, title) {
  
  # Check if model is NULL
  if(is.null(model)) {
    return(ggplot() + 
             ggtitle(paste(title, "- Model not available")) + 
             theme_bw())
  }
  
  tryCatch({
    # Get predictions
    model_vars <- names(model$coefficients)
    complete_rows <- complete.cases(data[, model_vars])
    data_plot <- data[complete_rows, ]
    
    lp <- predict(model, newdata = data_plot, type = "lp")
    
    basehaz <- survfit(model)
    idx_180 <- which.min(abs(basehaz$time - 180))
    baseline_surv <- ifelse(length(idx_180) > 0, basehaz$surv[idx_180], 0.8)
    
    data_plot$pred <- 1 - baseline_surv^exp(lp)
    data_plot$pred <- pmax(pmin(data_plot$pred, 0.999), 0.001)
    
    # Create risk groups
    n_groups <- min(10, floor(nrow(data_plot)/20))
    if(n_groups < 2) n_groups <- 2
    
    data_plot$risk_group <- cut(data_plot$pred,
                                breaks = quantile(data_plot$pred, 
                                                  probs = seq(0, 1, length.out = n_groups + 1)),
                                include.lowest = TRUE)
    
    # Calculate observed vs expected
    cal_data <- data_plot %>%
      group_by(risk_group) %>%
      summarise(
        n = n(),
        predicted = mean(pred),
        observed = mean(death_180),
        .groups = "drop"
      )
    
    # Add confidence intervals
    for(i in 1:nrow(cal_data)) {
      ci <- binom.test(round(cal_data$observed[i] * cal_data$n[i]), 
                       cal_data$n[i])$conf.int
      cal_data$lower[i] <- ci[1]
      cal_data$upper[i] <- ci[2]
    }
    
    # Create plot
    p <- ggplot(cal_data, aes(x = predicted, y = observed)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01, alpha = 0.5) +
      geom_point(aes(size = n), color = "darkblue", alpha = 0.7) +
      geom_smooth(method = "loess", se = TRUE, color = "red", size = 1) +
      scale_x_continuous(limits = c(0, max(cal_data$predicted) * 1.1),
                         labels = scales::percent) +
      scale_y_continuous(limits = c(0, max(cal_data$upper) * 1.1),
                         labels = scales::percent) +
      labs(title = title,
           x = "Predicted Mortality",
           y = "Observed Mortality",
           size = "N patients") +
      theme_bw(base_size = 12) +
      coord_equal()
    
    return(p)
    
  }, error = function(e) {
    return(ggplot() + 
             ggtitle(paste(title, "- Error creating plot")) + 
             theme_bw())
  })
}

# Plot original and best model
p1 <- plot_calibration(model1, data_model, "Original Model Calibration")
p_best <- plot_calibration(best_model, data_model, paste("Best Model:", best_model_name))

print(p1)
print(p_best)

# Save plots
ggsave("calibration_original.png", plot = p1, width = 8, height = 8, dpi = 300)
ggsave("calibration_best.png", plot = p_best, width = 8, height = 8, dpi = 300)

# ============================================================
# STEP 5: FEATURE IMPORTANCE
# ============================================================

cat("\n=== Feature Importance (Best Model) ===\n")

if(!is.null(best_model)) {
  # Get coefficients
  tryCatch({
    coef_summary <- summary(best_model)$coefficients
    importance <- data.frame(
      Variable = rownames(coef_summary),
      HR = round(exp(coef_summary[,1]), 2),
      Coef = round(coef_summary[,1], 3),
      SE = round(coef_summary[,3], 3),
      Z = round(coef_summary[,4], 2),
      P = round(coef_summary[,5], 4)
    )
    
    # Sort by absolute Z-score
    importance <- importance[order(abs(importance$Z), decreasing = TRUE), ]
    
    cat("\n--- Top 10 Most Important Predictors ---\n")
    print(head(importance, 10))
    
    # Identify key findings
    if("cre_age_int" %in% importance$Variable) {
      cre_age_row <- importance[importance$Variable == "cre_age_int", ]
      if(nrow(cre_age_row) > 0 && cre_age_row$P[1] < 0.05) {
        cat("\n* CRE × Age interaction is significant (p = ", cre_age_row$P[1], ")\n", sep="")
        cat("  This means CRE effect varies by patient age\n")
      }
    }
    
    if("cre_meld_int" %in% importance$Variable) {
      cre_meld_row <- importance[importance$Variable == "cre_meld_int", ]
      if(nrow(cre_meld_row) > 0 && cre_meld_row$P[1] < 0.05) {
        cat("\n* CRE × MELD interaction is significant (p = ", cre_meld_row$P[1], ")\n", sep="")
        cat("  This means CRE effect varies by disease severity\n")
      }
    }
    
    if("cre_compli_int" %in% importance$Variable) {
      cre_compli_row <- importance[importance$Variable == "cre_compli_int", ]
      if(nrow(cre_compli_row) > 0 && cre_compli_row$P[1] < 0.05) {
        cat("\n* CRE × Complications interaction is significant (p = ", cre_compli_row$P[1], ")\n", sep="")
        cat("  This means CRE effect varies by complication burden\n")
      }
    }
    
  }, error = function(e) {
    cat("Could not extract feature importance\n")
    importance <- NULL
  })
} else {
  cat("No valid model for feature importance\n")
  importance <- NULL
}

# ============================================================
# STEP 6: PREDICTIONS DISTRIBUTION
# ============================================================

cat("\n=== Prediction Distribution Analysis ===\n")

# Compare prediction distributions
get_predictions <- function(model, data) {
  if(is.null(model)) {
    return(rep(NA, nrow(data)))
  }
  
  tryCatch({
    model_vars <- names(model$coefficients)
    complete_rows <- complete.cases(data[, model_vars])
    data_pred <- data[complete_rows, ]
    
    lp <- predict(model, newdata = data_pred, type = "lp")
    basehaz <- survfit(model)
    idx_180 <- which.min(abs(basehaz$time - 180))
    baseline_surv <- ifelse(length(idx_180) > 0, basehaz$surv[idx_180], 0.8)
    
    pred <- 1 - baseline_surv^exp(lp)
    return(pmax(pmin(pred, 0.999), 0.001))
  }, error = function(e) {
    return(rep(NA, nrow(data)))
  })
}

pred_orig <- get_predictions(model1, data_model)
pred_best <- get_predictions(best_model, data_model)

# Remove NAs for comparison
pred_orig <- pred_orig[!is.na(pred_orig)]
pred_best <- pred_best[!is.na(pred_best)]

if(length(pred_orig) > 0) {
  cat("\nOriginal model:\n")
  cat("  Prediction range: ", round(min(pred_orig), 3), " - ", round(max(pred_orig), 3), "\n", sep="")
  cat("  Prediction SD: ", round(sd(pred_orig), 4), "\n", sep="")
} else {
  cat("\nOriginal model: No valid predictions\n")
}

if(length(pred_best) > 0) {
  cat("\nBest model:\n")
  cat("  Prediction range: ", round(min(pred_best), 3), " - ", round(max(pred_best), 3), "\n", sep="")
  cat("  Prediction SD: ", round(sd(pred_best), 4), "\n", sep="")
  
  if(length(pred_orig) > 0 && sd(pred_orig) > 0) {
    cat("  SD improvement: ", round(sd(pred_best)/sd(pred_orig), 2), "x\n", sep="")
  }
} else {
  cat("\nBest model: No valid predictions\n")
}

# Create distribution plot
if(length(pred_orig) > 0 & length(pred_best) > 0) {
  dist_data <- data.frame(
    Prediction = c(pred_orig, pred_best),
    Model = factor(rep(c("Original", "Best"), c(length(pred_orig), length(pred_best))))
  )
  
  dist_plot <- ggplot(dist_data, aes(x = Prediction, fill = Model)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    scale_x_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Original" = "red", "Best" = "blue")) +
    labs(title = "Prediction Distribution: Original vs Best Model",
         x = "Predicted 180-day Mortality",
         y = "Count") +
    theme_bw(base_size = 12)
  
  print(dist_plot)
  ggsave("prediction_distributions.png", plot = dist_plot, width = 10, height = 6, dpi = 300)
} else {
  cat("Insufficient predictions for distribution plot\n")
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("MODEL IMPROVEMENT SUMMARY\n")
cat("============================================\n")

orig_cal <- cal_results[1, ]
best_cal <- cal_results[best_idx, ]

cat("\nCalibration Improvement:\n")
if(!is.na(orig_cal$Cal_Slope)) {
  cat("  Original slope: ", orig_cal$Cal_Slope, " (target = 1.0)\n", sep="")
} else {
  cat("  Original slope: NA\n")
}

if(!is.na(best_cal$Cal_Slope)) {
  cat("  Best slope: ", best_cal$Cal_Slope, " (target = 1.0)\n", sep="")
} else {
  cat("  Best slope: NA\n")
}

if(!is.na(orig_cal$Cal_Slope) & !is.na(best_cal$Cal_Slope)) {
  orig_error <- abs(orig_cal$Cal_Slope - 1)
  best_error <- abs(best_cal$Cal_Slope - 1)
  if(orig_error > 0) {
    improvement <- (orig_error - best_error) / orig_error * 100
    cat("  Improvement: ", round(improvement), "% closer to perfect calibration\n", sep="")
  }
}

cat("\nDiscrimination:\n")
if(!is.na(orig_cal$AUC)) {
  cat("  Original C-statistic: ", orig_cal$AUC, "\n", sep="")
} else {
  cat("  Original C-statistic: NA\n")
}

if(!is.na(best_cal$AUC)) {
  cat("  Best C-statistic: ", best_cal$AUC, "\n", sep="")
} else {
  cat("  Best C-statistic: NA\n")
}

cat("\nPrediction Spread:\n")
if(!is.na(orig_cal$Pred_SD)) {
  cat("  Original SD: ", orig_cal$Pred_SD, "\n", sep="")
} else {
  cat("  Original SD: NA\n")
}

if(!is.na(best_cal$Pred_SD)) {
  cat("  Best SD: ", best_cal$Pred_SD, "\n", sep="")
} else {
  cat("  Best SD: NA\n")
}

if(!is.na(orig_cal$Pred_SD) & !is.na(best_cal$Pred_SD) & orig_cal$Pred_SD > 0) {
  cat("  Improvement: ", round(best_cal$Pred_SD / orig_cal$Pred_SD, 1), "x wider spread\n", sep="")
}

cat("\nKey Improvements:\n")
if(grepl("Complications", best_model_name)) {
  cat("  ✓ Complication burden scores helped\n")
}
if(grepl("Interactions", best_model_name)) {
  cat("  ✓ Interaction terms were crucial\n")
}
if(grepl("Non-linear", best_model_name)) {
  cat("  ✓ Non-linear relationships captured\n")
}

cat("\nRECOMMENDATION:\n")
cat("Use ", best_model_name, " for best natural calibration\n", sep="")

# Save everything
models_list <- list()
if(!is.null(model1)) models_list$model1 <- model1
if(!is.null(model2)) models_list$model2 <- model2
if(!is.null(model3)) models_list$model3 <- model3
if(!is.null(model4)) models_list$model4 <- model4
if(!is.null(best_model)) models_list$best_model <- best_model
models_list$cal_results <- cal_results
if(!is.null(importance)) models_list$importance <- importance

save(models_list, file = "improved_models.RData")

write.csv(cal_results, "calibration_comparison.csv", row.names = FALSE)

if(!is.null(importance)) {
  write.csv(importance, "feature_importance.csv", row.names = FALSE)
}

cat("\nFiles saved:\n")
cat("  - improved_models.RData\n")
cat("  - calibration_comparison.csv\n")
if(!is.null(importance)) {
  cat("  - feature_importance.csv\n")
}
if(file.exists("calibration_original.png")) {
  cat("  - calibration_original.png\n")
}
if(file.exists("calibration_best.png")) {
  cat("  - calibration_best.png\n")
}
if(file.exists("prediction_distributions.png")) {
  cat("  - prediction_distributions.png\n")
}


#############################################
## CRE MANAGEMENT STRATEGY SURVIVAL ANALYSIS
## Comparing three management approaches
## 1 = No strategy, 2 = Targeted prophylaxis, 3 = Pre-emptive
#############################################

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

cat("============================================\n")
cat("CRE MANAGEMENT STRATEGY ANALYSIS\n")
cat("============================================\n")

# Load data
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}

data_analysis <- baseline_data

# ============================================================
# STEP 1: CREATE MANAGEMENT STRATEGY VARIABLE
# ============================================================

cat("\n=== Creating Management Strategy Variable ===\n")

# Check the checkbox variables
cat("\nCheckbox variables found:\n")
cat("cre_colonization___1 (No strategy):\n")
print(table(data_analysis$cre_colonization___1, useNA = "ifany"))

cat("\ncre_colonization___2 (Targeted prophylaxis):\n")
print(table(data_analysis$cre_colonization___2, useNA = "ifany"))

cat("\ncre_colonization___3 (Pre-emptive strategy):\n")
print(table(data_analysis$cre_colonization___3, useNA = "ifany"))

# Create composite management strategy variable
# Assuming these are mutually exclusive (only one strategy per patient)
data_analysis$cre_management_strategy <- NA

# Assign based on which checkbox is checked
data_analysis$cre_management_strategy[data_analysis$cre_colonization___1 == 1] <- 1
data_analysis$cre_management_strategy[data_analysis$cre_colonization___2 == 1] <- 2
data_analysis$cre_management_strategy[data_analysis$cre_colonization___3 == 1] <- 3

# Check for patients with multiple strategies (data quality check)
multiple_strategies <- rowSums(data_analysis[, c("cre_colonization___1", 
                                                 "cre_colonization___2", 
                                                 "cre_colonization___3")], na.rm = TRUE)
n_multiple <- sum(multiple_strategies > 1, na.rm = TRUE)
if(n_multiple > 0) {
  cat("\nWARNING: ", n_multiple, " patients have multiple strategies recorded\n", sep="")
  cat("These patients will use the highest strategy number\n")
}

# Create factor with meaningful labels
data_analysis$cre_strategy <- factor(
  data_analysis$cre_management_strategy,
  levels = c(1, 2, 3),
  labels = c("No Strategy", "Targeted Prophylaxis", "Pre-emptive Strategy")
)

cat("\n--- Management Strategy Distribution ---\n")
print(table(data_analysis$cre_strategy, useNA = "ifany"))

# ============================================================
# STEP 2: BASELINE CHARACTERISTICS BY STRATEGY
# ============================================================

cat("\n=== Baseline Characteristics by Management Strategy ===\n")

baseline_table <- data_analysis %>%
  group_by(cre_strategy) %>%
  summarise(
    n = n(),
    age_mean = round(mean(age_clean, na.rm = TRUE), 1),
    age_sd = round(sd(age_clean, na.rm = TRUE), 1),
    meld_mean = round(mean(meld_clean, na.rm = TRUE), 1),
    meld_sd = round(sd(meld_clean, na.rm = TRUE), 1),
    cre_infections = sum(cre_180),
    infection_rate = round(mean(cre_180) * 100, 1),
    deaths = sum(death_180),
    mortality_rate = round(mean(death_180) * 100, 1),
    .groups = "drop"
  )

cat("\n--- Patient Characteristics ---\n")
print(as.data.frame(baseline_table))

# Statistical comparisons
cat("\n--- Statistical Comparisons ---\n")

# Age comparison
if(sum(!is.na(data_analysis$cre_strategy)) > 2) {
  age_anova <- aov(age_clean ~ cre_strategy, data = data_analysis)
  cat("Age difference p-value: ", format.pval(summary(age_anova)[[1]][["Pr(>F)"]][1]), "\n", sep="")
  
  # MELD comparison
  meld_anova <- aov(meld_clean ~ cre_strategy, data = data_analysis)
  cat("MELD difference p-value: ", format.pval(summary(meld_anova)[[1]][["Pr(>F)"]][1]), "\n", sep="")
  
  # CRE infection rate comparison
  cre_chi <- chisq.test(table(data_analysis$cre_strategy, data_analysis$cre_180))
  cat("CRE infection rate difference p-value: ", format.pval(cre_chi$p.value), "\n", sep="")
}

# ============================================================
# STEP 3: UNADJUSTED SURVIVAL ANALYSIS
# ============================================================

cat("\n=== Unadjusted Survival Analysis ===\n")

# Remove missing strategy for survival analysis
data_complete <- data_analysis[!is.na(data_analysis$cre_strategy), ]

# Create survival object
surv_obj <- Surv(time = data_complete$time_simple, 
                 event = data_complete$death_180)

# Kaplan-Meier curves
km_fit <- survfit(surv_obj ~ cre_strategy, data = data_complete)

# Print survival at key time points
cat("\n--- Survival Probability at Key Time Points ---\n")
print(summary(km_fit, times = c(30, 60, 90, 120, 150, 180)))

# Log-rank test
survdiff_result <- survdiff(surv_obj ~ cre_strategy, data = data_complete)
log_rank_p <- 1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)

cat("\nLog-rank test p-value: ", format.pval(log_rank_p), "\n", sep="")

# Pairwise comparisons if significant
if(log_rank_p < 0.05) {
  cat("\n--- Pairwise Log-rank Tests ---\n")
  strategies <- levels(data_complete$cre_strategy)
  for(i in 1:(length(strategies)-1)) {
    for(j in (i+1):length(strategies)) {
      subset_data <- data_complete[data_complete$cre_strategy %in% c(strategies[i], strategies[j]), ]
      subset_surv <- Surv(time = subset_data$time_simple, event = subset_data$death_180)
      pairwise_test <- survdiff(subset_surv ~ cre_strategy, data = subset_data)
      p_val <- 1 - pchisq(pairwise_test$chisq, 1)
      cat(strategies[i], " vs ", strategies[j], ": p = ", format.pval(p_val), "\n", sep="")
    }
  }
}

# Create KM plot
km_plot <- ggsurvplot(
  km_fit,
  data = data_complete,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  xlab = "Days from Transplant",
  ylab = "Survival Probability",
  title = "Survival by CRE Management Strategy (Unadjusted)",
  legend.title = "Management Strategy",
  legend.labs = levels(data_complete$cre_strategy),
  palette = c("#E41A1C", "#377EB8", "#4DAF4A"),  # Red, Blue, Green
  risk.table = TRUE,
  risk.table.title = "Number at Risk",
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.3,
  ggtheme = theme_bw(base_size = 12),
  break.time.by = 30,
  xlim = c(0, 180)
)

print(km_plot)

# Save plot
ggsave("km_cre_management_strategy.png", 
       plot = km_plot$plot, 
       width = 10, height = 8, dpi = 300)

# ============================================================
# STEP 4: STRATIFIED BY CRE INFECTION
# ============================================================

cat("\n=== Stratified Analysis by CRE Infection ===\n")

# Among patients who developed CRE infection
data_cre_pos <- data_complete[data_complete$cre_180 == 1, ]

if(nrow(data_cre_pos) > 20) {
  cat("\n--- Among Patients with CRE Infection (n=", nrow(data_cre_pos), ") ---\n", sep="")
  
  strategy_cre_table <- table(data_cre_pos$cre_strategy)
  print(strategy_cre_table)
  
  if(length(unique(data_cre_pos$cre_strategy)) > 1) {
    surv_cre_pos <- Surv(time = data_cre_pos$time_simple, 
                         event = data_cre_pos$death_180)
    km_cre_pos <- survfit(surv_cre_pos ~ cre_strategy, data = data_cre_pos)
    
    survdiff_cre_pos <- survdiff(surv_cre_pos ~ cre_strategy, data = data_cre_pos)
    p_cre_pos <- 1 - pchisq(survdiff_cre_pos$chisq, length(survdiff_cre_pos$n) - 1)
    
    cat("Log-rank p-value among CRE+: ", format.pval(p_cre_pos), "\n", sep="")
    
    # Mortality rates
    mort_cre_pos <- data_cre_pos %>%
      group_by(cre_strategy) %>%
      summarise(
        n = n(),
        deaths = sum(death_180),
        mortality = round(mean(death_180) * 100, 1)
      )
    print(as.data.frame(mort_cre_pos))
  }
}

# ============================================================
# STEP 5: COX REGRESSION ANALYSIS
# ============================================================

cat("\n=== Cox Regression Analysis ===\n")

# Model 1: Strategy alone (unadjusted)
cat("\n--- Model 1: Management Strategy Only ---\n")
cox_strategy <- coxph(surv_obj ~ cre_strategy, data = data_complete)
print(summary(cox_strategy))

# Model 2: Strategy + basic adjustments
cat("\n--- Model 2: Adjusted for Age, MELD ---\n")
cox_adjusted <- coxph(surv_obj ~ cre_strategy + age_clean + meld_clean, 
                      data = data_complete)
print(summary(cox_adjusted))

# Model 3: Strategy + CRE infection
cat("\n--- Model 3: Strategy + CRE Infection ---\n")
cox_with_cre <- coxph(surv_obj ~ cre_strategy + cre_180, 
                      data = data_complete)
print(summary(cox_with_cre))

# Model 4: With interaction
cat("\n--- Model 4: Strategy × CRE Infection Interaction ---\n")
cox_interaction <- coxph(surv_obj ~ cre_strategy * cre_180, 
                         data = data_complete)
print(summary(cox_interaction))

# Compare models
cat("\n--- Model Comparison (AIC) ---\n")
cat("Model 1 (Strategy only): ", AIC(cox_strategy), "\n", sep="")
cat("Model 2 (+ Age, MELD): ", AIC(cox_adjusted), "\n", sep="")
cat("Model 3 (+ CRE infection): ", AIC(cox_with_cre), "\n", sep="")
cat("Model 4 (+ Interaction): ", AIC(cox_interaction), "\n", sep="")

# ============================================================
# STEP 6: INTEGRATE WITH BEST MODEL (Model 3 + Interactions)
# ============================================================

cat("\n=== Integration with Best Model (All Interactions) ===\n")

# Create necessary variables for Model 3
# Complication scores
major_compli_vars <- c("compli_1", "compli_3", "compli_5")
data_complete$major_complications <- 0
for(v in major_compli_vars) {
  if(v %in% names(data_complete)) {
    data_complete$major_complications <- data_complete$major_complications + data_complete[[v]]
  }
}

all_compli_vars <- names(data_complete)[grep("^compli_[0-9]", names(data_complete))]
data_complete$total_complications <- 0
for(v in all_compli_vars) {
  data_complete$total_complications <- data_complete$total_complications + data_complete[[v]]
}

# Create interaction terms
if(sum(!is.na(data_complete$age_clean)) > 30) {
  age_centered <- data_complete$age_clean - mean(data_complete$age_clean, na.rm=T)
  data_complete$cre_age_int <- data_complete$cre_180 * age_centered
}

if(sum(!is.na(data_complete$meld_clean)) > 30) {
  meld_centered <- data_complete$meld_clean - mean(data_complete$meld_clean, na.rm=T)
  data_complete$cre_meld_int <- data_complete$cre_180 * meld_centered
}

data_complete$cre_compli_int <- data_complete$cre_180 * data_complete$major_complications

if(sum(!is.na(data_complete$age_clean)) > 30 & sum(!is.na(data_complete$meld_clean)) > 30) {
  age_std <- scale(data_complete$age_clean)[,1]
  meld_std <- scale(data_complete$meld_clean)[,1]
  data_complete$age_meld_int <- age_std * meld_std
}

# Best model formula (Model 3)
best_formula <- surv_obj ~ cre_180 + age_clean + meld_clean + 
  major_complications + total_complications +
  cre_age_int + cre_meld_int + cre_compli_int + age_meld_int

# Fit best model without strategy
best_model_base <- coxph(best_formula, data = data_complete)

# Fit best model WITH strategy
best_model_strategy <- coxph(update(best_formula, ~ . + cre_strategy), 
                             data = data_complete)

# Compare
cat("\n--- Best Model Comparison ---\n")
cat("Best Model without strategy AIC: ", AIC(best_model_base), "\n", sep="")
cat("Best Model with strategy AIC: ", AIC(best_model_strategy), "\n", sep="")

# Likelihood ratio test
lr_test <- anova(best_model_base, best_model_strategy)
cat("Likelihood ratio test p-value: ", format.pval(lr_test$`P(>|Chi|)`[2]), "\n", sep="")

cat("\n--- Best Model with Management Strategy ---\n")
print(summary(best_model_strategy))

# ============================================================
# STEP 7: ADJUSTED SURVIVAL CURVES
# ============================================================

cat("\n=== Creating Adjusted Survival Curves ===\n")

# First, check what variables are in the model
model_vars <- names(best_model_strategy$coefficients)
cat("Model variables needed:\n")
print(model_vars)

# Create reference patient with ALL required variables
ref_patient <- data.frame(
  age_clean = mean(data_complete$age_clean, na.rm = TRUE),
  meld_clean = mean(data_complete$meld_clean, na.rm = TRUE),
  cre_180 = 0,  # No CRE infection
  major_complications = 0,
  total_complications = median(data_complete$total_complications, na.rm = TRUE)
)

# Add interaction terms based on the reference values
ref_patient$cre_age_int <- ref_patient$cre_180 * (ref_patient$age_clean - mean(data_complete$age_clean, na.rm=T))
ref_patient$cre_meld_int <- ref_patient$cre_180 * (ref_patient$meld_clean - mean(data_complete$meld_clean, na.rm=T))
ref_patient$cre_compli_int <- ref_patient$cre_180 * ref_patient$major_complications

# Calculate age_meld_int properly
age_mean <- mean(data_complete$age_clean, na.rm=T)
age_sd <- sd(data_complete$age_clean, na.rm=T)
meld_mean <- mean(data_complete$meld_clean, na.rm=T)
meld_sd <- sd(data_complete$meld_clean, na.rm=T)
ref_patient$age_meld_int <- ((ref_patient$age_clean - age_mean)/age_sd) * ((ref_patient$meld_clean - meld_mean)/meld_sd)

# Create three scenarios (one for each strategy)
no_strategy <- ref_patient
no_strategy$cre_strategy <- factor("No Strategy", levels = levels(data_complete$cre_strategy))

targeted <- ref_patient
targeted$cre_strategy <- factor("Targeted Prophylaxis", levels = levels(data_complete$cre_strategy))

preemptive <- ref_patient
preemptive$cre_strategy <- factor("Pre-emptive Strategy", levels = levels(data_complete$cre_strategy))

# Check that all variables are present
cat("\nChecking newdata completeness:\n")
cat("Variables in no_strategy: ", names(no_strategy), "\n")
cat("Any NAs in no_strategy: ", any(is.na(no_strategy)), "\n")

# Get survival curves with error handling
tryCatch({
  surv_no_strategy <- survfit(best_model_strategy, newdata = no_strategy)
  surv_targeted <- survfit(best_model_strategy, newdata = targeted)
  surv_preemptive <- survfit(best_model_strategy, newdata = preemptive)
  
  # Create plot data
  plot_data <- data.frame(
    time = c(surv_no_strategy$time, surv_targeted$time, surv_preemptive$time),
    survival = c(surv_no_strategy$surv, surv_targeted$surv, surv_preemptive$surv),
    lower = c(surv_no_strategy$lower, surv_targeted$lower, surv_preemptive$lower),
    upper = c(surv_no_strategy$upper, surv_targeted$upper, surv_preemptive$upper),
    strategy = factor(c(rep("No Strategy", length(surv_no_strategy$time)),
                        rep("Targeted Prophylaxis", length(surv_targeted$time)),
                        rep("Pre-emptive Strategy", length(surv_preemptive$time))),
                      levels = c("No Strategy", "Targeted Prophylaxis", "Pre-emptive Strategy"))
  )
  
  # Create adjusted plot
  adj_plot <- ggplot(plot_data, aes(x = time, y = survival, color = strategy)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = strategy), alpha = 0.2) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    labs(title = "Adjusted Survival Curves by CRE Management Strategy",
         subtitle = "Adjusted for age, MELD, complications, CRE infection, and all interactions",
         x = "Days from Transplant",
         y = "Survival Probability",
         color = "Strategy",
         fill = "Strategy") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, 30)) +
    scale_y_continuous(limits = c(0.5, 1))
  
  print(adj_plot)
  
  ggsave("adjusted_survival_strategies.png", plot = adj_plot, 
         width = 10, height = 8, dpi = 300)
  
}, error = function(e) {
  cat("\nError creating adjusted survival curves: ", e$message, "\n")
  cat("Skipping adjusted curves visualization\n")
})

# ============================================================
# STEP 8: RISK PREDICTIONS
# ============================================================

cat("\n=== Risk Predictions by Patient Profile ===\n")

# Get means and SDs for standardization
age_mean <- mean(data_complete$age_clean, na.rm=T)
age_sd <- sd(data_complete$age_clean, na.rm=T)
meld_mean <- mean(data_complete$meld_clean, na.rm=T)
meld_sd <- sd(data_complete$meld_clean, na.rm=T)

# Create scenarios
scenarios <- expand.grid(
  age = c(45, 60),
  meld = c(15, 30),
  cre_infection = c(0, 1),
  strategy = levels(data_complete$cre_strategy)
)

# Calculate 180-day mortality risk
scenarios$mortality_risk <- NA

for(i in 1:nrow(scenarios)) {
  patient <- data.frame(
    age_clean = scenarios$age[i],
    meld_clean = scenarios$meld[i],
    cre_180 = scenarios$cre_infection[i],
    major_complications = 0,
    total_complications = 0,
    cre_strategy = factor(scenarios$strategy[i], levels = levels(data_complete$cre_strategy))
  )
  
  # Calculate interactions using the same centering as in model building
  patient$cre_age_int <- patient$cre_180 * (patient$age_clean - age_mean)
  patient$cre_meld_int <- patient$cre_180 * (patient$meld_clean - meld_mean)
  patient$cre_compli_int <- patient$cre_180 * patient$major_complications
  patient$age_meld_int <- ((patient$age_clean - age_mean)/age_sd) * ((patient$meld_clean - meld_mean)/meld_sd)
  
  # Get survival at 180 days with error handling
  tryCatch({
    surv_pred <- survfit(best_model_strategy, newdata = patient)
    idx <- which.min(abs(surv_pred$time - 180))
    if(length(idx) > 0) {
      scenarios$mortality_risk[i] <- round((1 - surv_pred$surv[idx]) * 100, 1)
    }
  }, error = function(e) {
    cat("Warning: Could not predict for scenario ", i, "\n")
  })
}

# Display results if predictions succeeded
if(sum(!is.na(scenarios$mortality_risk)) > 0) {
  cat("\n--- 180-Day Mortality Risk (%) ---\n")
  
  # Without CRE
  cat("\nWithout CRE Infection:\n")
  no_cre <- scenarios[scenarios$cre_infection == 0, ]
  no_cre_wide <- tidyr::pivot_wider(no_cre, 
                                    id_cols = c(age, meld),
                                    names_from = strategy,
                                    values_from = mortality_risk)
  print(as.data.frame(no_cre_wide))
  
  # With CRE
  cat("\nWith CRE Infection:\n")
  with_cre <- scenarios[scenarios$cre_infection == 1, ]
  with_cre_wide <- tidyr::pivot_wider(with_cre,
                                      id_cols = c(age, meld),
                                      names_from = strategy,
                                      values_from = mortality_risk)
  print(as.data.frame(with_cre_wide))
  
  # Calculate absolute risk reductions
  cat("\n--- Absolute Risk Reduction (%) ---\n")
  
  # For each scenario, calculate ARR compared to no strategy
  for(age_val in unique(scenarios$age)) {
    for(meld_val in unique(scenarios$meld)) {
      for(cre_val in unique(scenarios$cre_infection)) {
        no_strat <- scenarios$mortality_risk[scenarios$age == age_val & 
                                               scenarios$meld == meld_val & 
                                               scenarios$cre_infection == cre_val & 
                                               scenarios$strategy == "No Strategy"]
        targeted <- scenarios$mortality_risk[scenarios$age == age_val & 
                                               scenarios$meld == meld_val & 
                                               scenarios$cre_infection == cre_val & 
                                               scenarios$strategy == "Targeted Prophylaxis"]
        preempt <- scenarios$mortality_risk[scenarios$age == age_val & 
                                              scenarios$meld == meld_val & 
                                              scenarios$cre_infection == cre_val & 
                                              scenarios$strategy == "Pre-emptive Strategy"]
        
        if(!is.na(no_strat) && !is.na(targeted) && !is.na(preempt)) {
          cat("\nAge ", age_val, ", MELD ", meld_val, ", CRE ", 
              ifelse(cre_val == 1, "Yes", "No"), ":\n", sep="")
          cat("  ARR with Targeted: ", round(no_strat - targeted, 1), "%\n", sep="")
          cat("  ARR with Pre-emptive: ", round(no_strat - preempt, 1), "%\n", sep="")
          
          # Calculate NNT if ARR > 0
          if((no_strat - targeted) > 0) {
            nnt_targeted <- round(100 / (no_strat - targeted))
            cat("  NNT (Targeted): ", nnt_targeted, "\n", sep="")
          }
          if((no_strat - preempt) > 0) {
            nnt_preempt <- round(100 / (no_strat - preempt))
            cat("  NNT (Pre-emptive): ", nnt_preempt, "\n", sep="")
          }
        }
      }
    }
  }
} else {
  cat("\nWarning: Risk predictions could not be calculated\n")
  cat("This may be due to missing values or model convergence issues\n")
}

# ============================================================
# STEP 9: SUMMARY AND CLINICAL IMPLICATIONS
# ============================================================

cat("\n============================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("============================================\n")

# Extract key hazard ratios from best model
coef_table <- summary(best_model_strategy)$coefficients

cat("\n--- Key Hazard Ratios (from Best Model) ---\n")

# Targeted prophylaxis vs No strategy
if("cre_strategyTargeted Prophylaxis" %in% rownames(coef_table)) {
  hr_targeted <- exp(coef_table["cre_strategyTargeted Prophylaxis", "coef"])
  ci_targeted_low <- exp(coef_table["cre_strategyTargeted Prophylaxis", "coef"] - 
                           1.96 * coef_table["cre_strategyTargeted Prophylaxis", "se(coef)"])
  ci_targeted_high <- exp(coef_table["cre_strategyTargeted Prophylaxis", "coef"] + 
                            1.96 * coef_table["cre_strategyTargeted Prophylaxis", "se(coef)"])
  p_targeted <- coef_table["cre_strategyTargeted Prophylaxis", "Pr(>|z|)"]
  
  cat("Targeted Prophylaxis vs No Strategy:\n")
  cat("  HR = ", round(hr_targeted, 2), " (95% CI: ", round(ci_targeted_low, 2), 
      "-", round(ci_targeted_high, 2), "), p = ", format.pval(p_targeted), "\n", sep="")
}

# Pre-emptive vs No strategy
if("cre_strategyPre-emptive Strategy" %in% rownames(coef_table)) {
  hr_preemptive <- exp(coef_table["cre_strategyPre-emptive Strategy", "coef"])
  ci_preemptive_low <- exp(coef_table["cre_strategyPre-emptive Strategy", "coef"] - 
                             1.96 * coef_table["cre_strategyPre-emptive Strategy", "se(coef)"])
  ci_preemptive_high <- exp(coef_table["cre_strategyPre-emptive Strategy", "coef"] + 
                              1.96 * coef_table["cre_strategyPre-emptive Strategy", "se(coef)"])
  p_preemptive <- coef_table["cre_strategyPre-emptive Strategy", "Pr(>|z|)"]
  
  cat("Pre-emptive Strategy vs No Strategy:\n")
  cat("  HR = ", round(hr_preemptive, 2), " (95% CI: ", round(ci_preemptive_low, 2), 
      "-", round(ci_preemptive_high, 2), "), p = ", format.pval(p_preemptive), "\n", sep="")
}

cat("\n--- Clinical Implications ---\n")
cat("1. Management strategies show ", 
    ifelse(log_rank_p < 0.05, "significant", "no significant"), 
    " differences in survival (p = ", format.pval(log_rank_p), ")\n", sep="")

cat("2. After adjusting for patient characteristics and CRE infection:\n")
if(exists("hr_targeted") && hr_targeted < 1) {
  cat("   - Targeted prophylaxis reduces mortality by ", 
      round((1 - hr_targeted) * 100), "%\n", sep="")
}
if(exists("hr_preemptive") && hr_preemptive < 1) {
  cat("   - Pre-emptive strategy reduces mortality by ", 
      round((1 - hr_preemptive) * 100), "%\n", sep="")
}

cat("\n3. Strategy effectiveness varies by patient risk profile\n")
cat("   See risk prediction tables above for personalized estimates\n")

# Save all results
save(data_complete, 
     cox_strategy, cox_adjusted, cox_with_cre, cox_interaction,
     best_model_base, best_model_strategy,
     scenarios,
     file = "cre_management_strategy_results.RData")

cat("\nAnalysis complete. Results saved to 'cre_management_strategy_results.RData'\n")

#############################################
## CRE MANAGEMENT STRATEGY SURVIVAL ANALYSIS
## STRATIFIED BY COLLECTION METHOD
## Variable: retro_or_pros
## 
## FULLY ROBUST VERSION
## - All variables properly initialized
## - No conditional statements on potentially undefined variables
## - Manual likelihood ratio test calculation
## - Comprehensive error handling
##
## EXPECTED BEHAVIOR:
## - Will attempt stratified analysis if retro_or_pros exists
## - Falls back to overall analysis if stratification fails
## - Continues even if some models don't converge
## - Saves whatever results were successfully generated
#############################################

# Load libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tidyr)

cat("============================================\n")
cat("CRE MANAGEMENT STRATEGY ANALYSIS\n")
cat("STRATIFIED BY DATA COLLECTION METHOD\n")
cat("Version: Fully Robust with Error Handling\n")
cat("============================================\n")

# Load data
if(!exists("baseline_data")) {
  load("baseline_data_final.RData")
}
data_analysis <- baseline_data

# Initialize flags and variables
skip_stratified <- FALSE
p_interaction <- NA  # Always initialize
p_retro <- NA  # Always initialize
p_prosp <- NA  # Always initialize
data_retro <- NULL  # Initialize
data_prosp <- NULL  # Initialize

# ============================================================
# STEP 1: CREATE MANAGEMENT STRATEGY VARIABLE
# ============================================================

cat("\n=== Creating Management Strategy Variable ===\n")

# Check the checkbox variables
cat("\nCheckbox variables found:\n")
cat("cre_colonization___1 (No strategy):\n")
print(table(data_analysis$cre_colonization___1, useNA = "ifany"))
cat("\ncre_colonization___2 (Targeted prophylaxis):\n")
print(table(data_analysis$cre_colonization___2, useNA = "ifany"))
cat("\ncre_colonization___3 (Pre-emptive strategy):\n")
print(table(data_analysis$cre_colonization___3, useNA = "ifany"))

# Create composite management strategy variable
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
  cat("\nWARNING: ", n_multiple, " patients have multiple strategies recorded\n", sep="")
}

# Create factor with meaningful labels
data_analysis$cre_strategy <- factor(
  data_analysis$cre_management_strategy,
  levels = c(1, 2, 3),
  labels = c("No Strategy", "Targeted Prophylaxis", "Pre-emptive Strategy")
)

cat("\n--- Management Strategy Distribution ---\n")
print(table(data_analysis$cre_strategy, useNA = "ifany"))

# ============================================================
# STEP 2: CHECK AND PROCESS COLLECTION METHOD VARIABLE
# ============================================================

cat("\n=== Data Collection Method ===\n")

# Check if retro_or_pros exists
if(!"retro_or_pros" %in% names(data_analysis)) {
  cat("\nWARNING: 'retro_or_pros' variable not found\n")
  cat("Proceeding with unstratified analysis only\n")
  skip_stratified <- TRUE
  data_analysis$collection_method <- NA
} else {
  # Check the variable
  cat("\n--- Collection Method Distribution (raw values) ---\n")
  print(table(data_analysis$retro_or_pros, useNA = "ifany"))
  
  unique_vals <- unique(data_analysis$retro_or_pros)
  cat("\nUnique values in retro_or_pros: ", paste(unique_vals, collapse=", "), "\n")
  
  if(all(is.na(data_analysis$retro_or_pros))) {
    cat("\nAll values are NA - cannot stratify\n")
    skip_stratified <- TRUE
    data_analysis$collection_method <- NA
  } else {
    # Process the variable
    valid_vals <- unique_vals[!is.na(unique_vals)]
    
    # Convert to character for safer comparison
    valid_vals_char <- as.character(valid_vals)
    
    # Try to determine coding
    if(all(tolower(valid_vals_char) %in% c("retro", "pros"))) {
      cat("Detected text coding (retro/pros)\n")
      retro_val <- valid_vals[tolower(valid_vals_char) == "retro"][1]
      pros_val <- valid_vals[tolower(valid_vals_char) == "pros"][1]
      if(!is.na(retro_val) && !is.na(pros_val)) {
        data_analysis$collection_method <- factor(
          data_analysis$retro_or_pros,
          levels = c(retro_val, pros_val),
          labels = c("Retrospective", "Prospective")
        )
      } else {
        cat("Error: Could not properly map retro/pros values\n")
        skip_stratified <- TRUE
        data_analysis$collection_method <- NA
      }
    } else if(all(valid_vals_char %in% c("0", "1"))) {
      cat("Detected binary coding (0/1)\n")
      data_analysis$collection_method <- factor(
        data_analysis$retro_or_pros,
        levels = c(0, 1),
        labels = c("Retrospective", "Prospective")
      )
    } else if(all(valid_vals_char %in% c("1", "2"))) {
      cat("Detected binary coding (1/2)\n")
      data_analysis$collection_method <- factor(
        data_analysis$retro_or_pros,
        levels = c(1, 2),
        labels = c("Retrospective", "Prospective")
      )
    } else {
      cat("Unknown coding - using as is\n")
      cat("Values found:", paste(valid_vals_char, collapse=", "), "\n")
      data_analysis$collection_method <- factor(data_analysis$retro_or_pros)
    }
  }
}

cat("\n--- Collection Method Distribution (cleaned) ---\n")
if(!skip_stratified) {
  print(table(data_analysis$collection_method, useNA = "ifany"))
} else {
  cat("Collection method not available\n")
}

# ============================================================
# STEP 3: BASELINE CHARACTERISTICS
# ============================================================

cat("\n=== Baseline Characteristics ===\n")

if(!skip_stratified) {
  # Stratified summary
  baseline_stratified <- data_analysis %>%
    group_by(collection_method, cre_strategy) %>%
    summarise(
      n = n(),
      age_mean = round(mean(age_clean, na.rm = TRUE), 1),
      meld_mean = round(mean(meld_clean, na.rm = TRUE), 1),
      cre_infections = sum(cre_180, na.rm = TRUE),
      infection_rate = round(mean(cre_180, na.rm = TRUE) * 100, 1),
      deaths = sum(death_180, na.rm = TRUE),
      mortality_rate = round(mean(death_180, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
  cat("\n--- Patient Characteristics (Stratified) ---\n")
  print(as.data.frame(baseline_stratified))
} else {
  # Unstratified summary
  baseline_summary <- data_analysis %>%
    group_by(cre_strategy) %>%
    summarise(
      n = n(),
      age_mean = round(mean(age_clean, na.rm = TRUE), 1),
      meld_mean = round(mean(meld_clean, na.rm = TRUE), 1),
      cre_infections = sum(cre_180, na.rm = TRUE),
      infection_rate = round(mean(cre_180, na.rm = TRUE) * 100, 1),
      deaths = sum(death_180, na.rm = TRUE),
      mortality_rate = round(mean(death_180, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
  cat("\n--- Patient Characteristics (Unstratified) ---\n")
  print(as.data.frame(baseline_summary))
}

# ============================================================
# STEP 4: SURVIVAL ANALYSIS
# ============================================================

# Remove missing strategy
data_complete <- data_analysis[!is.na(data_analysis$cre_strategy), ]

if(nrow(data_complete) < 10) {
  cat("\n=== ERROR: Insufficient data for analysis ===\n")
  cat("Only", nrow(data_complete), "patients with valid strategy data\n")
  cat("Analysis cannot proceed.\n")
  stop("Insufficient data")
}

# Overall survival analysis
cat("\n=== Overall Survival Analysis ===\n")
cat("Analyzing", nrow(data_complete), "patients\n")
surv_obj <- Surv(time = data_complete$time_simple, event = data_complete$death_180)
km_fit <- survfit(surv_obj ~ cre_strategy, data = data_complete)

# Log-rank test
log_rank_p <- NA  # Initialize
tryCatch({
  survdiff_result <- survdiff(surv_obj ~ cre_strategy, data = data_complete)
  log_rank_p <- 1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)
  cat("Overall log-rank test p-value: ", format.pval(log_rank_p), "\n")
}, error = function(e) {
  cat("Error calculating log-rank test:", e$message, "\n")
})

# Create overall KM plot
tryCatch({
  km_plot_overall <- ggsurvplot(
    km_fit,
    data = data_complete,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Days from Transplant",
    ylab = "Survival Probability",
    title = "Overall Survival by CRE Management Strategy",
    legend.title = "Strategy",
    palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
    risk.table = TRUE,
    ggtheme = theme_bw(),
    xlim = c(0, 180)
  )
  print(km_plot_overall)
  ggsave("km_overall_strategy.png", plot = km_plot_overall$plot, 
         width = 10, height = 8, dpi = 300)
}, error = function(e) {
  cat("Error creating overall KM plot:", e$message, "\n")
})

if(!skip_stratified) {
  # RETROSPECTIVE ANALYSIS
  cat("\n=== Retrospective Data Analysis ===\n")
  data_retro <- data_complete[!is.na(data_complete$collection_method) & 
                                data_complete$collection_method == "Retrospective", ]
  
  if(nrow(data_retro) > 20 && length(unique(data_retro$cre_strategy)) > 1) {
    cat("N =", nrow(data_retro), "\n")
    cat("Strategies present:", paste(unique(data_retro$cre_strategy), collapse=", "), "\n")
    
    surv_retro <- Surv(time = data_retro$time_simple, event = data_retro$death_180)
    km_retro <- survfit(surv_retro ~ cre_strategy, data = data_retro)
    
    tryCatch({
      survdiff_retro <- survdiff(surv_retro ~ cre_strategy, data = data_retro)
      p_retro <- 1 - pchisq(survdiff_retro$chisq, length(survdiff_retro$n) - 1)
      cat("Retrospective log-rank p-value: ", format.pval(p_retro), "\n")
    }, error = function(e) {
      cat("Error calculating retrospective log-rank test:", e$message, "\n")
      p_retro <- NA
    })
    
    tryCatch({
      km_plot_retro <- ggsurvplot(
        km_retro,
        data = data_retro,
        pval = TRUE,
        conf.int = TRUE,
        xlab = "Days from Transplant",
        ylab = "Survival Probability",
        title = "Retrospective: Survival by CRE Management Strategy",
        legend.title = "Strategy",
        palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
        risk.table = TRUE,
        ggtheme = theme_bw(),
        xlim = c(0, 180)
      )
      print(km_plot_retro)
      ggsave("km_retrospective.png", plot = km_plot_retro$plot,
             width = 10, height = 8, dpi = 300)
    }, error = function(e) {
      cat("Error creating retrospective KM plot:", e$message, "\n")
    })
  } else {
    if(nrow(data_retro) <= 20) {
      cat("Insufficient retrospective data (n=", nrow(data_retro), ")\n")
    } else {
      cat("Insufficient variation in strategies for retrospective analysis\n")
    }
  }
  
  # PROSPECTIVE ANALYSIS
  cat("\n=== Prospective Data Analysis ===\n")
  data_prosp <- data_complete[!is.na(data_complete$collection_method) & 
                                data_complete$collection_method == "Prospective", ]
  
  if(nrow(data_prosp) > 20 && length(unique(data_prosp$cre_strategy)) > 1) {
    cat("N =", nrow(data_prosp), "\n")
    cat("Strategies present:", paste(unique(data_prosp$cre_strategy), collapse=", "), "\n")
    
    surv_prosp <- Surv(time = data_prosp$time_simple, event = data_prosp$death_180)
    km_prosp <- survfit(surv_prosp ~ cre_strategy, data = data_prosp)
    
    tryCatch({
      survdiff_prosp <- survdiff(surv_prosp ~ cre_strategy, data = data_prosp)
      p_prosp <- 1 - pchisq(survdiff_prosp$chisq, length(survdiff_prosp$n) - 1)
      cat("Prospective log-rank p-value: ", format.pval(p_prosp), "\n")
    }, error = function(e) {
      cat("Error calculating prospective log-rank test:", e$message, "\n")
      p_prosp <- NA
    })
    
    tryCatch({
      km_plot_prosp <- ggsurvplot(
        km_prosp,
        data = data_prosp,
        pval = TRUE,
        conf.int = TRUE,
        xlab = "Days from Transplant",
        ylab = "Survival Probability",
        title = "Prospective: Survival by CRE Management Strategy",
        legend.title = "Strategy",
        palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
        risk.table = TRUE,
        ggtheme = theme_bw(),
        xlim = c(0, 180)
      )
      print(km_plot_prosp)
      ggsave("km_prospective.png", plot = km_plot_prosp$plot,
             width = 10, height = 8, dpi = 300)
    }, error = function(e) {
      cat("Error creating prospective KM plot:", e$message, "\n")
    })
  } else {
    if(nrow(data_prosp) <= 20) {
      cat("Insufficient prospective data (n=", nrow(data_prosp), ")\n")
    } else {
      cat("Insufficient variation in strategies for prospective analysis\n")
    }
  }
}

# ============================================================
# STEP 5: COX REGRESSION MODELS
# ============================================================

cat("\n=== Cox Regression Models ===\n")

# Overall Cox model
cat("\n--- Overall Model (Adjusted) ---\n")
tryCatch({
  cox_overall <- coxph(surv_obj ~ cre_strategy + age_clean + meld_clean + cre_180, 
                       data = data_complete)
  print(summary(cox_overall))
}, error = function(e) {
  cat("Error fitting overall Cox model:", e$message, "\n")
  cat("Trying simpler model...\n")
  tryCatch({
    cox_overall <- coxph(surv_obj ~ cre_strategy, data = data_complete)
    print(summary(cox_overall))
  }, error = function(e2) {
    cat("Error fitting simple Cox model:", e2$message, "\n")
  })
})

if(!skip_stratified) {
  # Test for interaction
  cat("\n--- Testing Strategy x Collection Method Interaction ---\n")
  
  # Initialize p_interaction outside tryCatch to ensure it exists
  p_interaction <- NA
  
  if(sum(!is.na(data_complete$collection_method)) > 30) {
    # Try to fit models and test interaction
    test_succeeded <- FALSE
    
    tryCatch({
      cox_no_int <- coxph(surv_obj ~ cre_strategy + collection_method + 
                            age_clean + meld_clean + cre_180, 
                          data = data_complete)
      
      cox_with_int <- coxph(surv_obj ~ cre_strategy * collection_method + 
                              age_clean + meld_clean + cre_180, 
                            data = data_complete)
      
      # Use logLik to manually calculate LR test
      ll_no_int <- as.numeric(logLik(cox_no_int))
      ll_with_int <- as.numeric(logLik(cox_with_int))
      lr_stat <- 2 * (ll_with_int - ll_no_int)
      df_diff <- length(coef(cox_with_int)) - length(coef(cox_no_int))
      
      # Calculate p-value
      if(df_diff > 0 && lr_stat >= 0 && !is.na(lr_stat) && !is.na(df_diff)) {
        p_interaction <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
        test_succeeded <- TRUE
        
        cat("Likelihood ratio statistic:", round(lr_stat, 3), "\n")
        cat("Degrees of freedom:", df_diff, "\n")
        cat("Interaction p-value: ", format.pval(p_interaction), "\n")
        
        # Only check p-value if it's valid
        if(!is.na(p_interaction) && is.numeric(p_interaction) && length(p_interaction) == 1) {
          if(p_interaction < 0.05) {
            cat("SIGNIFICANT INTERACTION detected\n")
            cat("Strategy effectiveness differs by collection method\n")
          } else {
            cat("No significant interaction\n")
            cat("Strategy effectiveness similar across collection methods\n")
          }
        }
      } else {
        cat("Could not calculate valid test statistics\n")
        p_interaction <- NA
      }
    }, error = function(e) {
      cat("Error testing interaction:", e$message, "\n")
      p_interaction <- NA
      test_succeeded <- FALSE
    })
    
    if(!test_succeeded) {
      cat("Interaction test did not complete successfully\n")
      p_interaction <- NA
    }
  } else {
    cat("Insufficient data to test interaction\n")
    p_interaction <- NA
  }
  
  # Stratified Cox models if enough data
  if(exists("data_retro") && nrow(data_retro) > 30 && 
     length(unique(data_retro$cre_strategy)) > 1) {
    cat("\n--- Retrospective Cox Model ---\n")
    tryCatch({
      surv_retro_cox <- Surv(time = data_retro$time_simple, event = data_retro$death_180)
      cox_retro <- coxph(surv_retro_cox ~ cre_strategy + age_clean + meld_clean, 
                         data = data_retro)
      print(summary(cox_retro))
    }, error = function(e) {
      cat("Error fitting retrospective Cox model:", e$message, "\n")
    })
  }
  
  if(exists("data_prosp") && nrow(data_prosp) > 30 && 
     length(unique(data_prosp$cre_strategy)) > 1) {
    cat("\n--- Prospective Cox Model ---\n")
    tryCatch({
      surv_prosp_cox <- Surv(time = data_prosp$time_simple, event = data_prosp$death_180)
      cox_prosp <- coxph(surv_prosp_cox ~ cre_strategy + age_clean + meld_clean, 
                         data = data_prosp)
      print(summary(cox_prosp))
    }, error = function(e) {
      cat("Error fitting prospective Cox model:", e$message, "\n")
    })
  }
} else {
  cat("\n--- Stratified Cox models not performed ---\n")
  p_interaction <- NA  # Ensure it's defined even if skip_stratified is TRUE
}

# ============================================================
# STEP 6: SUMMARY
# ============================================================

cat("\n============================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("============================================\n")

cat("\n1. Overall Analysis:\n")
cat("   - Total N: ", nrow(data_complete), "\n")
if(!is.na(log_rank_p)) {
  cat("   - Overall log-rank p: ", format.pval(log_rank_p), "\n")
} else {
  cat("   - Overall log-rank p: Could not calculate\n")
}

if(!skip_stratified) {
  cat("\n2. Stratified Analysis:\n")
  if(exists("data_retro") && !is.null(data_retro)) {
    cat("   - Retrospective N: ", nrow(data_retro), "\n")
    if(exists("p_retro") && !is.na(p_retro)) {
      cat("   - Retrospective p: ", format.pval(p_retro), "\n")
    }
  }
  if(exists("data_prosp") && !is.null(data_prosp)) {
    cat("   - Prospective N: ", nrow(data_prosp), "\n")
    if(exists("p_prosp") && !is.na(p_prosp)) {
      cat("   - Prospective p: ", format.pval(p_prosp), "\n")
    }
  }
  if(exists("p_interaction") && !is.na(p_interaction)) {
    cat("\n3. Interaction Test:\n")
    cat("   - p-value: ", format.pval(p_interaction), "\n")
  }
} else {
  cat("\n2. Stratified analysis not performed (collection method unavailable)\n")
}

cat("\nAnalysis complete.\n")
cat("Plots saved:\n")
if(file.exists("km_overall_strategy.png")) {
  cat("  - km_overall_strategy.png\n")
}
if(!skip_stratified) {
  if(file.exists("km_retrospective.png")) {
    cat("  - km_retrospective.png\n")
  }
  if(file.exists("km_prospective.png")) {
    cat("  - km_prospective.png\n")
  }
}

# Save results
tryCatch({
  save.image(file = "cre_stratified_analysis_results.RData")
  cat("\nResults saved to cre_stratified_analysis_results.RData\n")
}, error = function(e) {
  cat("\nError saving results:", e$message, "\n")
})

