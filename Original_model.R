# Step 1: Check your data structure first
str(CRE3)
summary(CRE3)

# Step 2: Specifically check the variables used in your model
# Check if time variable is numeric
class(CRE3$time)
summary(CRE3$time)

# Check if status variable is appropriate
class(CRE3$status)
table(CRE3$status)

# Check all predictor variables
predictors <- c("reint", "mv", "arf", "crepre_60", "crepost_60", "multipost")
for(var in predictors) {
  cat("Variable:", var, "\n")
  cat("Class:", class(CRE3[[var]]), "\n")
  print(table(CRE3[[var]], useNA = "always"))
  cat("\n")
}

# Step 3: Clean and prepare your data
library(dplyr)

# Create a cleaned dataset
CRE3_clean <- CRE3 %>%
  mutate(
    # Ensure time is numeric
    time = as.numeric(time),
    
    # Ensure status is numeric (0/1 or 1/2)
    status = as.numeric(status),
    
    # Ensure all predictors are numeric (assuming they're binary 0/1)
    reint = as.numeric(reint),
    mv = as.numeric(mv),
    arf = as.numeric(arf),
    crepre_60 = as.numeric(crepre_60),
    crepost_60 = as.numeric(crepost_60),
    multipost = as.numeric(multipost)
  ) %>%
  # Remove rows with missing values in key variables
  filter(
    !is.na(time),
    !is.na(status),
    time > 0  # Ensure positive survival times
  )

# Step 4: Check for any remaining issues
summary(CRE3_clean[, c("time", "status", predictors)])

# Step 5: Rebuild your model with cleaned data
library(survival)
library(regplot)
library(DynNom)
library(shiny)

# Fit the Cox model
f.csc <- coxph(Surv(time, status==1) ~ reint + mv + arf + crepre_60 + crepost_60 + multipost, 
               data = CRE3_clean)

# Check model summary
summary(f.csc)

# Step 6: Create the dynamic nomogram with error handling
tryCatch({
  # First try with original function
  DynNom(f.csc, 
         data = CRE3_clean, 
         m.summary = "formatted", 
         covariate = "slider", 
         DNtitle = "Move slider to (0) when risk factor is absent; Move slider to (1) when risk factor is present",
         DNxlab = "Predicted CRE infection probability +/-95% CI", 
         DNlimits = c(0, 180), 
         DNylab = "Group",
         KMtitle = "CRE infection probability", 
         KMxlab = "Days post transplant", 
         KMylab = "Probability of CRE infection", 
         ptype = "1-st")
}, error = function(e) {
  cat("Error with DynNom:", e$message, "\n")
  
  # Try alternative approach with different parameters
  cat("Trying alternative parameters...\n")
  
  DynNom(f.csc, 
         data = CRE3_clean, 
         covariate = "slider",
         DNtitle = "CRE Infection Risk Calculator",
         ptype = "1-st")
})

# Step 7: Build the app
tryCatch({
  DNbuilder(f.csc, 
            data = CRE3_clean, 
            m.summary = "formatted", 
            covariate = "slider", 
            DNtitle = "Move slider to (0) when risk factor is absent; Move slider to (1) when risk factor is present",
            DNxlab = "Predicted CRE infection probability +/-95% CI", 
            DNlimits = c(0, 180), 
            DNylab = "Group",
            KMtitle = "CRE infection probability", 
            KMxlab = "Days post transplant", 
            KMylab = "Probability of CRE infection", 
            ptype = "1-st")
}, error = function(e) {
  cat("Error with DNbuilder:", e$message, "\n")
  
  # Simple version
  DNbuilder(f.csc, 
            data = CRE3_clean,
            covariate = "slider")
})

# Step 8: Alternative approach if DynNom still fails
# Sometimes you need to ensure factor variables are properly coded
CRE3_alt <- CRE3_clean %>%
  mutate(
    # Convert to factors with explicit levels if they're binary
    reint = factor(reint, levels = c(0, 1)),
    mv = factor(mv, levels = c(0, 1)),
    arf = factor(arf, levels = c(0, 1)),
    crepre_60 = factor(crepre_60, levels = c(0, 1)),
    crepost_60 = factor(crepost_60, levels = c(0, 1)),
    multipost = factor(multipost, levels = c(0, 1))
  )

# Try with factor variables
f.csc_alt <- coxph(Surv(time, status==1) ~ reint + mv + arf + crepre_60 + crepost_60 + multipost, 
                   data = CRE3_alt)

# Test this version
tryCatch({
  DynNom(f.csc_alt, data = CRE3_alt, covariate = "slider")
}, error = function(e) {
  cat("Still having issues. Here are some diagnostic steps:\n")
  cat("1. Check if your time variable has any negative or zero values\n")
  cat("2. Ensure status is coded as 0/1 or 1/2\n")
  cat("3. Make sure all predictor variables are either numeric (0/1) or factors\n")
  cat("4. Check for any infinite or NaN values\n")
  
  # Show problematic data
  cat("\nTime variable summary:\n")
  print(summary(CRE3_clean$time))
  cat("\nStatus variable summary:\n")
  print(table(CRE3_clean$status))
})

# Final deployment code (only run after successful testing)
# library(rsconnect)
# rsconnect::setAccountInfo(name='idbologna',
#                          token='DA368A4CBFC7435F67911388C6A3E9B4',
#                          secret='N4GLNTKq0RMD4HqdByj+CB6MD05VL1XqFnkRx8cp')
# runApp("DynNomapp")
# deployApp()