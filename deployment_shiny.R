rsconnect::setAccountInfo(
  name='idbologna',
  token='DA368A4CBFC7435F67911388C6A3E9B4',  # CHANGE THIS AFTER DEPLOYMENT!
  secret='N4GLNTKq0RMD4HqdByj+CB6MD05VL1XqFnkRx8cp'  # CHANGE THIS AFTER DEPLOYMENT!
)


# deploy_script.R - DO NOT INCLUDE IN APP DEPLOYMENT
library(rsconnect)

# This file is for running deployment, not to be deployed itself
rsconnect::deployApp(
  appFiles = c("app.R", "CRE3.csv"),
  appName = "CRE-Risk-Calculator",
  account = "idbologna",
  forceUpdate = TRUE
)
