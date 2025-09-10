library(shiny)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# Check your current working directory
getwd()

# List files in current directory
list.files()

# Check if CRE3.csv exists
file.exists("CRE3.csv")

# Define UI
ui <- fluidPage(
  # CSS styling
  tags$head(
    tags$style(HTML("
      .well {
        background-color: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .checkbox-container {
        padding: 10px;
        margin: 5px 0;
        background-color: white;
        border-radius: 5px;
        border: 1px solid #e0e0e0;
      }
      h3 {
        color: #2c3e50;
        font-weight: 600;
      }
      .risk-text {
        font-size: 18px;
        font-weight: 500;
        padding: 15px;
        border-radius: 5px;
        margin-top: 20px;
      }
      .low-risk { background-color: #d4edda; color: #155724; }
      .medium-risk { background-color: #fff3cd; color: #856404; }
      .high-risk { background-color: #f8d7da; color: #721c24; }
      .model-info {
        background-color: #e7f3ff;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 20px;
        border-left: 4px solid #007bff;
      }
    "))
  ),
  
  titlePanel("CRE Bloodstream Infection Risk Predictor After Liver Transplantation"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      # Model information at top
      div(class = "model-info",
          h4(icon("info-circle"), "How to Use"),
          p("Select which risk factors are present in your patient, then click 'Calculate Risk' to see the predicted probability of CRE infection."),
          p(strong("Model trained on:"), textOutput("dataInfo", inline = TRUE))
      ),
      
      h3("Clinical Risk Factors"),
      p("Check all factors that apply to the patient:"),
      br(),
      
      div(class = "checkbox-container",
          checkboxInput("reint", 
                        "Reintubation", 
                        value = FALSE),
          em("Need for reintubation after initial extubation")
      ),
      
      div(class = "checkbox-container",
          checkboxInput("mv", 
                        "Mechanical Ventilation", 
                        value = FALSE),
          em("Prolonged mechanical ventilation requirement")
      ),
      
      div(class = "checkbox-container",
          checkboxInput("arf", 
                        "Acute Renal Failure", 
                        value = FALSE),
          em("Development of acute kidney injury")
      ),
      
      div(class = "checkbox-container",
          checkboxInput("crepre_60", 
                        "CRE Pre-transplant (within 60 days)", 
                        value = FALSE),
          em("CRE colonization/infection before transplant")
      ),
      
      div(class = "checkbox-container",
          checkboxInput("crepost_60", 
                        "CRE Post-transplant (within 60 days)", 
                        value = FALSE),
          em("CRE colonization/infection after transplant")
      ),
      
      div(class = "checkbox-container",
          checkboxInput("multipost", 
                        "Multi-drug Resistant Organism Post-transplant", 
                        value = FALSE),
          em("MDRO infection after transplant")
      ),
      
      br(),
      actionButton("predict", "Calculate Risk", 
                   class = "btn-primary btn-lg btn-block",
                   style = "background-color: #007bff; color: white; font-weight: bold;"),
      
      actionButton("reset", "Reset All", 
                   class = "btn-secondary btn-lg btn-block",
                   style = "margin-top: 10px;"),
      
      br(),
      p(strong("Disclaimer:"), "This tool is for educational and research purposes only. Clinical decisions should be based on comprehensive patient assessment and clinical judgment.")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Risk Prediction", 
                 br(),
                 conditionalPanel(
                   condition = "input.predict == 0",
                   div(style = "text-align: center; padding: 50px; background-color: #f8f9fa; border-radius: 10px;",
                       icon("calculator", "fa-4x", style = "color: #007bff;"),
                       h3("Ready to Calculate Risk"),
                       p("Select the relevant risk factors from the panel on the left,"),
                       p("then click 'Calculate Risk' to see the prediction.")
                   )
                 ),
                 conditionalPanel(
                   condition = "input.predict > 0",
                   h4("Predicted Risk of CRE Infection"),
                   uiOutput("riskSummary"),
                   br(),
                   plotOutput("survivalPlot", height = "500px"),
                   br(),
                   h4("Risk Estimates at Key Time Points"),
                   tableOutput("riskTable"),
                   br(),
                   div(style = "background-color: #f0f0f0; padding: 15px; border-radius: 5px;",
                       h5("Clinical Interpretation:"),
                       uiOutput("clinicalInterpretation")
                   )
                 )
        ),
        
        tabPanel("Model Information",
                 br(),
                 h4("Cox Proportional Hazards Model"),
                 p("Formula: Surv(time, status==1) ~ reint + mv + arf + crepre_60 + crepost_60 + multipost"),
                 br(),
                 verbatimTextOutput("modelSummary"),
                 br(),
                 h4("Hazard Ratios and Confidence Intervals"),
                 tableOutput("hazardTable"),
                 br(),
                 h4("Model Performance Metrics"),
                 verbatimTextOutput("modelMetrics")
        ),
        
        tabPanel("Dataset Overview",
                 br(),
                 h4("Training Dataset Summary"),
                 verbatimTextOutput("dataSummary"),
                 br(),
                 h4("Risk Factor Distribution in Training Data"),
                 plotOutput("riskFactorPlot", height = "400px"),
                 br(),
                 h4("Kaplan-Meier Survival Curve (Overall)"),
                 plotOutput("kmPlot", height = "400px")
        ),
        
        tabPanel("Help & About",
                 br(),
                 h4("About This Tool"),
                 p("This risk calculator predicts the probability of carbapenem-resistant Enterobacteriaceae (CRE) bloodstream infections following liver transplantation using a Cox proportional hazards model."),
                 br(),
                 
                 h4("Understanding the Risk Factors"),
                 tags$dl(
                   tags$dt("Reintubation"),
                   tags$dd("Requirement for reintubation after initial extubation attempt, indicating respiratory complications"),
                   br(),
                   tags$dt("Mechanical Ventilation"),
                   tags$dd("Prolonged need for mechanical ventilatory support post-transplant"),
                   br(),
                   tags$dt("Acute Renal Failure"),
                   tags$dd("Development of acute kidney injury requiring medical intervention"),
                   br(),
                   tags$dt("CRE Pre-transplant"),
                   tags$dd("Documented CRE colonization or infection within 60 days before transplantation"),
                   br(),
                   tags$dt("CRE Post-transplant"),
                   tags$dd("Documented CRE colonization or infection within 60 days after transplantation"),
                   br(),
                   tags$dt("MDRO Post-transplant"),
                   tags$dd("Multi-drug resistant organism infection following transplantation")
                 ),
                 br(),
                 
                 h4("Interpreting the Results"),
                 tags$ul(
                   tags$li(strong("Low Risk (< 10% at 180 days):"), " Standard post-transplant monitoring may be appropriate"),
                   tags$li(strong("Moderate Risk (10-30% at 180 days):"), " Consider enhanced surveillance and preventive strategies"),
                   tags$li(strong("High Risk (> 30% at 180 days):"), " Intensive monitoring and aggressive preventive measures recommended")
                 ),
                 br(),
                 
                 h4("Technical Details"),
                 p("The model uses Cox regression to estimate the hazard of CRE infection based on the presence or absence of six key risk factors. The survival curve shows the probability of remaining CRE-free over time, with 95% confidence intervals."),
                 br(),
                 
                 h4("References"),
                 p("This calculator is based on clinical research in liver transplant recipients. For detailed methodology and validation, please consult the original research publication."),
                 br(),
                 
                 h4("Contact"),
                 p("For questions or technical support regarding this calculator, please contact your institution's transplant team or clinical research department.")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Load and prepare data at startup
  cre_data <- reactive({
    tryCatch({
      # Try to read the CSV file
      df <- read.csv("CRE3.csv", stringsAsFactors = FALSE)
      
      # Ensure variables are properly coded
      df$reint <- as.numeric(df$reint)
      df$mv <- as.numeric(df$mv)
      df$arf <- as.numeric(df$arf)
      df$crepre_60 <- as.numeric(df$crepre_60)
      df$crepost_60 <- as.numeric(df$crepost_60)
      df$multipost <- as.numeric(df$multipost)
      df$time <- as.numeric(df$time)
      df$status <- as.numeric(df$status)
      
      return(df)
    }, error = function(e) {
      # If file not found, show error message
      showModal(modalDialog(
        title = "Data File Not Found",
        HTML(paste("Unable to load CRE3.csv. Please ensure the file is in the app directory.<br><br>",
                   "Error details:", e$message)),
        footer = NULL,
        easyClose = FALSE
      ))
      return(NULL)
    })
  })
  
  # Fit Cox model
  cox_model <- reactive({
    df <- cre_data()
    if(is.null(df)) return(NULL)
    
    model <- coxph(Surv(time, status == 1) ~ reint + mv + arf + crepre_60 + crepost_60 + multipost, 
                   data = df)
    return(model)
  })
  
  # Data info for display
  output$dataInfo <- renderText({
    df <- cre_data()
    if(is.null(df)) return("Data not loaded")
    paste(nrow(df), "patients")
  })
  
  # Reset button functionality
  observeEvent(input$reset, {
    updateCheckboxInput(session, "reint", value = FALSE)
    updateCheckboxInput(session, "mv", value = FALSE)
    updateCheckboxInput(session, "arf", value = FALSE)
    updateCheckboxInput(session, "crepre_60", value = FALSE)
    updateCheckboxInput(session, "crepost_60", value = FALSE)
    updateCheckboxInput(session, "multipost", value = FALSE)
  })
  
  # Calculate predictions when button is clicked
  predictions <- eventReactive(input$predict, {
    model <- cox_model()
    if(is.null(model)) return(NULL)
    
    # Create new data frame with user inputs
    new_data <- data.frame(
      reint = as.numeric(input$reint),
      mv = as.numeric(input$mv),
      arf = as.numeric(input$arf),
      crepre_60 = as.numeric(input$crepre_60),
      crepost_60 = as.numeric(input$crepost_60),
      multipost = as.numeric(input$multipost)
    )
    
    # Calculate survival probabilities
    surv_fit <- survfit(model, newdata = new_data)
    
    # Create time points from 0 to 180 days
    time_points <- seq(0, 180, by = 1)
    
    # Get survival probabilities at these time points
    summary_surv <- summary(surv_fit, times = time_points)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      time = summary_surv$time,
      surv_prob = summary_surv$surv,
      risk_prob = 1 - summary_surv$surv,
      lower = summary_surv$lower,
      upper = summary_surv$upper
    )
    
    # Calculate risk at specific time points
    risk_30 <- 1 - summary(surv_fit, times = 30)$surv
    risk_60 <- 1 - summary(surv_fit, times = 60)$surv
    risk_90 <- 1 - summary(surv_fit, times = 90)$surv
    risk_180 <- 1 - summary(surv_fit, times = 180)$surv
    
    return(list(
      plot_data = plot_data,
      risk_30 = risk_30,
      risk_60 = risk_60,
      risk_90 = risk_90,
      risk_180 = risk_180,
      new_data = new_data
    ))
  })
  
  # Create survival plot
  output$survivalPlot <- renderPlot({
    pred <- predictions()
    if(is.null(pred)) return(NULL)
    
    ggplot(pred$plot_data, aes(x = time, y = risk_prob)) +
      geom_line(size = 2, color = "#007bff") +
      geom_ribbon(aes(ymin = 1 - upper, ymax = 1 - lower), 
                  alpha = 0.2, fill = "#007bff") +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 180, 30)) +
      labs(
        title = "Predicted Probability of CRE Bloodstream Infection",
        subtitle = "With 95% Confidence Interval",
        x = "Days After Liver Transplantation",
        y = "Cumulative Risk of CRE Infection"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA)
      ) +
      geom_hline(yintercept = c(0.1, 0.3), linetype = "dashed", alpha = 0.5) +
      annotate("text", x = 170, y = 0.12, label = "10% threshold", size = 3) +
      annotate("text", x = 170, y = 0.32, label = "30% threshold", size = 3)
  })
  
  # Create risk summary
  output$riskSummary <- renderUI({
    pred <- predictions()
    if(is.null(pred)) return(NULL)
    
    risk_180 <- pred$risk_180 * 100
    
    # Determine risk category
    if (risk_180 < 10) {
      risk_class <- "low-risk"
      risk_text <- "LOW RISK"
    } else if (risk_180 <= 30) {
      risk_class <- "medium-risk"
      risk_text <- "MODERATE RISK"
    } else {
      risk_class <- "high-risk"
      risk_text <- "HIGH RISK"
    }
    
    # Count number of risk factors
    n_factors <- sum(as.numeric(c(input$reint, input$mv, input$arf, 
                                  input$crepre_60, input$crepost_60, input$multipost)))
    
    div(
      class = paste("risk-text", risk_class),
      h3(sprintf("%s: %.1f%% probability of CRE infection at 180 days", risk_text, risk_180)),
      p(sprintf("Number of risk factors present: %d out of 6", n_factors))
    )
  })
  
  # Clinical interpretation
  output$clinicalInterpretation <- renderUI({
    pred <- predictions()
    if(is.null(pred)) return(NULL)
    
    risk_180 <- pred$risk_180 * 100
    
    if (risk_180 < 10) {
      tags$ul(
        tags$li("Continue standard post-transplant infection prevention protocols"),
        tags$li("Routine surveillance cultures as per institutional guidelines"),
        tags$li("Monitor for changes in clinical status that may increase risk")
      )
    } else if (risk_180 <= 30) {
      tags$ul(
        tags$li("Consider enhanced surveillance with more frequent cultures"),
        tags$li("Review and optimize infection prevention strategies"),
        tags$li("Consider prophylactic measures based on local resistance patterns"),
        tags$li("Close monitoring for early signs of infection")
      )
    } else {
      tags$ul(
        tags$li("Implement intensive surveillance protocols"),
        tags$li("Consider targeted prophylaxis if appropriate"),
        tags$li("Multidisciplinary team review recommended"),
        tags$li("Aggressive source control and early intervention strategies"),
        tags$li("Consider infectious disease consultation")
      )
    }
  })
  
  # Create risk table
  output$riskTable <- renderTable({
    pred <- predictions()
    if(is.null(pred)) return(NULL)
    
    risk_table <- data.frame(
      "Time Point" = c("30 days", "60 days", "90 days", "180 days"),
      "Risk of CRE Infection" = sprintf("%.1f%%", 
                                        c(pred$risk_30, pred$risk_60, 
                                          pred$risk_90, pred$risk_180) * 100),
      "Infection-Free Probability" = sprintf("%.1f%%", 
                                             c(1-pred$risk_30, 1-pred$risk_60, 
                                               1-pred$risk_90, 1-pred$risk_180) * 100),
      "Number Needed to Screen" = round(100 / c(pred$risk_30, pred$risk_60, 
                                                pred$risk_90, pred$risk_180) / 100)
    )
    
    return(risk_table)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Model summary output
  output$modelSummary <- renderPrint({
    model <- cox_model()
    if(is.null(model)) return(cat("Model not loaded"))
    summary(model)
  })
  
  # Model metrics
  output$modelMetrics <- renderPrint({
    model <- cox_model()
    if(is.null(model)) return(cat("Model not loaded"))
    
    cat("Model Fit Statistics:\n")
    cat("=====================================\n")
    cat("Concordance:", round(summary(model)$concordance[1], 3), "\n")
    cat("Likelihood ratio test p-value:", format.pval(summary(model)$logtest[3], digits = 3), "\n")
    cat("Wald test p-value:", format.pval(summary(model)$waldtest[3], digits = 3), "\n")
    cat("Score test p-value:", format.pval(summary(model)$sctest[3], digits = 3), "\n")
  })
  
  # Hazard ratio table
  output$hazardTable <- renderTable({
    model <- cox_model()
    if(is.null(model)) return(NULL)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)
    hr_table <- data.frame(
      "Risk Factor" = c("Reintubation", "Mechanical Ventilation", 
                        "Acute Renal Failure", "CRE Pre-transplant", 
                        "CRE Post-transplant", "MDRO Post-transplant"),
      "Hazard Ratio" = round(coef_summary$conf.int[, 1], 2),
      "95% CI Lower" = round(coef_summary$conf.int[, 3], 2),
      "95% CI Upper" = round(coef_summary$conf.int[, 4], 2),
      "P-value" = format.pval(coef_summary$coefficients[, 5], digits = 3)
    )
    
    return(hr_table)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Data summary
  output$dataSummary <- renderPrint({
    df <- cre_data()
    if(is.null(df)) return(cat("Data not loaded"))
    
    cat("Dataset Information:\n")
    cat("=====================================\n")
    cat("Total patients:", nrow(df), "\n")
    cat("CRE infections:", sum(df$status == 1), "(", round(mean(df$status == 1) * 100, 1), "%)\n")
    cat("Censored observations:", sum(df$status == 0), "(", round(mean(df$status == 0) * 100, 1), "%)\n\n")
    
    cat("Follow-up Time:\n")
    cat("-------------------------------------\n")
    cat("Median:", round(median(df$time), 1), "days\n")
    cat("Mean:", round(mean(df$time), 1), "days\n")
    cat("Range:", min(df$time), "-", max(df$time), "days\n\n")
    
    cat("Risk Factor Prevalence:\n")
    cat("-------------------------------------\n")
    cat("Reintubation:", sum(df$reint == 1), "(", round(mean(df$reint == 1) * 100, 1), "%)\n")
    cat("Mechanical Ventilation:", sum(df$mv == 1), "(", round(mean(df$mv == 1) * 100, 1), "%)\n")
    cat("Acute Renal Failure:", sum(df$arf == 1), "(", round(mean(df$arf == 1) * 100, 1), "%)\n")
    cat("CRE Pre-transplant:", sum(df$crepre_60 == 1), "(", round(mean(df$crepre_60 == 1) * 100, 1), "%)\n")
    cat("CRE Post-transplant:", sum(df$crepost_60 == 1), "(", round(mean(df$crepost_60 == 1) * 100, 1), "%)\n")
    cat("MDRO Post-transplant:", sum(df$multipost == 1), "(", round(mean(df$multipost == 1) * 100, 1), "%)\n")
  })
  
  # Risk factor distribution plot
  output$riskFactorPlot <- renderPlot({
    df <- cre_data()
    if(is.null(df)) return(NULL)
    
    # Calculate prevalence
    prevalence_data <- data.frame(
      Risk_Factor = c("Reintubation", "Mechanical\nVentilation", "Acute Renal\nFailure", 
                      "CRE\nPre-transplant", "CRE\nPost-transplant", "MDRO\nPost-transplant"),
      Prevalence = c(mean(df$reint == 1), mean(df$mv == 1), mean(df$arf == 1),
                     mean(df$crepre_60 == 1), mean(df$crepost_60 == 1), mean(df$multipost == 1)) * 100
    )
    
    ggplot(prevalence_data, aes(x = reorder(Risk_Factor, Prevalence), y = Prevalence)) +
      geom_bar(stat = "identity", fill = "#007bff", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", Prevalence)), 
                hjust = -0.2, size = 4) +
      coord_flip() +
      labs(title = "Risk Factor Prevalence in Training Dataset",
           x = "",
           y = "Prevalence (%)") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 11)
      ) +
      scale_y_continuous(limits = c(0, max(prevalence_data$Prevalence) * 1.2))
  })
  
  # Kaplan-Meier plot
  output$kmPlot <- renderPlot({
    df <- cre_data()
    if(is.null(df)) return(NULL)
    
    # Create overall survival object
    km_fit <- survfit(Surv(time, status == 1) ~ 1, data = df)
    
    # Plot
    plot(km_fit, 
         xlab = "Days After Liver Transplantation",
         ylab = "Probability of CRE Infection-Free Survival",
         main = "Kaplan-Meier Curve: Overall CRE-Free Survival",
         xlim = c(0, 180),
         col = "#007bff",
         lwd = 2)
    
    # Add confidence intervals
    lines(km_fit, conf.int = TRUE, col = "#007bff", lty = 2)
    
    # Add grid
    grid()
    
    # Add legend
    legend("bottomleft", 
           legend = c("Survival Probability", "95% CI"),
           col = "#007bff",
           lty = c(1, 2),
           lwd = 2)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

rsconnect::setAccountInfo(
  name='idbologna',
  token='DA368A4CBFC7435F67911388C6A3E9B4',  # CHANGE THIS AFTER DEPLOYMENT!
  secret='N4GLNTKq0RMD4HqdByj+CB6MD05VL1XqFnkRx8cp'  # CHANGE THIS AFTER DEPLOYMENT!
)

# Step 5: Deploy the app
rsconnect::deployApp(
  appDir = ".",  # Current directory
  appName = "CRE-Risk-Calculator",  # Name for your app on shinyapps.io
  appTitle = "CRE Risk After Liver Transplant",  # Title shown on shinyapps.io
  forceUpdate = TRUE,
  launch.browser = TRUE
)
