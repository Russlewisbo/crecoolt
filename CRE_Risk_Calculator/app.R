library(shiny)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

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
      .upload-box {
        padding: 20px;
        border: 2px dashed #007bff;
        border-radius: 10px;
        background-color: #f0f8ff;
        text-align: center;
        margin-bottom: 20px;
      }
    "))
  ),
  
  titlePanel("CRE Bloodstream Infection Risk Predictor After Liver Transplantation"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      # File upload section
      div(class = "upload-box",
        h4("Step 1: Upload Data"),
        fileInput("datafile", 
                 "Choose CRE3.csv file",
                 accept = c(".csv"),
                 buttonLabel = "Browse...",
                 placeholder = "No file selected"),
        uiOutput("uploadStatus")
      ),
      
      # Risk factors section (initially hidden)
      conditionalPanel(
        condition = "output.fileUploaded",
        h3("Step 2: Select Clinical Risk Factors"),
        p("Select the presence of each risk factor:"),
        br(),
        
        div(class = "checkbox-container",
          checkboxInput("reint", 
                       "Reintubation", 
                       value = FALSE)
        ),
        
        div(class = "checkbox-container",
          checkboxInput("mv", 
                       "Mechanical Ventilation", 
                       value = FALSE)
        ),
        
        div(class = "checkbox-container",
          checkboxInput("arf", 
                       "Acute Renal Failure", 
                       value = FALSE)
        ),
        
        div(class = "checkbox-container",
          checkboxInput("crepre_60", 
                       "CRE Pre-transplant (within 60 days)", 
                       value = FALSE)
        ),
        
        div(class = "checkbox-container",
          checkboxInput("crepost_60", 
                       "CRE Post-transplant (within 60 days)", 
                       value = FALSE)
        ),
        
        div(class = "checkbox-container",
          checkboxInput("multipost", 
                       "Multi-drug Resistant Organism Post-transplant", 
                       value = FALSE)
        ),
        
        br(),
        actionButton("predict", "Calculate Risk", 
                    class = "btn-primary btn-lg btn-block",
                    style = "background-color: #007bff; color: white; font-weight: bold;")
      ),
      
      br(),
      h4("Model Information"),
      p("This model uses Cox proportional hazards regression to predict the probability of CRE bloodstream infection within 180 days after liver transplantation."),
      p(strong("Note:"), "This tool is for educational purposes only and should not replace clinical judgment.")
    ),
    
    mainPanel(
      width = 8,
      conditionalPanel(
        condition = "!output.fileUploaded",
        br(),
        div(style = "text-align: center; padding: 50px;",
          icon("upload", "fa-5x", style = "color: #007bff;"),
          h3("Please upload your CRE3.csv file to begin"),
          p("The file should contain the following columns:"),
          p(code("time, status, reint, mv, arf, crepre_60, crepost_60, multipost")),
          br(),
          p("Once uploaded, you'll be able to select risk factors and generate predictions.")
        )
      ),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        tabsetPanel(
          tabPanel("Prediction", 
                  br(),
                  h4("Predicted Risk of CRE Infection"),
                  uiOutput("riskSummary"),
                  br(),
                  plotOutput("survivalPlot", height = "500px"),
                  br(),
                  h4("Risk Estimates at Key Time Points"),
                  tableOutput("riskTable")
          ),
          
          tabPanel("Model Details",
                  br(),
                  h4("Cox Regression Model"),
                  verbatimTextOutput("modelSummary"),
                  br(),
                  h4("Hazard Ratios"),
                  tableOutput("hazardTable"),
                  br(),
                  p("The model formula: Surv(time, status==1) ~ reint + mv + arf + crepre_60 + crepost_60 + multipost")
          ),
          
          tabPanel("Data Summary",
                  br(),
                  h4("Uploaded Data Summary"),
                  verbatimTextOutput("dataSummary"),
                  br(),
                  h4("First Few Rows of Data"),
                  tableOutput("dataPreview")
          ),
          
          tabPanel("About",
                  br(),
                  h4("About This Tool"),
                  p("This application predicts the risk of carbapenem-resistant Enterobacteriaceae (CRE) bloodstream infections following liver transplantation."),
                  br(),
                  h4("Risk Factors Explained:"),
                  tags$ul(
                    tags$li(strong("Reintubation:"), " Need for reintubation after initial extubation"),
                    tags$li(strong("Mechanical Ventilation:"), " Prolonged mechanical ventilation requirement"),
                    tags$li(strong("Acute Renal Failure:"), " Development of acute kidney injury requiring intervention"),
                    tags$li(strong("CRE Pre-transplant:"), " CRE colonization or infection within 60 days before transplant"),
                    tags$li(strong("CRE Post-transplant:"), " CRE colonization or infection within 60 days after transplant"),
                    tags$li(strong("MDRO Post-transplant:"), " Multi-drug resistant organism infection after transplant")
                  ),
                  br(),
                  h4("Interpretation:"),
                  p("The survival curve shows the probability of remaining CRE infection-free over time. Lower curves indicate higher risk of infection."),
                  p("Risk categories:"),
                  tags$ul(
                    tags$li("Low risk: < 10% probability at 180 days"),
                    tags$li("Moderate risk: 10-30% probability at 180 days"),
                    tags$li("High risk: > 30% probability at 180 days")
                  )
          )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive value to store upload status
  values <- reactiveValues(
    data = NULL,
    model = NULL,
    fileUploaded = FALSE
  )
  
  # Handle file upload
  observeEvent(input$datafile, {
    req(input$datafile)
    
    tryCatch({
      # Read the CSV file
      df <- read.csv(input$datafile$datapath, stringsAsFactors = FALSE)
      
      # Check for required columns
      required_cols <- c("time", "status", "reint", "mv", "arf", 
                        "crepre_60", "crepost_60", "multipost")
      missing_cols <- required_cols[!required_cols %in% names(df)]
      
      if(length(missing_cols) > 0) {
        showModal(modalDialog(
          title = "Missing Columns",
          paste("The following required columns are missing:", 
                paste(missing_cols, collapse = ", ")),
          easyClose = TRUE
        ))
        return()
      }
      
      # Ensure variables are properly coded
      df$reint <- as.numeric(df$reint)
      df$mv <- as.numeric(df$mv)
      df$arf <- as.numeric(df$arf)
      df$crepre_60 <- as.numeric(df$crepre_60)
      df$crepost_60 <- as.numeric(df$crepost_60)
      df$multipost <- as.numeric(df$multipost)
      df$time <- as.numeric(df$time)
      df$status <- as.numeric(df$status)
      
      # Store data
      values$data <- df
      values$fileUploaded <- TRUE
      
      # Fit Cox model
      values$model <- coxph(Surv(time, status == 1) ~ reint + mv + arf + 
                           crepre_60 + crepost_60 + multipost, 
                           data = df)
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Error reading file:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Output for file upload status
  output$fileUploaded <- reactive({
    values$fileUploaded
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  # Upload status message
  output$uploadStatus <- renderUI({
    if(values$fileUploaded) {
      tags$div(style = "color: green; font-weight: bold;",
               icon("check-circle"),
               " Data loaded successfully!")
    }
  })
  
  # Data summary
  output$dataSummary <- renderPrint({
    req(values$data)
    cat("Dataset dimensions:", nrow(values$data), "rows x", ncol(values$data), "columns\n\n")
    cat("Number of events (CRE infections):", sum(values$data$status == 1), "\n")
    cat("Number of censored observations:", sum(values$data$status == 0), "\n")
    cat("Median follow-up time:", median(values$data$time), "days\n")
    cat("Max follow-up time:", max(values$data$time), "days\n\n")
    
    cat("Risk factor prevalence:\n")
    cat("- Reintubation:", mean(values$data$reint == 1) * 100, "%\n")
    cat("- Mechanical Ventilation:", mean(values$data$mv == 1) * 100, "%\n")
    cat("- Acute Renal Failure:", mean(values$data$arf == 1) * 100, "%\n")
    cat("- CRE Pre-transplant:", mean(values$data$crepre_60 == 1) * 100, "%\n")
    cat("- CRE Post-transplant:", mean(values$data$crepost_60 == 1) * 100, "%\n")
    cat("- MDRO Post-transplant:", mean(values$data$multipost == 1) * 100, "%\n")
  })
  
  # Data preview
  output$dataPreview <- renderTable({
    req(values$data)
    head(values$data[, c("time", "status", "reint", "mv", "arf", 
                         "crepre_60", "crepost_60", "multipost")], 10)
  })
  
  # Calculate predictions when button is clicked
  predictions <- eventReactive(input$predict, {
    req(values$model)
    
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
    surv_fit <- survfit(values$model, newdata = new_data)
    
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
    req(input$predict)
    
    pred <- predictions()
    
    ggplot(pred$plot_data, aes(x = time, y = risk_prob)) +
      geom_line(size = 2, color = "#007bff") +
      geom_ribbon(aes(ymin = 1 - upper, ymax = 1 - lower), 
                  alpha = 0.2, fill = "#007bff") +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(breaks = seq(0, 180, 30)) +
      labs(
        title = "Predicted Probability of CRE Bloodstream Infection",
        x = "Days After Liver Transplantation",
        y = "Cumulative Risk of CRE Infection"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA)
      ) +
      geom_hline(yintercept = c(0.1, 0.3), linetype = "dashed", alpha = 0.5) +
      annotate("text", x = 170, y = 0.12, label = "10% risk", size = 3) +
      annotate("text", x = 170, y = 0.32, label = "30% risk", size = 3)
  })
  
  # Create risk summary
  output$riskSummary <- renderUI({
    req(input$predict)
    
    pred <- predictions()
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
  
  # Create risk table
  output$riskTable <- renderTable({
    req(input$predict)
    
    pred <- predictions()
    
    risk_table <- data.frame(
      "Time Point" = c("30 days", "60 days", "90 days", "180 days"),
      "Risk of CRE Infection" = sprintf("%.1f%%", 
                                        c(pred$risk_30, pred$risk_60, 
                                          pred$risk_90, pred$risk_180) * 100),
      "Infection-Free Probability" = sprintf("%.1f%%", 
                                             c(1-pred$risk_30, 1-pred$risk_60, 
                                               1-pred$risk_90, 1-pred$risk_180) * 100)
    )
    
    return(risk_table)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Model summary output
  output$modelSummary <- renderPrint({
    req(values$model)
    summary(values$model)
  })
  
  # Hazard ratio table
  output$hazardTable <- renderTable({
    req(values$model)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(values$model)
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
}

# Run the application
shinyApp(ui = ui, server = server)