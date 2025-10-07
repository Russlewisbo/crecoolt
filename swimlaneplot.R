# Load required libraries
library(tidyverse)
library(ggplot2)
library(lubridate)
library(patchwork)
library(viridis)
library(scales)

# Read the data
df <- read.csv("CRECOOLT_overall.csv", stringsAsFactors = FALSE)

# Data preprocessing
df_processed <- df %>%
  mutate(
    # Convert dates to Date format
    date_of_transplant = as.Date(date_of_transplant, format = "%Y-%m-%d"),
    cre_infection_date = as.Date(cre_infection_date, format = "%Y-%m-%d"),
    death_date = as.Date(death_date, format = "%Y-%m-%d"),
    discharge_date = as.Date(discharge_date, format = "%Y-%m-%d"),
    
    # Calculate days from transplant to events
    days_to_infection = as.numeric(cre_infection_date - date_of_transplant),
    days_to_death = as.numeric(death_date - date_of_transplant),
    days_to_discharge = as.numeric(discharge_date - date_of_transplant),
    
    # Create treatment strategy from REDCap checkbox variables
    # Assuming cre_colonization___1, ___2, ___3 represent different strategies
    # Note: Adjust these labels based on your actual treatment protocols
    # If a patient has multiple strategies checked, the first matching condition will be used
    treatment_strategy = case_when(
      cre_colonization___1 == 1 ~ "No Prophylaxis",
      cre_colonization___2 == 1 ~ "Prophylaxis",
      cre_colonization___3 == 1 ~ "Pre-emptive Therapy",
      # Check if patient is colonized but no specific strategy
      (is.na(cre_colonization___1) | cre_colonization___1 == 0) & 
        (is.na(cre_colonization___2) | cre_colonization___2 == 0) & 
        (is.na(cre_colonization___3) | cre_colonization___3 == 0) ~ "No CRE Colonization",
      TRUE ~ "No CRE Colonization"
    ),
    treatment_strategy = factor(treatment_strategy,
                                levels = c("No CRE Colonization", 
                                           "No Prophylaxis",
                                           "Prophylaxis", 
                                           "Pre-emptive Therapy")),
    
    # Determine primary outcome for each patient
    primary_event = case_when(
      !is.na(cre_infection_date) & cre_infection == 1 ~ "CRE Infection",
      !is.na(death_date) & death == 1 ~ "Death",
      !is.na(discharge_date) ~ "Discharge",
      TRUE ~ "Censored"
    ),
    
    # Calculate time to primary event
    time_to_event = case_when(
      primary_event == "CRE Infection" ~ days_to_infection,
      primary_event == "Death" ~ days_to_death,
      primary_event == "Discharge" ~ days_to_discharge,
      TRUE ~ 180  # Censoring at 180 days
    ),
    
    # Cap at 180 days for visualization
    time_to_event_capped = pmin(time_to_event, 180, na.rm = TRUE),
    
    # Calculate 180-day mortality
    mortality_180d = ifelse(!is.na(days_to_death) & days_to_death <= 180, 1, 0)
  ) %>%
  filter(!is.na(treatment_strategy)) %>%
  arrange(treatment_strategy, time_to_event)

# Assign patient IDs for y-axis
df_processed <- df_processed %>%
  group_by(treatment_strategy) %>%
  mutate(patient_order = row_number()) %>%
  ungroup() %>%
  mutate(patient_id = paste0(treatment_strategy, "_", patient_order))

# Create the main swimlane plot
swimlane_plot <- ggplot(df_processed) +
  # Draw horizontal lines for each patient (swimlanes)
  geom_segment(aes(x = 0, xend = time_to_event_capped, 
                   y = patient_order, yend = patient_order,
                   color = treatment_strategy),
               size = 0.3, alpha = 0.5) +
  
  # Add event markers
  geom_point(aes(x = time_to_event_capped, y = patient_order, 
                 shape = primary_event, fill = primary_event),
             size = 1.5, alpha = 0.8) +
  
  # Facet by treatment strategy - 2 columns for better space usage
  facet_wrap(~ treatment_strategy, scales = "free_y", ncol = 2) +
  
  # Customize scales
  scale_x_continuous(breaks = seq(0, 180, 30),
                     limits = c(0, 185),
                     expand = c(0.02, 0)) +
  
  # Define colors and shapes
  scale_color_manual(values = c("No CRE Colonization" = "#2E7D32",
                                "No Prophylaxis" = "#FF9800",
                                "Prophylaxis" = "#1976D2",
                                "Pre-emptive Therapy" = "#9C27B0")) +
  
  scale_shape_manual(values = c("CRE Infection" = 21,
                                "Death" = 24,
                                "Discharge" = 22,
                                "Censored" = 3)) +
  
  scale_fill_manual(values = c("CRE Infection" = "#FF6B6B",
                               "Death" = "#4E4E4E",
                               "Discharge" = "#4ECDC4",
                               "Censored" = "#95E1D3")) +
  
  # Labels and theme
  labs(title = "Swimlane Plot: CRE Infection with Competing Risks",
       subtitle = "Stratified by CRE Colonization Treatment Strategy",
       x = "Days from Transplant",
       y = "Patients",
       shape = "Event Type",
       fill = "Event Type",
       color = "Treatment Strategy") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",  # Arrange legends horizontally
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey90", color = NA)
  ) +
  
  guides(color = "none")  # Hide treatment strategy legend as it's shown in facets

# Create summary statistics table
summary_stats <- df_processed %>%
  group_by(treatment_strategy) %>%
  summarise(
    n_patients = n(),
    n_cre_infections = sum(primary_event == "CRE Infection", na.rm = TRUE),
    n_deaths = sum(primary_event == "Death", na.rm = TRUE),
    n_discharges = sum(primary_event == "Discharge", na.rm = TRUE),
    mortality_180d_pct = mean(mortality_180d, na.rm = TRUE) * 100,
    median_time_to_event = median(time_to_event_capped, na.rm = TRUE)
  ) %>%
  mutate(
    cre_infection_rate = n_cre_infections / n_patients * 100,
    death_rate = n_deaths / n_patients * 100
  )

# Create heatmap for 180-day mortality
heatmap_data <- df_processed %>%
  group_by(treatment_strategy) %>%
  summarise(
    mortality_rate = mean(mortality_180d, na.rm = TRUE) * 100,
    n = n()
  )

heatmap_plot <- ggplot(heatmap_data, aes(x = treatment_strategy, y = 1)) +
  geom_tile(aes(fill = mortality_rate), color = "grey40", size = 1) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", mortality_rate, n)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#D4F1D4", mid = "#FFF8DC", high = "#FFCCCB",
                       midpoint = 15,
                       limits = c(0, max(heatmap_data$mortality_rate) * 1.1),
                       name = "180-Day\nMortality (%)") +
  labs(title = "180-Day Mortality by Treatment Strategy",
       x = "",
       y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    legend.position = "none",  # Remove legend to avoid overlap
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed(ratio = 3)

# Create cumulative incidence plot
cumulative_data <- df_processed %>%
  filter(!is.na(time_to_event_capped)) %>%
  arrange(treatment_strategy, time_to_event_capped) %>%
  group_by(treatment_strategy) %>%
  mutate(
    cumulative_infection = cumsum(primary_event == "CRE Infection") / n() * 100,
    cumulative_death = cumsum(primary_event == "Death") / n() * 100
  )

cumulative_plot <- ggplot(cumulative_data) +
  geom_step(aes(x = time_to_event_capped, y = cumulative_infection, 
                color = treatment_strategy, linetype = "CRE Infection"),
            size = 0.8) +
  geom_step(aes(x = time_to_event_capped, y = cumulative_death, 
                color = treatment_strategy, linetype = "Death"),
            size = 0.8) +
  scale_color_manual(values = c("No CRE Colonization" = "#2E7D32",
                                "No Prophylaxis" = "#FF9800",
                                "Prophylaxis" = "#1976D2",
                                "Pre-emptive Therapy" = "#9C27B0"),
                     name = "Treatment Strategy") +
  scale_linetype_manual(values = c("CRE Infection" = "solid", "Death" = "dashed"),
                        name = "Event Type") +
  labs(title = "Cumulative Incidence of Events",
       x = "Days from Transplant",
       y = "Cumulative Incidence (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "right",  # Move legend to right side
    legend.direction = "vertical",  # Stack legends vertically
    legend.box = "vertical",  # Arrange multiple legends vertically
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  xlim(0, 180)

# Combine all plots using patchwork
# Option 1: Vertical layout with better spacing
final_plot <- swimlane_plot / (heatmap_plot + cumulative_plot) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "CRE Infection Outcomes Analysis in Liver Transplant Recipients",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# You can choose which option to display by commenting/uncommenting

# Option 1: Display combined plot with adjusted layout
print(final_plot)

# Option 2: Display just the simple version (swimlane + heatmap)
# print(simple_plot)

# Option 3: Display plots separately to completely avoid overlap
# print(swimlane_plot)
# print(heatmap_plot) 
# print(cumulative_plot)

# For presentations/publications, you might want to use the separate files:
# - swimlane_only.png (main swimlane visualization)
# - heatmap_only.png (mortality heatmap with black text)
# - cumulative_only.png (cumulative incidence curves)

# Print summary statistics
cat("\n===== SUMMARY STATISTICS =====\n\n")
print(summary_stats)

# Save the combined plot with larger dimensions
ggsave("cre_swimlane_analysis.png", final_plot, width = 16, height = 12, dpi = 300)

# Alternative: Create separate plots for better readability
# Save individual components with optimized sizes
ggsave("swimlane_main.png", swimlane_plot, width = 14, height = 10, dpi = 300)
ggsave("heatmap_mortality.png", heatmap_plot, width = 10, height = 4, dpi = 300)
ggsave("cumulative_incidence.png", cumulative_plot, width = 10, height = 6, dpi = 300)

# Option 2: Create a simpler combined plot with just swimlane and heatmap
simple_plot <- swimlane_plot / heatmap_plot +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "CRE Infection Outcomes with 180-Day Mortality",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave("cre_analysis_simple.png", simple_plot, width = 14, height = 10, dpi = 300)

# Option 3: Side-by-side layout for bottom plots
alternative_plot <- swimlane_plot / 
  (heatmap_plot | plot_spacer() | cumulative_plot) +
  plot_layout(heights = c(3, 1), widths = c(1, 0.1, 1)) +
  plot_annotation(
    title = "CRE Infection Outcomes Analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave("cre_analysis_alternative.png", alternative_plot, width = 16, height = 10, dpi = 300)