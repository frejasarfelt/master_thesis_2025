# Import libraries
library(ggplot2)
source("functions.R")

# Load and preprocess data
data <- read.csv("Data/Tessy_data_cleaned.csv", header = TRUE, sep = ",")

# Filter countries with ≥20 outbreaks >4
# Save the start dates
countries <- unique(data$ReportingCountry)
outbreak_dates <- list()
for (ctry in countries) {
  outbreak_sizes <- data[data$ReportingCountry == ctry, c("Count", "MinDate")]
  
  if(sum(outbreak_sizes$Count>=5) >= 20) {
    outbreak_dates[[ctry]] <- outbreak_sizes
  }
}


# Combine all valid countries into one data frame
plot_data <- do.call(rbind, lapply(names(outbreak_dates), function(ctry) {
  df <- outbreak_dates[[ctry]]
  df$Country <- ctry
  return(df)
}))

# Convert MinDate to actual date (e.g., first of month)
plot_data$MinDate <- as.Date(paste0(plot_data$MinDate, "-01"))


#### Removed largest point for ES:
plot_data <- plot_data[!(plot_data$Country == "ES" & plot_data$Count == max(plot_data$Count[plot_data$Country == "ES"])), ]

country_names <- c(
  AT = "Austria",
  ES = "Spain",
  IE = "Ireland",
  PL = "Poland",
  UK = "United Kingdom"
)

p_outbreaks <- ggplot(plot_data, aes(x = MinDate, y = Count)) +
  geom_point(size = 2, alpha = 0.8, colour = "firebrick3") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor.x    = element_blank(),
    panel.grid.minor.y    = element_blank(),
    legend.position       = "none" 
  ) +
  labs(
    x = "Outbreak Start Date",
    y = "Outbreak Size"
  ) +
  facet_wrap(
    ~ Country,
    ncol     = 1,
    scales   = "free_y",
    labeller = labeller(Country = country_names)
  )



ggsave(
  filename = "figures/outbreaks_over_time_by_country.png",
  plot     = p_outbreaks,
  width    = 8,     # adjust for desired width (in inches)
  height   = 1.8 * length(unique(plot_data$Country)),  
  # e.g. 2" per panel vertically
  units    = "in",
  dpi      = 300
)
