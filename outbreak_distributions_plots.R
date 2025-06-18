# Import libraries
library(zoo)
library(ggplot2)
library(scales)
library(here)
source("functions.R")


################### Load data as data frame ################### 
# load cleaned data from csv file - comma separated
data <- read.csv("Data/Tessy_data_cleaned.csv", header = TRUE, sep = ",")
cumulative_data_all <- prepare_cumulative_data(data)



################# PLOT MULTIPLE COUNTRIES WITH Q(X): cumulative distribution ################# 
################################### Excluding x = 1 on axis ################################## 
my_colors <- c("#00B4D8", "#FF007F", "green3", "#FF6B35", "#6A4C93")
my_colors <- c("firebrick4","hotpink", "orange2","firebrick1", "gold1")
my_colors <- c(
  "#4477AA",
  "#EE6677",
  "#228833",
  "#CCBB44",
  "#AA3377" 
  )


# Create theoretical adjusted cumulative probability data for outbreak sizes 2:1000
q_data_censored <- do.call(rbind, lapply(seq(0.1, 0.9, by = 0.1), function(R) {
  data.frame(
    x = 2:1000, 
    y = sapply(2:1000, function(n) compute_cumulative_logtransform_q(n, R)) /
      (1 - compute_q(1, R)),
    R = R
  )
}))


# Plot cumulative distribution for all countries with x > 1
p_censored_all <- ggplot(cumulative_data_all, 
                         aes(x = sorted_outbreak_sizes, y = cumulative_prob, color = ReportingCountry)) +
  geom_step(size = 1, direction = "vh") +
  scale_y_continuous(
    trans = "log10", 
    breaks = 10^(-seq(0, 3, by = 1)), 
    labels = scales::math_format(.x),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    trans = "log10", 
    breaks = c(1, 2, 5, 10, 100, 1000), 
    labels = c("1", "2", "5", "10", "100", "1000"),
    minor_breaks = NULL
  ) +
  coord_cartesian(xlim = c(2, 1000), ylim = c(0.001, 1), clip = "off") +
  labs(
    x = "Outbreak Size", y = "Cumulative Probability",
    color = "Country"
  ) +
  scale_color_manual(values = my_colors,labels = c("AT" = "Austria", "ES" = "Spain",  # Map country codes to full names
                                                       "IE" = "Ireland", "PL" = "Poland",
                                                       "UK" = "United Kingdom") )+
  geom_vline(xintercept = 2, linetype = "solid", color = "grey90") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
  theme_minimal()+
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.position = "right",
    panel.grid.minor.x = element_blank(),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14, hjust = 0.5),
    legend.title.align = 0
    
  ) +
  # Theoretical curves layer with a constant linetype mapping 
  geom_line(
    aes(x = x, y = y, group = factor(R), linetype = "h"), 
    data = subset(q_data_censored, y >= 0.001 & x >= 2),  
    color = "gray50", size = 0.2
  ) +
  scale_linetype_manual(name = NULL, 
                        values = c("h" = "dashed"),
                        labels = c("h" = expression(q(x)~ " for "  ~R[e] %in%~"{0.1, 0.2, ..., 0.9}"))) + 
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

# Print and save plot
print(p_censored_all)
save_plot(p_censored_all, "cum_prop_all_countries_q(x)", width = 9, height = 5)




####################### ONLY ONE country, x > 1 ####################### 
#### CHOSE COUNTRY 
Country = "ES"
country_prop_cum_censored <- cumulative_data_all[cumulative_data_all$ReportingCountry == Country, ]

ggplot(country_prop_cum_censored, aes(x = sorted_outbreak_sizes, y = cumulative_prob)) +
  geom_step(size = 1, direction = "vh", color = "black") +
  scale_y_continuous(
    trans = "log10",
    breaks = 10^(-seq(0, 3)),
    labels = scales::math_format(.x)
  ) +
  coord_cartesian(ylim = c(0.001, 1), xlim = c(2, 1000)) + # Use coord_cartesian instead of limits +
  scale_x_continuous(
    trans = "log10",
    limits = c(2, 1000),    # Force x axis from 2 to 1000
    breaks = c( 5, 10, 100, 1000),
    labels = c("5", "10", "100", "1000")
  ) +
  labs(
    x = "Outbreak Size", y = "Cumulative Probability"
  ) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.minor.x = element_blank()
  ) +
  # Plot the adjusted cumulative theoretical curves for R values from 0.1 to 0.9
  geom_line(
    aes(x = x, y = y, color = factor(R)), 
    data = q_data_censored %>% 
      filter(y >= 0.001, x >= 2) %>%  # Match plot limits
      mutate(R = factor(R)),           # Convert R to factor explicitly
    size = 0.5, 
    linetype = "dashed"
  )


######## AGE:

ggplot(data, aes(x = AgeGroup, y = Year, fill = OutbreakSize)) +  
  geom_tile() +  
  facet_wrap(~ReportingCountry)  






############################# ABOUT DATA YEARS AND COUNTS ############################# 

# Extract year from MinDate
data$Year <- substr(data$MinDate, 1, 4)

# Filter only valid countries
valid_country_names <- names(valid_countries)
data_valid <- subset(data, ReportingCountry %in% valid_country_names)

# Create summary list per country
summary_list <- lapply(valid_country_names, function(country) {
  subset_data <- subset(data_valid, ReportingCountry == country)
  years <- sort(unique(subset_data$Year))
  total_outbreaks <- nrow(subset_data)
  large_outbreaks <- sum(subset_data$Count >= 5)
  
  list(
    Years = years,
    TotalOutbreaks = total_outbreaks,
    LargeOutbreaks = large_outbreaks
  )
})

# Assign names to the list
names(summary_list) <- valid_country_names

# View the result
print(summary_list)




