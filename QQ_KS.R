# Import libraries
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
source("functions.R")


#############
# For model validation we will first create a QQ-plot with theoretical and empirical data
# and then a KS-test to compare the two distributions.

# An R value for each country has been computed
R_estimates <- data.frame(
  CountryName = c("Ireland", "Austria", "Poland", "United Kingdom", "Spain"),
  ReportingCountry = c("IE", "AT", "PL", "UK", "ES"),
  R = c(0.66, 0.68, 0.78, 0.88, 0.92)
)


### Collect data + clean
data <- read.csv("Data/Tessy_data_cleaned.csv", stringsAsFactors = FALSE)
data <- data[data$ReportingCountry %in% R_estimates$ReportingCountry, ]


### FUNCTION to compute quantiles for each country - both theoretically and empirically
get_qq_data <- function(data, country_code, R){
  
  # Calculate theoretical q(x) for the given country
  q_data <- data.frame(x = 2:2000) %>%
    mutate(qx = sapply(x, function(xi) compute_q(xi, R)/(1-compute_q(1, R)))) %>%
    mutate(qx = ifelse(is.nan(qx) | is.infinite(qx), 0, qx)) %>%
    filter(qx > 0 & x >= 5)
  
  # Sample from the theoretical distribution
  set.seed(175)
  theory_sample <- sample(q_data$x, size = 1000, replace = TRUE, prob = q_data$qx)
  
  # Empirical data
  emp_sample <- data %>%
    filter(ReportingCountry == country_code, Count >= 5) %>%
    pull(Count)
  
  # Quantiles
  probs <- seq(0, 1, by = 1/length(emp_sample))
  theory_q <- quantile(theory_sample, probs = probs, type = 8)
  emp_q <- quantile(emp_sample, probs = probs, type = 8)
  
  # Combine into a data frame
  data.frame(
    ReportingCountry = country_code,
    Probability = probs,
    Theoretical = as.numeric(theory_q),
    Empirical = as.numeric(emp_q)
  )
  
}


# Apply get_qq_data to all countries in R_estimates
qq_data_all <- purrr::map2_dfr(
  R_estimates$ReportingCountry,
  R_estimates$R,
  ~get_qq_data(data, country_code = .x, R = .y)
)

# Add country names to the corresponding country codes
qq_data_all <- qq_data_all %>%
  left_join(R_estimates, by = "ReportingCountry") %>%
  select(ReportingCountry, CountryName, Probability, Theoretical, Empirical)



### QQ-PLOT ###
p <- ggplot(qq_data_all, aes(x = Theoretical, y = Empirical)) +
  geom_point(alpha = 0.7, size = 1.5, color = "firebrick3") +
  geom_abline(slope = 1, intercept = 0, color = "firebrick4") +
  facet_wrap(~CountryName, scales = "free") +
  labs(
    x = "Theoretical Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
  )

# Print and save plot
print(p)
save_plot(p, "QQ_plot")







### Chi-squared test for goodness of fit ###
chisq_test <- function(data, country_code, R){
  # Collect data for a specific country
  emp_sample <- data %>%
    filter(ReportingCountry == country_code, Count >= 5) %>%
    pull(Count)
  
  # Observed frequencies
  obs_counts <- table(emp_sample)
  x_vals <- as.integer(names(obs_counts))
  
  # Expected probabilities using q(x)
  pmf <- sapply(x_vals, function(x) compute_q(x, R) / (1 - compute_q(1, R)))
  
  # Chi-squared test
  res <- chisq.test(x = as.numeric(obs_counts), p = pmf, rescale.p = TRUE, simulate.p.value = TRUE, B = 10000)
  
  return(tibble(
    p_value = res$p.value,
    statistic = res$statistic
  ))
}

# Run for all countries
results <- R_estimates %>%
  mutate(test_result = map2(ReportingCountry, R, ~chisq_test(data, .x, .y))) %>%
  unnest(test_result)



print(results)






