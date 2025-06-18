source("functions.R")
library(dplyr)
library(tidyr)
library(purrr)

set.seed("2025")
# Load and preprocess data
data <- read.csv("Data/Tessy_data_cleaned.csv", header = TRUE, sep = ",")

# Filter countries with ≥20 outbreaks >4
countries <- unique(data$ReportingCountry)
valid_countries <- list()
for (ctry in countries) {
  outbreak_sizes <- data[data$ReportingCountry == ctry & data$Count > 4, "Count"]
  if(length(outbreak_sizes) >= 20) {
    valid_countries[[ctry]] <- outbreak_sizes
  }
}

### Load MLE and CI results from file Cencored_MLE_R.R
results_df <- read.csv("R_estimates_countries.csv", 
                         stringsAsFactors = FALSE)



################################## TEST COUNTRIES AGAINST EACH OTHER ##########################################
######## PAIRWISE TESTING

# Generate country pairs
target_countries <- unique(results_df$Country)
country_pairs <- combn(target_countries, 2, simplify = FALSE)

# Initialize storage
test_results <- data.frame()

# Main pairwise comparison loop
for (pair in country_pairs) {
  country1 <- pair[1]
  country2 <- pair[2]
  
  data1 <- valid_countries[[country1]]
  data2 <- valid_countries[[country2]]
  
  # Initialize empty results
  llrt_result <- data.frame(LLRT_Statistic = NA, LLRT_CriticalValue = NA, LLRT_Reject = NA)
  wald_result <- data.frame(Wald_Z = NA, Wald_p_value = NA, Wald_Reject = NA)
  bootstrap_result <- data.frame(Bootstrap_p = NA, Bootstrap_Reject = NA)
  
  try({
    # Likelihood Ratio Test
    llrt <- likelihood_ratio_censored(data1, data2, c = 4)
    llrt_result <- data.frame(
      LLRT_Statistic = llrt$Lambda,
      LLRT_CriticalValue = llrt$CriticalValue,
      LLRT_Reject = llrt$Lambda > llrt$CriticalValue
    )
    
    # Wald Test
    R1 <- estimate_R_censored(data1, c = 4)
    R2 <- estimate_R_censored(data2, c = 4)
    
    var1 <- 1/(-hessian_at_estimate(R1, data1, c = 4))
    var2 <- 1/(-hessian_at_estimate(R2, data2, c = 4))
    
    Z <- (R1 - R2)/sqrt(var1 + var2)
    p_value <- 2 * pnorm(-abs(Z))
    
    wald_result <- data.frame(
      Wald_Z = Z,
      Wald_p_value = p_value,
      Wald_Reject = p_value < 0.05
    )
    
    # Bootstrap Test
    bootstrap_result <- bootstrap_distribution_test(data1, data2)
  })
  
  # Combine results
  pair_results <- data.frame(
    Country1 = country1,
    Country2 = country2,
    N1 = length(data1),
    N2 = length(data2),
    R1 = if(exists("R1")) R1 else NA,
    R2 = if(exists("R2")) R2 else NA
  )
  
  pair_results <- cbind(pair_results, llrt_result, wald_result, bootstrap_result)
  test_results <- rbind(test_results, pair_results)
}

# Final formatting
formatted_results <- test_results %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  select(
    Country1, Country2, N1, N2, R1, R2,
    LLRT_Statistic, LLRT_CriticalValue, LLRT_Reject,
    Wald_Z, Wald_p_value, Wald_Reject,
    Bootstrap_p, Bootstrap_Reject
  )

print(formatted_results)
# write.csv(formatted_results, "pairwise_test_results.csv", row.names = FALSE)










################################ TEST FROM EARLIER ################################ 
####################### Censored LLRT between two datasets ####################### 
censored_data <- read.csv("Data/Tessy_data_cleaned.csv", header = TRUE, sep = ",")

################### Subset UK Data and Divide into Periods ###################
uk_data <- censored_data[censored_data$ReportingCountry == "AT", ]
uk_data$year <- as.numeric(substr(uk_data$MinDate, 1, 4))

# Define period boundaries (modify these as needed)
period1_start <- 2009; period1_end <- 2013
period2_start <- 2016; period2_end <- 2019

# Subset UK data into two periods
uk_data_p1 <- uk_data[uk_data$year >= period1_start & uk_data$year <= period1_end, ]$Count
uk_data_p2 <- uk_data[uk_data$year >= period2_start & uk_data$year <= period2_end, ]$Count

R_hat1 <- estimate_R_censored(uk_data_p1$Count, c = 4)
R_hat2 <- estimate_R_censored(uk_data_p2$Count, c = 4)

# Perform LLRT
result <- likelihood_ratio_censored(uk_data_p1, uk_data_p2, c = 4)
cat("Likelihood Ratio Statistic:", result$Lambda, "\n")
cat("Critical Value (alpha=0.05):", result$CriticalValue, "\n")


######################## LLRT: UK vs IE ####################### 
IE_data <- censored_data[censored_data$ReportingCountry == "IE", ]
IE_data$year <- as.numeric(substr(IE_data$MinDate, 1, 4))

UK_data = uk_data$Count
IE_data = IE_data$Count

# Perform LLRT
result <- likelihood_ratio_censored(UK_data, IE_data, c = 4)
cat("Likelihood Ratio Statistic:", result$Lambda, "\n")
cat("Critical Value (alpha=0.05):", result$CriticalValue, "\n")


########################  WALD test: UK vs IE ########################  

# Example for Dataset 1 (censored at c=4)
R1_hat <- estimate_R_censored(IE_data, c = 4)
fisher_info1 <- -hessian_at_estimate(R1_hat, IE_data, c = 4)
var1 <- 1 / fisher_info1

# Repeat for Dataset 2
R2_hat <- estimate_R_censored(UK_data, c = 4)
fisher_info2 <- -hessian_at_estimate(R2_hat, UK_data, c = 4)
var2 <- 1 / fisher_info2

# Compute Z-statistic
Z <- (R1_hat - R2_hat) / sqrt(var1 + var2)
# Two-tailed p-value
p_value <- 2 * pnorm(-abs(Z))
cat("Z =", Z, "\np-value =", p_value, "\n")
if (p_value < 0.05) {
  print("Reject H0: Distributions differ significantly.")
} else {
  print("Fail to reject H0: No significant difference.")
}





