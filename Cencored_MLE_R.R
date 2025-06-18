# Import libraries
source("functions.R")
library(dplyr)
library(tidyr)

set.seed("2025") # the bootstrap CIs change alot per run..

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

# Initialize results storage
results <- list()

# Main analysis loop
for (country in names(valid_countries)) {
  # Get outbreak sizes for current country
  outbreak_sizes <- valid_countries[[country]]
  
  tryCatch({
    # MLE Estimation
    R_hat <- estimate_R_censored(outbreak_sizes, c = 4)
    
    # Fisher Information CI
    hessian <- hessian_at_estimate(R_hat, outbreak_sizes, c = 4)
    observed_fisher <- -hessian
    if(observed_fisher > 0) {
      se <- sqrt(1 / observed_fisher)
      fisher_ci <- R_hat + c(-1.96, 1.96) * se
    } else {
      fisher_ci <- c(NA, NA)
      warning("Non-positive Fisher information for ", country)
    }
    
    # Bootstrap CI
    bootstrap_estimates <- bootstrap_R_censored(outbreak_sizes, c = 4, n_bootstraps = 1000)
    if(length(bootstrap_estimates) >= 100) {  # Ensure sufficient bootstrap samples
      bootstrap_ci <- quantile(bootstrap_estimates, probs = c(0.025, 0.975), na.rm = TRUE)
    } else {
      bootstrap_ci <- c(NA, NA)
      warning("Insufficient bootstrap samples for ", country)
    }
    
    # Store results
    results[[country]] <- data.frame(
      Country = country,
      R_hat = R_hat,
      Fisher_CI_low = fisher_ci[1],
      Fisher_CI_high = fisher_ci[2],
      Bootstrap_CI_low = bootstrap_ci[1],
      Bootstrap_CI_high = bootstrap_ci[2],
      N_Outbreaks = length(outbreak_sizes),
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    message("\nFailed for ", country, ": ", e$message)
  })
}

# Combine results into dataframe
results_df <- do.call(rbind, results)

# Save results
write.csv(results_df, "R_estimates_countries.csv", row.names = FALSE)

# Print formatted results
print(results_df)




#################################### SIMPLE CALCULATIONS PREVIOUS ############################################ 

################### Load and Clean Data ###################
censored_data <- read.csv("Data/Tessy_data_cleaned.csv", header = TRUE, sep = ",")


#### Choose country
Country = "ES"
censored_data_Country <- censored_data[censored_data$ReportingCountry == Country, ]
#### Only counts are saved
outbreak_sizes <- censored_data_Country$Count


######################## Validatation of estimate_R_censored function ######################### 
############### with c = 1 closed-form solution ###############
censored_data_c1 <- outbreak_sizes[outbreak_sizes > 1]
m <- mean(censored_data_c1)
R_hat_closed <- 1 - 2 / m
R_hat_numerical <- estimate_R_censored(outbreak_sizes, c = 1)

cat("Closed-form R_hat (c=1):", R_hat_closed, "\n")
cat("Numerical R_hat (c=1):", R_hat_numerical, "\n")

#################  Estimate R for c = 4 ##############  
R_hat <- estimate_R_censored(outbreak_sizes, c = 4)
cat("MLE for R (censored at c):", R_hat, "\n")

##### Validate with bootstrap confidence interval:
bootstrap_estimates <- bootstrap_R_censored(outbreak_sizes, n_bootstraps = 1000,c = 4)
R_ci <- quantile(bootstrap_estimates, probs = c(0.025, 0.975))
cat("95% CI for R:", R_ci, "\n")

####################### Fisher Information Method of uncertainty estimation #################### 
# Example usage
R_hat <- estimate_R_censored(outbreak_sizes, c = 4)
observed_fisher_info <- -hessian_at_estimate(R_hat, outbreak_sizes, c = 4)
var_R_hat <- 1 / observed_fisher_info

# 95% Confidence interval
se_R_hat <- sqrt(var_R_hat)
ci <- R_hat + c(-1.96, 1.96) * se_R_hat
cat("95% CI for R:", ci, "\n")






