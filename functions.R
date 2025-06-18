# This R file should contain all relevant functions for the project.
# The functions should be called from other R files.

##### Library
library(hypergeo)


####################### Simulate outbreaks ####################### 
simulate_outbreak <- function(R) {
  current_infected <- 1  
  event <- 0
  p_up <- R / (1 + R)  # Probability of increasing I
  
  # Simulation loop until the outbreak dies out
  while (current_infected > 0) {
    event <- event + 1
    # Determine if the number of infected individuals increases or decreases
    if (runif(1) < p_up) {
      current_infected <- current_infected + 1  # Increase number of infected individuals
    } else {
      current_infected <- current_infected - 1  # Decrease number of infected individuals
    }
  }
  outbreak_size = 1 + (event-1)/2
  return(outbreak_size)
}


####################### Simulate ABM ####################### 
simulate_run <- function(run_id, input_df) {
  library(data.table)
  
  # Inputs
  gruppe <- input_df$group
  gruppe_population <- input_df$N
  coverage <- input_df$coverage
  
  # Calculate group measures
  N <- sum(gruppe_population)
  overall_coverage <- mean(coverage)
  population_size <- length(gruppe_population)
  group_proportions <- gruppe_population / N
  
  # Parameters
  gamma <- 1/5
  k <- 2
  t0 <- 1
  R0 <- 15
  
  # Beta matrix (based on factor matrix with within and between group transmission factors)
  epsilon <- 0.1  # Adjust epsilon to control between-group transmission (0.1 = 10%)
  factor_matrix <- matrix(epsilon / (population_size - 1), nrow = population_size, ncol = population_size)
  diag(factor_matrix) <- 1 - epsilon
  beta_matrix <- t(t(factor_matrix) * (R0 * gamma / gruppe_population))
  
  # Initial SIR-values
  I0 <- 1
  R_per_group <- setNames(round(coverage * gruppe_population), gruppe)  # Vaccinated individuals per group
  S0 <- N - I0 - sum(R_per_group)
  
  
  
  
  
  ### Create population with group assignments ###
  population <- data.table(
    ID = 1:N,
    group = sample(rep(gruppe, gruppe_population))
  )
  
  # Initialize all as S. Then assign R per group and randomly assign the initial I
  population[, tilstand := "S"]
  population[, tilstand := fifelse(.I %in% sample(.I, R_per_group[group]), "R", tilstand), by = group]
  population[tilstand == "S", tilstand := fifelse(.I %in% sample(.I, I0), "I", tilstand)]
  
  # Assign tid_tilbage for I (countdown to recovery)
  population[, tid_tilbage := NA_real_]
  population[tilstand == "I", tid_tilbage := round(rgamma(.N, shape = k, rate = k * gamma))]
  
  
  
  ### Outbreak knowledge before simulation ###
  # Where did the outbreak originate?
  outbreak_origin <- population[tilstand == "I", .(group = group[1])]
  
  # Effective reproductive number for each group
  S <- population[tilstand=="S", .N, by=group][order(group)]$N
  effective_R0 <- sapply(1:population_size, function(i) {
    sum(beta_matrix[,i] * S / gamma)
  })
  
  # Compute overall effective R using the next-generation matrix
  K <- (diag(S) %*% beta_matrix) / gamma
  eigen_values <- eigen(K)$values
  effective_R_total <- max(Re(eigen_values))  # Take the real part of the largest eigenvalue
  
  
  
  ### Simulation ###
  for (t in seq(0, 1000, by = t0)) {
    # Update tid_tilbage for I, move to R if tid_tilbage <= 0
    population[tilstand == "I", tid_tilbage := tid_tilbage - t0]
    population[tilstand == "I" & tid_tilbage <= 0, tilstand := "R"]
    
    # Break loop if no infected individuals remain
    if (population[tilstand == "I", .N] == 0) break
    
    
    # Calculate infectious counts per group
    infectious_counts <- population[tilstand == "I", .N, by = group][CJ(group = gruppe, unique = TRUE), on = "group"][is.na(N), N := 0]
    # Add infection pressure per group
    infectious_counts[, pressure := sum(infectious_counts$N * beta_matrix[.GRP, ]), by = group]
    # Add random_r and P_j to susceptible individuals
    population[tilstand == "S", `:=`(
      random_r = runif(.N),
      P_j = 1 - exp(-infectious_counts[.BY$group, pressure, on = "group"])
    ), by = group]
    # Infect susceptible individuals
    population[tilstand == "S" & random_r < P_j, tilstand := "I"]
    population[tilstand == "I" & is.na(tid_tilbage), tid_tilbage := round(rgamma(.N, shape = k, rate = k * gamma))]
    
  }
  
  
  ### Final outbreak size and groups infected ###
  outbreak_size <- (S0 + I0) - population[tilstand == "S", .N]
  
  group_infected <- population[!is.na(tid_tilbage), .N, by = group]
  group_infected <- merge(data.table(group = gruppe), group_infected, by = "group", all.x = TRUE)
  group_infected[is.na(N), N := 0]
  group_infected[, outbreak_origin := outbreak_origin$group]   # add outbreak origin
  
  
  ### Return results ###
  return(list(outbreak_size = outbreak_size, 
              group_infected = group_infected,
              outbreak_origin = outbreak_origin,
              effective_R0 = effective_R0,
              effective_R_total = effective_R_total))
}

##### Run simulations in parallel #####
parallel_run_simulations <- function(input_df, N_simulations) {
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  
  # Set reproducible RNG stream for each worker
  parallel::clusterSetRNGStream(cl, iseed = 1234)
  
  # Create a custom environment to hold input_df
  env <- new.env()
  env$input_df <- input_df
  env$simulate_run <- simulate_run
  
  # Export simulate_run and input_df from env
  parallel::clusterExport(cl, varlist = c("simulate_run", "input_df"), envir = env)
  parallel::clusterEvalQ(cl, library(data.table))
  
  results <- parallel::parLapply(cl, 1:N_simulations, function(run_id) {simulate_run(run_id, input_df)})
  
  parallel::stopCluster(cl)
  
  ## How it was run outside function:
  # cl <- makeCluster(detectCores() - 1)
  # clusterExport(cl, c("simulate_run", "input_df"))
  # clusterEvalQ(cl, library(data.table))
  # results <- parLapply(cl, 1:1000, function(run_id) simulate_run(run_id, input_df))
  # stopCluster(cl)
  
  return(results)
}

##### Test overall coverage as function #####
overall_coverage <- function(input_df) {
  return( sum(input_df$N * input_df$coverage) / sum(input_df$N) )
}


############### Calculate cumulative probability of data #################
cumulative_prob_func <- function(outbreak_sizes) {
  # Count the frequency of each outbreak size
  outbreak_size_freq <- table(outbreak_sizes)
  # Sort the outbreak sizes in ascending order
  sorted_outbreak_sizes <- sort(as.numeric(names(outbreak_size_freq)), decreasing = FALSE)
  
  # Calculate cumulative probability (proportion of outbreaks >= outbreak_size)
  reverse_cumsum <- rev(cumsum(rev(outbreak_size_freq)))
  cumulative_prob <- reverse_cumsum / length(outbreak_sizes)
  
  return(list(sorted_outbreak_sizes = sorted_outbreak_sizes, cumulative_prob = cumulative_prob))
}


####################### Calculate q(x) #######################
compute_q <- function(x, R) {
  # Compute the log of q(x) using lgamma for the factorial parts.
  log_q <- (x-1)*log(R) - (2*x-1)*log(R+1) +
    lgamma(2*x-1) - (lgamma(x+1) + lgamma(x))
  q <- exp(log_q)
  if (is.nan(q)) return(0)
  return(q)
}


####################### Calculate cumulative q(x) #######################
compute_cumulative_logtransform_q <- function(n, R) {
  if (R == 1) {
    # Handle R = 1 separately to avoid division by zero
    log_term1 <- (n - 1) * log(4 * R) - (2 * n - 1) * log(1 + R)
    log_term2 <- lgamma(n - 0.5) - (log(sqrt(pi)) + lgamma(n + 1))
    log_total <- log_term1 + log_term2
    hyper_term <- Re(hypergeo(1, n - 0.5, n + 1, 4 * R / (1 + R)^2))
    cum_q <- exp(log_total) * hyper_term
  } else {
    cum_q = max(0,1-(1/R)) + 
      exp((log(4*R)*(n-1))+lgamma(n-(1/2)) - ((log(1+R)*(2*n-1))+log(sqrt(pi))+lgamma(n+1))) * Re(hypergeo(1, n-(1/2), n+1, 4*R/((1+R)**(2))))
  }
  return(cum_q)
}

##########################################################################################
#################################### NON CENSORED DATA ###################################
##########################################################################################

####################### Calculate R #######################
calculate_R <- function(outbreak_sizes, min_outbreak_size = 1) {
  R_estimate <- 1 - (min_outbreak_size / mean(outbreak_sizes))
  
  return(R_estimate)
}

####################### Logarithm of Total Likelihood Function ####################### 
# is the product of the likelihoods for each individual outbreak
# Closed form solution
log_likelihood <- function(data) {
  n <- length(data)
  m <- mean(data)
  Rc = 1 - 1/m
  term1 <- (m * n - n) * log(Rc)
  term2 <- -(2 * m * n - n) * log(Rc + 1)
  term3 <- sum(lfactorial(2 * data - 2) - lfactorial(data) - lfactorial(data - 1))
  logL <- term1 + term2 + term3
  return(logL)
}


####################### Likelihood Ratio between two datasets ####################### 
likelihood_ratio <- function(data1,data2, alpha = 0.05){
  # log likelihood
  log_likelihood_1 = log_likelihood(data1)
  log_likelihood_2 = log_likelihood(data2)
  
  # pooled data
  pooled = c(data1,data2)
  log_likelihood_pooled = log_likelihood(pooled)
  
  # Likelihood ratio Lambda
  Lambda = -2 * (log_likelihood_pooled - log_likelihood_1 - log_likelihood_2)
  # Critical value for chi-squared distribution (alpha = 0.05)
  critical_value <- qchisq(1-alpha, df = 1)
  
  if (Lambda > critical_value) {
    print("Reject null hypothesis: Rc is different between the two datasets.")
  } else {
    print("Fail to reject null hypothesis: No significant difference in Rc.")
  }
  
  return(cat("Lambda = ",Lambda,"\nCritical value = ", critical_value))
}


####################### Bootstrap R (not censored) #######################
bootstrap_R <- function(outbreak_sizes, n_bootstraps = 1000) {
  R_estimates <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    bootstrap_sample <- sample(outbreak_sizes, size = length(outbreak_sizes), replace = TRUE)
    R_estimates[i] <- 1 - (1/mean(bootstrap_sample))
  }
  return(R_estimates)
}

####################### Calculate p-values of bootstrap #######################
calculate_p_value <- function(R_observed, R_bootstrap_estimates, R) {
  # Calculate the proportion of bootstrapped estimates as extreme as or more extreme than the observed estimate
  extreme_cases <- sum(abs(R_bootstrap_estimates - R) >= abs(R_observed - R))
  p_value <- extreme_cases / length(R_bootstrap_estimates)
  return(p_value)
}


####################### Fisher Variance #######################
# Fisher variance from second derivative of log likelihood function. 
# The formula stems from the L(R) and has validated through calculations in Maple.
Fisher_Variance <- function(R_hat,n,m){
  d2_log_likelihood = (-n*(m-1)/(R_hat**2)) + n*(2*m-1)/((R_hat+1)**2)
  1/(-d2_log_likelihood) 
}




##########################################################################################
#################################### CENSORED DATA #######################################
##########################################################################################

####################### Censoring in log likelihood ####################### 
# Nummerical solution, not closed form solution
censored_log_likelihood <- function(R, data, c) {
  if (R >= 1) return(-Inf)  # Constraint: R < 1
  
  # Compute normalization factor S = 1 - sum(q(1:c))
  sum_q_censored <- sum(sapply(1:c, function(x) compute_q(x, R)))
  S <- 1 - sum_q_censored
  if (S <= 0) return(-Inf)  # Guard against invalid S
  
  # Log-likelihood components
  term1 <- sum(sapply(data, function(x) log(compute_q(x, R))))
  term2 <- length(data) * log(S)
  return(term1 - term2)
}


####################### MLE of R censored ####################### 
# data must be a list of outbreak sizes, not a whole dataframe of informations
# c is a threshold, x>c is included in the distribution.
# uses censored_log_likelihood function above
estimate_R_censored <- function(data, c) {
  # Filter data to x > c
  censored_data <- data[data > c]
  if (length(censored_data) == 0) {
    stop("No data above censoring threshold c = ", c)
  }
  # Optimize the log-likelihood
  result <- optim(
    par = 0.5,  # Initial guess for R
    fn = function(R) -censored_log_likelihood(R, censored_data, c),  # Minimize negative log-likelihood
    method = "Brent",
    lower = 0.001,
    upper = 0.999
  )
  return(result$par)
}

####################### Fisher Information Method of uncertainty estimation #################### 
######### Compute the Hessian (second derivative) at R_hat
hessian_at_estimate <- function(R_hat, data, c, eps = 1e-5) {
  # Central difference approximation
  f <- function(R) censored_log_likelihood(R, data, c)
  (f(R_hat + eps) - 2 * f(R_hat) + f(R_hat - eps)) / (eps^2)
}


####################### Bootstrap CENSORED MLE #######################
bootstrap_R_censored <- function(outbreak_sizes, n_bootstraps = 1000, c) {
  R_estimates <- numeric(n_bootstraps)
  
  for (i in 1:n_bootstraps) {
    # Resample WITH replacement from the original censored data
    bootstrap_sample <- sample(outbreak_sizes, size = length(outbreak_sizes), replace = TRUE)
    
    # Ensure the bootstrap sample contains valid data (x > c)
    valid_bootstrap_sample <- bootstrap_sample[bootstrap_sample > c]
    
    # Skip if no valid data (unlikely if original data is large enough)
    if (length(valid_bootstrap_sample) == 0) {
      R_estimates[i] <- NA
      next
    }
    # Estimate R for the bootstrap sample
    R_estimates[i] <- estimate_R_censored(valid_bootstrap_sample, c)
  }
  # Remove NA values (failed bootstrap iterations)
  R_estimates <- na.omit(R_estimates)
  
  return(R_estimates)
}



####################### Censored Likelihood Ratio between two datasets ####################### 
likelihood_ratio_censored <- function(data1, data2, c, alpha = 0.05) {
  # Filter data to include only x > c (ensure censoring)
  data1_censored <- data1[data1 > c]
  data2_censored <- data2[data2 > c]
  
  # Estimate MLEs under separate models
  R1 <- estimate_R_censored(data1_censored, c)
  R2 <- estimate_R_censored(data2_censored, c)
  
  # Log-likelihoods for separate models
  logL1 <- censored_log_likelihood(R1, data1_censored, c)
  logL2 <- censored_log_likelihood(R2, data2_censored, c)
  logL_separate <- logL1 + logL2
  
  # Pool data and estimate MLE under the null (same R)
  pooled_data <- c(data1_censored, data2_censored)
  R_pooled <- estimate_R_censored(pooled_data, c)
  logL_pooled <- censored_log_likelihood(R_pooled, pooled_data, c)
  
  # Likelihood ratio statistic (Lambda)
  Lambda <- -2 * (logL_pooled - logL_separate)
  
  # Critical value (chi-squared, df = 1)
  critical_value <- qchisq(1 - alpha, df = 1)
  
  # Hypothesis test
  if (Lambda > critical_value) {
    print("Reject H0: R differs between datasets (censored).")
  } else {
    print("Fail to reject H0: No evidence of difference in R (censored).")
  }
  
  return(list(Lambda = Lambda, CriticalValue = critical_value))
}

####################### Save plots ####################### 
save_plot <- function(p, filename, width = 6, height = 4, dpi = 300) {
  ggsave(paste0("Figures/", filename, ".png"), plot = p, width = width, height = height, dpi = dpi)
}






###############################################################################
################################# DATA HANDELING ##############################
###############################################################################


##################### Prepare cumulative data ##################### 
# from raw data
## Removes data with at less than 20 outbreaks

prepare_cumulative_data <- function(data, min_valid_counts = 20, min_outbreak_size = 5) {
  countries <- unique(data$ReportingCountry)
  cumulative_results <- list()
  
  for (ctry in countries) {
    # Get raw outbreak counts (including <5 coded as 1)
    outbreak_sizes <- data[data$ReportingCountry == ctry, "Count"]
    
    # Calculate cumulative probability before filtering
    result <- cumulative_prob_func(outbreak_sizes)
    
    # Remove the "<5" category if it exists
    if(length(result$cumulative_prob) > 0) {
      result$cumulative_prob <- result$cumulative_prob[-1]
      result$sorted_outbreak_sizes <- result$sorted_outbreak_sizes[-1]
    }
    
    # Count actual outbreaks ≥5 (original counts before cumulative_prob_func)
    num_valid_outbreaks <- sum(outbreak_sizes >= min_outbreak_size)
    
    if(num_valid_outbreaks >= min_valid_counts) {
      cumulative_results[[ctry]] <- result
    } else {
      message("Skipping ", ctry, 
              ": only ", num_valid_outbreaks, 
              " outbreaks ≥", min_outbreak_size)
    }
  }
  
  # Combine results
  cumulative_data_all <- do.call(rbind, lapply(names(cumulative_results), function(ctry) {
    res <- cumulative_results[[ctry]]
    data.frame(
      ReportingCountry = rep(ctry, length(res$cumulative_prob)),
      sorted_outbreak_sizes = res$sorted_outbreak_sizes,
      cumulative_prob = res$cumulative_prob
    )
  }))
  
  return(cumulative_data_all)
}





############################### Bootstrap Test H0: R1 != R2 ############################### 
# Is the difference in two R estimates larger than the differences when drawn from the pooled data?

# Define bootstrap test function
bootstrap_distribution_test <- function(data1, data2, c = 4, n_boot = 1000) {
  tryCatch({
    pooled_data <- c(data1, data2)
    n1 <- length(data1)
    n2 <- length(data2)
    
    R_diff_obs <- abs(estimate_R_censored(data1, c) - estimate_R_censored(data2, c))
    
    boot_stats <- replicate(n_boot, {
      boot1 <- sample(pooled_data, n1, replace = TRUE)
      boot2 <- sample(pooled_data, n2, replace = TRUE)
      abs(estimate_R_censored(boot1, c) - estimate_R_censored(boot2, c))
    })
    
    p_value <- mean(boot_stats >= R_diff_obs, na.rm = TRUE)
    data.frame(Bootstrap_p = p_value, Bootstrap_Reject = p_value < 0.05)
    
  }, error = function(e) {
    data.frame(Bootstrap_p = NA, Bootstrap_Reject = NA)
  })
}