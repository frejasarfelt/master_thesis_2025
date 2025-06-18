
#########################################################################   
######################## Plot of q(x) for different R ######################## 
######### q(x): the probability that the final size of the outbreak is x

library(ggplot2)
library(dplyr)
library(tidyr)
source("functions.R")



# Define R values
R_values <- seq(0.1, 0.9, by = 0.1) # R from 0.1 to 0.9


# Compute q(x) for each R
q_data <- expand.grid(x = 1:1000, R = R_values) %>%
  mutate(qx = mapply(compute_q, x, R)) %>%
  mutate(qx = ifelse(is.nan(qx), 0, qx)) %>%
  mutate(qx = ifelse(is.infinite(qx), 0, qx)) %>%
  filter(qx > 0)



# Generate simulation data
simulation_data <- do.call(rbind, lapply(R_values, function(R) {
  outbreak_sizes <- replicate(100000, simulate_outbreak(R))
  data.frame(R = R, outbreak_size = outbreak_sizes)
}))

# Convert simulation data into probabilities
simulated_probs <- simulation_data %>%
  group_by(R, outbreak_size) %>%
  summarise(prob = n()/100000, .groups = "drop")



# Convert R to factor for color consistency
q_data$R <- factor(q_data$R)
simulated_probs$R <- factor(simulated_probs$R)



### PLOT ###
# color palette
red_scale <- colorRampPalette(c("coral", "firebrick1", "firebrick4", "tan3", "goldenrod1"))(length(R_values))

# Plot with ggplot2
p <- ggplot() +
  geom_line(data = q_data, aes(x = x, y = qx, color = R), size = 1)+
  geom_point(data = simulated_probs, aes(x = outbreak_size, y = prob, color = R), size = 1.5) +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(-seq(0, 3, by = 1)),
                     minor_breaks = NULL,             
                     # labels = scales::trans_format("log10", scales::math_format(10^.x))
                     labels = scales::math_format(.x)
                     ) +
  #xlim(1, 100) +
  scale_x_continuous(
    limits     = c(0, 100),
    expand     = c(0, 0),                # no padding below 1
    breaks     = c(1,seq(20, 100, by = 20)),   # only these ticks
    minor_breaks = NULL                  # no minor x–gridlines
  ) +
  coord_cartesian(ylim = c(0.0001, 1)) +
  labs(x = "Outbreak Size", y = "Probability", color = "R Value") +
  theme_minimal() +
  theme(
    # drop vertical major grid lines entirely
    panel.grid.major.x   = element_blank(),
    # keep only horizontal majors (already there), no minors anywhere
    panel.grid.minor     = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 12, hjust = 0.5)
  ) +
  scale_color_manual(values = red_scale, name = expression(R[e]))
  # scale_color_viridis_d(option = "inferno")

print(p)


save_plot(p, "q(x)_simulations")









########################### Uncertaincy depending on sample size ########################


library(parallel)
library(ggplot2)
library(tidyr)
library(dplyr)
source("functions.R")

# Parameters
true_R       <- 0.5
sample_sizes <- c(10, 20, 50, 100, 200, 300, 500)
sim_count    <- 10000
boot_B       <- 500

# Simulate outbreak pool
pool_size <- 100000
outbreak_pool <- replicate(pool_size, simulate_outbreak(true_R))

# Prepare results data frame
results <- data.frame(
  n               = sample_sizes,
  avg_R_hat       = NA,
  empirical_SE    = NA,
  avg_fisher_SE   = NA,
  avg_boot_SE     = NA,
  fisher_coverage = NA,
  boot_coverage   = NA
)

# Main loop
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  cat("Processing n =", n, "\n")
  
  stats <- mclapply(1:sim_count, function(j) {
    tryCatch({
      xs <- sample(outbreak_pool, n)
      m_s <- mean(xs)
      R_hat <- 1 - 1 / m_s 
      
      # Fisher SE
      fisher_SE = sqrt(Fisher_Variance(R_hat,n,m_s))
      
      # Bootstrap
      idx_mat <- matrix(sample.int(n, n * boot_B, replace = TRUE), nrow = boot_B)
      xb_means <- rowMeans(matrix(xs[idx_mat], nrow = boot_B, ncol = n))
      Rb <- 1 - 1 / xb_means
      boot_SE <- sd(Rb, na.rm = TRUE)
      
      # CIs and coverage
      ci_fisher <- R_hat + c(-1.96, 1.96) * fisher_SE
      ci_boot <- quantile(Rb, c(0.025, 0.975), na.rm = TRUE)
      
      fisher_cov <- as.numeric(ci_fisher[1] <= true_R && true_R <= ci_fisher[2])
      boot_cov <- as.numeric(ci_boot[1] <= true_R && true_R <= ci_boot[2])
      
      
      c(R_hat = R_hat, fisher_SE = fisher_SE, boot_SE = boot_SE,
        fisher_cov = fisher_cov, boot_cov = boot_cov)
    }, error = function(e) {
      return(rep(NA, 5))
    })
  }, mc.cores = detectCores() - 1)
  
  stats <- do.call(rbind, stats)
  stats <- stats[complete.cases(stats), , drop = FALSE]
  
  results[i, ] <- c(
    n,
    mean(stats[, "R_hat"]),
    sd(stats[, "R_hat"]),
    sqrt(mean(stats[, "fisher_SE"]^2)),
    sqrt(mean(stats[, "boot_SE"]^2)),
    mean(stats[, "fisher_cov"]),
    mean(stats[, "boot_cov"])
  )
}

print(results)

# Plotting
se_df <- results %>%
  select(n, avg_fisher_SE, avg_boot_SE) %>%
  pivot_longer(-n, names_to = "method", values_to = "SE")

p1 <- ggplot(se_df, aes(x = n, y = SE, color = method)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  #  geom_abline(intercept = log10(0.6), slope = -0.5, linetype = "dashed", color = "#1f77b4") +
  scale_x_log10() +
  scale_y_log10() +
  scale_y_log10(limits = c(0.02, 0.3)) +
  labs(
    x = "Sample Size (n)",
    y = "Standard Error",
    color = "Estimator"
  ) +
  scale_color_manual(
    values = c(
      #"empirical_SE" = "black",
      "avg_fisher_SE" = "#1f77b4",
      "avg_boot_SE" = "#ff7f0e"
    ),
    labels = c(
      #"empirical_SE" = "Empirical",
      "avg_fisher_SE" = "Fisher",
      "avg_boot_SE" = "Bootstrap"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p1)
save_plot(p1, "SE_vs_R",width = 5, height = 4)

########### SLOPE OF FISHER ########### 

# extract just the Fisher SE and sample sizes
fisher_df <- results %>% 
  select(n, avg_fisher_SE)

# Fit linear model on log10‐scale
lm_fit <- lm(log10(avg_fisher_SE) ~ log10(n), data = fisher_df)

# Look at the slope estimate and its 95% CI
summary(lm_fit)$coefficients





# Coverage probability vs. sample size
cov_df <- results %>% 
  select(n, fisher_coverage, boot_coverage) %>%
  tidyr::pivot_longer(-n, names_to="method", values_to="coverage")

p2 <- ggplot(cov_df, aes(x = n, y = coverage, color = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3) +
  labs(
    #title = "Coverage of 95% CIs vs Sample Size",
    x = "Sample Size (n)", 
    y = "Coverage Probability",
    color = NULL
  ) +
  scale_color_manual(
    values = c("fisher_coverage" = "#1f77b4", "boot_coverage" = "#ff7f0e"),
    labels = c("fisher_coverage" = "Fisher",
               "boot_coverage" = "Bootstrap")
  ) +
  ylim(0.7, 1) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )


print(p2)
save_plot(p2, "Coverage_CI_vs_n",width = 7, height = 4)









