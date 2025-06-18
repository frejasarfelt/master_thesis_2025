# Import libraries
library(ggplot2)
library(scales)
source("functions.R")



####################### Function simulate outbreak ####################### 

# Simulate multiple outbreaks for each R value from 0.1 to 0.9
R_values <- seq(0.1, 0.9, by = 0.1)
outbreak_sizes_list <- list()
for (R in R_values) {
  outbreak_sizes <- replicate(10000, simulate_outbreak(R))
  outbreak_sizes_list[[paste("R =", R)]] <- outbreak_sizes
}


# Precompute the cumulative probability data for each R value
cumulative_data_list <- lapply(R_values, function(R) {
  outbreak_sizes <- outbreak_sizes_list[[paste("R =", R)]]
  result <- cumulative_prob_func(outbreak_sizes)
  data.frame(sorted_outbreak_sizes = result$sorted_outbreak_sizes, 
             cumulative_prob = result$cumulative_prob, R = R)
})

# Combine all the data into a single data frame
cumulative_data <- do.call(rbind, cumulative_data_list)


# Theoretical cumulative probability data
q_data <- do.call(rbind, lapply(R_values, function(R) {
  data.frame(
    x = 1:1000, 
    y = sapply(1:1000, function(n) compute_cumulative_logtransform_q(n, R)),
    R = R
  )
}))



### PLOT ###

# color palette
red_scale <- colorRampPalette(c("coral", "firebrick1", "firebrick4", "tan3", "goldenrod1"))(length(R_values))

p <- ggplot() +
  geom_step(aes(x = sorted_outbreak_sizes, y = cumulative_prob, color = factor(R)), 
            data = cumulative_data, size = 1, direction = "vh")+
  scale_y_continuous(trans = "log10", 
                     breaks = 10^(-seq(0, 3, by = 1)), 
                     labels = scales::math_format(.x), 
                     minor_breaks = NULL,
                     ) +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(0, 3), 
                     # labels = c("1", "10", "100", "1000")
                     minor_breaks = NULL,
                     ) +
  coord_cartesian(
    xlim = c(1, 1000),
    ylim = c(0.0001, 1)) +
  labs(x = "Outbreak Size", y = "Cumulative Probability", color = "R") +
  scale_color_manual(values = red_scale, name = expression(R[e])) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +   # or whatever family
  theme(
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.text            = element_text(size = 10),
    axis.title.x         = element_text(margin = margin(t = 15)),
    axis.title.y         = element_text(margin = margin(r = 15)),
    legend.title         = element_text(size = 12, hjust = 0.5)
  )+

  # Power-law curve for R = 1 as a dashed line (from 1-10 with cum(q(x)) and 10-100 with power law)
  stat_function(fun = function(n) compute_cumulative_logtransform_q(n, R=1),
                color = "black",
                size = 0.7,
                linetype = "dashed",
                xlim = c(0,1)) +
  stat_function(fun = function(x) x^(-1/2) / sqrt(pi),
                color = "black",
                size = 0.7,
                linetype = "dashed",
                xlim = c(1, 10)) +
  annotate("text", x = 100, y = 0.07, label = expression(R[e] %~~% 1), color = "black", size = 4, hjust = 0) +
  
  # Cumulative q(x) for R = 0.1 to 0.9 (theoretical)
  geom_line(aes(x = x, y = y, color = factor(R)), 
            data = q_data, size = 0.5, linetype = "dashed")

# Print and save plot  
print(p)
save_plot(p, "q(x)_cumulative")



# %-chance of having an outbreak >=100 for R=0.9
q_data[q_data$x == 100 & q_data$R == 0.9,]$y * 100
