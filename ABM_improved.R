library(data.table)
library(parallel)
library(reshape2)
library(ggplot2)
library(dplyr)
source("functions.R")



### INPUT SCENARIOS ###
scenario_1 <- data.frame(
  group = c("A", "B", "C", "D"),
  N = c(2500,2500,2500, 2500),
  coverage = c(0.95, 0.95, 0.95, 0.95),
  scenario = 1)

scenario_2 <- data.frame(
  group = c("A", "B", "C", "D"),
  N = c(2500,2500,2500, 2500),
  coverage = c(0.99, 0.97, 0.93, 0.91),
  scenario = 2)

scenario_3 <- data.frame(
  group = c("A", "B"),
  N = c(9000, 1000),
  coverage = c(0.972, 0.752),
  scenario = 3)


# Observe the overall vaccination coverage for each scenario
overall_coverage(scenario_1)
overall_coverage(scenario_2)
overall_coverage(scenario_3)



### RUN SIMULATIONS FOR ALL SCENARIOS ###
N_sim <- 1000
results_scenario_1 <- parallel_run_simulations(scenario_1, N_sim)
results_scenario_2 <- parallel_run_simulations(scenario_2, N_sim)
results_scenario_3 <- parallel_run_simulations(scenario_3, N_sim)

# Add a tag to results that indicates the scenario
results_scenario_1 <- lapply(results_scenario_1, function(x) {x$scenario <- 1; return(x)})
results_scenario_2 <- lapply(results_scenario_2, function(x) {x$scenario <- 2; return(x)})
results_scenario_3 <- lapply(results_scenario_3, function(x) {x$scenario <- 3; return(x)})

# Gather the results into one data frame
results <- c(results_scenario_1, results_scenario_2, results_scenario_3)


### Extract results ###
outbreak_info_dt <- data.table(outbreak_size = sapply(results, `[[`, "outbreak_size"),
                                outbreak_origin = sapply(results, `[[`, "outbreak_origin"),
                               scenario = sapply(results, `[[`, "scenario"))

group_infected_dt <- rbindlist(lapply(seq_along(results), function(i) {
  dt <- results[[i]]$group_infected
  dt$run_id <- i
  dt$scenario <- results[[i]]$scenario
  return(dt)}) )

### Subpopulation R_e
R_eff_1 <- data.table(
  group = scenario_1$group,
  effective_R0 = rowMeans(do.call(cbind, lapply(results_scenario_1, `[[`, "effective_R0"))),
  scenario = 1)
R_eff_2 <- data.table(
  group = scenario_2$group,
  effective_R0 = rowMeans(do.call(cbind, lapply(results_scenario_2, `[[`, "effective_R0"))),
  scenario = 2)
R_eff_3 <- data.table(
  group = scenario_3$group,
  effective_R0 = rowMeans(do.call(cbind, lapply(results_scenario_3, `[[`, "effective_R0"))),
  scenario = 3)
R_eff <- rbindlist(list(R_eff_1, R_eff_2, R_eff_3))

### Total R_e
R_eff_total_1 <- data.table(
  effective_R_total = rowMeans(do.call(cbind, lapply(results_scenario_1, `[[`, "effective_R_total"))),
  scenario = 1)
R_eff_total_2 <- data.table(
  effective_R_total = rowMeans(do.call(cbind, lapply(results_scenario_2, `[[`, "effective_R_total"))),
  scenario = 2)
R_eff_total_3 <- data.table(
  effective_R_total = rowMeans(do.call(cbind, lapply(results_scenario_3, `[[`, "effective_R_total"))),
  scenario = 3)
R_eff_total <- rbindlist(list(R_eff_total_1, R_eff_total_2, R_eff_total_3))






############# PLOTTING #############
# Plot only scenario x
scenario_x <- 3


# Define a red color scale based on the number of unique groups
red_scale <- colorRampPalette(c("coral", "firebrick1", "firebrick4"))(length(unique(group_infected_dt$group)))

# Histogram of outbreak sizes with y-axis scaled from 0 to 1
p <- ggplot(outbreak_info_dt[scenario == scenario_x], aes(x = outbreak_size, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "firebrick3", color = "firebrick4", alpha = 0.7) +
  labs(
    x = "Outbreak Size",
    y = "Density") +
  theme_minimal() +
  theme(
        text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        )


### SAVE DIFFERENT SCENARIOS
# save_plot(p, "ABM_homogenitet")
# save_plot(p, "ABM_heterogenitet_4grupper")
# save_plot(p, "ABM_heterogenitet_2grupper")

print(p)




#############
## Plot the step function of the cumulative distribution function of outbreak sizes
cdf_by_scenario <- rbindlist(lapply(unique(outbreak_info_dt$scenario), function(scn) {
  result <- cumulative_prob_func(outbreak_info_dt[scenario == scn]$outbreak_size)
  data.table(
    sorted_outbreak_sizes = result$sorted_outbreak_sizes,
    cumulative_prob = result$cumulative_prob,
    scenario = as.factor(scn)
  )
}))

p <- ggplot() +
  geom_step(aes(x = sorted_outbreak_sizes, y = cumulative_prob, color = scenario), 
            data = cdf_by_scenario, size = 1, direction = "vh") +
  scale_y_continuous(trans = "log10", 
                     breaks = 10^(-seq(0, 3, by = 1)), 
                     labels = scales::math_format(.x)) +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(0, 2)) +
  coord_cartesian(xlim = c(1, 1000),
                  ylim = c(0.001, 1)) +
  labs(x = "Outbreak Size", 
       y = "Cumulative Probability",
       color = "Scenario") +
  scale_color_manual(values = red_scale, name = "Scenario") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p)

### SAVE STEP FUNCTION PLOT
# save_plot(p, "ABM_cumulative_step_plot")


# Find the outbreak sizes where 90% of the cumulative probability is reached
head(cdf_by_scenario[scenario == 1],10) # 0.116 -> 5
head(cdf_by_scenario[scenario == 2],20) # 0.099 -> 14
cdf_by_scenario[scenario == 3] # 0.108 -> 313

# Find the cumulative_prob for sorted_outbreak_sizes = 1 in scenario = 1 in the cdf_by_scenario
cdf_by_scenario[scenario == 1 & sorted_outbreak_sizes == 1, cumulative_prob] - cdf_by_scenario[scenario == 1 & sorted_outbreak_sizes == 2, cumulative_prob] # =0.613
cdf_by_scenario[scenario == 2 & sorted_outbreak_sizes == 1, cumulative_prob] - cdf_by_scenario[scenario == 2 & sorted_outbreak_sizes == 2, cumulative_prob] # =0.543
cdf_by_scenario[scenario == 3 & sorted_outbreak_sizes == 1, cumulative_prob] - cdf_by_scenario[scenario == 3 & sorted_outbreak_sizes == 2, cumulative_prob] # =0.494



######################## Plot PDF: ########################  
pmf_all <- outbreak_info_dt %>%
  # Step 1: count per (scenario, outbreak_size)
  count(scenario, outbreak_size) %>%          # creates columns: scenario, outbreak_size, n
  # Step 2: for each scenario, compute prob = n / sum(n)
  group_by(scenario) %>%
  mutate(prob = n / sum(n)) %>%
  ungroup()


# 2) Re-create your red_scale (one color per scenario)
#    (or pick any palette you like)
n_scenarios <- length(unique(pmf_all$scenario))
red_scale    <- colorRampPalette(c("coral", "firebrick1", "firebrick4"))(n_scenarios)


ggplot(pmf_all, aes(x = outbreak_size, y = prob, color = factor(scenario))) +
  #geom_line(size = 1) +
  geom_point(size = 2, alpha = 0.8) +
  scale_y_log10(
    breaks = 10^(-seq(0, 3, by = 1)),
    labels = scales::math_format(.x),
    limits = c(min(pmf_all$prob[pmf_all$prob > 0]), 1)
  ) +
  scale_x_log10(
    breaks = 10^seq(0, 2),
    labels = scales::comma
  ) +
  coord_cartesian(xlim = c(1, 1000)) +
  scale_color_manual(values = red_scale, name = "Scenario") +
  labs(
    x = "Outbreak Size",
    y = "Probability (log₁₀ scale)"
  ) +
  theme_minimal() +
  theme(
    text         = element_text(size = 12),
    axis.text    = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(hjust = 0.5)
  )









###### OTHER PLOTS TO VISUALIZE THE DIFFERENCE BETWEEN GROUPS INFECTIONS ######
## Boxplot of number of infected in each group
p <- ggplot(group_infected_dt[scenario == scenario_x], aes(x = group, y = N)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Fordeling af smittede i hver gruppe",
       x = "Gruppe",
       y = "Antal smittede") +
  theme_minimal()


## Histogram of number of infected in each group
p <- ggplot(group_infected_dt[scenario == scenario_x], aes(x = N, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Histogram over antal smittede i hver gruppe",
       x = "Antal smittede",
       y = "Frekvens") +
  # scale_fill_manual(values = c("red", "blue")) +
  theme_minimal()




## Facetted scatter plot of number of infected in each group
group_long <- merge(group_infected_dt[scenario == scenario_x], 
                    group_infected_dt[scenario == scenario_x], 
                    by = "run_id", 
                    allow.cartesian = TRUE)
group_long <- group_long[group.x != group.y]  # Fjern sammenligning af samme gruppe
group_long <- group_long[group.x < group.y]     # Beholder kun unikke gruppepar


p <- ggplot(group_long, 
            aes(x = N.x, y = N.y)) +
            # aes(x = perc.x, y = perc.y)) +
  geom_point(alpha = 0.5, fill = "firebrick3", color = "firebrick4", shape = 21, size = 3) + 
  facet_grid(group.y ~ group.x, scales = "free", switch = "both") +  
  labs(x = "Subpopulation Infections",
       y = "Subpopulation Infections") +
  # labs(x = "Procent smittede i gruppe",
  #      y = "Procent smittede i gruppe") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(size = 12, face = "bold"),  
    strip.text.y = element_text(size = 12, face = "bold"),  
    strip.placement = "outside",
    panel.grid.minor = element_blank()
  )

print(p)

### SAVE PLOT OF FACETTED SCATTERPLOT
# save_plot(p, "ABM_facetted_scatterplot_scn3")


# the max value of small outbreaks
max(group_long[group_long$N.x < 25]$N.x)
max(group_long[group_long$N.y < 50]$N.y)




##### INSIGHTS ON OUTBREAK ORIGIN #####
# Create all possible combinations of group and outbreak_origin
all_combinations <- CJ(group = unique(group_infected_dt[scenario == scenario_x]$group), 
                       outbreak_origin = unique(group_infected_dt[scenario == scenario_x]$outbreak_origin))
group_spread_count <- group_infected_dt[N > 0 & scenario == scenario_x, .N, by = .(group, outbreak_origin)]
total_counts <- group_spread_count[group==outbreak_origin, .(outbreak_origin, N)]
group_spread_count <- merge(all_combinations, group_spread_count, by = c("group", "outbreak_origin"), all.x = TRUE)
group_spread_count[is.na(N), N := 0]
group_spread_count <- group_spread_count[outbreak_origin != group]
group_spread_count <- merge(group_spread_count, total_counts, by = "outbreak_origin", all.x = TRUE)
group_spread_count[, fraction := N.x / N.y]

# Define a red color scale based on the number of unique groups
red_scale <- colorRampPalette(c("coral", "firebrick1", "firebrick4"))(length(unique(group_spread_count$group)))


p <- ggplot(group_spread_count, aes(x = outbreak_origin, y = fraction, fill = as.factor(group))) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.9) +
  scale_fill_manual(values = red_scale, name = "Infected group") +
  geom_text(aes(label = group),  # Puts text slightly below the bars
            position = position_dodge(width = 0.9),
            vjust = -1,
            size = 3
            ) +
  labs(x = "Origin group", y = "Fraction of times infected") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 10, hjust = 0.5)
  )


print(p)





# Calculate the fraction of infected individuals in each group
infected_by_origin <- group_infected_dt[N > 0 & scenario == scenario_x, .(total_infected = sum(N)), by = .(outbreak_origin, group)]
total_infected_origin <- infected_by_origin[, .(total = sum(total_infected))]
infected_by_origin[, fraction := total_infected / total_infected_origin$total]

red_scale <- colorRampPalette(c("coral", "firebrick1", "firebrick4"))(length(unique(group_infected_dt$group)))


p <- ggplot(infected_by_origin, aes(x = outbreak_origin, y = fraction, fill = group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.9, width = 0.5) +
  scale_fill_manual(values = red_scale, name = "Group") +
  labs(
    x = "Outbreak origin",
    y = "Fraction of infections",
  ) +
  theme_minimal() +
  theme( #standing alone in Latex
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) #+
  # theme( #standing side by side in Latex
  #   text = element_text(size = 20),
  #   axis.text = element_text(size = 16),
  #   axis.title.x = element_text(margin = margin(t = 15)),
  #   axis.title.y = element_text(margin = margin(r = 15)),
  #   legend.title = element_text(size = 20),
  #   panel.grid.minor = element_blank(),
  #   panel.grid.major.x = element_blank()
  # ) +
  # guides(fill = "none")
  # ylab(NULL)

print(p)

### SAVE PLOT OF INFECTED BY ORIGIN
# save_plot(p, "ABM_infected_by_origin_scn1")
# save_plot(p, "ABM_infected_by_origin_scn2")
# save_plot(p, "ABM_infected_by_origin_scn3")


# Scenario 2: Find the relative fraction of infections that originated from group D and did not end in group D
with(infected_by_origin, sum(total_infected[outbreak_origin == "D" & group != "D"]) / sum(total_infected[outbreak_origin == "D"]))

# Scenario 3: Find the fraction of infections that originated from group A vs B
with(infected_by_origin, sum(fraction[outbreak_origin == "A"]))
with(infected_by_origin, sum(fraction[outbreak_origin == "B"]))

# Scenario 3: Originating in A --> infections in B (and vice versa)
with(infected_by_origin, sum(total_infected[outbreak_origin == "A" & group == "B"]) / sum(total_infected[outbreak_origin == "A"]))
with(infected_by_origin, sum(total_infected[outbreak_origin == "B" & group == "A"]) / sum(total_infected[outbreak_origin == "B"]))



# # TEST AT BETA_MATRIX GÅR OP MED R0
# # Compute the next-generation matrix (NGM)
# G <- beta_matrix * matrix(gruppe_population, nrow = m, ncol = m, byrow = TRUE) / gamma
# 
# # Calculate eigenvalues and extract the largest one (spectral radius)
# eigenvalues <- eigen(G, only.values = TRUE)$values
# spectral_radius <- max(Re(eigenvalues))
