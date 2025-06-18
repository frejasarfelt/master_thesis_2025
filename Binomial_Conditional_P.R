library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)
library(binom)
library(dplyr)


######### Data 

age_groups <- c("X0.4y","X5.14y","X15.25y","X26.50y","X50.")
age_labels <- c(
  "X0.4y"   = "0–4",
  "X5.14y"  = "5–14",
  "X15.25y" = "15–25",
  "X26.50y" = "26–50",
  "X50."    = "50+"
)

widths     <- c(X0.4y=5, X5.14y=10, X15.25y=11, X26.50y=25, X50.=30)

data <- read_csv("Data/Tessy_data_cleaned.csv") %>% 
  mutate(across(all_of(age_groups),
                ~ case_when(. == "Y" ~ 1L,
                            . == "N" ~ 0L,
                            TRUE      ~ NA_integer_))) %>% 
  filter(ReportingCountry %in% c("ES", "PL", "UK", "AT", "IE"))

# Conditional probabilities: Use ALL data (Not depending on N_i)
cond_data = data  # No filtering

# Binomial model: Exclude Count <5 (excact N_i are needed)
data <- data %>% filter(Count >= 5)

# Define the mapping
country_names <- c(
  "ES" = "Spain",
  "PL" = "Poland",
  "UK" = "United Kingdom",
  "AT" = "Austria",
  "IE" = "Ireland"
)

########## Binomial model

# by country
country_nest <- data %>%
  group_by(ReportingCountry) %>%
  nest()

# Fit MLEs inside each country
country_p <- country_nest %>%
  mutate(
    p_df = map(data, ~ {
      df <- .
      map_dfr(age_groups, function(ag) {
        # Negative log-likelihood function
        fn <- function(p) {
          ll <- ifelse(df[[ag]] == 1,
                       log(1 - (1 - p)^df$Count),
                       log((1 - p)^df$Count))
          -sum(ll, na.rm = TRUE) #negative log-likelihood to be minimized by solving for p
        }
        
        # Find MLE
        opt <- optim(0.5, fn, method = "Brent", lower = 1e-9, upper = 1-1e-9)
        p_mle <- opt$par
        min_nll <- opt$value
        
        # Likelihood ratio CI (95% confidence)
        lr_threshold <- min_nll + qchisq(0.95, 1)/2
        
        # Find lower bound
        lower <- tryCatch({
          uniroot(function(p) fn(p) - lr_threshold,
                  interval = c(1e-9, p_mle))$root
        }, error = function(e) NA)
        
        # Find upper bound
        upper <- tryCatch({
          uniroot(function(p) fn(p) - lr_threshold,
                  interval = c(p_mle, 1-1e-9))$root
        }, error = function(e) NA)
        
        # Return results
        tibble(
          AgeBand = ag,
          p_j      = p_mle,
          lower    = lower,
          upper    = upper
        ) %>%
          mutate(
            Width      = widths[ag],
            PerYear   = p_j / Width,
            LowerYear = lower / Width,
            UpperYear = upper / Width
          )
      })
    })
  ) %>%
  unnest(p_df) %>%
  mutate(
    AgeBand = factor(AgeBand, 
                     levels = age_groups,
                     labels = age_labels)  
  )

# Investigate sums of p_j across countries:
country_p_sums <- country_p %>%
  group_by(ReportingCountry) %>%
  summarize(
    sum_pj = sum(p_j, na.rm = TRUE)
  )

country_p_sums




# Calculate the percentage of 0-4 year olds for each country
percent_0_4 <- country_p %>%
  group_by(ReportingCountry) %>%
  summarise(
    percent_0_4 = 100 * (PerYear[AgeBand == "0–4"] / sum(PerYear)),
    .groups = "drop"
  )

print(percent_0_4)





###### Conditional probabilities

country_cond <- cond_data %>%
  group_by(ReportingCountry) %>%
  nest() %>%
  mutate(
    cond_df = map(data, ~{ # 'data' is the nested data table from country_cond
      df <- . #. is the current element in the map iteration, df is the nested dataframe for each country
      total_outbreaks <- nrow(df)
      
      # new columns Age_J and Age_K in tibble with pairs
      expand_grid(Age_J = age_groups, Age_K = age_groups) %>%
        filter(Age_J != Age_K) %>%
        rowwise() %>%
        mutate(
          n_j       = sum(df[[Age_J]] == 1, na.rm = TRUE), # How many of Age_J in that country
          n_k       = sum(df[[Age_K]] == 1, na.rm = TRUE), # How many of Age_K in that country
          n_both    = sum(df[[Age_J]] == 1 & df[[Age_K]] == 1, na.rm = TRUE), # How many of both in that country
          n_neither = total_outbreaks - n_j - n_k + n_both,
          # Chi2 test with simulated p value:
          test_result = list(
            chisq.test(
              matrix(
                c(n_both,
                  n_j - n_both,
                  n_k - n_both,
                  n_neither),
                2, 2
              ),
              simulate.p.value = TRUE,
              B = 2000
            )
          )
        ) %>%
        ungroup() %>%
        mutate(
          p_value = map_dbl(test_result, "p.value"),
          p_adj   = p.adjust(p_value, method = "bonferroni"),
          sig     = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ ""
          ),
          
          # actual conditional probabilies and difference Delta:
          P_KgJ  = n_both / n_j,
          P_K    = n_k / total_outbreaks,
          Enrich = P_KgJ - P_K
        ) %>%
        select(Age_J, Age_K, P_KgJ, P_K, Enrich, p_value, p_adj, sig)
    })
  ) %>%
  select(-data) %>%     # drop the original nested data
  unnest(cond_df) %>%   # bring your Age_J, Age_K, Enrich, etc. into the main flat tibble
  mutate(
    Age_J = factor(Age_J, levels = age_groups),
    Age_K = factor(Age_K, levels = age_groups),
    ReportingCountry = recode(ReportingCountry, !!!country_names)
  )






####### Plotting


# Common y-axis limits across all plots
y_limit <- range(country_p$UpperYear, na.rm = TRUE)

# Replace codes with names
country_p <- country_p %>%
  mutate(ReportingCountry = recode(ReportingCountry, !!!country_names))

# Country order for reproducibility
countries <- unique(country_p$ReportingCountry)
plots <- map(seq_along(countries), function(i) {
  ctry <- countries[i]
  df <- filter(country_p, ReportingCountry == ctry)
  
  ggplot(df, aes(x = AgeBand, y = PerYear)) +
    geom_col(fill = "grey",alpha = 1) +
    geom_errorbar(aes(ymin = LowerYear, ymax = UpperYear),
                  width = 0.2,color = "red3") +
    scale_x_discrete(labels = age_labels) +
    coord_cartesian(ylim = y_limit) +
    labs(
      title = ctry,
      x     = NULL,
      y     = if (i %% 3 == 1) "Probability (per year)" else NULL  # y-axis only on left
    ) +
    theme_minimal() +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "lightgrey"),
      axis.ticks        = element_line(color = "lightgrey"),
      plot.title        = element_text(hjust = 0.5, size = 12),
      axis.title.y      = if (i %% 3 == 1) element_text() else element_blank(),
      axis.text.y       = if (i %% 3 == 1) element_text() else element_blank()
    )
})

# Assemble in 3+2 layout
binom_p <-(plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )
binom_p




######## Enrichment heatmaps by country

# Replace codes with names
country_cond <- country_cond %>%
  mutate(ReportingCountry = recode(ReportingCountry, !!!country_names))

ggplot(country_cond, aes(x = Age_K, y = Age_J, fill = Enrich)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sig),
            size = 3,            # adjust point size as you like
            color = "grey38",     # or white if your fill is dark
            vjust = 0.5,         # center vertically
            hjust = 0.5) +       # center horizontally
  #geom_text(aes(label = sig), size = 3.5, vjust = 0.8) +
  scale_x_discrete(breaks = age_groups, labels = age_labels, drop = FALSE) +
  scale_y_discrete(breaks = age_groups, labels = age_labels, drop = FALSE) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = bquote("  " * Delta)
  ) +
  facet_wrap(~ ReportingCountry, scales = "free", drop = FALSE) +
  labs(
    x = "Target age band (K)",
    y = "Given age band (J)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.grid = element_blank()
  )



