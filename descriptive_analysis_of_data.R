# Import libraries
library(ggplot2)
library(dplyr)
source("functions.R")

### Collect data
data <- read.csv("Data/Tessy_data.csv", stringsAsFactors = FALSE)
data$MinDate <- as.Date(data$MinDate, format = "%Y-%m-%d")
data$MaxDate <- as.Date(data$MaxDate, format = "%Y-%m-%d")
data_UK <- data[data$ReportingCountry == "UK", ]

n <- nrow(data)


##### Numbers on outbreaks #####
# How many outbreaks per country
total_outbreaks_countries <- data %>%
  group_by(ReportingCountry) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  arrange(desc(Count)) %>%
  print(n = Inf)

#How many outbreaks in total
total_outbreaks <- total_outbreaks_countries %>%
  summarise(total_count = sum(Count)) %>%
  print(n = Inf)

# How many outbreaks per country, but exclude counts with less than 5
total_outbreaks_countries_filtered <- data %>%
  group_by(ReportingCountry) %>%
  filter(Count != "<5") %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  arrange(desc(Count)) %>%
  print(n = Inf)

# How many outbreaks in total, but exclude counts with less than 5
total_outbreaks_filtered <- total_outbreaks_countries_filtered %>%
  summarise(total_count = sum(Count)) %>%
  print(n = Inf)

# Which countries only have outbreaks with size <5
setdiff(unique(data$ReportingCountry),
        data %>%
          group_by(ReportingCountry) %>%
          filter(Count != "<5") %>%
          summarise(Count = n()) %>%
          ungroup() %>%
          arrange(desc(Count)) %>%
          select(ReportingCountry) %>%
          # summarise(total_count = sum(Count)) %>%
          pull(ReportingCountry)
)

# How long did the outbreaks last
outbreak_duration <- as.numeric(data$MaxDate - data$MinDate)
min_duration <- min(outbreak_duration, na.rm = TRUE)
max_duration <- max(outbreak_duration, na.rm = TRUE)

min_indexes <- which(outbreak_duration == min_duration)
max_indexes <- which(outbreak_duration == max_duration)

data[min_indexes, ] %>%
  select(ReportingCountry, Count) %>%
  distinct() %>%
  arrange(Count)
data[max_indexes, ] %>%
  select(ReportingCountry, MinDate, MaxDate, Count)

# Find second longest outbreak duration
sort(outbreak_duration, decreasing = TRUE) #2 --> 690 days
which(outbreak_duration == 690)            #  --> 255 index
data[255, ] %>%
  select(ReportingCountry, MinDate, MaxDate, Count)

# Median duration of outbreaks
median_duration <- median(outbreak_duration, na.rm = TRUE)




##### Numbers on age #####
# Normalized age group counts
age_groups <- c("X0.4y", "X5.14y", "X15.25y", "X26.50y", "X50.")
age_labels <- c("0-4", "5-14", "15-25", "26-50", "50+")
age_counts <- colSums(data[age_groups] == "Y")
age_perc <- round(age_counts / n * 100,2)
age_spans <- c(5, 10, 11, 25, 30)  # adjust 30 if you want a different upper bound
age_density <- age_counts / age_spans
age_density_perc <- age_perc / age_spans

# Create a data frame
age_df <- data.frame(
  AgeGroup = gsub("X", "", age_groups),
  Count = as.integer(age_counts),
  Span = age_spans,
  PerYearDensity = age_density,
  PerYearDensityPerc = age_density_perc)
age_df$AgeGroup <- factor(age_df$AgeGroup, levels = gsub("X", "", age_groups), labels = age_labels)


# Plot: normalized count per year
p <- ggplot(age_df, aes(x = AgeGroup, y = PerYearDensityPerc)) +
  geom_bar(stat = "identity", fill = "firebrick3", width = 0.7) +
  labs(x = "Age Group",
       y = "Prevalence % (normalized)") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor = element_blank()     # Remove all minor grid lines
  )
print(p)
save_plot(p, "Agegroup_prevalence", height = 3)









##### Numbers on vaccination status for UK #####
# Vaccination group counts
vaccination_groups <- c("X0_dose", "X1_dose", "X..2_dose", "UNKDOSE", "Unknown")
vaccination_labels <- c("0 dose", "1 dose", "2 doses", "Unknown dose", "Unknown")
vaccination_counts <- colSums(data[vaccination_groups] == "Y")

# Create a data frame
vaccination_df <- data.frame(
  VaccinationGroup = gsub("X", "", vaccination_groups),
  Count = as.integer(vaccination_counts))
vaccination_df$VaccinationGroup <- factor(vaccination_df$VaccinationGroup, levels = gsub("X", "", vaccination_groups), labels = vaccination_labels)


# Plot: count per year
p <- ggplot(vaccination_df, aes(x = VaccinationGroup, y = Count)) +
  geom_bar(stat = "identity", fill = "firebrick3", width = 0.7) +
  labs(x = "Vaccination Group",
       y = "Prevalence") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15))
  )
print(p)
save_plot(p, "vaccinationgroup_prevalence", height = 2)




# How many people were vaccinated (1 dose, 2 doses or more, unknown number of doses)
total_vaccinated <- vaccination_counts["X1_dose"] +
  vaccination_counts["X..2_dose"] +
  vaccination_counts["UNKDOSE"]
not_vaccinated <- vaccination_counts["X0_dose"]
unknown_vaccinated <- vaccination_counts["Unknown"]
unknown_vaccinated <- unknown_vaccinated + (data %>% filter(if_all(all_of(vaccination_groups), ~ . == "N")) %>% nrow())


# Calculate the percentage of vaccinated individuals
vaccinated_percentage <- (total_vaccinated / n) * 100
not_vaccinated_percentage <- (not_vaccinated / n) * 100
unknown_vacc_status <- unknown_vaccinated / n * 100

cat(paste0("Percentage of vaccinated individuals: ", round(vaccinated_percentage, 2), "%\n",
          "Percentage of not vaccinated individuals: ", round(not_vaccinated_percentage, 2), "%\n",
          "Percentage of unknown vaccinated individuals: ", round(unknown_vacc_status, 2), "%"))


# Calculate the percentage of vaccinated individuals in <5 outbreaks
total_vaccinated_small_outbreaks <- nrow(data[data[vaccination_groups[2:4]] == "Y" & data$Count == "<5",])
n_small_outbreaks <- nrow(data[data$Count == "<5",])

total_vaccinated_small_outbreaks / n_small_outbreaks * 100








################### AGE PER COUNTRY ###################

# the list of countries you care about
countries <- c("UK", "PL", "ES", "AT", "IE")

# your age-group definitions
age_groups <- c("X0.4y", "X5.14y", "X15.25y", "X26.50y", "X50.")
age_labels <- c("0-4", "5-14", "15-25", "26-50", "50+")
age_spans  <- c(5, 10, 11, 25, 30)

# this will accumulate each country’s age_df in a list
age_dfs_list <- lapply(countries, function(ctry) {
  
  # filter to that country
  data_ctry <- subset(data, ReportingCountry == ctry)
  n <- nrow(data_ctry)
  if(n == 0) return(NULL)
  
  # same calculations you had
  age_counts       <- colSums(data_ctry[ , age_groups] == "Y")
  age_perc         <- round(age_counts / n * 100, 2)
  age_density      <- age_counts / age_spans
  age_density_perc <- age_perc   / age_spans
  
  # build the per-country data frame
  df <- data.frame(
    Country            = ctry,
    AgeGroup           = factor(gsub("X", "", age_groups),
                                levels = gsub("X", "", age_groups),
                                labels = age_labels),
    Count              = as.integer(age_counts),
    Span               = age_spans,
    PerYearDensity     = age_density,
    PerYearDensityPerc = age_density_perc
  )
  
  # also assign it globally if you really want data_UK, data_PL, etc.
  assign(paste0("age_df_", ctry), df, envir = .GlobalEnv)
  
  return(df)
})

# combine into one big DF (optional)
age_df_all <- do.call(rbind, age_dfs_list)

age_df_all[age_df_all$AgeGroup=="0-4",]


