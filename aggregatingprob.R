
p_month <- 0.0005

# Number of years
years <- 20

# Probability at least one outbreak in 20 years
p_20years <- 1 - (1 - p_month)^(12 * years)

p_20years


# Calculate 20-year probability for each district
district_data$p_20years <- 1 - (1 - district_data$p_month)^(12 * years)

# Show results
district_data

split_df <- split(pred_df, pred_df$year)



df <- split_df %>%
  group_by(district_country) %>% 
  summarise(
    # Product of probabilities of no outbreak across all months
    p_no_outbreak_20yrs = prod(1 - fitted, na.rm = TRUE)
  ) %>%
  mutate(
    # Probability of at least one outbreak in 20 years
    p_20years = 1 - p_no_outbreak_20yrs
  )
yearly_probs <- pred_df %>%
  group_by(district_country, year) %>%
  summarise(
    # Probability of no outbreak in all months of this year
    p_no_outbreak_year = prod(1 - fitted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(



    yearly_probs <- pred_df %>%
  group_by(district_country, year) %>%
  summarise(
    # Probability of no outbreak in all months of this year
    p_no_outbreak_year = prod(1 - fitted, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Probability of at least one outbreak in the year
    p_year = 1 - p_no_outbreak_year
  )

yearly_probs
cumulative_probs <- yearly_probs %>%
  group_by(district_country) %>%
  summarise(
    p_no_outbreak_20yrs = prod(1 - p_year, na.rm = TRUE),
    p_20years = 1 - p_no_outbreak_20yrs,
    .groups = "drop"
  )

cumulative_probs
# Show result
df
# ==========================
    # Probability of at least one outbreak in the year
    p_year = 1 - p_no_outbreak_year
  )
