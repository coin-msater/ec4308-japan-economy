library(readr)
library(dplyr)

# evaluation functions
RMSE <- function(pred, actual) {
  sqrt(mean((actual - pred)^2))
}

MAE <- function(pred, actual) {
  mean(abs(actual - pred))
}

# evaluation periods
validation_end <- as.Date("1992-08-01")
test_start <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")
covid_end <- as.Date("2021-12-01")

horizons <- c(1, 3, 6, 12)

# empty table for results
evaluation_results <- data.frame()

for (h in horizons) {
  results <- read_csv(
    paste0("lasso/results/lasso_h", h, ".csv"),
    show_col_types = FALSE
  )
  
  # define evaluation periods
  results <- results %>%
    mutate(
      period = case_when(
        forecast_origin <= validation_end ~ "Validation",
        forecast_origin >= test_start &
          forecast_origin <= pre_covid_end ~ "Pre-COVID",
        forecast_origin >= as.Date("2020-01-01") &
          forecast_origin <= covid_end ~ "COVID period"
      )
    )
  
  # calculate performance for each period
  metrics <- results %>%
    group_by(period) %>%
    summarise(
      N = n(),
      RMSE = RMSE(predicted, actual),
      MAE = MAE(predicted, actual),
      .groups = "drop"
    ) %>%
    mutate(horizon = h)
  
  evaluation_results <- bind_rows(
    evaluation_results,
    metrics
  )
}

print(evaluation_results)