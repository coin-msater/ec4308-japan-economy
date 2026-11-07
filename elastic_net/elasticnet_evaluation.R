# EVALUATION
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

# results folder
output_dir <- "elastic_net/results/periodic_10block_aligned"

# empty table for results
evaluation_results <- data.frame()

for (h in horizons) {
  
  # test forecasts
  test_results <- read_csv(
    paste0(output_dir, "/elastic_net_h", h, ".csv"),
    show_col_types = FALSE
  )
  
  # validation forecasts (saved separately from the test forecasts)
  validation_results <- read_csv(
    paste0(output_dir, "/elastic_net_validation_h", h, ".csv"),
    show_col_types = FALSE
  )
  
  # combine validation and test forecasts
  results <- bind_rows(
    validation_results %>% select(forecast_origin, target_date, actual, predicted),
    test_results %>% select(forecast_origin, target_date, actual, predicted)
  ) %>%
    mutate(target_date = as.Date(target_date))
  
  # define evaluation periods
  results <- results %>%
    mutate(
      period = case_when(
        target_date <= validation_end ~ "Validation",
        target_date >= test_start & target_date <= pre_covid_end ~ "Pre-COVID",
        target_date >= as.Date("2020-01-01") &
          target_date <= covid_end ~ "COVID period"
      )
    )
  
  # validation forecasts from later retunes overlap the test period,
  # so only the test forecasts are used for the test periods
  results <- results %>%
    filter(
      period == "Validation" |
        paste(forecast_origin, target_date) %in%
        paste(test_results$forecast_origin, test_results$target_date)
    )
  
  full_test <- results %>%
    filter(
      target_date >= test_start &
        target_date <= covid_end
    ) %>%
    summarise(
      period = "Full test",
      N = n(),
      RMSE = RMSE(predicted, actual),
      MAE = MAE(predicted, actual)
    ) %>%
    mutate(horizon = h)
  
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
    full_test,
    metrics
  )
}

# order periods within each horizon
evaluation_results <- evaluation_results %>%
  mutate(
    period = factor(
      period,
      levels = c("Validation", "Full test", "Pre-COVID", "COVID period")
    )
  ) %>%
  arrange(horizon, period) %>%
  select(horizon, period, N, RMSE, MAE)

print(evaluation_results)

# save evaluation results
write_csv(
  evaluation_results,
  paste0(output_dir, "/elastic_net_evaluation.csv")
)


# PLOT ACTUAL VS PREDICTED

# Set 2x2 plotting layout BEFORE the loop
# (smaller margins so 4 plots fit in the Plots pane)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

for (h in horizons) {
  
  results <- read_csv(
    paste0(output_dir, "/elastic_net_h", h, ".csv"),
    show_col_types = FALSE
  )
  
  results$target_date <- as.Date(results$target_date)
  results <- results[order(results$target_date), ]
  
  # Actual vs Predicted
  plot(
    results$target_date,
    results$actual,
    type = "l",
    col = "black",
    lwd = 1.5,
    xlab = "Target date",
    ylab = "Industrial production",
    main = paste0("Elastic Net h", h, ": Actual vs Predicted")
  )
  
  lines(
    results$target_date,
    results$predicted,
    col = "red",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("Actual", "Elastic Net forecast"),
    col = c("black", "red"),
    lty = 1,
    cex = 0.7
  )
}

# Reset layout AFTER all 4 plots
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))