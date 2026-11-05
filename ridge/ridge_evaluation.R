library(readr)
library(dplyr)

RMSE <- function(pred, actual) {
  sqrt(mean((actual - pred)^2))
}

MAE <- function(pred, actual) {
  mean(abs(actual - pred))
}

horizons <- c(1, 3, 6, 12)

comparison <- data.frame()

for (h in horizons) {
  
  # Original monthly-tuning Ridge
  old <- read_csv(
    paste0("ridge/results/ridge_h", h, ".csv"),
    show_col_types = FALSE
  )
  
  # New periodic-tuning Ridge
  new <- read_csv(
    paste0(
      "ridge/results/periodic_10block/ridge_h",
      h, ".csv"
    ),
    show_col_types = FALSE
  )
  
  # Match forecasts by origin AND target date
  matched <- inner_join(
    old %>%
      select(forecast_origin, target_date,
             actual_old = actual,
             predicted_old = predicted),
    new %>%
      select(forecast_origin, target_date,
             actual_new = actual,
             predicted_new = predicted),
    by = c("forecast_origin", "target_date")
  )
  
  # Verify identical realized targets
  stopifnot(
    all(abs(matched$actual_old -
              matched$actual_new) < 1e-8)
  )
  
  # Define evaluation periods
  matched <- matched %>%
    mutate(
      period = case_when(
        forecast_origin <= as.Date("2019-12-01") ~ "Pre-COVID",
        TRUE ~ "COVID period"
      )
    )
  
  # Compare forecasting errors
  metrics <- matched %>%
    group_by(period) %>%
    summarise(
      N = n(),
      RMSE_old = RMSE(predicted_old, actual_old),
      RMSE_new = RMSE(predicted_new, actual_new),
      MAE_old = MAE(predicted_old, actual_old),
      MAE_new = MAE(predicted_new, actual_new),
      .groups = "drop"
    ) %>%
    mutate(
      horizon = h,
      RMSE_change_pct =
        100 * (RMSE_new / RMSE_old - 1)
    )
  
  comparison <- bind_rows(comparison, metrics)
}

print(comparison)