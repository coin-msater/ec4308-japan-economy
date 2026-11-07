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

# empty table for results
evaluation_results <- data.frame()

for (h in horizons) {
  results <- read_csv(
    paste0("lasso/results/periodic_10block_aligned/lasso_h", h, ".csv"),
    show_col_types = FALSE
  )
  
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

print(evaluation_results)
# save evaluation results
write_csv(
  evaluation_results,
  "lasso/results/periodic_10block_aligned/lasso_evaluation.csv"
)

par(mfrow = c(2, 2))
# PLOT ACTUAL VS PREDICTED
library(readr)

horizons <- c(1, 3, 6, 12)

# Set 2x2 plotting layout BEFORE the loop
par(mfrow = c(2, 2))

for (h in horizons) {
  
  results <- read_csv(
    paste0("lasso/results/periodic_10block/lasso_h", h, ".csv"),
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
    main = paste0("lasso h", h, ": Actual vs Predicted")
  )
  
  lines(
    results$target_date,
    results$predicted,
    col = "red",
    lwd = 1.5
  )
  
  legend(
    "topright",
    legend = c("Actual", "lasso forecast"),
    col = c("black", "red"),
    lty = 1,
    cex = 0.7
  )
}

# Reset layout AFTER all 4 plots
par(mfrow = c(1, 1))