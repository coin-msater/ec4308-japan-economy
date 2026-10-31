library(glmnet)
library(readr)

# forecast horizons
horizons <- c(1, 3, 6, 12)

# sample split
window_size <- 120

validation_start <- as.Date("1984-05-01")
validation_end   <- as.Date("1992-08-01")

test_start <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")
covid_end     <- as.Date("2021-12-01")

# create folder to save forecasts
dir.create("lasso/results",recursive = TRUE,showWarnings = FALSE)

# FORECAST HORIZONS LOOP
for (h in horizons) {
  # load horizon-specific dataset
  df <- read_csv(paste0("data/processed/h", h, ".csv"),show_col_types = FALSE)
  
  # prepare predictors and target
  X <- as.matrix(df[, !names(df) %in% c("month", "target_date", "target")])
  y <- df$target
  
  # forecast origins
  forecast_dates <- df$month[
    df$month >= validation_start &
      df$month <= covid_end
  ]
  
  # store results separately for each horizon
  lasso_results <- data.frame(
    forecast_origin = as.Date(character()),
    target_date = as.Date(character()),
    actual = numeric(),
    predicted = numeric(),
    lambda = numeric()
  )
  
  # ROLLING WINDOW LOOP
  for (forecast_date in forecast_dates) {
    forecast_date <- as.Date(forecast_date,origin = "1970-01-01")
    # only use observations with available targets
    eligible <- which(
      df$month < forecast_date &
        df$target_date <= forecast_date)
    
    # require 120 eligible training observations
    if (length(eligible) < window_size) {
      next
    }
    
    # latest 120 eligible observations
    train_idx <- tail(eligible, window_size)
    
    # update training data
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    
    # generate candidate lambda values by fitting initial lasso model
    lasso.mod <- glmnet(
      X_train,
      y_train,
      alpha = 1,
      standardize = TRUE)
    lambda_grid <- lasso.mod$lambda
    
    # RECURSIVE CROSS-VALIDATION LOOP (inner loop)
    n <- nrow(X_train)
    initial_size <- 60
    
    cv_mse <- matrix(
      NA_real_,
      nrow = n - initial_size,
      ncol = length(lambda_grid))
    
    for (i in (initial_size + 1):n) {
      # training outcomes must be available by the validation forecast origin
      inner_train <- which(
        seq_len(n) < i &
          df$target_date[train_idx] <=
          df$month[train_idx[i]])
      inner_val <- i
      if (length(inner_train) < 30) {
        next
      }
      # fit lasso on inner training observations
      lasso_cv <- glmnet(
        X_train[inner_train, , drop = FALSE],
        y_train[inner_train],
        alpha = 1,
        lambda = lambda_grid,
        standardize = TRUE
      )
      
      # predict validation observation
      cv_pred <- predict(
        lasso_cv,
        newx = X_train[inner_val, , drop = FALSE],
        s = lambda_grid
      )
      
      # store squared validation errors
      cv_mse[i - initial_size, ] <-
        (as.numeric(cv_pred) - y_train[inner_val])^2
    }
    
    # selecting optimal lambda
    mean_cv_mse <- colMeans(cv_mse, na.rm = TRUE)
    if (!any(is.finite(mean_cv_mse))) {
      stop("No valid CV results at ", as.character(forecast_date))
    }
    
    best_idx <- which.min(mean_cv_mse)
    best_lambda <- lambda_grid[best_idx]
    
    # refitting lasso and generate forecast
    lasso_final <- glmnet(
      X_train,
      y_train,
      alpha = 1,
      lambda = best_lambda,
      standardize = TRUE
    )
    
    forecast_idx <- which(df$month == forecast_date)
    
    if (length(forecast_idx) != 1) {
      stop("Forecast date not found")
    }
    
    lasso_forecast <- predict(
      lasso_final,
      newx = X[forecast_idx, , drop = FALSE]
    )
    
    # save forecast
    lasso_results <- rbind(
      lasso_results,
      data.frame(
        forecast_origin = forecast_date,
        target_date = df$target_date[forecast_idx],
        actual = y[forecast_idx],
        predicted = as.numeric(lasso_forecast),
        lambda = best_lambda
      )
    )
  }
  
  # saving results for each horzion
  write_csv(
    lasso_results,
    paste0("lasso/results/lasso_h", h, ".csv")
  )
  # status check
  cat(
    "Finished horizon:", h,
    "| Forecasts:", nrow(lasso_results),
    "\n"
  )
}
