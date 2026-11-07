library(glmnet)
library(readr)

# setup
horizons <- c(1, 3, 6, 12)

window_size <- 120
min_train_rows <- 50
validation_size <- 100
retune_interval <- 100
n_folds <- 10

validation_start <- as.Date("1984-05-01")
validation_end   <- as.Date("1992-08-01")

test_start <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")
covid_end <- as.Date("2021-12-01")

# manually specified lambda grid
lambda_grid <- 10^seq(5, -4, length.out = 200)

# elastic net mixing grid
# alpha = 0 is ridge and alpha = 1 is lasso, so both are excluded
alpha_grid <- seq(0.05, 0.95, by = 0.05)

# save results
output_dir <- "elastic_net/results/periodic_10block_aligned"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# helper function : add months
# all dates are the first day of a month
add_months <- function(date, n) {
  date <- as.Date(date)
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))
  total_months <- year * 12 + (month - 1) + n
  new_year <- total_months %/% 12
  new_month <- total_months %% 12 + 1
  as.Date(sprintf("%04d-%02d-01", new_year, new_month))
}

# helper function : select latest 120 eligible observations for rolling window
get_training_idx <- function(origin, df, window_size = 120) {
  eligible <- which(
    df$month < origin & df$target_date <= origin
  )
  
  tail(eligible, window_size)
}

# build the periodic tuning schedule
test_dates <- seq(
  from = test_start,
  to = covid_end,
  by = "month"
)

block_id <- ceiling(seq_along(test_dates) / retune_interval)
schedule <- do.call(
  rbind,
  lapply(split(test_dates, block_id), function(dates) {
    first_date <- min(dates)
    data.frame(
      tuning_date = first_date,
      validation_start = add_months(first_date, -validation_size),
      validation_end = add_months(first_date, -1),
      test_start = first_date,
      test_end = max(dates)
    )
  })
)
rownames(schedule) <- NULL

# verify the first validation period
stopifnot(
  schedule$validation_start[1] == validation_start,
  schedule$validation_end[1] == validation_end
)

# function : 10-BLOCK CHRONOLOGICAL CV OVER ALPHA AND LAMBDA
tune_enet <- function(df, X, y, validation_idx, tuning_date,
                      lambda_grid, alpha_grid, window_size = 120,
                      n_folds = 10) {
  
  # historical validation observations
  validation_idx <- validation_idx[
    order(df$target_date[validation_idx])
  ]
  
  stopifnot(
    length(validation_idx) >= n_folds,
    all(df$target_date[validation_idx] <= tuning_date)
  )
  
  # divide consecutive observations into chronological folds
  fold_id <- as.integer(
    cut(
      seq_along(validation_idx),
      breaks = n_folds,
      labels = FALSE
    )
  )
  
  folds <- split(validation_idx, fold_id)
  
  # squared errors: validation observation x lambda x alpha
  cv_errors <- array(
    NA_real_,
    dim = c(length(validation_idx), length(lambda_grid), length(alpha_grid)))
  
  # validation forecasts: validation observation x lambda x alpha
  cv_predictions <- array(
    NA_real_,
    dim = c(length(validation_idx), length(lambda_grid), length(alpha_grid)))
  
  for (k in seq_along(folds)) {
    fold_idx <- folds[[k]]
    
    # fit model at the first origin of each fold
    fold_start <- df$month[fold_idx[1]]
    train_idx <- get_training_idx(fold_start, df, window_size)
    
    # stop if insufficient historical data
    if (length(train_idx) < min_train_rows) {
      stop("Insufficient training observations at ", as.character(fold_start))
    }
    
    # score only targets observable at the first test forecast origin
    score_idx <- fold_idx
    
    if (length(score_idx) == 0) next
    
    positions <- match(score_idx, validation_idx)
    
    # fit elastic net once per alpha across the common lambda grid
    for (a in seq_along(alpha_grid)) {
      fit <- glmnet(
        X[train_idx, , drop = FALSE],
        y[train_idx],
        alpha = alpha_grid[a],
        lambda = lambda_grid,
        standardize = TRUE
      )
      
      # predict every scoreable origin in this fold
      pred <- predict(
        fit,
        newx = X[score_idx, , drop = FALSE],
        s = lambda_grid
      )
      
      # squared errors for each candidate lambda
      cv_errors[positions, , a] <- sweep(pred, 1, y[score_idx], FUN = "-")^2
      cv_predictions[positions, , a] <- pred
    }
  }
  
  # keep only fully scored validation rows
  valid_rows <- apply(is.finite(cv_errors), 1, all)
  
  n_scored <- sum(valid_rows)
  
  if (n_scored == 0) {
    stop("No valid CV forecasts at ", as.character(tuning_date))
  }
  
  # mean CV MSE as a lambda x alpha matrix
  mean_cv_mse <- apply(cv_errors[valid_rows, , , drop = FALSE], c(2, 3), mean)
  
  cv_rmse <- sqrt(mean_cv_mse)
  best <- arrayInd(which.min(cv_rmse), dim(cv_rmse))
  best_idx <- best[1, 1]
  best_alpha_idx <- best[1, 2]
  
  # save validation predictions using optimal alpha and lambda
  validation_forecasts <- data.frame(
    forecast_origin = df$month[validation_idx],
    target_date = df$target_date[validation_idx],
    actual = y[validation_idx],
    predicted = cv_predictions[, best_idx, best_alpha_idx])
  
  # remove validation observations that were not scored
  validation_forecasts <- validation_forecasts[
    valid_rows, ,
    drop = FALSE]
  
  list(
    alpha = alpha_grid[best_alpha_idx],
    lambda = lambda_grid[best_idx],
    cv_mse = mean_cv_mse[best_idx, best_alpha_idx],
    mean_cv_mse = mean_cv_mse,
    n_scored = n_scored,
    best_idx = best_idx,
    best_alpha_idx = best_alpha_idx,
    validation_forecasts = validation_forecasts,
    n_folds_used = sum(vapply(
      folds,
      function(idx) any(valid_rows[match(idx, validation_idx)]),
      logical(1)
    ))
  )
}

# function : MONTHLY TEST FORECASTING
forecast_block <- function(df, X, y, start_date, end_date, alpha, lambda,
                           lambda_grid, window_size = 120) {
  
  # select test observations using target dates
  forecast_idx <- which(
    df$target_date >= start_date &
      df$target_date <= end_date
  )
  
  results <- data.frame()
  
  for (idx in forecast_idx) {
    origin <- df$month[idx]
    train_idx <- get_training_idx(origin, df, window_size)
    
    if (length(train_idx) < min_train_rows) {
      stop("Insufficient training observations at ", as.character(origin))
    }
    
    # fit the complete path for numerical stability
    fit <- glmnet(
      X[train_idx, , drop = FALSE],
      y[train_idx],
      alpha = alpha,
      lambda = lambda_grid,
      standardize = TRUE
    )
    
    pred <- as.numeric(
      predict(fit, newx = X[idx, , drop = FALSE], s = lambda)
    )
    
    results <- rbind(
      results,
      data.frame(
        forecast_origin = origin,
        target_date = df$target_date[idx],
        actual = y[idx],
        predicted = pred,
        alpha = alpha,
        lambda = lambda
      )
    )
  }
  
  results
}

# run for all 4 horizons
for (h in horizons) {
  cat("\n====================================\n")
  cat("Running elastic net horizon:", h, "\n")
  cat("====================================\n")
  
  df <- read_csv(
    paste0("data/processed/h", h, ".csv"),
    show_col_types = FALSE
  )
  
  df$month <- as.Date(df$month)
  df$target_date <- as.Date(df$target_date)
  
  # sort observations chronologically
  df <- df[order(df$month), ]
  stopifnot(!anyDuplicated(df$month))
  
  X <- as.matrix(df[, !names(df) %in% c("month", "target_date", "target")])
  y <- df$target
  
  stopifnot(
    nrow(X) == length(y),
    nrow(X) == nrow(df),
    all(is.finite(X)),
    all(is.finite(y))
  )
  
  all_forecasts <- data.frame()
  all_tuning <- data.frame()
  all_validation <- data.frame()
  
  for (b in seq_len(nrow(schedule))) {
    s <- schedule[b, ]
    
    # first test forecast origin
    first_test_idx <- which(df$target_date == s$test_start)
    
    if (length(first_test_idx) != 1) {
      stop("First test target date not found")
    }
    
    first_test_origin <- df$month[first_test_idx]
    
    # select validation observations
    validation_idx <- which(
      df$target_date >= s$validation_start &
        df$target_date <= s$validation_end &
        df$target_date <= first_test_origin
    )
    
    # tune alpha and lambda
    tuning <- tune_enet(
      df = df,
      X = X,
      y = y,
      validation_idx = validation_idx,
      tuning_date = first_test_origin,
      lambda_grid = lambda_grid,
      alpha_grid = alpha_grid,
      window_size = window_size,
      n_folds = n_folds)
    
    # save validation forecasts for ensemble
    validation_h <- tuning$validation_forecasts
    validation_h$horizon <- h
    validation_h$retune_target_date <- s$test_start
    validation_h$information_cutoff <- first_test_origin
    
    all_validation <- rbind(
      all_validation,
      validation_h)
    
    forecasts <- forecast_block(
      df = df,
      X = X,
      y = y,
      start_date = s$test_start,
      end_date = s$test_end,
      alpha = tuning$alpha,
      lambda = tuning$lambda,
      lambda_grid = lambda_grid,
      window_size = window_size
    )
    
    all_forecasts <- rbind(
      all_forecasts,
      forecasts
    )
    
    all_tuning <- rbind(
      all_tuning,
      data.frame(
        horizon = h,
        retune_target_date = s$test_start,
        information_cutoff = first_test_origin,
        alpha = tuning$alpha,
        lambda = tuning$lambda,
        cv_mse = tuning$cv_mse,
        n_scored = tuning$n_scored,
        n_folds_used = tuning$n_folds_used,
        best_idx = tuning$best_idx,
        best_alpha_idx = tuning$best_alpha_idx
      )
    )
    
    cat(
      "Retune target:", as.character(s$test_start),
      "| Information cutoff:", as.character(first_test_origin),
      "| Alpha:", tuning$alpha,
      "| Lambda:", round(tuning$lambda, 4),
      "| CV forecasts:", tuning$n_scored,
      "| Folds:", tuning$n_folds_used,
      "\n"
    )
  }
  
  # verify target date alignment
  stopifnot(
    nrow(all_forecasts) == length(test_dates),
    min(all_forecasts$target_date) == test_start,
    max(all_forecasts$target_date) == covid_end
  )
  
  # Save this horizon's results
  write_csv(
    all_forecasts,
    file.path(output_dir, paste0("elastic_net_h", h, ".csv"))
  )
  
  write_csv(
    all_tuning,
    file.path(output_dir, paste0("elastic_net_h", h, "_tuning.csv"))
  )
  
  # save forecasts for ensemble models
  ensemble_results <- data.frame(
    horizon = h,
    forecast_origin = all_forecasts$forecast_origin,
    target_date = all_forecasts$target_date,
    actual = all_forecasts$actual,
    predicted = all_forecasts$predicted
  )
  
  write_csv(
    ensemble_results,
    file.path(output_dir, paste0("elastic_net_ensemble_h", h, ".csv"))
  )
  
  # save validation forecasts for ensemble
  write_csv(
    all_validation,
    file.path(output_dir, paste0("elastic_net_validation_h", h, ".csv"))
  )
  
  cat(
    "Finished horizon:", h,
    "| Test forecasts:", nrow(all_forecasts),
    "\n"
  )
}

cat("\nAll elastic net horizons completed.\n")