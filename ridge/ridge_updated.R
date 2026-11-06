library(glmnet)
library(readr)

# setup
horizons <- c(1, 3, 6, 12)

window_size <- 120
validation_size <- 100
retune_interval <- 100
n_folds <- 10

validation_start <- as.Date("1984-05-01")
validation_end   <- as.Date("1992-08-01")

test_start <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")
covid_end <- as.Date("2021-12-01")

# manually specified lambda grid 
lambda_grid <- 10^seq(5, -4,length.out = 200)
# save results
output_dir <- "ridge/results/periodic_10block"
dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)

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

# helper function : select latest 120 eligibile observations for rolling window 
get_training_idx <- function(origin, df, window_size = 120) {
  eligible <- which(
    df$month < origin & df$target_date <= origin)
  if (length(eligible) < window_size) {
    return(integer(0))
  }
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
stopifnot(schedule$validation_start[1] == validation_start,schedule$validation_end[1] == validation_end)

# function : 10-BLOCK CHRONOLOGICAL CV
tune_lambda <- function(df, X, y, validation_start, validation_end, tuning_date,lambda_grid,
    window_size = 120,n_folds = 10, alpha = 0) {
  # historical validation origins
  validation_idx <- which(df$month >= validation_start & df$month <= validation_end)
  if (length(validation_idx) == 0) {
    stop("No validation observations")
  }
  # divide consecutive origins into chronological folds.
  # with 100 observations, this gives 10 folds of 10.
  fold_id <- ceiling(
    seq_along(validation_idx) * n_folds /
      length(validation_idx))
  folds <- split(validation_idx, fold_id)
  cv_errors <- matrix(
    NA_real_,
    nrow = length(validation_idx),
    ncol = length(lambda_grid))
  
  for (k in seq_along(folds)) {
    fold_idx <- folds[[k]]
    # fit model at the first origin of each fold
    fold_start <- df$month[fold_idx[1]]
    train_idx <- get_training_idx(fold_start, df, window_size)
    # skip fold if insufficient historical data
    if (length(train_idx) == 0) next
    # score only targets observable at tuning date
    score_idx <- fold_idx[df$target_date[fold_idx] <= tuning_date]
    if (length(score_idx) == 0) next
    # fit Ridge once per fold across the common grid
    fit <- glmnet(X[train_idx, , drop = FALSE],y[train_idx],alpha = alpha,
      lambda = lambda_grid,standardize = TRUE)
    # predict every scoreable origin in this fold
    pred <- predict(fit,newx = X[score_idx, , drop = FALSE],s = lambda_grid)
    # squared errors for each candidate lambda
    errors <- sweep(pred,1,y[score_idx],FUN = "-")^2
    positions <- match(score_idx, validation_idx)
    cv_errors[positions, ] <- errors
  }
  
  # keep only fully scored validation rows
  valid_rows <- rowSums(is.finite(cv_errors)) ==
    length(lambda_grid)
  n_scored <- sum(valid_rows)
  if (n_scored == 0) {
    stop("No valid CV forecasts at ",as.character(tuning_date))
  }
  
  mean_cv_mse <- colMeans(cv_errors[valid_rows, , drop = FALSE])
  
  best_idx <- which.min(mean_cv_mse)
  
  list(
    lambda = lambda_grid[best_idx],
    cv_mse = mean_cv_mse[best_idx],
    mean_cv_mse = mean_cv_mse,
    n_scored = n_scored,
    best_idx = best_idx,
    n_folds_used = sum(vapply(
      folds,
      function(idx) any(valid_rows[match(idx, validation_idx)]),
      logical(1)
    ))
  )
}

# function : MONTHLY TEST FORECASTING
forecast_block <- function(df, X, y, start_date,end_date,lambda, lambda_grid,
    window_size = 120,alpha = 0) {
  forecast_idx <- which(df$month >= start_date & df$month <= end_date)
  results <- data.frame()
  
  for (idx in forecast_idx) {
    origin <- df$month[idx]
    train_idx <- get_training_idx(origin, df, window_size)
    
    if (length(train_idx) == 0) next
    # fit the complete path for numerical stability
    fit <- glmnet(
      X[train_idx, , drop = FALSE],
      y[train_idx],
      alpha = alpha,
      lambda = lambda_grid,
      standardize = TRUE)
    
    pred <- as.numeric(predict(fit,newx = X[idx, , drop = FALSE],s = lambda))
    
    results <- rbind(
      results,
      data.frame(
        forecast_origin = origin,
        target_date = df$target_date[idx],
        actual = y[idx],
        predicted = pred,
        lambda = lambda
      )
    )
  }
  results
}

# run for all 4 horzions 
for (h in horizons) {
  cat("\n====================================\n")
  cat("Running Ridge horizon:", h, "\n")
  cat("====================================\n")
  df <- read_csv(paste0("data/processed/h", h, ".csv"),show_col_types = FALSE)
  
  df$month <- as.Date(df$month)
  df$target_date <- as.Date(df$target_date)
  
  X <- as.matrix(df[, !names(df) %in% c("month", "target_date", "target")])
  y <- df$target
  
  stopifnot(nrow(X) == length(y),nrow(X) == nrow(df),all(is.finite(X)),all(is.finite(y)))
  
  all_forecasts <- data.frame()
  all_tuning <- data.frame()
  for (b in seq_len(nrow(schedule))) {
    s <- schedule[b, ]
    tuning <- tune_lambda(
      df = df,
      X = X,
      y = y,
      validation_start = s$validation_start,
      validation_end = s$validation_end,
      tuning_date = s$tuning_date,
      lambda_grid = lambda_grid,
      window_size = window_size,
      n_folds = n_folds,
      alpha = 0
    )
  
    forecasts <- forecast_block(
      df = df,
      X = X,
      y = y,
      start_date = s$test_start,
      end_date = s$test_end,
      lambda = tuning$lambda,
      lambda_grid = lambda_grid,
      window_size = window_size,
      alpha = 0
    )
    
    all_forecasts <- rbind(
      all_forecasts,
      forecasts
    )
    
    all_tuning <- rbind(
      all_tuning,
      data.frame(
        horizon = h,
        tuning_date = s$tuning_date,
        lambda = tuning$lambda,
        cv_mse = tuning$cv_mse,
        n_scored = tuning$n_scored,
        n_folds_used = tuning$n_folds_used,
        best_idx = tuning$best_idx
      )
    )
    
    cat(
      "Tuning date:", as.character(s$tuning_date),
      "| Lambda:", round(tuning$lambda, 4),
      "| CV forecasts:", tuning$n_scored,
      "| Folds:", tuning$n_folds_used,
      "\n"
    )
  }
  
  # Save this horizon's results
  write_csv(
    all_forecasts,
    file.path(output_dir,paste0("ridge_h", h, ".csv"))
  )
  
  write_csv(
    all_tuning,
    file.path(output_dir,paste0("ridge_h", h, "_tuning.csv")
    )
  )
  # save Ridge forecasts for ensemble models
  ensemble_results <- data.frame(
    horizon = h,
    forecast_origin = all_forecasts$forecast_origin,
    target_date = all_forecasts$target_date,
    actual = all_forecasts$actual,
    predicted = all_forecasts$predicted)
  
  write_csv(
    ensemble_results,
    file.path(output_dir, paste0("ridge_ensemble_h", h, ".csv")))
  
  cat(
    "Finished horizon:", h,
    "| Test forecasts:", nrow(all_forecasts),
    "\n"
  )
}
cat("\nAll Ridge horizons completed.\n")
