############################################
# Elastic Net
############################################

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

lambda_grid <- 10^seq(5, -4, length.out = 200)
alpha_grid <- seq(0.05, 0.95, by = 0.05)

# save results
output_dir <- "elastic_net/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################################
# Helpers
############################################

# add months; all dates are the first day of a month
add_months <- function(date, n) {
  date <- as.Date(date)
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))
  total_months <- year * 12 + (month - 1) + n
  new_year <- total_months %/% 12
  new_month <- total_months %% 12 + 1
  as.Date(sprintf("%04d-%02d-01", new_year, new_month))
}

# latest 120 eligible observations for the rolling window:
# origin before the forecast origin AND target already observed
get_training_idx <- function(origin, df, window_size = 120) {
  eligible <- which(df$month < origin & df$target_date <= origin)
  if (length(eligible) < window_size) {
    return(integer(0))
  }
  tail(eligible, window_size)
}


############################################
# Periodic tuning schedule
############################################

test_dates <- seq(from = test_start, to = covid_end, by = "month")

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


############################################
# 10-block chronological CV over alpha AND lambda, scored by RMSE
#
# Inside each fold the model
# is fitted once per alpha across the full lambda grid.
############################################

tune_enet <- function(df, X, y, validation_start, validation_end, tuning_date,
                      lambda_grid, alpha_grid,
                      window_size = 120, n_folds = 10) {
  
  # historical validation origins
  validation_idx <- which(df$month >= validation_start & df$month <= validation_end)
  if (length(validation_idx) == 0) {
    stop("No validation observations")
  }
  
  # consecutive origins into chronological folds (100 origins -> 10 folds of 10)
  fold_id <- ceiling(seq_along(validation_idx) * n_folds / length(validation_idx))
  folds <- split(validation_idx, fold_id)
  
  # squared errors: validation origin x lambda x alpha
  cv_errors <- array(
    NA_real_,
    dim = c(length(validation_idx), length(lambda_grid), length(alpha_grid))
  )
  
  for (k in seq_along(folds)) {
    fold_idx <- folds[[k]]
    
    # fit at the first origin of each fold
    fold_start <- df$month[fold_idx[1]]
    train_idx <- get_training_idx(fold_start, df, window_size)
    
    # skip fold if insufficient historical data
    if (length(train_idx) == 0) next
    
    # score only targets observable at the tuning date
    score_idx <- fold_idx[df$target_date[fold_idx] <= tuning_date]
    if (length(score_idx) == 0) next
    
    positions <- match(score_idx, validation_idx)
    
    for (a in seq_along(alpha_grid)) {
      fit <- glmnet(
        X[train_idx, , drop = FALSE],
        y[train_idx],
        alpha = alpha_grid[a],
        lambda = lambda_grid,
        standardize = TRUE
      )
      
      pred <- predict(fit, newx = X[score_idx, , drop = FALSE], s = lambda_grid)
      
      cv_errors[positions, , a] <- sweep(pred, 1, y[score_idx], FUN = "-")^2
    }
  }
  
  # keep only fully scored validation rows
  valid_rows <- apply(is.finite(cv_errors), 1, all)
  n_scored <- sum(valid_rows)
  if (n_scored == 0) {
    stop("No valid CV forecasts at ", as.character(tuning_date))
  }
  
  # CV RMSE as a lambda x alpha matrix
  cv_rmse <- sqrt(apply(cv_errors[valid_rows, , , drop = FALSE], c(2, 3), mean))
  
  best <- arrayInd(which.min(cv_rmse), dim(cv_rmse))
  best_lambda_idx <- best[1, 1]
  best_alpha_idx  <- best[1, 2]
  
  list(
    alpha = alpha_grid[best_alpha_idx],
    lambda = lambda_grid[best_lambda_idx],
    cv_rmse = cv_rmse[best_lambda_idx, best_alpha_idx],
    cv_rmse_grid = cv_rmse,
    best_curve = cv_rmse[, best_alpha_idx],
    n_scored = n_scored,
    best_idx = best_lambda_idx,
    best_alpha_idx = best_alpha_idx,
    lambda_at_edge = best_lambda_idx %in% c(1, length(lambda_grid)),
    n_folds_used = sum(vapply(
      folds,
      function(idx) any(valid_rows[match(idx, validation_idx)]),
      logical(1)
    ))
  )
}


############################################
# Monthly test forecasting
############################################

forecast_block <- function(df, X, y, start_date, end_date, alpha, lambda,
                           lambda_grid, window_size = 120) {
  
  forecast_idx <- which(df$month >= start_date & df$month <= end_date)
  rows <- vector("list", length(forecast_idx))
  
  for (j in seq_along(forecast_idx)) {
    idx <- forecast_idx[j]
    origin <- df$month[idx]
    train_idx <- get_training_idx(origin, df, window_size)
    
    if (length(train_idx) == 0) next
    
    # fit the complete path for numerical stability
    fit <- glmnet(
      X[train_idx, , drop = FALSE],
      y[train_idx],
      alpha = alpha,
      lambda = lambda_grid,
      standardize = TRUE
    )
    
    pred <- as.numeric(predict(fit, newx = X[idx, , drop = FALSE], s = lambda))
    beta <- as.matrix(coef(fit, s = lambda))
    
    rows[[j]] <- data.frame(
      forecast_origin = origin,
      target_date = df$target_date[idx],
      forecast = pred,
      actual = y[idx],
      alpha = alpha,
      lambda = lambda,
      n_nonzero = sum(beta[-1, 1] != 0)
    )
  }
  
  do.call(rbind, rows)
}


############################################
# Run all 4 horizons
############################################

all_results <- list()   # full forecast detail per horizon
cv_curves <- list()     # CV RMSE curve (at the chosen alpha) per horizon x tuning date

for (h in horizons) {
  cat("\n====================================\n")
  cat("Running Elastic Net horizon:", h, "\n")
  cat("====================================\n")
  
  df <- read_csv(paste0("data/processed/h", h, ".csv"), show_col_types = FALSE)
  
  df$month <- as.Date(df$month)
  df$target_date <- as.Date(df$target_date)
  
  X <- as.matrix(df[, !names(df) %in% c("month", "target_date", "target")])
  y <- df$target
  
  stopifnot(
    nrow(X) == length(y),
    nrow(X) == nrow(df),
    all(is.finite(X)),
    all(is.finite(y))
  )
  
  forecast_list <- list()
  tuning_list <- list()
  grid_list <- list()
  
  for (b in seq_len(nrow(schedule))) {
    s <- schedule[b, ]
    
    tuning <- tune_enet(
      df = df,
      X = X,
      y = y,
      validation_start = s$validation_start,
      validation_end = s$validation_end,
      tuning_date = s$tuning_date,
      lambda_grid = lambda_grid,
      alpha_grid = alpha_grid,
      window_size = window_size,
      n_folds = n_folds
    )
    
    forecast_list[[b]] <- forecast_block(
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
    
    tuning_list[[b]] <- data.frame(
      horizon = h,
      tuning_date = s$tuning_date,
      alpha = tuning$alpha,
      lambda = tuning$lambda,
      cv_rmse = tuning$cv_rmse,
      n_scored = tuning$n_scored,
      n_folds_used = tuning$n_folds_used,
      best_idx = tuning$best_idx,
      lambda_at_edge = tuning$lambda_at_edge
    )
    
    # full alpha x lambda CV RMSE surface (expand.grid order matches the matrix)
    grid_list[[b]] <- data.frame(
      horizon = h,
      tuning_date = s$tuning_date,
      expand.grid(lambda = lambda_grid, alpha = alpha_grid),
      cv_rmse = as.vector(tuning$cv_rmse_grid)
    )
    
    cv_curves[[paste0("h", h, " ", format(s$tuning_date, "%Y-%m"))]] <- tuning$best_curve
    
    cat(
      "Tuning date:", as.character(s$tuning_date),
      "| Alpha:", tuning$alpha,
      "| Lambda:", signif(tuning$lambda, 4),
      "| CV RMSE:", signif(tuning$cv_rmse, 4),
      "| CV forecasts:", tuning$n_scored,
      "| Folds:", tuning$n_folds_used,
      if (tuning$lambda_at_edge) "| WARNING: lambda at grid edge" else "",
      "\n"
    )
  }
  
  details <- do.call(rbind, forecast_list)
  all_results[[paste0("h", h)]] <- details
  
  # standardized team output: horizon, target_date, forecast, actual
  standard <- data.frame(
    horizon = h,
    target_date = details$target_date,
    forecast = details$forecast,
    actual = details$actual
  )
  write_csv(standard, file.path(output_dir, paste0("elastic_net_h", h, ".csv")))
  
  # extra Elastic Net diagnostics, kept separate from the standard file
  write_csv(details, file.path(output_dir, paste0("elastic_net_h", h, "_details.csv")))
  write_csv(do.call(rbind, tuning_list), file.path(output_dir, paste0("elastic_net_h", h, "_tuning.csv")))
  write_csv(do.call(rbind, grid_list),   file.path(output_dir, paste0("elastic_net_h", h, "_cv_grid.csv")))
  
  cat("Finished horizon:", h, "| Test forecasts:", nrow(details), "\n")
}

# all horizons stacked in the standard format
standard_all <- do.call(rbind, lapply(horizons, function(h) {
  d <- all_results[[paste0("h", h)]]
  data.frame(horizon = h, target_date = d$target_date, forecast = d$forecast, actual = d$actual)
}))
write_csv(standard_all, file.path(output_dir, "elastic_net_all_horizons.csv"))

cat("\nAll Elastic Net horizons completed.\n")


############################################
# CV RMSE curves at the chosen alpha (16 panels: 4 horizons x 4 tuning dates)
############################################

plot_cv_curves <- function() {
  old_par <- par(mfrow = c(4, 4), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0),
                 cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
  on.exit(par(old_par))
  for (nm in names(cv_curves)) {
    curve <- cv_curves[[nm]]
    plot(log10(lambda_grid), curve, type = "l",
         xlab = "log10(lambda)", ylab = "CV RMSE", main = nm)
    abline(v = log10(lambda_grid[which.min(curve)]), lty = 2)
  }
}

png(file.path(output_dir, "cv_curves.png"), width = 1600, height = 1400, res = 150)
plot_cv_curves()
dev.off()

try(plot_cv_curves())

length(cv_curves)   # should be 16

# best CV RMSE relative to the heaviest penalty (close to intercept-only).
# values near 1 mean the predictors add little over the mean.
sapply(cv_curves, function(curve) min(curve) / curve[1])


############################################
# Evaluation: RMSE and MAE
############################################

# accepts the standard format (forecast) or the older one (predicted)
read_forecasts <- function(path) {
  f <- as.data.frame(read_csv(path, show_col_types = FALSE))
  if (!"forecast" %in% names(f) && "predicted" %in% names(f)) {
    names(f)[names(f) == "predicted"] <- "forecast"
  }
  f$target_date <- as.Date(f$target_date)
  f
}

evaluate <- function(f) {
  pre <- f$target_date <= pre_covid_end
  data.frame(
    n = nrow(f),
    rmse = sqrt(mean((f$actual - f$forecast)^2)),
    mae = mean(abs(f$actual - f$forecast)),
    n_pre = sum(pre),
    rmse_pre = sqrt(mean((f$actual[pre] - f$forecast[pre])^2)),
    mae_pre = mean(abs(f$actual[pre] - f$forecast[pre]))
  )
}

enet_errors <- do.call(rbind, lapply(horizons, function(h) {
  cbind(horizon = h, evaluate(all_results[[paste0("h", h)]]))
}))
print(enet_errors, row.names = FALSE)
write_csv(enet_errors, file.path(output_dir, "elastic_net_errors.csv"))

# intercept-only forecasts per horizon (every predictor shrunk to zero)
sapply(all_results, function(f) sum(f$n_nonzero == 0))