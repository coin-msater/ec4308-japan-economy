library(readr)

# setup 
horizons <- c(1, 3, 6, 12)


window_size <- 120
min_train_rows <- 50
validation_size <- 100
retune_interval <- 100

validation_start <- as.Date("1984-05-01")
validation_end   <- as.Date("1992-08-01")

test_start <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")
covid_end <- as.Date("2021-12-01")

# AR(p) benchmark: BIC chooses p among AR(1) ... AR(max_p)
max_p <- 12
p_grid <- seq_len(max_p)

# save results
output_dir <- "AR/ar_results_block_cv"
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

# build the periodic retuning schedule 
# test_dates are TARGET dates
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

# read + validate one horizon file
read_horizon <- function(h) {
  f <- file.path(DATA_DIR, paste0("h", h, ".csv"))
  if (!file.exists(f)) stop("Missing input CSV: ", f)
  d <- as.data.frame(read_csv(f, show_col_types = FALSE))
  d$month <- as.Date(d$month)
  d$target_date <- as.Date(d$target_date)
  
  # sort observations chronologically
  d <- d[order(d$month), ]
  rownames(d) <- NULL
  n <- nrow(d)
  
  stopifnot(!anyDuplicated(d$month),
            all(is.finite(d$x1_lag0)), all(is.finite(d$target)))
  # months are consecutive
  stopifnot(identical(as.character(d$month),
                      as.character(seq(d$month[1], by = "month", length.out = n))))
  # target_date = month + h
  expected_target_dates <- seq(d$month[1], by = "month",
                               length.out = n + h)[seq_len(n) + h]
  stopifnot(identical(as.character(d$target_date),
                      as.character(expected_target_dates)))
  d
}

# build the full monthly x1 series
# position i of the series <-> d$month[i]; the last h positions come from target
make_series <- function(d, h) {
  n <- nrow(d)
  # sanity checks on the assumption target[i] = x1 at month[i] + h
  stopifnot(isTRUE(all.equal(d$target[seq_len(n - h)], d$x1_lag0[seq_len(n - h) + h])))
  if (all(c("x1_lag1", "x1_lag2", "x1_lag3") %in% names(d))) {
    for (l in 1:3) {
      col <- d[[paste0("x1_lag", l)]]
      stopifnot(isTRUE(all.equal(col[(l + 1):n], d$x1_lag0[1:(n - l)])))
    }
  }
  c(d$x1_lag0, tail(d$target, h))
}

# row j = (x[j], x[j-1], ..., x[j-max_p+1]) -> regressors known at origin month[j]
lag_matrix <- function(x, n, max_p) {
  out <- matrix(NA_real_, n, max_p, dimnames = list(NULL, paste0("L", 0:(max_p - 1))))
  for (l in 0:(max_p - 1)) {
    idx <- seq_len(n) - l
    ok <- idx >= 1
    out[ok, l + 1] <- x[idx[ok]]
  }
  out
}

# training rows for an origin & check that all max_p lags exist 
get_train <- function(origin, df, L) {
  idx <- get_training_idx(origin, df, window_size)
  if (length(idx) < min_train_rows) {
    stop("Insufficient training observations at ", as.character(origin))
  }
  if (anyNA(L[idx, , drop = FALSE])) {
    stop("Missing lags in training window at ", as.character(origin))
  }
  idx
}

# OLS fit of the direct h-step AR(p) and prediction from its coefficients
fit_ar <- function(L, y, train_idx, p) {
  X <- L[train_idx, 1:p, drop = FALSE]
  yy <- y[train_idx]
  coef(lm(yy ~ X))
}

predict_ar <- function(coefs, Lrows, p) {
  as.numeric(cbind(1, Lrows[, 1:p, drop = FALSE]) %*% coefs)
}

# function : BIC LAG-ORDER SELECTION
# AR(1) ... AR(max_p) are all fit on the same latest-120 eligible window
# as of the block's first forecast origin; the lowest BIC wins
select_p_bic <- function(origin, df, L, y) {
  train_idx <- get_train(origin, df, L)
  yy <- y[train_idx]
  X <- L[train_idx, , drop = FALSE]
  
  bics <- sapply(p_grid, function(p) BIC(lm(yy ~ X[, 1:p, drop = FALSE])))
  
  list(p = p_grid[which.min(bics)], bic = min(bics))
}

# function : MONTHLY TEST FORECASTING
forecast_block <- function(df, L, y, start_date, end_date, p) {
  
  # select test observations using target dates
  forecast_idx <- which(
    df$target_date >= start_date &
      df$target_date <= end_date
  )
  
  results <- data.frame()
  
  for (idx in forecast_idx) {
    origin <- df$month[idx]
    train_idx <- get_train(origin, df, L)
    
    coefs <- fit_ar(L, y, train_idx, p)
    pred <- predict_ar(coefs, L[idx, , drop = FALSE], p)
    
    results <- rbind(
      results,
      data.frame(
        forecast_origin = origin,
        target_date = df$target_date[idx],
        actual = y[idx],
        predicted = pred,
        p = p
      )
    )
  }
  
  results
}

# run one horizon
run_horizon <- function(h) {
  cat("\n====================================\n")
  cat("Running AR horizon:", h, "\n")
  cat("====================================\n")
  
  df <- read_horizon(h)
  n <- nrow(df)
  x <- make_series(df, h)
  L <- lag_matrix(x, n, max_p)
  y <- df$target
  
  all_forecasts <- data.frame()
  all_tuning <- data.frame()
  
  for (b in seq_len(nrow(schedule))) {
    s <- schedule[b, ]
    
    # first test forecast origin
    first_test_idx <- which(df$target_date == s$test_start)
    
    if (length(first_test_idx) != 1) {
      stop("First test target date not found")
    }
    
    first_test_origin <- df$month[first_test_idx]
    
    # choose p by BIC at the block's first forecast origin
    sel <- select_p_bic(first_test_origin, df, L, y)
    
    # test forecasts for the block: monthly refit, p fixed within the block
    forecasts <- forecast_block(
      df = df,
      L = L,
      y = y,
      start_date = s$test_start,
      end_date = s$test_end,
      p = sel$p
    )
    
    all_forecasts <- rbind(all_forecasts, forecasts)
    
    all_tuning <- rbind(
      all_tuning,
      data.frame(
        horizon = h,
        retune_target_date = s$test_start,
        information_cutoff = first_test_origin,
        p = sel$p,
        bic = sel$bic
      )
    )
    
    cat(
      "Retune target:", as.character(s$test_start),
      "| Information cutoff:", as.character(first_test_origin),
      "| p (BIC):", sel$p,
      "\n"
    )
  }
  
  # verify target date alignment
  stopifnot(
    nrow(all_forecasts) == length(test_dates),
    identical(as.character(all_forecasts$target_date), as.character(test_dates))
  )
  
  # save this horizon's results
  ar_results <- data.frame(
    horizon = h,
    forecast_origin = all_forecasts$forecast_origin,
    target_date = all_forecasts$target_date,
    actual = all_forecasts$actual,
    predicted = all_forecasts$predicted
  )
  write_csv(ar_results, file.path(output_dir, paste0("ar_h", h, ".csv")))
  write_csv(all_tuning, file.path(output_dir, paste0("ar_h", h, "_tuning.csv")))
  
  cat("Finished horizon:", h, "| Test forecasts:", nrow(all_forecasts), "\n")
  
  ar_results
}

# run all horizons
res <- lapply(horizons, run_horizon)

# summary (pre-COVID and through Dec 2021, by target date)
summary_tbl <- do.call(rbind, lapply(res, function(f) {
  ends <- list(pre_covid = pre_covid_end, through_dec_2021 = covid_end)
  do.call(rbind, lapply(names(ends), function(nm) {
    z <- f[f$target_date <= ends[[nm]], ]
    e <- z$actual - z$predicted
    data.frame(horizon = f$horizon[1], period = nm, test_months = nrow(z),
               RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)))
  }))
}))
rownames(summary_tbl) <- NULL
write_csv(summary_tbl, file.path(output_dir, "ar_summary.csv"))

cat("\nAll AR horizons completed. Files saved under: ", output_dir, "\n", sep = "")