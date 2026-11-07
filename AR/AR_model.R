library(readr)

# setup
DATA_DIR   <- "data/processed"       # folder containing h1.csv, h3.csv, h6.csv, h12.csv
OUTPUT_DIR <- "ar_results_block_cv"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

horizons      <- c(1L, 3L, 6L, 12L)
TRAIN_N       <- 120L   # latest 120 eligible observations 
MAX_P         <- 12L    # BIC may choose AR(1) ... AR(MAX_P)
RETUNE_MONTHS <- 100L   # re-select p every 100 test months 

TEST_START     <- as.Date("1992-09-01")   # first TARGET date of the test period
PRE_COVID_END  <- as.Date("2019-12-01")
TEST_END       <- as.Date("2021-12-01")   # last TARGET date of the test period
RETUNE_TARGETS <- as.Date(c("1992-09-01", "2001-01-01",
                            "2009-05-01", "2017-09-01"))

# read + validate one horizon file
read_horizon <- function(h) {
  f <- file.path(DATA_DIR, paste0("h", h, ".csv"))
  if (!file.exists(f)) stop("Missing input CSV: ", f)
  d <- as.data.frame(read_csv(f, show_col_types = FALSE))
  d$month <- as.Date(d$month)
  d$target_date <- as.Date(d$target_date)
  n <- nrow(d)
  
  stopifnot(!is.unsorted(d$month), !anyDuplicated(d$month),
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

# row j = (x[j], x[j-1], ..., x[j-MAX_P+1]) -> regressors known at origin month[j]
lag_matrix <- function(x, n, max_p) {
  out <- matrix(NA_real_, n, max_p, dimnames = list(NULL, paste0("L", 0:(max_p - 1))))
  for (l in 0:(max_p - 1)) {
    idx <- seq_len(n) - l
    ok <- idx >= 1
    out[ok, l + 1] <- x[idx[ok]]
  }
  out
}

# same eligibility rule as Lasso / RF: outcome must be known at the origin
training_rows <- function(d, origin) {
  eligible <- which(d$month < origin & d$target_date <= origin)
  tail(eligible, TRAIN_N)
}

get_train <- function(d, L, origin) {
  idx <- training_rows(d, origin)
  if (length(idx) < TRAIN_N)
    stop("Only ", length(idx), " eligible training rows at ", as.character(origin))
  if (anyNA(L[idx, , drop = FALSE]))
    stop("Missing lags in training window at ", as.character(origin))
  idx
}

# BIC lag-order selection on the 120-obs window at the retune origin
select_p <- function(d, L, y, origin) {
  idx <- get_train(d, L, origin)
  yy <- y[idx]
  X <- L[idx, , drop = FALSE]
  bics <- sapply(seq_len(MAX_P), function(p) BIC(lm(yy ~ X[, 1:p, drop = FALSE])))
  which.min(bics)
}

# direct h-step AR(p) forecast made at the origin of row `row`
forecast_ar <- function(d, L, y, row, p) {
  idx <- get_train(d, L, d$month[row])
  yy <- y[idx]
  X <- L[idx, 1:p, drop = FALSE]
  m <- lm(yy ~ X)
  sum(c(1, L[row, 1:p]) * coef(m))
}

run_horizon <- function(h) {
  d <- read_horizon(h)
  n <- nrow(d)
  x <- make_series(d, h)
  L <- lag_matrix(x, n, MAX_P)
  y <- d$target
  
  # test set defined by target date 
  test_rows <- which(d$target_date >= TEST_START & d$target_date <= TEST_END)
  stopifnot(identical(as.character(d$target_date[test_rows]),
                      as.character(seq(TEST_START, TEST_END, by = "month"))))
  block_id <- (seq_along(test_rows) - 1L) %/% RETUNE_MONTHS + 1L
  first_positions <- which(!duplicated(block_id))
  stopifnot(max(block_id) == length(RETUNE_TARGETS),
            identical(as.character(d$target_date[test_rows[first_positions]]),
                      as.character(RETUNE_TARGETS)))
  
  preds <- numeric(length(test_rows))
  sel_p <- integer(length(test_rows))
  
  for (b in seq_along(first_positions)) {
    retune_origin <- d$month[test_rows[first_positions[b]]]
    p <- select_p(d, L, y, retune_origin)
    cat("H", h, " block ", b, ": retune origin=", as.character(retune_origin),
        ", selected p=", p, "\n", sep = "")
    
    # monthly refit on the latest 120 eligible rows, p fixed within the block
    for (k in which(block_id == b)) {
      preds[k] <- forecast_ar(d, L, y, test_rows[k], p)
      sel_p[k] <- p
    }
  }
  
  out <- data.frame(horizon = h,
                    forecast_origin = d$month[test_rows],
                    target_date = d$target_date[test_rows],
                    actual = d$target[test_rows],
                    predicted = preds)
  write_csv(out, file.path(OUTPUT_DIR, paste0("ar_h", h, ".csv")))
  
  list(forecast = out, p = sel_p)
}

# run all horizons
res <- lapply(horizons, run_horizon)
names(res) <- paste0("h", horizons)

# summary 
summary_tbl <- do.call(rbind, lapply(res, function(r) {
  f <- r$forecast
  ends <- list(pre_covid = PRE_COVID_END, through_dec_2021 = TEST_END)
  do.call(rbind, lapply(names(ends), function(nm) {
    z <- f[f$target_date <= ends[[nm]], ]
    e <- z$actual - z$predicted
    data.frame(horizon = f$horizon[1], period = nm, test_months = nrow(z),
               RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)))
  }))
}))
rownames(summary_tbl) <- NULL
print(summary_tbl, row.names = FALSE)
cat("Completed. Files saved under: ", OUTPUT_DIR, "\n", sep = "")