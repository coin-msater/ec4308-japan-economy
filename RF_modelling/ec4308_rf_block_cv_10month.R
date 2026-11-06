

library(randomForest)
library(future)
library(future.apply)

DATA_DIR <- "."        # Folder containing h1.csv, h3.csv, h6.csv, h12.csv
OUTPUT_DIR <- "rf_results_block_cv"  
if (dir.exists(OUTPUT_DIR) && length(list.files(OUTPUT_DIR)) > 0L) {
  stop("Output directory already has files: ", OUTPUT_DIR,
       ". Change OUTPUT_DIR to a new name to avoid overwriting them.")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAIN_N <- 120L          # Latest 120 eligible observations 
MIN_TRAIN_ROWS <- 50L
N_FOLDS <- 10L          # 10 chronological validation groups
NTREE <- 500L
N_WORKERS <- 4L
MTRY_GRID <- c(100L, 300L, 600L, 900L)
RETUNE_MONTHS <- 100L

VALID_START <- as.Date("1984-05-01")
VALID_END <- as.Date("1992-08-01")
TEST_START <- as.Date("1992-09-01")
PRE_COVID_END <- as.Date("2019-12-01")
TEST_END <- as.Date("2021-12-01")
RETUNE_TARGETS <- as.Date(c("1992-09-01", "2001-01-01",
                           "2009-05-01", "2017-09-01"))


CHECK_ONLY <- FALSE

months_before <- function(date, n) {
  tail(seq(as.Date(date), by = "-1 month", length.out = n + 1L), 1L)
}

read_horizon <- function(h) {
  f <- file.path(DATA_DIR, paste0("h", h, ".csv"))
  if (!file.exists(f)) stop("Missing input CSV: ", f)
  d <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("month", "target_date", "target") %in% names(d)))
  d$month <- as.Date(d$month)
  d$target_date <- as.Date(d$target_date)
  p <- grep("^x[0-9]+_lag[0-3]$", names(d), value = TRUE)
  stopifnot(length(p) == 900L, !anyDuplicated(d$month),
            !anyDuplicated(d$target_date), !is.unsorted(d$month),
            all(complete.cases(d[, c(p, "target")])))
  expected_target_dates <- seq(d$month[1], by = "month",
                               length.out = nrow(d) + h)[seq_len(nrow(d)) + h]
  stopifnot(identical(as.character(d$target_date),
                      as.character(expected_target_dates)))
  list(d = d, p = p)
}


training_rows <- function(d, origin) {
  eligible <- which(d$month < origin & d$target_date <= origin)
  tail(eligible, TRAIN_N)
}

fit_rf <- function(d, p, origin, mtry, h, seed_id) {
  idx <- training_rows(d, origin)
  if (length(idx) < MIN_TRAIN_ROWS)
    stop("H", h, ": only ", length(idx),
         " eligible training observations at ", as.character(origin))
  # consistent seed for a given validation fold/origin.
  set.seed(as.integer(430800L + 1000L * h + seed_id))
  randomForest(x = d[idx, p, drop = FALSE], y = d$target[idx],
               mtry = as.integer(mtry), ntree = NTREE, importance = FALSE)
}

initial_cv_rows <- function(d, first_test_origin) {
  which(d$month >= VALID_START & d$month <= VALID_END &
        d$target_date <= first_test_origin)
}

retune_cv_rows <- function(d, retune_origin) {
  start <- months_before(retune_origin, RETUNE_MONTHS - 1L)
  which(d$target_date >= start & d$target_date <= retune_origin &
        d$month < retune_origin)
}

# Fit exactly 1 rf per fold for each mtry. The same fitted model predicts every validation observation in its fold
cv_tune_blockwise <- function(d, p, rows, h, tuning_origin, label) {
  rows <- rows[order(d$target_date[rows])]
  stopifnot(length(rows) >= N_FOLDS,
            all(d$target_date[rows] <= tuning_origin))
  n <- length(rows)
  # Always 10 chronological folds
  fold_ids <- as.integer(cut(seq_len(n), breaks = N_FOLDS, labels = FALSE))
  folds <- split(rows, fold_ids)
  candidates <- unique(MTRY_GRID[MTRY_GRID >= 1 & MTRY_GRID <= length(p)])
  if (!length(candidates)) stop("No mtry candidates available")

  cat("H", h, " ", label, ": ", n, " validation forecasts, ",
      length(folds), " chronological blocks (ONE fit per block / mtry)\n", sep = "")

  scores <- lapply(candidates, function(m) {
    actual <- numeric(0)
    predicted <- numeric(0)
    for (k in seq_along(folds)) {
      fold_rows <- folds[[k]]
      
      origin <- min(d$month[fold_rows])
      model <- fit_rf(d, p, origin, mtry = m, h = h,
                      seed_id = match(origin, d$month))
      fold_pred <- as.numeric(predict(model,
                             newdata = d[fold_rows, p, drop = FALSE]))
      predicted <- c(predicted, fold_pred)
      actual <- c(actual, d$target[fold_rows])
    }
    data.frame(mtry = m, RMSE = sqrt(mean((actual - predicted)^2)),
               MAE = mean(abs(actual - predicted)), n_forecasts = n,
               n_model_fits = length(folds))
  })
  tab <- do.call(rbind, scores)
  tab <- tab[order(tab$RMSE, tab$MAE), , drop = FALSE]
  best <- as.integer(tab$mtry[1])
  cat("H", h, " ", label, ": selected mtry=", best,
      "; RMSE=", round(tab$RMSE[1], 4), "\n", sep = "")
  list(mtry = best, scores = tab, n_cv = n, n_folds = length(folds))
}

run_horizon <- function(h) {
  obj <- read_horizon(h)
  d <- obj$d
  p <- obj$p

  test_rows <- which(d$target_date >= TEST_START & d$target_date <= TEST_END)
  expected <- seq(TEST_START, TEST_END, by = "month")
  stopifnot(identical(as.character(d$target_date[test_rows]),
                      as.character(expected)))
  block_id <- (seq_along(test_rows) - 1L) %/% RETUNE_MONTHS + 1L
  stopifnot(max(block_id) == length(RETUNE_TARGETS))
  first_positions <- which(!duplicated(block_id))
  stopifnot(identical(as.character(d$target_date[test_rows[first_positions]]),
                      as.character(RETUNE_TARGETS)))

  cv_rows <- vector("list", length(RETUNE_TARGETS))
  cv_rows[[1]] <- initial_cv_rows(d, d$month[test_rows[1]])
  for (b in 2:length(cv_rows)) {
    retune_origin <- d$month[test_rows[first_positions[b]]]
    cv_rows[[b]] <- retune_cv_rows(d, retune_origin)
  }
  for (b in seq_along(cv_rows)) {
    origin <- d$month[test_rows[first_positions[b]]]
    stopifnot(length(cv_rows[[b]]) >= N_FOLDS,
              all(d$target_date[cv_rows[[b]]] <= origin))
    cat("H", h, " block ", b, ": retune origin=", as.character(origin),
        ", known validation targets=", length(cv_rows[[b]]), "\n", sep = "")
  }
  if (CHECK_ONLY) {
    return(data.frame(horizon = h,
                      initial_cv = length(cv_rows[[1]]),
                      retune2_cv = length(cv_rows[[2]]),
                      retune3_cv = length(cv_rows[[3]]),
                      retune4_cv = length(cv_rows[[4]])))
  }

  preds <- numeric(length(test_rows))
  selected <- integer(length(test_rows))
  for (b in seq_along(cv_rows)) {
    origin <- d$month[test_rows[first_positions[b]]]
    tune <- cv_tune_blockwise(d, p, cv_rows[[b]], h, origin,
                              paste0("retune_block_", b))
    m <- tune$mtry
    positions <- which(block_id == b)
    # KEEP MONTHLY TEST RETRAINING, matching the original RF forecasting rule.
    for (j in seq_along(positions)) {
      k <- positions[j]
      row <- test_rows[k]
      model <- fit_rf(d, p, d$month[row], mtry = m, h = h, seed_id = row)
      preds[k] <- as.numeric(predict(model,
                            newdata = d[row, p, drop = FALSE]))
      selected[k] <- m
      if (j %% 20L == 0L || j == length(positions)) {
        cat("H", h, " block ", b, " test ", j, "/",
            length(positions), "\n", sep = "")
      }
    }
  }

  forecast <- data.frame(horizon = h, month = d$month[test_rows],
                         target_date = d$target_date[test_rows],
                         actual = d$target[test_rows], prediction = preds,
                         error = d$target[test_rows] - preds,
                         selected_mtry = selected, retune_block = block_id)


  write.csv(forecast,
            file.path(OUTPUT_DIR, paste0("blockcv_forecasts_h", h, ".csv")),
            row.names = FALSE)
  ensemble_rf <- data.frame(horizon = h,
                            forecast_origin = forecast$month,
                            target_date = forecast$target_date,
                            actual = forecast$actual,
                            predicted = forecast$prediction)
  write.csv(ensemble_rf,
            file.path(OUTPUT_DIR, paste0("blockcv_rf_ensemble_h", h, ".csv")),
            row.names = FALSE)


  dates <- c(pre_covid = PRE_COVID_END, through_dec_2021 = TEST_END)
  do.call(rbind, lapply(names(dates), function(name) {
    z <- forecast[forecast$target_date <= dates[name], , drop = FALSE]
    data.frame(horizon = h, period = name, test_months = nrow(z),
               RMSE = sqrt(mean(z$error^2)), MAE = mean(abs(z$error)))
  }))
}

workers <- min(N_WORKERS, 4L, future::availableCores())
future::plan(future::multisession, workers = workers)
results_by_h <- tryCatch(
  future.apply::future_lapply(c(1L, 3L, 6L, 12L), run_horizon,
                              future.seed = TRUE),
  finally = future::plan(future::sequential)
)
results <- do.call(rbind, results_by_h)
if (CHECK_ONLY) {
  print(results, row.names = FALSE)
} else {
  write.csv(results,
            file.path(OUTPUT_DIR, "blockcv_rf_summary_8_results.csv"),
            row.names = FALSE)
  files <- file.path(OUTPUT_DIR,
                     paste0("blockcv_rf_ensemble_h", c(1, 3, 6, 12), ".csv"))
  combined <- do.call(rbind, lapply(files, read.csv))
  combined <- combined[order(combined$horizon, combined$target_date), ]
  rownames(combined) <- NULL
  write.csv(combined,
            file.path(OUTPUT_DIR, "blockcv_rf_ensemble_all_horizons.csv"),
            row.names = FALSE)
  print(results, row.names = FALSE)
  cat("Completed. Files saved under: ", OUTPUT_DIR, "\n", sep = "")
}
