
library(randomForest)
library(future)
library(future.apply)

DATA_DIR <- "."  # Put h1.csv, h3.csv, h6.csv and h12.csv
OUTPUT_DIR <- "rf_results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAIN_MONTHS <- 120L
VALID_MONTHS <- 100L
CV_FOLDS <- 10L
RETUNE_MONTHS <- 100L
NTREE <- 500L
N_WORKERS <- 4L
MIN_TRAIN_ROWS <- 50L  
MTRY_GRID <- c(100L, 300L, 600L, 900L)  

TRAIN_START <- as.Date("1974-05-01")
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
  if (!file.exists(f)) stop("Input file not found: ", f)
  d <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("month", "target_date", "target") %in% names(d)))
  d$month <- as.Date(d$month)
  d$target_date <- as.Date(d$target_date)
  p <- grep("^x[0-9]+_lag[0-3]$", names(d), value = TRUE)
  stopifnot(length(p) == 900L, !anyDuplicated(d$month),
            !anyDuplicated(d$target_date),
            !is.unsorted(d$month),
            all(complete.cases(d[, c(p, "target")])))
  stopifnot(all(as.character(d$target_date) ==
                as.character(seq(d$month[1], by = "month", length.out = nrow(d) + h)[
                  seq_len(nrow(d)) + h])))
  list(d = d, p = p)
}


training_rows <- function(d, row) {
  origin <- d$month[row]
  start <- months_before(origin, TRAIN_MONTHS)
  which(d$month >= start & d$month < origin &
        d$target_date <= origin)
}

fit_predict <- function(d, p, row, mtry, h, ntree = NTREE) {
  idx <- training_rows(d, row)
  if (length(idx) < MIN_TRAIN_ROWS) {
    stop("Too few known labels at ", as.character(d$month[row]),
         " for H", h, ": ", length(idx))
  }
  # Same RNG seed for different mtry candidates
  set.seed(430800L + 1000L * as.integer(h) + as.integer(row))
  fit <- randomForest(x = d[idx, p, drop = FALSE], y = d$target[idx],
                      mtry = as.integer(mtry), ntree = ntree,
                      importance = FALSE)
  as.numeric(predict(fit, newdata = d[row, p, drop = FALSE]))
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


cv_tune <- function(d, p, rows, h, label) {
  rows <- rows[order(d$target_date[rows])]
  n <- length(rows)
  if (n < CV_FOLDS) stop("Too few validation rows for 10 folds: ", label)
  folds <- ceiling(seq_len(n) * CV_FOLDS / n)
  candidates <- MTRY_GRID[MTRY_GRID >= 1L & MTRY_GRID <= length(p)]
  candidates <- unique(candidates)
  if (!length(candidates)) stop("No valid mtry candidates")

  cat("H", h, " ", label, ": ", n, " historical validation forecasts / ",
      length(unique(folds)), " chronological folds; targets ",
      as.character(min(d$target_date[rows])), " to ",
      as.character(max(d$target_date[rows])), "\n", sep = "")

  scores <- lapply(candidates, function(m) {
    pred <- numeric(n)
    for (j in seq_len(n)) {
      pred[j] <- fit_predict(d, p, rows[j], m, h)
    }
    errors <- d$target[rows] - pred
    fold_mse <- vapply(seq_len(CV_FOLDS), function(k) {
      mean(errors[folds == k]^2)
    }, numeric(1))
    data.frame(mtry = m,
               RMSE = sqrt(mean(errors^2)),
               MAE = mean(abs(errors)),
               min_fold_RMSE = sqrt(min(fold_mse)),
               max_fold_RMSE = sqrt(max(fold_mse)))
  })
  tab <- do.call(rbind, scores)
  tab <- tab[order(tab$RMSE, tab$MAE), ]
  best <- as.integer(tab$mtry[1])
  cat("H", h, " ", label, ": selected mtry=", best,
      " (CV RMSE=", round(tab$RMSE[1], 4), ")\n", sep = "")
  list(mtry = best, scores = tab, n_cv = n, n_folds = CV_FOLDS)
}

run_horizon <- function(h) {
  obj <- read_horizon(h)
  d <- obj$d
  p <- obj$p
  test_rows <- which(d$target_date >= TEST_START &
                       d$target_date <= TEST_END)
  expected_targets <- seq(TEST_START, TEST_END, by = "month")
  if (!identical(as.character(d$target_date[test_rows]),
                 as.character(expected_targets))) {
    stop("Test TARGET dates are incomplete or misaligned for H", h)
  }

  first_test_origin <- d$month[test_rows[1]]
  init_rows <- initial_cv_rows(d, first_test_origin)
  if (!length(init_rows)) stop("No usable initial validation for H", h)

  block_id <- (seq_along(test_rows) - 1L) %/% RETUNE_MONTHS + 1L
  nblocks <- max(block_id)
  stopifnot(nblocks == length(RETUNE_TARGETS))
  actual_boundaries <- d$target_date[test_rows[which(!duplicated(block_id))]]
  stopifnot(identical(as.character(actual_boundaries),
                      as.character(RETUNE_TARGETS)))

  retune_rows <- vector("list", nblocks)
  retune_rows[[1]] <- init_rows
  for (b in 2:nblocks) {
    first_row <- test_rows[which(block_id == b)[1]]
    origin <- d$month[first_row]
    retune_rows[[b]] <- retune_cv_rows(d, origin)
  }

  cat("\nH", h, ": first test origin=", as.character(first_test_origin),
      "; initial validation forecasts=", length(init_rows), "\n", sep = "")
  for (b in seq_len(nblocks)) {
    rr <- retune_rows[[b]]
    cat("  block ", b, " target start=", as.character(RETUNE_TARGETS[b]),
        ", tuning rows=", length(rr),
        ", last known target=", as.character(max(d$target_date[rr])), "\n", sep = "")
  }

  if (CHECK_ONLY) {
    return(data.frame(horizon = h,
                      first_test_origin = first_test_origin,
                      initial_validation_forecasts = length(init_rows),
                      retune2_forecasts = length(retune_rows[[2]]),
                      retune3_forecasts = length(retune_rows[[3]]),
                      retune4_forecasts = length(retune_rows[[4]])))
  }

  pred <- numeric(length(test_rows))
  selected_mtry <- integer(length(test_rows))
  cv_log <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    first_row <- test_rows[which(block_id == b)[1]]
    retune_origin <- d$month[first_row]
    stopifnot(all(d$target_date[retune_rows[[b]]] <= retune_origin))
    tune <- cv_tune(d, p, retune_rows[[b]], h,
                    label = paste0("retune_block_", b))
    m <- tune$mtry
    cv_log[[b]] <- data.frame(horizon = h, block = b,
                              first_test_target = as.character(d$target_date[first_row]),
                              first_test_origin = as.character(retune_origin),
                              cv_forecasts = tune$n_cv,
                              cv_folds = tune$n_folds,
                              selected_mtry = m,
                              cv_RMSE = tune$scores$RMSE[1])
    positions <- which(block_id == b)
    for (j in seq_along(positions)) {
      k <- positions[j]
      pred[k] <- fit_predict(d, p, test_rows[k], m, h)
      selected_mtry[k] <- m
      if (j %% 20L == 0L || j == length(positions)) {
        cat("H", h, " block ", b, " test forecast ", j, "/",
            length(positions), "\n", sep = "")
      }
    }
  }

  forecast <- data.frame(
    horizon = h, month = d$month[test_rows],
    target_date = d$target_date[test_rows],
    actual = d$target[test_rows], prediction = pred,
    error = d$target[test_rows] - pred,
    selected_mtry = selected_mtry, retune_block = block_id
  )
  write.csv(forecast, file.path(OUTPUT_DIR,
                                 paste0("forecasts_h", h, ".csv")),
            row.names = FALSE)
  ensemble_rf <- data.frame(
    horizon = forecast$horizon,
    forecast_origin = forecast$month,
    target_date = forecast$target_date,
    actual = forecast$actual,
    predicted = forecast$prediction
  )
  write.csv(ensemble_rf,
            file.path(OUTPUT_DIR, paste0("rf_ensemble_h", h, ".csv")),
            row.names = FALSE)

  dates <- c(pre_covid = PRE_COVID_END, through_dec_2021 = TEST_END)
  out <- lapply(names(dates), function(name) {
    z <- forecast[forecast$target_date <= dates[name], ]
    data.frame(horizon = h, period = name,
               RMSE = sqrt(mean(z$error^2)),
               MAE = mean(abs(z$error)))
  })
  do.call(rbind, out)
}

# Run each horizon in a separate process when N_WORKERS > 1.
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
  write.csv(results, file.path(OUTPUT_DIR, "rf_summary_8_results.csv"),
            row.names = FALSE)
  print(results, row.names = FALSE)
}

if (!CHECK_ONLY) {
  ensemble_files <- file.path(OUTPUT_DIR,
                              paste0("rf_ensemble_h", c(1, 3, 6, 12), ".csv"))
  rf_ensemble_all <- do.call(rbind, lapply(ensemble_files, read.csv))
  rf_ensemble_all <- rf_ensemble_all[
    order(rf_ensemble_all$horizon, rf_ensemble_all$target_date), ]
  rownames(rf_ensemble_all) <- NULL
  write.csv(rf_ensemble_all,
            file.path(OUTPUT_DIR, "rf_ensemble_all_horizons.csv"),
            row.names = FALSE)
  cat("Saved ensemble-ready RF forecasts: ",
      file.path(OUTPUT_DIR, "rf_ensemble_all_horizons.csv"), "\n", sep = "")
}
