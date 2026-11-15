
#load packages for random forest and parallel processing
library(randomForest)
library(future)
library(future.apply)

#set input and output directories
DATA_DIR <- "."        #folder containing h1.csv, h3.csv, h6.csv, h12.csv
OUTPUT_DIR <- "rf_results"  
if (dir.exists(OUTPUT_DIR) && length(list.files(OUTPUT_DIR)) > 0L) {
  stop("Output directory already has files: ", OUTPUT_DIR,
       ". Change OUTPUT_DIR to a new name to avoid overwriting them.")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

#set training window and random forest parameters
TRAIN_N <- 120L          #latest 120 eligible observations (initial sample starts May 1974) 
MIN_TRAIN_ROWS <- 50L #minimum number of observations required to fit a forest
NTREE <- 500L #number of trees per forest
N_WORKERS <- 4L #maximum number of parallel workers
TUNE_STEP <- 1.5 #minimum relative oob error improvement to continue the search
TUNE_IMPROVE <- 0.01
RETUNE_MONTHS <- 100L     #deliberately retune only every 100 test months

#define evaluation periods using target dates
TEST_START <- as.Date("1992-09-01")
PRE_COVID_END <- as.Date("2019-12-01")
TEST_END <- as.Date("2021-12-01")
#set the first target date of each retuning block
RETUNE_TARGETS <- as.Date(c("1992-09-01", "2001-01-01",
                           "2009-05-01", "2017-09-01"))


CHECK_ONLY <- FALSE

#load and validate the data for each forecast horizon
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

#select the latest eligible observations for the rolling training window
training_rows <- function(d, origin) {
  eligible <- which(d$month < origin & d$target_date <= origin)
  tail(eligible, TRAIN_N)
}

#fit a new random forest at the current forecast origin
fit_rf <- function(d, p, origin, mtry, h, seed_id) {
  idx <- training_rows(d, origin)
  if (length(idx) < MIN_TRAIN_ROWS)
    stop("H", h, ": only ", length(idx),
         " eligible training observations at ", as.character(origin))
  #consistent seed for a given forecast origin.
  set.seed(as.integer(430800L + 1000L * h + seed_id))
  #fit the forest using the selected mtry and enable variable importance
  randomForest(x = d[idx, p, drop = FALSE], y = d$target[idx],
               mtry = as.integer(mtry), ntree = NTREE, importance = TRUE)
}

#tune mtry on the latest eligible rolling-training observations at each
#retuning origin. OOB error only; no validation set or chronological CV.
tune_mtry <- function(d, p, rows, h, tuning_origin, block) {
  stopifnot(length(rows) == TRAIN_N,
            all(d$target_date[rows] <= tuning_origin),
            all(d$month[rows] < tuning_origin))
  set.seed(as.integer(430800L + 1000L * h + match(tuning_origin, d$month)))
  #search for mtry using out-of-bag mean squared error
  scores <- randomForest::tuneRF(
    x = d[rows, p, drop = FALSE], y = d$target[rows],
    mtryStart = max(1L, floor(length(p) / 3)),
    ntreeTry = NTREE, stepFactor = TUNE_STEP, improve = TUNE_IMPROVE,
    trace = TRUE, plot = FALSE, doBest = FALSE
  )
  #select the candidate with the lowest observed oob mse
  best <- which.min(scores[, 2L])
  m <- as.integer(scores[best, 1L])
  cat("H", h, " block ", block, ": selected mtry=", m,
      "; OOB MSE=", scores[best, 2L], "\n", sep = "")
  list(mtry = m, scores = scores)
}

#run the complete forecasting procedure for one horizon
run_horizon <- function(h) {
  obj <- read_horizon(h)
  d <- obj$d
  p <- obj$p
  #select test observations by target date and verify monthly coverage
  test_rows <- which(d$target_date >= TEST_START & d$target_date <= TEST_END)
  expected <- seq(TEST_START, TEST_END, by = "month")
  stopifnot(identical(as.character(d$target_date[test_rows]),
                      as.character(expected)))
  #divide the test period into consecutive blocks of 100 target months
  block_id <- (seq_along(test_rows) - 1L) %/% RETUNE_MONTHS + 1L
  stopifnot(max(block_id) == length(RETUNE_TARGETS))
  #identify and verify the first forecast position of each retuning block
  first_positions <- which(!duplicated(block_id))
  stopifnot(identical(as.character(d$target_date[test_rows[first_positions]]),
                      as.character(RETUNE_TARGETS)))

  #check the training sample at each retune origin without running the forest.
  retune_train_n <- vapply(seq_along(first_positions), function(b) {
    origin <- d$month[test_rows[first_positions[b]]]
    rows <- training_rows(d, origin)
    stopifnot(length(rows) >= MIN_TRAIN_ROWS,
              all(d$month[rows] < origin),
              all(d$target_date[rows] <= origin))
    cat("H", h, " block ", b, ": retune origin=", as.character(origin),
        ", eligible training rows=", length(rows), "\n", sep = "")
    length(rows)
  }, integer(1L))
  if (CHECK_ONLY) {
    return(data.frame(horizon = h,
                      retune_block = seq_along(first_positions),
                      retune_target = as.character(d$target_date[test_rows[first_positions]]),
                      retune_origin = as.character(d$month[test_rows[first_positions]]),
                      training_rows = retune_train_n))
  }

  preds <- numeric(length(test_rows))
  selected <- integer(length(test_rows))
  vi <- matrix(NA_real_, nrow = length(test_rows), ncol = length(p),
               dimnames = list(NULL, p))
  #retune mtry once at the beginning of each test block
  for (b in seq_along(first_positions)) {
    origin <- d$month[test_rows[first_positions[b]]]
    rows <- training_rows(d, origin)
    tune <- tune_mtry(d, p, rows, h, origin, b)
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
      # OOB permutation importance for this monthly fitted forest.
      # Unscaled mean increase in MSE, in squared target units.
      imp <- randomForest::importance(model, type = 1, scale = FALSE)
      stopifnot(all(p %in% rownames(imp)))
      vi[k, ] <- imp[p, 1L]
      if (j %% 20L == 0L || j == length(positions)) {
        cat("H", h, " block ", b, " test ", j, "/",
            length(positions), "\n", sep = "")
      }
    }
  }

  #assemble forecasts, actual values, errors and tuning information
  forecast <- data.frame(horizon = h, month = d$month[test_rows],
                         target_date = d$target_date[test_rows],
                         actual = d$target[test_rows], prediction = preds,
                         error = d$target[test_rows] - preds,
                         selected_mtry = selected, retune_block = block_id)


  write.csv(forecast,
            file.path(OUTPUT_DIR, paste0("rf_forecasts_h", h, ".csv")),
            row.names = FALSE)
  ensemble_rf <- data.frame(horizon = h,
                            forecast_origin = forecast$month,
                            target_date = forecast$target_date,
                            actual = forecast$actual,
                            predicted = forecast$prediction)
  write.csv(ensemble_rf,
            file.path(OUTPUT_DIR, paste0("rf_ensemble_h", h, ".csv")),
            row.names = FALSE)


  #save one importance vector per forecast (wide format)
  stopifnot(all(is.finite(vi)))
  vi_monthly <- data.frame(horizon = h, forecast_origin = forecast$month,
                           target_date = forecast$target_date,
                           selected_mtry = selected, retune_block = block_id,
                           as.data.frame(vi), check.names = FALSE)
  write.csv(vi_monthly,
            file.path(OUTPUT_DIR, paste0("rf_importance_monthly_h", h, ".csv")),
            row.names = FALSE)

  #equal weight per monthly forecast, separately for each evaluation period.
  vi_periods <- list(pre_covid = which(forecast$target_date <= PRE_COVID_END),
                     through_dec_2021 = seq_len(nrow(forecast)))
  #calculate and rank average importance separately for each evaluation period
  vi_average <- do.call(rbind, lapply(names(vi_periods), function(name) {
    idx <- vi_periods[[name]]
    avg <- colMeans(vi[idx, , drop = FALSE])
    ord <- order(avg, decreasing = TRUE)
    data.frame(horizon = h, period = name, n_forecasts = length(idx),
               variable = p[ord], mean_importance = unname(avg[ord]),
               rank = seq_along(ord))
  }))
  write.csv(vi_average,
            file.path(OUTPUT_DIR, paste0("rf_importance_average_h", h, ".csv")),
            row.names = FALSE)

  #calculate rmse and mae for the pre-covid and full test periods
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
            file.path(OUTPUT_DIR, "rf_summary_8_results.csv"),
            row.names = FALSE)
  files <- file.path(OUTPUT_DIR,
                     paste0("rf_ensemble_h", c(1, 3, 6, 12), ".csv"))
  combined <- do.call(rbind, lapply(files, read.csv))
  combined <- combined[order(combined$horizon, combined$target_date), ]
  rownames(combined) <- NULL
  write.csv(combined,
            file.path(OUTPUT_DIR, "rf_ensemble_all_horizons.csv"),
            row.names = FALSE)
  print(results, row.names = FALSE)
  cat("Completed. Files saved under: ", OUTPUT_DIR, "\n", sep = "")
}
