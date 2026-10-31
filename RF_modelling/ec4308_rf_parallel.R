# EC4308: Rolling Random Forest forecasts for Japan industrial production (x1)
# Adapted from the professor's runrf() / rf.rolling.window() structure.
# Inputs: h1.csv, h3.csv, h6.csv, h12.csv in DATA_DIR.
# R packages: install.packages(c('randomForest', 'future', 'future.apply'))
# NOTE: This is computationally intensive: each candidate mtry is refitted monthly.

library(randomForest)
library(future)
library(future.apply)
set.seed(4308)
DATA_DIR <- "."             # change to the folder containing the four CSV files
OUTPUT_DIR <- "rf_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
TRAIN_N <- 120L
VALID_N <- 100L
NTREE_TUNE <- 500L        # 500 trees for initial tuneRF OOB search
NTREE_ROLL <- 500L        # trees for all rolling validation / test fits
N_WORKERS <- 4L          # independent horizon processes; lower if memory is limited
STEP_FACTOR <- 1.5
IMPROVE <- 0.01


read_horizon <- function(h) {
  f <- file.path(DATA_DIR, paste0("h", h, ".csv"))
  d <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("month", "target_date", "target") %in% names(d)))
  d$month <- as.Date(d$month)
  d$target_date <- as.Date(d$target_date)
  predictors <- grep("^x[0-9]+_lag[0-3]$", names(d), value = TRUE)
  stopifnot(length(predictors) == 900L)
  stopifnot(all(complete.cases(d[, c(predictors, "target")])))
  stopifnot(!is.unsorted(d$month))
  
  list(data = d, predictors = predictors)
}

score <- function(actual, prediction) {
  c(RMSE = sqrt(mean((actual - prediction)^2)),
    MAE = mean(abs(actual - prediction)))
}

# A forecast issued at 'month[i]' may only train on rows for which
# target_date <= month[i]. Use the LAST 120 such labeled observations.
# This additionally guards against look-ahead bias for horizons h > 1.
training_rows <- function(d, i) {
  eligible <- which(d$target_date <= d$month[i] & seq_len(nrow(d)) < i)
  if (length(eligible) < TRAIN_N) return(integer(0))
  tail(eligible, TRAIN_N)
}


runrf <- function(d, predictors, i, mtry, ntree = NTREE_ROLL) {
  idx <- training_rows(d, i)
  if (length(idx) != TRAIN_N) stop("Not enough observed labels at forecast origin")
  x <- d[idx, predictors, drop = FALSE]
  y <- d$target[idx]
  model <- randomForest(x = x, y = y,
                        mtry = mtry, ntree = ntree,
                        importance = FALSE)
  new_x <- d[i, predictors, drop = FALSE]
  pred <- as.numeric(predict(model, newdata = new_x))
  list(model = model, pred = pred)
}

# re-estimate at every forecast origin.
rf.rolling.window <- function(d, predictors, rows, mtry, ntree = NTREE_ROLL) {
  pred <- numeric(length(rows))
  for (j in seq_along(rows)) {
    pred[j] <- runrf(d, predictors, rows[j], mtry, ntree)$pred
    if (j %% 20L == 0L) cat("  mtry", mtry, "forecast", j, "/", length(rows), "\n")
  }
  data.frame(month = d$month[rows], target_date = d$target_date[rows],
             actual = d$target[rows], prediction = pred,
             error = d$target[rows] - pred)
}


# selection using rolling predictions over the 100 validation observations.
# Validation and test outcomes are NEVER included when fitting their forecasts.
run_horizon <- function(h) {
  set.seed(4308L + as.integer(h))  # reproducible and distinct per horizon
  obj <- read_horizon(h)
  d <- obj$data
  p <- obj$predictors
  available <- which(vapply(seq_len(nrow(d)), function(i)
    length(training_rows(d, i)) == TRAIN_N, logical(1)))
  stopifnot(length(available) > VALID_N)
  val_rows <- head(available, VALID_N)
  first_train <- training_rows(d, val_rows[1])
  cat("\nHORIZON", h, "first validation target:",
      as.character(d$target_date[val_rows[1]]), "\n")

 
  oob_model <- tuneRF(
    x = d[first_train, p, drop = FALSE], y = d$target[first_train],
    mtryStart = floor(length(p) / 3), ntreeTry = NTREE_TUNE,
    stepFactor = STEP_FACTOR, improve = IMPROVE,
    trace = TRUE, plot = FALSE, doBest = TRUE)
  oob_mtry <- as.integer(oob_model$mtry)

  candidate_mtry <- sort(unique(pmax(1L, pmin(length(p),
                        as.integer(round(oob_mtry * c(1 / STEP_FACTOR, 1, STEP_FACTOR)))))))
  cat("OOB best mtry:", oob_mtry, "Validation candidates:",
      paste(candidate_mtry, collapse = ", "), "\n")

  val_results <- vector("list", length(candidate_mtry))
  for (k in seq_along(candidate_mtry)) {
    m <- candidate_mtry[k]
    f <- rf.rolling.window(d, p, val_rows, m)
    sc <- score(f$actual, f$prediction)
    val_results[[k]] <- data.frame(h = h, mtry = m,
                                    validation_RMSE = sc["RMSE"],
                                    validation_MAE = sc["MAE"])
  }
  validation_scores <- do.call(rbind, val_results)
  # Select on validation RMSE; use validation MAE to break ties.
  validation_scores <- validation_scores[order(
    validation_scores$validation_RMSE, validation_scores$validation_MAE), ]
  best_mtry <- validation_scores$mtry[1]
  cat("Selected validation mtry:", best_mtry, "\n")


  test_rows <- available[seq_along(available) > VALID_N]
  test_rows <- test_rows[d$target_date[test_rows] <= as.Date("2021-12-01")]
  stopifnot(
    max(d$target_date[test_rows]) == as.Date("2021-12-01")
  )
  stopifnot(length(test_rows) > 0)
  test_forecasts <- rf.rolling.window(d, p, test_rows, best_mtry)

  cutoffs <- c(pre_covid = "2019-12-01",
               through_dec_2021 = "2021-12-01")
  summary_rows <- lapply(names(cutoffs), function(period) {
    last_date <- as.Date(cutoffs[[period]])
    z <- test_forecasts[test_forecasts$target_date <= last_date, ]
    sc <- score(z$actual, z$prediction)
    data.frame(
      horizon = h,
      period = period,
      selected_mtry = best_mtry,
      RMSE = unname(sc["RMSE"]),
      MAE = unname(sc["MAE"])
    )
  })
  do.call(rbind, summary_rows)
}


workers <- min(N_WORKERS, 4L, future::availableCores())
future::plan(future::multisession, workers = workers)
results_by_h <- tryCatch(
  future.apply::future_lapply(c(1L, 3L, 6L, 12L), run_horizon,
                              future.seed = TRUE),
  finally = future::plan(future::sequential)
)
results <- do.call(rbind, results_by_h)
write.csv(results, file.path(OUTPUT_DIR, "rf_summary_results.csv"), row.names = FALSE)
print(results, row.names = FALSE)
