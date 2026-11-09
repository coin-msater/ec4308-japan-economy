#load the forecast files and compare the models
models <- c("ar", "lasso", "ridge", "elastic_net", "rf", "ensemble")
horizons <- c(1, 3, 6, 12)

#use the same test periods as the model evaluations
periods <- list(
  "Pre-COVID" = c("1992-09-01", "2019-12-01"),
  "Full test" = c("1992-09-01", "2021-12-01")
)

#set up the model pairs for the dm tests
comparisons <- list(
  c("ar", "lasso"),
  c("ar", "ridge"),
  c("ar", "elastic_net"),
  c("ar", "rf"),
  c("ar", "ensemble"),
  c("rf", "ridge"),
  c("rf", "ensemble"),
  c("ridge", "ensemble")
)

#find each model's csv inside its folder
get_path <- function(model, h) {
  folder <- if (model == "elastic_net") "en" else model
  filename <- if (model %in% c("ar", "ensemble")) {
    sprintf("%s_h%d.csv", model, h)
  } else {
    sprintf("%s_ensemble_h%d.csv", model, h)
  }
  path <- file.path(folder, filename)
  if (!file.exists(path)) stop("Missing file: ", path)
  path
}

#convert dates to month numbers for the horizon check
month_number <- function(x) {
  as.integer(format(x, "%Y")) * 12L + as.integer(format(x, "%m"))
}

#run a two sided dm test using squared errors
#bartlett weights account for overlapping multi step forecast errors
dm_test <- function(error_1, error_2, h) {
  d <- error_1^2 - error_2^2
  n <- length(d)
  centered <- d - mean(d)
  gamma <- vapply(0:(h - 1), function(k) {
    if (k == 0) return(sum(centered^2) / n)
    sum(centered[(k + 1):n] * centered[1:(n - k)]) / n
  }, numeric(1))
  weights <- 1 - (0:(h - 1)) / h
  long_run_var <- gamma[1] + 2 * sum(weights[-1] * gamma[-1])
  if (!is.finite(long_run_var) || long_run_var <= 0) {
    stop("Could not calculate the DM variance")
  }
  dm <- mean(d) / sqrt(long_run_var / n)

  #apply the small sample correction
  correction <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  dm <- dm * correction
  p <- 2 * pt(-abs(dm), df = n - 1)
  c(dm_statistic = dm, p_value_two_sided = p)
}

results <- list()

for (h in horizons) {
  #read the six models for this horizon
  forecasts <- setNames(lapply(models, function(model) {
    path <- get_path(model, h)
    df <- read.csv(path, stringsAsFactors = FALSE)
    df$forecast_origin <- as.Date(df$forecast_origin)
    df$target_date <- as.Date(df$target_date)
    df <- df[order(df$target_date), ]
    if (anyNA(df[, c("forecast_origin", "target_date", "actual", "predicted")])) {
      stop("Missing values in ", path)
    }
    if (anyDuplicated(df$target_date)) stop("Duplicate target dates in ", path)
    if (!all(month_number(df$target_date) - month_number(df$forecast_origin) == h)) {
      stop("Incorrect horizon in ", path)
    }
    df
  }), models)

  #check that all models have matching dates and actual values
  ar <- forecasts$ar
  for (model in models[-1]) {
    other <- forecasts[[model]]
    if (!identical(ar$target_date, other$target_date) ||
        !identical(ar$forecast_origin, other$forecast_origin) ||
        !isTRUE(all.equal(ar$actual, other$actual, tolerance = 1e-9))) {
      stop("Forecasts do not align for ", model, " H", h)
    }
  }

  #check that the ensemble is the average of the four ml models
  ensemble_check <- (forecasts$lasso$predicted + forecasts$ridge$predicted +
                     forecasts$elastic_net$predicted + forecasts$rf$predicted) / 4
  if (!isTRUE(all.equal(ensemble_check, forecasts$ensemble$predicted,
                        tolerance = 1e-10))) {
    stop("Ensemble predictions do not match for H", h)
  }

  for (period in names(periods)) {
    #keep only the target dates in this test period
    dates <- as.Date(periods[[period]])
    keep <- ar$target_date >= dates[1] & ar$target_date <= dates[2]
    if (!any(keep)) stop("No forecasts in ", period)

    for (pair in comparisons) {
      #calculate errors for the two models being compared
      model_1 <- pair[1]
      model_2 <- pair[2]
      first <- forecasts[[model_1]]
      second <- forecasts[[model_2]]
      error_1 <- first$actual[keep] - first$predicted[keep]
      error_2 <- second$actual[keep] - second$predicted[keep]

      #run the dm test and calculate rmse
      test <- dm_test(error_1, error_2, h)
      rmse_1 <- sqrt(mean(error_1^2))
      rmse_2 <- sqrt(mean(error_2^2))

      results[[length(results) + 1]] <- data.frame(
        period = period,
        horizon = h,
        comparison = paste(toupper(model_1), "vs", toupper(model_2)),
        model_1 = model_1,
        model_2 = model_2,
        rmse_1 = rmse_1,
        rmse_2 = rmse_2,
        mean_loss_diff_first_minus_second = mean(error_1^2 - error_2^2),
        dm_statistic = unname(test["dm_statistic"]),
        p_value_two_sided = unname(test["p_value_two_sided"]),
        significant_5pct = unname(test["p_value_two_sided"] < 0.05),
        better_rmse = if (rmse_1 < rmse_2) model_1 else if (rmse_2 < rmse_1) model_2 else "tie"
      )
    }
  }
}

#save the results without mae or sample size
output <- do.call(rbind, results)
write.csv(output, "dm_results.csv", row.names = FALSE)
print(output[, c("period", "horizon", "comparison", "rmse_1", "rmse_2",
                 "dm_statistic", "p_value_two_sided")], row.names = FALSE)
