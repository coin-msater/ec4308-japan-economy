############################################
# Sanity checks for Elastic Net output
# Run from the project root after elastic_net.R.
############################################

library(readr)

horizons <- c(1, 3, 6, 12)
window_size <- 120
validation_size <- 100
n_folds <- 10
pre_covid_end <- as.Date("2019-12-01")
expected_n <- length(seq(as.Date("1992-09-01"), as.Date("2021-12-01"), by = "month"))  # 352
standard_cols <- c("horizon", "target_date", "forecast", "actual")
results_dir <- "elastic_net/results"

res_path <- function(h, suffix = "") {
  sprintf("%s/elastic_net_h%d%s.csv", results_dir, h, suffix)
}

read_res <- function(h, suffix = "") {
  f <- as.data.frame(read_csv(res_path(h, suffix), show_col_types = FALSE))
  for (v in intersect(c("forecast_origin", "target_date", "tuning_date"), names(f))) {
    f[[v]] <- as.Date(f[[v]])
  }
  f
}

read_data <- function(h) {
  df <- as.data.frame(read_csv(sprintf("data/processed/h%d.csv", h), show_col_types = FALSE))
  df$month <- as.Date(df$month)
  df$target_date <- as.Date(df$target_date)
  df
}

mdiff <- function(a, b) {
  12 * (as.integer(format(b, "%Y")) - as.integer(format(a, "%Y"))) +
    as.integer(format(b, "%m")) - as.integer(format(a, "%m"))
}

get_training_idx <- function(origin, df, window_size = 120) {
  eligible <- which(df$month < origin & df$target_date <= origin)
  if (length(eligible) < window_size) return(integer(0))
  tail(eligible, window_size)
}


############################################
# 1. Structure of the output files
############################################

cat("\n== 1. Structure of output ==\n")
for (h in horizons) {
  f <- read_res(h)
  d <- read_res(h, "_details")
  gap <- mdiff(d$forecast_origin, d$target_date)
  cat(sprintf(
    "h=%2d | standard columns: %s | n=%d (expect %d) | NA forecasts=%d | dup targets=%d | target = origin + h: %s | targets %s to %s\n",
    h, identical(names(f), standard_cols), nrow(f), expected_n,
    sum(is.na(f$forecast)), sum(duplicated(f$target_date)), all(gap == h),
    format(min(f$target_date)), format(max(f$target_date))
  ))
}


############################################
# 2. CV design: counts vs what the design implies
#    Validation origins are the 100 months before the tuning date.
#    A forecast is scored only if its target is observed by the
#    tuning date, so 101 - h forecasts should be scored,
#    in ceiling((101 - h) / 10) folds.
############################################

cat("\n== 2. CV design ==\n")
for (h in horizons) {
  t <- read_res(h, "_tuning")
  exp_scored <- validation_size + 1 - h
  exp_folds <- ceiling(exp_scored / (validation_size / n_folds))
  cat(sprintf("h = %d | expect %d scored forecasts in %d folds\n", h, exp_scored, exp_folds))
  print(data.frame(
    tuning_date = t$tuning_date,
    n_scored = t$n_scored,
    n_folds_used = t$n_folds_used,
    as_expected = t$n_scored == exp_scored & t$n_folds_used == exp_folds
  ), row.names = FALSE)
}


############################################
# 3. CV RMSE surface
#    best:               chosen (alpha, lambda)
#    vs_intercept:       best / heaviest penalty (all predictors out)
#                        below 1 = predictors help in CV
#    vs_alpha_0.05/0.95: best / best at the near-Ridge and
#                        near-LASSO ends of the alpha grid
#                        should be <= 1
############################################

cat("\n== 3. CV RMSE surface ==\n")
chk3 <- do.call(rbind, lapply(horizons, function(h) {
  t <- read_res(h, "_tuning")
  g <- read_res(h, "_cv_grid")
  do.call(rbind, lapply(seq_len(nrow(t)), function(i) {
    sub <- g[g$tuning_date == t$tuning_date[i], ]
    best <- min(sub$cv_rmse)
    heavy <- min(sub$cv_rmse[sub$lambda == max(sub$lambda)])
    lo <- min(sub$cv_rmse[sub$alpha == min(sub$alpha)])
    hi <- min(sub$cv_rmse[sub$alpha == max(sub$alpha)])
    data.frame(
      horizon = h,
      tuning_date = t$tuning_date[i],
      alpha = t$alpha[i],
      lambda = signif(t$lambda[i], 4),
      cv_rmse = best,
      matches_tuning_file = isTRUE(all.equal(best, t$cv_rmse[i])),
      vs_intercept = best / heavy,
      vs_alpha_0.05 = best / lo,
      vs_alpha_0.95 = best / hi
    )
  }))
}))
print(chk3, digits = 4, row.names = FALSE)


############################################
# 4. Grid edges and sparsity
############################################

cat("\n== 4. Grid edges and sparsity ==\n")
for (h in horizons) {
  t <- read_res(h, "_tuning")
  d <- read_res(h, "_details")
  cat(sprintf(
    "h=%2d | lambda at edge: %d of %d | alpha at 0.05 or 0.95: %d | nonzero predictors median %g (range %d to %d) | intercept-only forecasts: %d\n",
    h, sum(t$lambda_at_edge), nrow(t), sum(t$alpha %in% c(0.05, 0.95)),
    median(d$n_nonzero), as.integer(min(d$n_nonzero)), as.integer(max(d$n_nonzero)),
    sum(d$n_nonzero == 0)
  ))
}


############################################
# 5. Accuracy vs simple benchmarks (pre-COVID targets)
#    hist_mean: mean of the same 120 training targets
#    last_obs:  most recent observed target (no-change forecast)
#    rel_*:     RMSE / hist_mean RMSE  (below 1 = beats the mean)
############################################

cat("\n== 5. Pre-COVID RMSE relative to the historical mean ==\n")
chk5 <- do.call(rbind, lapply(horizons, function(h) {
  df <- read_data(h)
  d <- read_res(h, "_details")
  bench <- t(sapply(d$forecast_origin, function(o) {
    tr <- get_training_idx(o, df, window_size)
    c(hist_mean = mean(df$target[tr]), last_obs = df$target[tail(tr, 1)])
  }))
  pre <- d$target_date <= pre_covid_end
  rmse <- function(p) sqrt(mean((d$actual[pre] - p[pre])^2))
  data.frame(
    horizon = h,
    rmse_enet = rmse(d$forecast),
    rmse_hist_mean = rmse(bench[, "hist_mean"]),
    rel_enet = rmse(d$forecast) / rmse(bench[, "hist_mean"]),
    rel_last_obs = rmse(bench[, "last_obs"]) / rmse(bench[, "hist_mean"])
  )
}))
print(chk5, digits = 3, row.names = FALSE)


############################################
# 6. Forecast variability (pre-COVID)
############################################

cat("\n== 6. Forecast variability ==\n")
chk6 <- do.call(rbind, lapply(horizons, function(h) {
  f <- read_res(h)
  d <- read_res(h, "_details")
  pre <- f$target_date <= pre_covid_end
  data.frame(
    horizon = h,
    files_agree = isTRUE(all.equal(f$forecast, d$forecast)) &&
      isTRUE(all.equal(f$actual, d$actual)),
    sd_ratio = sd(f$forecast[pre]) / sd(f$actual[pre]),
    cor_actual = cor(f$forecast[pre], f$actual[pre])
  )
}))
print(chk6, digits = 3, row.names = FALSE)


############################################
# 7. Leakage audit of every training window
############################################

cat("\n== 7. Leakage audit ==\n")
for (h in horizons) {
  df <- read_data(h)
  d <- read_res(h, "_details")
  ok <- sapply(d$forecast_origin, function(o) {
    tr <- get_training_idx(o, df, window_size)
    length(tr) == window_size &&
      max(df$target_date[tr]) <= o &&           # every training target already observed
      max(df$month[tr]) < o &&
      mdiff(max(df$month[tr]), o) == h &&       # newest row is exactly h months old
      mdiff(min(df$month[tr]), max(df$month[tr])) == window_size - 1  # 120 consecutive months
  })
  cat(sprintf("h=%2d | windows passing all checks: %d of %d\n", h, sum(ok), length(ok)))
}


############################################
# 8. Actual vs forecast plot
############################################

plot_path <- file.path(results_dir, "actual_vs_forecast.png")
png(plot_path, width = 1600, height = 1200, res = 150)
old_par <- par(mfrow = c(2, 2), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0))
for (h in horizons) {
  f <- read_res(h)
  plot(f$target_date, f$actual, type = "l", col = "grey50",
       xlab = "target date", ylab = "target", main = paste("h =", h))
  lines(f$target_date, f$forecast, col = "red")
  abline(v = pre_covid_end, lty = 2)
  if (h == horizons[1]) legend("bottomleft", c("actual", "Elastic Net"),
                               col = c("grey50", "red"), lty = 1, bty = "n", cex = 0.8)
}
par(old_par)
dev.off()
cat("\nPlot saved to", plot_path, "\n")