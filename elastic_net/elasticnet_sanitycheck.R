############################################
# Sanity checks for Elastic Net (and Ridge/LASSO) output
# Run from the project root after all three model scripts.
#
# Works with the standard format (horizon, target_date, forecast, actual)
# and with the older Ridge/LASSO format (... predicted ...).
# Models are matched by target_date, which is unique within a horizon.
############################################

library(readr)

horizons <- c(1, 3, 6, 12)
models <- c("ridge", "lasso", "elastic_net")
window_size <- 120
pre_covid_end <- as.Date("2019-12-01")
expected_n <- length(seq(as.Date("1992-09-01"), as.Date("2021-12-01"), by = "month"))  # 352
standard_cols <- c("horizon", "target_date", "forecast", "actual")

res_path <- function(model, h, suffix = "") {
  sprintf("%s/results/periodic_10block/%s_h%d%s.csv", model, model, h, suffix)
}

read_res <- function(model, h, suffix = "") {
  f <- as.data.frame(read_csv(res_path(model, h, suffix), show_col_types = FALSE))
  if (!"forecast" %in% names(f) && "predicted" %in% names(f)) {
    names(f)[names(f) == "predicted"] <- "forecast"
  }
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
# 1. Structure of the Elastic Net output
############################################

cat("\n== 1. Structure of Elastic Net output ==\n")
for (h in horizons) {
  f <- read_res("elastic_net", h)
  d <- read_res("elastic_net", h, "_details")
  gap <- mdiff(d$forecast_origin, d$target_date)
  cat(sprintf(
    "h=%2d | standard columns: %s | n=%d (expect %d) | NA forecasts=%d | dup targets=%d | target = origin + h: %s | targets %s to %s\n",
    h, identical(names(f), standard_cols), nrow(f), expected_n,
    sum(is.na(f$forecast)), sum(duplicated(f$target_date)), all(gap == h),
    format(min(f$target_date)), format(max(f$target_date))
  ))
}


############################################
# 2. Same CV design across the three models
############################################

cat("\n== 2. Validation forecasts per tuning date (should match across models) ==\n")
for (h in horizons) {
  tl <- lapply(models, function(m) read_res(m, h, "_tuning"))
  tab <- sapply(tl, function(t) t$n_scored)
  colnames(tab) <- models
  rownames(tab) <- format(tl[[1]]$tuning_date)
  folds_same <- all(sapply(tl, function(t) identical(as.numeric(t$n_folds_used),
                                                     as.numeric(tl[[1]]$n_folds_used))))
  cat("h =", h, "| folds identical:", folds_same, "\n")
  print(tab)
}


############################################
# 3. CV RMSE: Elastic Net vs Ridge and LASSO (same folds)
############################################

cv_rmse_of <- function(t) {
  if ("cv_rmse" %in% names(t)) t$cv_rmse else sqrt(t$cv_mse)
}

cat("\n== 3. CV RMSE at each tuning date ==\n")
chk3 <- do.call(rbind, lapply(horizons, function(h) {
  tl <- lapply(models, function(m) read_res(m, h, "_tuning"))
  names(tl) <- models
  data.frame(
    horizon = h,
    tuning_date = tl$elastic_net$tuning_date,
    alpha = tl$elastic_net$alpha,
    ridge = cv_rmse_of(tl$ridge),
    lasso = cv_rmse_of(tl$lasso),
    enet = cv_rmse_of(tl$elastic_net),
    enet_vs_best = cv_rmse_of(tl$elastic_net) / pmin(cv_rmse_of(tl$ridge), cv_rmse_of(tl$lasso))
  )
}))
print(chk3, digits = 4, row.names = FALSE)


############################################
# 4. Grid edges and sparsity
############################################

cat("\n== 4. Grid edges and sparsity ==\n")
for (h in horizons) {
  t <- read_res("elastic_net", h, "_tuning")
  d <- read_res("elastic_net", h, "_details")
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
  base <- read_res("elastic_net", h, "_details")
  bench <- t(sapply(base$forecast_origin, function(o) {
    tr <- get_training_idx(o, df, window_size)
    c(hist_mean = mean(df$target[tr]), last_obs = df$target[tail(tr, 1)])
  }))
  pre <- base$target_date <= pre_covid_end
  rmse <- function(p) sqrt(mean((base$actual[pre] - p[pre])^2))

  out <- data.frame(
    horizon = h,
    rmse_hist_mean = rmse(bench[, "hist_mean"]),
    rel_last_obs = rmse(bench[, "last_obs"]) / rmse(bench[, "hist_mean"])
  )
  for (m in models) {
    f <- read_res(m, h)
    p <- f$forecast[match(base$target_date, f$target_date)]
    out[[paste0("rel_", m)]] <- rmse(p) / out$rmse_hist_mean
  }
  out
}))
print(chk5, digits = 3, row.names = FALSE)


############################################
# 6. Forecast variability and agreement (pre-COVID)
############################################

cat("\n== 6. Variability and agreement ==\n")
chk6 <- do.call(rbind, lapply(horizons, function(h) {
  fs <- lapply(models, function(m) read_res(m, h))
  names(fs) <- models
  base <- fs$elastic_net
  pre <- base$target_date <= pre_covid_end
  pick <- function(m, col) fs[[m]][[col]][match(base$target_date, fs[[m]]$target_date)]
  data.frame(
    horizon = h,
    same_actuals = all(sapply(models, function(m) isTRUE(all.equal(pick(m, "actual"), base$actual)))),
    sd_ratio = sd(base$forecast[pre]) / sd(base$actual[pre]),
    cor_actual = cor(base$forecast[pre], base$actual[pre]),
    cor_with_ridge = cor(base$forecast, pick("ridge", "forecast"), use = "complete.obs"),
    cor_with_lasso = cor(base$forecast, pick("lasso", "forecast"), use = "complete.obs")
  )
}))
print(chk6, digits = 3, row.names = FALSE)


############################################
# 7. Leakage audit of every training window
############################################

cat("\n== 7. Leakage audit ==\n")
for (h in horizons) {
  df <- read_data(h)
  d <- read_res("elastic_net", h, "_details")
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

png("elastic_net/results/periodic_10block/actual_vs_forecast.png",
    width = 1600, height = 1200, res = 150)
old_par <- par(mfrow = c(2, 2), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0))
for (h in horizons) {
  f <- read_res("elastic_net", h)
  plot(f$target_date, f$actual, type = "l", col = "grey50",
       xlab = "target date", ylab = "target", main = paste("h =", h))
  lines(f$target_date, f$forecast, col = "red")
  abline(v = pre_covid_end, lty = 2)
  if (h == horizons[1]) legend("bottomleft", c("actual", "Elastic Net"),
                               col = c("grey50", "red"), lty = 1, bty = "n", cex = 0.8)
}
par(old_par)
dev.off()
cat("\nPlot saved to elastic_net/results/periodic_10block/actual_vs_forecast.png\n")