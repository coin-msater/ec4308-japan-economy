library(readr)

# setup
horizons      <- c(1L, 3L, 6L, 12L)
PRE_COVID_END <- as.Date("2019-12-01")   # pre-COVID cutoff, on target date
TEST_END      <- as.Date("2021-12-01")

OUTPUT_DIR <- "ensemble/ensemble_results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# where each model's forecast files live; {h} is replaced by the horizon
model_files <- list(
  RF         = "RF_modelling/rf_results_block_cv/blockcv_rf_ensemble_h{h}.csv",
  Lasso      = "lasso/results/periodic_10block_aligned/lasso_ensemble_h{h}.csv",
  Ridge      = "ridge/results/periodic_10block_aligned/ridge_ensemble_h{h}.csv",
  ElasticNet = "elastic_net/results/periodic_10block_aligned/elastic_net_ensemble_h{h}.csv"
)

# read one model's forecasts for one horizon
read_model <- function(model, h) {
  path <- sub("{h}", h, model_files[[model]], fixed = TRUE)
  if (!file.exists(path)) stop("Missing file for ", model, ", h=", h, ": ", path)
  d <- as.data.frame(read_csv(path, show_col_types = FALSE))
  need <- c("forecast_origin", "target_date", "actual", "predicted")
  if (!all(need %in% names(d)))
    stop(model, " h=", h, ": file must contain columns ", paste(need, collapse = ", "))
  d <- d[, need]
  d$forecast_origin <- as.Date(d$forecast_origin)
  d$target_date     <- as.Date(d$target_date)
  stopifnot(all(is.finite(d$predicted)), all(is.finite(d$actual)))
  d
}

score <- function(actual, pred) {
  c(RMSE = sqrt(mean((actual - pred)^2)), MAE = mean(abs(actual - pred)))
}

run_horizon <- function(h) {
  models <- names(model_files)
  dl <- setNames(lapply(models, read_model, h = h), models)
  base <- dl[[1]]
  
  # all models must cover exactly the same origins, target dates and actuals
  for (m in models[-1]) {
    d <- dl[[m]]
    if (!identical(as.character(d$forecast_origin), as.character(base$forecast_origin)) ||
        !identical(as.character(d$target_date), as.character(base$target_date)))
      stop("h=", h, ": dates of ", m, " do not match ", models[1],
           " (", nrow(d), " vs ", nrow(base), " rows). ",
           "Make sure every model uses the same test period.")
    if (!isTRUE(all.equal(d$actual, base$actual)))
      stop("h=", h, ": 'actual' of ", m, " does not match ", models[1])
  }
  
  # simple average of the forecasts (equal weights)
  P <- do.call(cbind, lapply(dl, function(d) d$predicted))
  colnames(P) <- models
  ens <- rowMeans(P)
  
  out <- data.frame(horizon = h,
                    forecast_origin = base$forecast_origin,
                    target_date = base$target_date,
                    actual = base$actual,
                    predicted = ens)
  write_csv(out, file.path(OUTPUT_DIR, paste0("ensemble_h", h, ".csv")))
  
  # RMSE / MAE for the ensemble and each individual model
  P_all <- cbind(P, Ensemble = ens)
  periods <- list(pre_covid = PRE_COVID_END, through_dec_2021 = TEST_END)
  do.call(rbind, lapply(names(periods), function(pn) {
    keep <- base$target_date <= periods[[pn]]
    do.call(rbind, lapply(colnames(P_all), function(m) {
      s <- score(base$actual[keep], P_all[keep, m])
      data.frame(horizon = h, period = pn, model = m, test_months = sum(keep),
                 RMSE = s[["RMSE"]], MAE = s[["MAE"]])
    }))
  }))
}

summary_tbl <- do.call(rbind, lapply(horizons, run_horizon))
rownames(summary_tbl) <- NULL
write_csv(summary_tbl, file.path(OUTPUT_DIR, "ensemble_summary.csv"))
print(summary_tbl, row.names = FALSE)
cat("Completed. Files saved under: ", OUTPUT_DIR, "\n", sep = "")