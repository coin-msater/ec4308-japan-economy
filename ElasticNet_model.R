# Elastic Net rolling-window forecasts

library(tidyverse)
library(lubridate)
library(glmnet)


# ################################
# # Prof Denis' Elastic Net Code
# ################################
# aa=0.5 #set alpha - useful for experimenting without changing the main code line
# #Same glmnet syntax, but now alpha is between 0 and 1:
# elnet.mod <- glmnet(x[train,], y[train], alpha = aa, lambda = grid)
# plot(elnet.mod) #see how shrinkage goes
# 
# set.seed(1789424)
# cv.out10el = cv.glmnet(x[train,], y[train], alpha = aa, lambda=grid)
# plot(cv.out10el$lambda,cv.out10el$cvm, main="10-fold CV user-defined grid", xlab="Lambda", ylab="CV MSE")
# bestlam10el = cv.out10el$lambda.min
# bestlam10el
# 
# #What is the MSE associated with the value of lambda chosen by cross-validation?
# elnet.pred <- predict(elnet.mod, s = bestlam10el, newx = x[test,])
# MSE(elnet.pred, y[test])


############################################
# Elastic Net: direct forecasts of Japan IP
# CV logic aligned with the Ridge/LASSO team stats
#
# Tuning blocks (100 months each, by forecast ORIGIN):
#   used at Sep 1992: May 1984 - Aug 1992
#   used at Jan 2001: Sep 1992 - Dec 2000
#   used at May 2009: Jan 2001 - Apr 2009
#   used at Sep 2017: May 2009 - Aug 2017
#
# A validation forecast is kept only if:
#   (1) its origin is inside the block
#   (2) its target_date is inside the block (<= block end)
#   (3) it has 120 eligible training rows
#       (month < origin and target_date <= origin)
# Kept forecasts are split into consecutive folds of 10
# (last fold may be smaller).
#
# Every fit (CV and test) uses the latest 120 eligible rows,
# so there is no look-ahead for any horizon.
############################################

library(tidyverse)
library(glmnet)
library(lubridate)


############################################
# 1. Settings
############################################

horizons <- c(1, 3, 6, 12)

window_size <- 120
block_size  <- 100
fold_size   <- 10

test_start    <- as.Date("1992-09-01")
pre_covid_end <- as.Date("2019-12-01")   # by TARGET date
covid_end     <- as.Date("2021-12-01")

# How to end the test sample:
#   "target": keep forecasts whose target_date <= covid_end
#   "origin": keep forecasts whose origin <= covid_end
#             (same as the LASSO script; h = 12 then reaches Dec 2022)
# Use whatever rule the whole team agrees on.
end_by <- "target"

# Elastic Net only: exclude alpha = 0 (Ridge) and alpha = 1 (LASSO)
alpha.grid.default <- seq(0.05, 0.95, by = 0.05)

# Decreasing, because glmnet sorts lambda that way internally.
# Predictions always use s = lambda.grid so columns line up.
lambda.grid.default <- exp(seq(log(20), log(1e-4), length.out = 100))


############################################
# 2. Prepare X, y and dates (same column rule as LASSO script)
############################################

make_xy <- function(df_h) {
  
  df_h <- df_h %>%
    mutate(
      month = as.Date(month),
      target_date = as.Date(target_date)
    ) %>%
    arrange(month)
  
  predictor_names <- setdiff(names(df_h), c("month", "target_date", "target"))
  
  if (length(predictor_names) == 0) {
    stop("No predictor columns found.")
  }
  
  odd <- predictor_names[!grepl("^x", predictor_names)]
  if (length(odd) > 0) {
    warning(
      "Predictor columns not starting with 'x': ",
      paste(odd, collapse = ", "),
      ". Check these are genuine predictors."
    )
  }
  
  X <- as.matrix(df_h[, predictor_names])
  y <- df_h$target
  
  if (anyNA(X) || anyNA(y)) {
    stop("Missing values detected in predictors or target.")
  }
  
  if (anyDuplicated(df_h$month) > 0) {
    stop("Duplicate forecast origins (month) found.")
  }
  
  list(
    X = X,
    y = y,
    month = df_h$month,
    target_date = df_h$target_date,
    predictor_names = predictor_names
  )
}


############################################
# 3. Eligible training rows at an origin
############################################

eligible_rows <- function(month, target_date, origin) {
  which(month < origin & target_date <= origin)
}


############################################
# 4. Build folds for one 100-month block
############################################

make_block_folds <- function(month,
                             target_date,
                             block_start,
                             block_end,
                             train.size = window_size,
                             fold.size = fold_size) {
  
  cand <- which(
    month >= block_start &
      month <= block_end &
      target_date <= block_end
  )
  
  has_window <- vapply(
    cand,
    function(r) length(eligible_rows(month, target_date, month[r])) >= train.size,
    logical(1)
  )
  
  val.rows <- cand[has_window]
  
  if (length(val.rows) == 0) {
    stop("No usable validation forecasts in block starting ", block_start)
  }
  
  fold.id <- ceiling(seq_along(val.rows) / fold.size)
  
  lapply(split(val.rows, fold.id), function(v) {
    fold.origin <- month[v[1]]
    list(
      train = tail(eligible_rows(month, target_date, fold.origin), train.size),
      val = v
    )
  })
}


############################################
# 5. Count-only summary (no model fitting)
#    Compare this output with the Ridge/LASSO stats
############################################

retune_origins_for <- function(last_origin) {
  r <- seq(test_start, by = paste(block_size, "months"), length.out = 20)
  r[r <= last_origin]
}

cv_block_summary <- function(df_h) {
  
  dat <- make_xy(df_h)
  origins <- retune_origins_for(max(dat$month))
  
  map_dfr(origins, function(r) {
    b.start <- r %m-% months(block_size)
    b.end   <- r %m-% months(1)
    folds <- make_block_folds(dat$month, dat$target_date, b.start, b.end)
    tibble(
      used_from = r,
      block = paste(format(b.start, "%b %Y"), "-", format(b.end, "%b %Y")),
      forecasts = sum(sapply(folds, function(f) length(f$val))),
      folds = length(folds)
    )
  })
}


############################################
# 6. Tune Elastic Net on given folds
############################################

tuneENET <- function(X,
                     y,
                     folds,
                     alpha.grid = alpha.grid.default,
                     lambda.grid = lambda.grid.default) {
  
  best.mse <- Inf
  best.alpha <- NA_real_
  best.lambda <- NA_real_
  
  tuning.results <- tibble(
    alpha = numeric(),
    lambda = numeric(),
    cv_mse = numeric()
  )
  
  for (aa in alpha.grid) {
    
    sse <- rep(0, length(lambda.grid))
    n.validation <- 0
    
    for (f in folds) {
      
      fit <- glmnet(
        x = X[f$train, , drop = FALSE],
        y = y[f$train],
        alpha = aa,
        lambda = lambda.grid,
        standardize = TRUE
      )
      
      pred <- predict(
        fit,
        newx = X[f$val, , drop = FALSE],
        s = lambda.grid
      )
      
      sse <- sse + colSums((pred - y[f$val])^2)
      n.validation <- n.validation + length(f$val)
    }
    
    cv.mse <- sse / n.validation
    idx <- which.min(cv.mse)
    
    tuning.results <- bind_rows(
      tuning.results,
      tibble(alpha = aa, lambda = lambda.grid[idx], cv_mse = cv.mse[idx])
    )
    
    if (cv.mse[idx] < best.mse) {
      best.mse <- cv.mse[idx]
      best.alpha <- aa
      best.lambda <- lambda.grid[idx]
    }
  }
  
  if (is.na(best.alpha) || is.na(best.lambda)) {
    stop("Unable to select Elastic Net hyperparameters.")
  }
  
  if (best.lambda %in% range(lambda.grid)) {
    warning(
      "Selected lambda ", signif(best.lambda, 4),
      " is at the edge of the grid. Consider widening it."
    )
  }
  
  list(
    alpha = best.alpha,
    lambda = best.lambda,
    cv_mse = best.mse,
    tuning_results = tuning.results
  )
}


############################################
# 7. Count nonzero coefficients
############################################

count_nonzero <- function(coef_matrix) {
  sum(coef_matrix[-1, 1] != 0)
}


############################################
# 8. Rolling-window forecasts
############################################

enet.rolling.window <- function(df_h,
                                alpha.grid = alpha.grid.default,
                                lambda.grid = lambda.grid.default) {
  
  dat <- make_xy(df_h)
  
  X <- dat$X
  y <- dat$y
  month <- dat$month
  target_date <- dat$target_date
  
  
  ############################################
  # Test rows
  ############################################
  
  if (end_by == "target") {
    test.rows <- which(month >= test_start & target_date <= covid_end)
  } else {
    test.rows <- which(month >= test_start & month <= covid_end)
  }
  
  ntest <- length(test.rows)
  
  if (ntest == 0) {
    stop("No observations found in the test period.")
  }
  
  retune.origins <- retune_origins_for(max(month[test.rows]))
  
  
  ############################################
  # Storage
  ############################################
  
  save.pred    <- rep(NA_real_, ntest)
  save.alpha   <- rep(NA_real_, ntest)
  save.lambda  <- rep(NA_real_, ntest)
  save.cv.mse  <- rep(NA_real_, ntest)
  save.nonzero <- rep(NA_integer_, ntest)
  save.coef    <- vector("list", ntest)
  
  retuning.history <- tibble(
    forecast_number = integer(),
    forecast_origin = as.Date(character()),
    block = character(),
    val_forecasts = integer(),
    n_folds = integer(),
    alpha = numeric(),
    lambda = numeric(),
    cv_mse = numeric()
  )
  
  tuning.history <- list()
  
  current.alpha <- NA_real_
  current.lambda <- NA_real_
  current.cv.mse <- NA_real_
  
  
  ############################################
  # Loop
  ############################################
  
  for (i in seq_along(test.rows)) {
    
    forecast.row <- test.rows[i]
    forecast.origin <- month[forecast.row]
    
    
    ############################################
    # Retune at Sep 1992, Jan 2001, May 2009, Sep 2017
    ############################################
    
    if (forecast.origin %in% retune.origins) {
      
      b.start <- forecast.origin %m-% months(block_size)
      b.end   <- forecast.origin %m-% months(1)
      
      folds <- make_block_folds(month, target_date, b.start, b.end)
      
      tuning <- tuneENET(X, y, folds, alpha.grid, lambda.grid)
      
      current.alpha  <- tuning$alpha
      current.lambda <- tuning$lambda
      current.cv.mse <- tuning$cv_mse
      
      n.val <- sum(sapply(folds, function(f) length(f$val)))
      block.label <- paste(format(b.start, "%b %Y"), "-", format(b.end, "%b %Y"))
      
      retuning.history <- bind_rows(
        retuning.history,
        tibble(
          forecast_number = i,
          forecast_origin = forecast.origin,
          block = block.label,
          val_forecasts = n.val,
          n_folds = length(folds),
          alpha = current.alpha,
          lambda = current.lambda,
          cv_mse = current.cv.mse
        )
      )
      
      tuning.history[[length(tuning.history) + 1]] <-
        tuning$tuning_results %>%
        mutate(forecast_number = i, forecast_origin = forecast.origin)
      
      cat(
        "\nRetune at", as.character(forecast.origin),
        "| block:", block.label,
        "| forecasts:", n.val,
        "| folds:", length(folds),
        "| alpha:", round(current.alpha, 3),
        "| lambda:", signif(current.lambda, 4),
        "\n"
      )
    }
    
    if (is.na(current.alpha)) {
      stop("No hyperparameters available at ", forecast.origin)
    }
    
    
    ############################################
    # Latest 120 eligible rows
    ############################################
    
    train.rows <- tail(
      eligible_rows(month, target_date, forecast.origin),
      window_size
    )
    
    if (length(train.rows) < window_size) {
      stop("Fewer than 120 eligible rows at ", forecast.origin)
    }
    
    final.model <- glmnet(
      x = X[train.rows, , drop = FALSE],
      y = y[train.rows],
      alpha = current.alpha,
      lambda = lambda.grid,
      standardize = TRUE
    )
    
    X.out <- X[forecast.row, , drop = FALSE]
    
    pred <- predict(final.model, newx = X.out, s = current.lambda)
    coef.current <- as.matrix(coef(final.model, s = current.lambda))
    
    save.pred[i]    <- as.numeric(pred)
    save.alpha[i]   <- current.alpha
    save.lambda[i]  <- current.lambda
    save.cv.mse[i]  <- current.cv.mse
    save.coef[[i]]  <- coef.current
    save.nonzero[i] <- count_nonzero(coef.current)
    
    if (i %% 50 == 0) {
      cat("Completed", i, "of", ntest, "forecasts\n")
    }
  }
  
  
  ############################################
  # Evaluation (pre-COVID split by TARGET date)
  ############################################
  
  actual <- y[test.rows]
  test.targets <- target_date[test.rows]
  
  pre.idx <- which(test.targets <= pre_covid_end)
  
  errors <- c(
    rmse = sqrt(mean((actual - save.pred)^2)),
    mae = mean(abs(actual - save.pred)),
    n = ntest,
    rmse.pre = sqrt(mean((actual[pre.idx] - save.pred[pre.idx])^2)),
    mae.pre = mean(abs(actual[pre.idx] - save.pred[pre.idx])),
    n.pre = length(pre.idx)
  )
  
  forecasts <- tibble(
    forecast_origin = month[test.rows],
    target_date = test.targets,
    precovid = test.targets <= pre_covid_end,
    actual = actual,
    predicted = save.pred,
    alpha = save.alpha,
    lambda = save.lambda,
    cv_mse = save.cv.mse,
    nonzero_predictors = save.nonzero
  )
  
  list(
    forecasts = forecasts,
    nonzero = save.nonzero,
    coef = save.coef,
    errors = errors,
    retuning_history = retuning.history,
    tuning_history = bind_rows(tuning.history)
  )
}


############################################
# 9. Load data and check block counts first
############################################

dfs <- setNames(
  lapply(horizons, function(h) {
    read_csv(paste0("data/processed/h", h, ".csv"), show_col_types = FALSE)
  }),
  paste0("h", horizons)
)

# Should match the Ridge/LASSO stats
for (nm in names(dfs)) {
  cat("\n==", nm, "==\n")
  print(cv_block_summary(dfs[[nm]]))
}


############################################
# 10. Run all horizons
############################################

enet_all <- lapply(dfs, enet.rolling.window)

enet_errors <- do.call(rbind, lapply(enet_all, `[[`, "errors"))
print(enet_errors)

lapply(enet_all, `[[`, "retuning_history")

# Nonzero predictors right after each retune
lapply(enet_all, function(res) {
  res$nonzero[res$retuning_history$forecast_number]
})

# Intercept-only forecasts
sapply(enet_all, function(res) sum(res$nonzero == 0))


############################################
# 11. Prediction variability
############################################

prediction_diagnostics <- tibble(
  horizon = names(enet_all),
  actual_sd = sapply(enet_all, function(r) sd(r$forecasts$actual)),
  forecast_sd = sapply(enet_all, function(r) sd(r$forecasts$predicted)),
  correlation = sapply(enet_all, function(r) {
    cor(r$forecasts$actual, r$forecasts$predicted)
  })
)

print(prediction_diagnostics)


############################################
# 12. Save results
############################################

dir.create("elastic_net/results", recursive = TRUE, showWarnings = FALSE)

for (nm in names(enet_all)) {
  write_csv(
    enet_all[[nm]]$forecasts,
    paste0("elastic_net/results/elastic_net_", nm, ".csv")
  )
}

stack_results <- function(field) {
  bind_rows(lapply(names(enet_all), function(nm) {
    enet_all[[nm]][[field]] %>% mutate(horizon = as.integer(sub("h", "", nm)))
  }))
}

write_csv(stack_results("retuning_history"), "elastic_net/results/elastic_net_retuning.csv")
write_csv(stack_results("tuning_history"),   "elastic_net/results/elastic_net_tuning_grid.csv")
write_csv(stack_results("forecasts"),        "elastic_net/results/elastic_net_predictions.csv")
write_csv(prediction_diagnostics,            "elastic_net/results/elastic_net_diagnostics.csv")