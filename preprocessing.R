library(tidyverse)

df <- read_csv("data/data_raw.csv")

# 221 - CPIAUCSL
# 222 - GS1
# 223 - GS10
# 224 - M1SL
# 225 - M2SL
# 226 - UNRATE

df <- df %>%
  mutate(month = dmy(month)) %>%
  arrange(month) %>%
  filter(month >= as.Date("1989-10-01") &
           month <= as.Date("2020-12-01"))

colSums(is.na(df))
df %>% count(month) %>% filter(n > 1) # No duplicates and no missing values

#################

target_var <- 'x1'
predictors <- setdiff(names(df), c("month"))

# Lagged predictors
df_lagged <- df %>%
  mutate(
    across(
      all_of(predictors),
      list(
        lag0 = identity,
        lag1 = ~ lag(.x, 1),
        lag2 = ~ lag(.x, 2),
        lag3 = ~ lag(.x, 3)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  select(month, ends_with(c("_lag0", "_lag1", "_lag2", "_lag3")))

# Target variables in different horizons
df_targets <- df %>%
  mutate(
    y_h1  = lead(.data[[target_var]], 1),
    y_h3  = lead(.data[[target_var]], 3),
    y_h6  = lead(.data[[target_var]], 6),
    y_h12 = lead(.data[[target_var]], 12)
  ) %>%
  select(month, starts_with("y_h"))

# Combine predictors and targets
df_model <- df_lagged %>%
  left_join(df_targets, by = "month") %>%
  filter(month <= as.Date("2019-12-01") & month >= as.Date("1990-01-01"))

colSums(is.na(df_model))
df_model %>% count(month) %>% filter(n > 1) # No duplicates and no missing values

#################

# Create datasets for each horizon
horizons <- c(1, 3, 6, 12)

datasets <- map(
  horizons,
  function(h) {
    df_model %>%
      mutate(
        target_date = month %m+% months(h),
        target = .data[[paste0("y_h", h)]]
      ) %>%
      select(
        month,
        target_date,
        starts_with(paste0(target_var, "_lag")),
        ends_with(c("_lag0", "_lag1", "_lag2", "_lag3")),
        target
      ) %>%
      select(-any_of(target_var)) %>%
      select(-matches("^y_h(1|3|6|12)$"))
  }
)

names(datasets) <- paste0("h", horizons)

df_h1  <- datasets$h1
df_h3  <- datasets$h3
df_h6  <- datasets$h6
df_h12 <- datasets$h12

#################

# Save as csv
write_csv(df_h1,  "data/processed/h1.csv")
write_csv(df_h3,  "data/processed/h3.csv")
write_csv(df_h6,  "data/processed/h6.csv")
write_csv(df_h12, "data/processed/h12.csv")
