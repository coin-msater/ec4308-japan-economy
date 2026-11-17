# EC4308: A Machine Learning Approach to Economic Forecasting in Japan

This repository contains the data, R code, model outputs, and forecast evaluation for our EC4308 group project on forecasting Japan's industrial production using machine learning methods.

## Project Overview

The project compares the out-of-sample forecasting performance of six approaches:

- Autoregressive model (AR)
- Ridge regression
- LASSO
- Elastic Net
- Random Forest
- Equal-weighted ensemble

Forecasts are produced for four horizons:

- 1 month
- 3 months
- 6 months
- 12 months

The machine learning models use a high-dimensional macroeconomic dataset containing 225 underlying variables. The contemporaneous value and first three lags of each variable are included, giving 900 predictors.

The forecasting target is Japan's industrial production series, `x1`.

## Repository Structure

```text
ec4308-japan-economy/
│
├── data/
│   ├── data_raw.csv
│   ├── processed/
│   └── ADF/
│
├── AR/
│   ├── AR_model.R
│   └── ar_results/
│
├── ridge/
│   ├── ridge_updated.R
│   ├── ridge_evaluation.R
│   └── results/
│
├── lasso/
│   ├── lasso_updated.R
│   ├── lasso_evaluation.R
│   └── results/
│
├── elastic_net/
│   ├── elasticnet_model.R
│   ├── elasticnet_evaluation.R
│   └── results/
│
├── RF_modelling/
│   ├── ec4308_rf.R
│   └── rf_results/
│
├── ensemble/
│   ├── ensemble_model.R
│   └── ensemble_results/
│
├── DM test/
│   ├── dm_tests.R
│   └── dm_results.csv
│
├── preprocessing.R
└── README.md
```

## Data

The raw dataset contains 225 macroeconomic variables.

For each variable, the modelling dataset includes:

```text
x_t, x_(t-1), x_(t-2), x_(t-3)
```

giving:

```text
225 × 4 = 900 predictors
```

The processed horizon-specific datasets are stored under `data/processed/`.

For each forecast horizon `h`, the target is:

```text
target = x1_(t+h)
```

where `h` is 1, 3, 6, or 12 months.

## Forecasting Framework

Models are evaluated using a rolling pseudo-out-of-sample forecasting framework.

- Rolling estimation window: 120 observations
- Forecast horizons: 1, 3, 6, and 12 months
- Initial validation period: 100 observations
- Regularised models use chronological cross-validation for hyperparameter tuning
- AR lag order is selected using BIC
- Random Forest tuning uses out-of-bag error
- The ensemble combines Ridge, LASSO, Elastic Net, and Random Forest forecasts using equal weights

Forecast performance is evaluated primarily using RMSE, with Diebold-Mariano tests used to compare predictive accuracy.

## Running the Project

Start by running:

```text
preprocessing.R
```

to generate the processed horizon-specific datasets.

The individual model scripts can then be run from their respective folders:

```text
AR/
ridge/
lasso/
elastic_net/
RF_modelling/
```

After generating the individual model forecasts, run the ensemble and Diebold-Mariano test scripts for the combined evaluation.

## Authors

EC4308 Group 5, National University of Singapore

- Lai Jing Xuan
- Ong Hai Liang
- Owen Lim Wen Xuan
- Wang Ruoyu
- Yang Xiuli