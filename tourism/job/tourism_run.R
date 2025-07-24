library(Matrix)
library(tidyr)
library(ggplot2)

library(vctrs)
library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(lubridate)
library(dplyr)

library(ReconCov)


# Load

visnights_full <- readRDS("tourism/data/visnights_full.rds")
S <- readRDS("tourism/data/S.rds")
row_names <- readRDS("tourism/data/row_names.rds")
seq_dates <- readRDS("tourism/job/seq_dates.rds")

run <- function(
    data,
    S,
    start, end
) {
  fit <- data |>
    filter(Month >= start & Month <= end) |>
    model(base = ARIMA(Nights))

  # base forecasts
  h <- 12
  fc <- fit %>%
    forecast(h = h)

  # Get data from fit
  fit_augment <- fit %>%
    augment() %>%
    filter(.model == "base") %>%
    left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) %>%
    select(Month, State, Region, Purpose, name, .fitted, Nights)

  y_hat <- fit_augment %>%
    as_tibble() %>%
    select(.fitted, name, Month) %>%
    pivot_wider(names_from = name, values_from = .fitted) %>%
    select(-Month) %>%
    as.matrix()

  y <- fit_augment %>%
    as_tibble() %>%
    select(Nights, name, Month) %>%
    pivot_wider(names_from = name, values_from = Nights) %>%
    select(-Month) %>%
    as.matrix()

  base_fc <- fc %>%
    as_tibble() %>%
    filter(.model == "base") %>%
    left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) %>%
    select(name, .mean, Month) %>%
    pivot_wider(names_from = name, values_from = .mean) %>%
    select(-Month) %>%
    as.matrix()

  actual <- visnights_full %>%
    filter(Month > end & Month <= end + 12) %>%
    as_tibble() %>%
    left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) %>%
    select(name, Nights, Month) %>%
    pivot_wider(names_from = name, values_from = Nights) %>%
    select(-Month) %>%
    as.matrix()

  actual <- actual[, rownames(S)]

  # Cov estimators
  W_shr <- shrinkage_est(
    y - y_hat
  )

  window <- round(dim(y)[1] * 0.7)
  W_n <- novelist_cv(
    y,
    y_hat,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    ensure_PD = TRUE,
    message = TRUE
  )

  # Reconcile
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)
  recon_mint_n <- reconcile_mint(base_fc, S, W_n$cov)

  sample_cov <- compute_cov_matrix(y - y_hat, zero_mean = T)
  if (any(eigen(sample_cov)$values < 1e-8)) {
    # cat("Sample covariance for mint_sample is singular, using nearPD\n")
    sample_cov <- as.matrix(nearPD(sample_cov)$mat)
  }
  recon_mint_sample <- reconcile_mint(base_fc, S, sample_cov)

  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix

  e2 <- list(
    base = ((actual - base_fc)^2),
    ols = ((actual - recon_ols)^2),
    mint_shr = ((actual - recon_mint_shr)^2),
    mint_n = ((actual - recon_mint_n)^2),
    mint_sample = ((actual - recon_mint_sample)^2)
  )

  return(list(
    e2 = e2,
    resid_base = actual - base_fc,
    W_shr = W_shr$lambda,
    W_n = c(W_n$lambda, W_n$delta)
  ))
}

args  <- commandArgs(trailingOnly = TRUE)
index <- as.integer(args[1])

start_month <- seq_dates$start[index]
end_month <- seq_dates$end[index]

result <- run(
  data = visnights_full,
  S = S,
  start = start_month,
  end = end_month
)
file <- paste("tourism/job/results/result_", index, ".rds")

saveRDS(result, file = file)
