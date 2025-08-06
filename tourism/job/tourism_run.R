library(Matrix)
library(tidyr)

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
    select(Month, State, Zone, Region, Purpose, name, .fitted, Nights)

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
  Ks <- c(0, 1, 2, 5, 10, 20)

  # shrinkage estimator
  W_shr <- shrinkage_est(
    y - y_hat
  )
  # shrinkage estimator with PC adjustment
  W_shr_pc <- lapply(Ks[-1], function(K) {
    shrinkage_pc_est(y - y_hat, K = K)
  })
  names(W_shr_pc) <- paste0("K", Ks[-1])


  window <- round(dim(y)[1] * 0.7)
  # NOVELIST estimator with and without PC adjustment
  W_n_list <- novelist_pc_cv(
    y,
    y_hat,
    Ks = Ks,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    ensure_PD = TRUE,
    message = FALSE
  )
  W_n <- W_n_list$cov[[1]]

  # Reconcile
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)  
  recon_mint_n <- reconcile_mint(base_fc, S, W_n)

  recon_mint_shr_pc <- lapply(W_shr_pc, function(W) {
    reconcile_mint(base_fc, S, W$cov)
  })
  recon_mint_n_pc <- lapply(W_n_list$cov[-1], function(W) {
    reconcile_mint(base_fc, S, W)
  })

  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix

  # list of out-of-sample errors
  e <- list(
    base = actual - base_fc,
    ols = actual - recon_ols,
    mint_shr = actual - recon_mint_shr,
    mint_n = actual - recon_mint_n,
    mint_shr_pc = lapply(recon_mint_shr_pc, function(x) actual - x), # sublist
    mint_n_pc = lapply(recon_mint_n_pc, function(x) actual - x) # sublist
  )

  # lambdas for shrinkage
  W_shr_lambdas <- c(
    K0 = W_shr$lambda,
    sapply(W_shr_pc, function(x) x$lambda)
  )
  
  return(list(
    e = e,
    resid = y - y_hat,
    base_fit = fit,
    W_shr = W_shr_lambdas,
    W_n = rbind(
      delta = W_n_list$delta,
      lambda = W_n_list$lambda
    )
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
file <- paste0("tourism/job/results/result_", index, ".rds")

saveRDS(result, file = file)
