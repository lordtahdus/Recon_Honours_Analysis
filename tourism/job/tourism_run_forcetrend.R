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
  
  # condition for top, state-, and purpose-level series
  data_high <- data |>
    filter(Month >= start & Month <= end) |>
    filter(is_aggregated(State) | (is_aggregated(Zone) & is_aggregated(Purpose)))
  data_low <- data |>
    filter(Month >= start & Month <= end) |>
    filter(
      !(is_aggregated(State) | (is_aggregated(Zone) & is_aggregated(Purpose)))
    )

  # constraint: constant + d=1 for high-level series to force trend
  fit_high <- data_high |>
    model(base = ARIMA(
      Nights, 
      order_constraint = (constant == 1) & (d == 1)
    ))

  fit_low <- data_low |>
    model(base = ARIMA(Nights))

  fit <- bind_rows(fit_high, fit_low)

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

  y <- y[, rownames(S)]
  y_hat <- y_hat[, rownames(S)]
  base_fc <- base_fc[, rownames(S)]
  actual <- actual[, rownames(S)]

  # shrinkage estimator
  W_shr <- shrinkage_est(
    y - y_hat
  )
  window_size <- round(dim(y)[1] * 0.7)
  # NOVELIST estimator
  W_n <- novelist_cv(
    y,
    y_hat,
    S,
    window_size = window_size,
    deltas = seq(0, 1, by = 0.05),
    ensure_PD = TRUE,
    message = FALSE
  )

  # Reconcile
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)  
  recon_mint_n <- reconcile_mint(base_fc, S, W_n$cov)
  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix

  # list of out-of-sample errors
  e <- list(
    base_tr = actual - base_fc,
    ols_tr = actual - recon_ols,
    mint_shr_tr = actual - recon_mint_shr,
    mint_n_tr = actual - recon_mint_n
  )
  
  return(list(
    e = e
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
file <- paste0("tourism/job/results_forcetrend/result_", index, ".rds")

saveRDS(result, file = file)
