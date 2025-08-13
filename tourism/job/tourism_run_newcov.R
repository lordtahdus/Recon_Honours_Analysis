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

args  <- commandArgs(trailingOnly = TRUE)
index <- as.integer(args[1])


# New Covs
# 1. _sv Scale Variance:
#     Cov construct by 1-step-ahead Cor and pre- and post-multiply with variances of h-step-ahead
# 2. _hcov Direct h-step Covariance:
#     Cov construct directly from h-step-ahead residuals

path <- "tourism/job/results/"
result <- readRDS(paste0(path, "result_", index, ".rds"))
h <- 12

# START ------------------------------------------
fit <- result$base_fit

fc <- fit |>
  forecast(h = h)

# Get 1-step-ahead fitted values and residuals
fit_augment_1 <- augment(fit, h = 1) |> 
  filter(.model == "base") |> 
  left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) |> 
  select(Month, State, Zone, Region, Purpose, name, .fitted, Nights)

y_hat_1 <- fit_augment_1 |>
  as_tibble() |>
  select(.fitted, name, Month) |>
  pivot_wider(names_from = name, values_from = .fitted) |>
  select(-Month) |>
  as.matrix()

y <- fit_augment_1 |>
  as_tibble() |>
  select(Nights, name, Month) |>
  pivot_wider(names_from = name, values_from = Nights) |>
  select(-Month) |>
  as.matrix()

base_fc <- fc |>
  as_tibble() |>
  filter(.model == "base") |>
  left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) |>
  select(name, .mean, Month) |>
  pivot_wider(names_from = name, values_from = .mean) |>
  select(-Month) |>
  as.matrix()

end <- seq_dates$end[index]
actual <- visnights_full |>
  filter(Month > end & Month <= end + 12) |>
  as_tibble() |>
  left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) |>
  select(name, Nights, Month) |>
  pivot_wider(names_from = name, values_from = Nights) |>
  select(-Month) |>
  as.matrix()

actual <- actual[, rownames(S)]

# get original shrinkage and novelist covariance matrices
shr_lambda <- result$W_shr[[1]]
n_lambda <- result$W_n["lambda", 1]
n_delta <- result$W_n["delta", 1]

W_shr_sv_1 <- shrinkage_est(
  y - y_hat_1,
  lambda = shr_lambda,
  zero_mean = TRUE
)
W_n_sv_1 <- novelist_est(
  y - y_hat_1,
  delta = n_delta,
  lambda = n_lambda,
  zero_mean = TRUE,
  ensure_PD = TRUE
)

# initialise the error matrices
e <- setNames(
  lapply(1:4, function(x) {
    matrix(0, h, length(rownames(S)), dimnames = list(1:h, rownames(S)))
  }),
  c(
    "mint_shr_sv",
    "mint_n_sv",
    "mint_shr_hcov",
    "mint_n_hcov"
  )
)

# 1-step-ahead error is the same with original for all 4 methods
e[["mint_shr_sv"]][1, ] <- e[["mint_shr_hcov"]][1, ] <- result$e$mint_shr[1, ]
e[["mint_n_sv"]][1, ] <- e[["mint_n_hcov"]][1, ] <- result$e$mint_n[1, ]

for (h_i in 2:h) {

  fit_augment_h <- augment(fit, h = h_i) |> 
    filter(.model == "base") |> 
    left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) |> 
    select(Month, State, Zone, Region, Purpose, name, .fitted, Nights)
  
  y_hat_h <- fit_augment_h |>
    as_tibble() |>
    select(.fitted, name, Month) |>
    pivot_wider(names_from = name, values_from = .fitted) |>
    select(-Month) |>
    as.matrix()

  # max number of NAs in each column
  na_length <- max(colSums(is.na(y_hat_h)))

  y_h <- y[(na_length + 1):nrow(y), , drop = FALSE] # remove the first `na_length` rows
  y_hat_h <- y_hat_h[(na_length + 1):nrow(y_hat_h), , drop = FALSE] # remove the first `na_length` rows

  base_fc_h <- base_fc[h_i, , drop = FALSE] # a row vector of the h-step-ahead forecasts
  actual_h <- actual[h_i, , drop = FALSE]

  # variance of h-step-ahead in-sample residuals
  D_half_h <- diag(sqrt(diag(
    compute_cov_matrix(y_h - y_hat_h, zero_mean = TRUE)
  ))) # diagonal matrix

  # 1) Scaled-correlation --------------------------

  # convert to correlation then pre-multiply and post-multiply with D_half_h
  W_shr_sv_h <- D_half_h %*% cov2cor(W_shr_sv_1$cov) %*% D_half_h
  W_n_sv_h <- D_half_h %*% cov2cor(W_n_sv_1$cov) %*% D_half_h

  # 2) Direct h-step covariance ----------------------
  W_shr_hcov_h <- shrinkage_est(
    y_h - y_hat_h,
    zero_mean = TRUE
  )
  
  window <- round(nrow(y_hat_h) * 0.7)
  # NOVELIST estimator with and without PC adjustment
  W_n_hcov_h <- novelist_cv(
    y_h, y_hat_h,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    h = h_i,
    ensure_PD = TRUE, message = FALSE
  )

  # reconcile
  recon_shr_sv_h <- reconcile_mint(base_fc_h, S, W_shr_sv_h)
  recon_n_sv_h <- reconcile_mint(base_fc_h, S, W_n_sv_h)
  recon_shr_hcov_h <- reconcile_mint(base_fc_h, S, W_shr_hcov_h$cov)
  recon_n_hcov_h <- reconcile_mint(base_fc_h, S, W_n_hcov_h$cov)

  # error 
  e[["mint_shr_sv"]][h_i, ] <- actual_h - recon_shr_sv_h
  e[["mint_n_sv"]][h_i, ] <- actual_h - recon_n_sv_h
  e[["mint_shr_hcov"]][h_i, ] <- actual_h - recon_shr_hcov_h
  e[["mint_n_hcov"]][h_i, ] <- actual_h - recon_n_hcov_h
}

# END ------------------------------------------

# append to error list in result
result$e <- c(result$e, e)

file <- paste0(path, "result_", index, ".rds")
saveRDS(result, file = file)