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

visnights_full <- readRDS("tourism/data/visnights_full.rds")
# only keep State level data
visnights <- visnights_full |>
  filter(is_aggregated(Region) & is_aggregated(Purpose) & is_aggregated(Zone) == TRUE) |> 
  select(Month, State, Nights)

# fabletools -----------------------------
fit <- visnights |> 
  model(arima = ARIMA(Nights))

fit_fable <- fit |>
  reconcile(mint = min_trace(arima, method = "mint_shrink"))

fc_fable <- fit_fable |> 
  forecast(h = "4 months")

# fc <- fc_fable |>
#   filter(.model == "arima")

# map(fc_fable, function(x) x[[distribution_var(x)]])
# fc_dist <- fc[[distribution_var(fc)]]
# as.matrix(invoke(cbind, map(fc_dist, mean)))


# ReconCov -----------------------------

key_data <- key_data(fit)
agg_data <- fabletools:::build_key_data_smat(key_data)
S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L

# rownames
row_names <- key_data |>
  distinct(State) |>
  mutate(
    name = case_when(
      State == "<aggregated>" ~ "Total",
      TRUE ~ as.character(State)
    )
  )
# colnames
colnames(S) <- row_names |> 
  filter(!is_aggregated(State)) |> 
  pull(name)
rownames(S) <- row_names$name |> 
  as.character()

# Get data from fit
res <- residuals(fit, type = "response") |> 
  left_join(row_names, by = c("State")) |> 
  as_tibble() |>
  select(.resid, name, Month) |>
  pivot_wider(names_from = name, values_from = .resid) |> 
  select(-Month) |> 
  as.matrix()

base_fc <- fc_fable |> 
  as_tibble() |> 
  filter(.model == "arima") |>
  left_join(row_names, by = c("State")) |>
  pivot_wider(id_cols = Month, names_from = name, values_from = .mean) |>
  select(-Month) |> 
  as.matrix()

W_shr <- shrinkage_est(res)
recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)
# pivotlonger
fc_ReconCov <- as_tibble(recon_mint_shr) |>
  mutate(h = row_number()) |>
  pivot_longer(-h, names_to = "State", values_to = ".mean")

fc_ReconCov <- fc_ReconCov |>
  mutate(.var = NA_real_)

# using h-step-ahead residual to get variance
for (h_i in 1:4) {

  fit_augment_h <- augment(fit, h = h_i) |> 
    filter(.model == "arima") |> 
    left_join(row_names, by = c("State")) |> 
    select(Month, State, name, .fitted, Nights)
  
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

  # variance of h-step-ahead in-sample residuals
  D_half_h <- diag(sqrt(diag(
    compute_cov_matrix(y_h - y_hat_h, zero_mean = TRUE)
  ))) # diagonal matrix

  # Scaled-correlation
  # convert to correlation then pre-multiply and post-multiply with D_half_h
  W_shr_h <- D_half_h %*% cov2cor(W_shr$cov) %*% D_half_h

  fc_ReconCov[fc_ReconCov$h==h_i, ".var"] <- diag(W_shr_h)
}



fc_compare <- fc_fable |> 
  filter(.model == "mint") |>
  select(State, Month, Nights) |>
  bind_cols(fc_ReconCov |> arrange(State == "Total", State) )

fc_compare |> 
  print(n=40)
