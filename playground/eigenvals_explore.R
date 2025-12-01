library(dplyr)
library(tidyr)
library(fabletools)
library(fable)
library(tsibble)
library(ReconCov)
library(ggplot2)
library(vctrs)


# ------------------------------------------------
# Tourism data -----------------------------------
# ------------------------------------------------

visnights_full <- readRDS("tourism/data/visnights_full.rds")

# Eigenvalues of raw data
y_mat <- visnights_full |>
  as_tibble() |>
  pivot_wider(id_cols = Month, names_from = c(State, Purpose, Zone, Region), values_from = Nights) |>
  select(-Month) |>
  as.matrix()
bottom_idx <- which(
  grepl("<aggregated>", colnames(y_mat), fixed = TRUE) == FALSE
)
bottom_eigen <- eigen(cor(y_mat[, bottom_idx]))$values
hts_eigen <- eigen(cor(y_mat))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_eigen, 20), color = "Bottom")) +
  geom_line(aes(x = 1:20, y = head(hts_eigen, 20), color = "HTS"))

# Eigenvalues of residuals from fitted model
fit <- readRDS("tourism/data/fit.rds")
# fit <- visnights_full |> 
#   model(base = ARIMA(Nights)) |> 
#   reconcile(mint = min_trace(base, method = "mint_shrink"))

key_data <- key_data(fit)
agg_data <- fabletools:::build_key_data_smat(key_data)
S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L

resid <- fit |>
  residuals() |>
  filter(.model == "base") %>%
  as_tibble() %>%
  select(Month, State, Zone, Region, Purpose, .resid) %>%
  pivot_wider(names_from = c(State, Purpose, Zone, Region), values_from = .resid) %>%
  select(-Month) %>%
  as.matrix()

bottom_idx <- which(!is_aggregated(key_data$State) & !is_aggregated(key_data$Purpose) & !is_aggregated(key_data$Zone) & !is_aggregated(key_data$Region))

bottom_eigen <- eigen(cor(resid[, bottom_idx]))$values
hts_eigen <- eigen(cor(resid))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_eigen, 20), color = "Bottom Residuals")) +
  geom_line(aes(x = 1:20, y = head(hts_eigen, 20), color = "HTS Residuals"))


# ------------------------------------------------
# Controlled settings ---------------------------
# ------------------------------------------------

params <- readRDS("sim/thesis_sim/params_S4_2_1.rds")

## 1. 2-level hierarchy ---------------------
bottom <- simulate_bottom_var(params$groups, 100, intercept = 100, A=params$A, Sig=params$Sigma)$Y
hts_mat <- bottom %*% t(params$S)

# Eigenvalues
bottom_eigen <- eigen(cor(bottom))$values
hts_eigen <- eigen(cor(hts_mat))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_eigen, 20), color = "Bottom")) +
  geom_line(aes(x = 1:20, y = head(hts_eigen, 20), color = "HTS"))

# Eigenvalues of residuals
hts <- hts_mat |>
  as_tibble() |>
  mutate(time = seq(1, nrow(hts_mat))) |>
  as_tsibble(index = time) |>
  pivot_longer(
    cols = -time,
    names_to = "series",
    values_to = "value"
  )
fit <- hts |>
  model(base = ARIMA(value))
resid <- fit |>
  augment() |>
  select(series, .resid) |>
  pivot_wider(names_from = series, values_from = .resid, names_sort = FALSE) |>
  as_tibble() |>
  select(-time) |>
  as.matrix()

bottom_resid <- resid[, colnames(params$S)]
bottom_res_eigen <- eigen(cor(bottom_resid))$values
hts_res_eigen <- eigen(cor(resid))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_res_eigen, 20), color = "Bottom Residuals")) +
  geom_line(aes(x = 1:20, y = head(hts_res_eigen, 20), color = "HTS Residuals"))


## 2. Simulate from a factor structure ---------------------
set.seed(1)

# params <- readRDS("sim/thesis_sim/params_S36_6_1.rds")
params <- readRDS("sim/thesis_sim/params_S100_10_3_1_dense.rds")

n_factors <- 5
n_bottom <- ncol(params$S)
# generate factors from 5 ARIMA
random_arima_sim <- function(n = 100,
                             max_p = 3, max_d = 0, max_q = 3) {
  # random orders
  p <- sample(0:max_p, 1)
  d <- sample(0:max_d, 1)
  q <- sample(0:max_q, 1)
  # random coefficients (kept moderate to reduce instability)
  cap <- 0.7 * 0.9^(p + q)
  ar <- if (p > 0) runif(p, -cap, cap) else NULL
  ma <- if (q > 0) runif(q, -cap, cap) else NULL
  x <- arima.sim(
    model = list(order = c(p, d, q), ar = ar, ma = ma),
    n = n + 50
  )
  as.numeric(tail(x, n))
  # list(order = c(p, d, q), ar = ar, ma = ma, series = x)
}
n <- 10000
factors <- sapply(1:n_factors, function(i) random_arima_sim(n = n))
# matplot(factors, type = "l", main = "Simulated Factors")
loadings <- matrix(rnorm(n_bottom * n_factors), nrow = n_bottom, ncol = n_factors)
# loadings <- loadings/10

innov <- matrix(rnorm(n * n_bottom, sd = 1), nrow = n, ncol = n_bottom)

bottom <- factors %*% t(loadings) + innov
hts_mat <- bottom %*% t(params$S)

# Eigenvalues
bottom_eigen <- eigen(cor(bottom))$values
hts_eigen <- eigen(cor(hts_mat))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_eigen, 20), color = "Bottom")) +
  geom_line(aes(x = 1:20, y = head(hts_eigen, 20), color = "HTS"))

# Eigenvalues of residuals
hts <- hts_mat |>
  as_tibble() |>
  mutate(time = seq(1, nrow(hts_mat))) |>
  as_tsibble(index = time) |>
  pivot_longer(
    cols = -time,
    names_to = "series",
    values_to = "value"
  )
fit <- hts |>
  model(base = ARIMA(value))
resid <- fit |>
  augment() |>
  select(series, .resid) |>
  pivot_wider(names_from = series, values_from = .resid, names_sort = FALSE) |>
  as_tibble() |>
  select(-time) |>
  as.matrix()
bottom_resid <- resid[, colnames(params$S)]
bottom_res_eigen <- eigen(cor(bottom_resid))$values
hts_res_eigen <- eigen(cor(resid))$values
ggplot() +
  geom_line(aes(x = 1:20, y = head(bottom_res_eigen, 20), color = "Bottom Residuals")) +
  geom_line(aes(x = 1:20, y = head(hts_res_eigen, 20), color = "HTS Residuals"))



## 3. Nearly uniform eigenvalues ---------------------
set.seed(1)

S <- readRDS("sim/thesis_sim/params_S100_10_3_1_dense.rds")$S
p <- 100
lambda <- runif(p, min = 1, max = 3)

# random orthogonal Q via QR
A <- matrix(rnorm(p * p), nrow = p, ncol = p)
qrA <- qr(A)
Q <- qr.Q(qrA)

# Sigma = Q Λ Q^T
Sigma <- Q %*% diag(lambda) %*% t(Q)
Sigma

eigen(cov2cor(Sigma))$values |> plot(ylim = c(0, 5))

# eigen(S %*% Sigma %*% t(S))$values[1:20] |> plot()
eigen(cov2cor(S %*% Sigma %*% t(S)))$values[1:114] |> plot(ylim = c(0, 5))


# testing with different S structure
structure <- list(
  rep(5,20),
  as.list(1:20),
  list(1:8, 9:12, 13:20),
  list(1:3)
)

S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
)
S |> plot_heatmap()


# # Simulate the data from Sigma
# set.seed(1)
# n <- 1000
# bottom <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
# hts_mat <- bottom %*% t(S)
# # Eigenvalues
# bottom_eigen <- eigen(cor(bottom))$values
# hts_eigen <- eigen(cor(hts_mat))$values
# ggplot() +
#   geom_point(aes(x = 1:50, y = head(bottom_eigen, 50), color = "Bottom")) +
#   geom_point(aes(x = 1:50, y = head(hts_eigen, 50), color = "HTS"))

# # fit ARIMA and get residuals
# hts <- hts_mat |>
#   as_tibble() |>
#   mutate(time = seq(1, nrow(hts_mat))) |>
#   as_tsibble(index = time) |>
#   pivot_longer(
#     cols = -time,
#     names_to = "series",
#     values_to = "value"
#   )
# fit <- hts |>
#   model(base = ARIMA(value))
# resid <- fit |>
#   augment() |>
#   select(series, .resid) |>
#   pivot_wider(names_from = series, values_from = .resid, names_sort = FALSE) |>
#   as_tibble() |>
#   select(-time) |>
#   as.matrix()
# bottom_resid <- resid[, colnames(S)]
# bottom_res_eigen <- eigen(cor(bottom_resid))$values
# hts_res_eigen <- eigen(cor(resid))$values
# ggplot() +
#   geom_point(aes(x = 1:50, y = head(bottom_res_eigen, 50), color = "Bottom Residuals")) +
#   geom_point(aes(x = 1:50, y = head(hts_res_eigen, 50), color = "HTS Residuals"))

