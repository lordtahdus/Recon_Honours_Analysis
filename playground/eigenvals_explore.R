library(dplyr)
library(tidyr)
library(fabletools)
library(fable)
library(tsibble)
library(ReconCov)
library(ggplot2)
library(vctrs)

structure <- list(
  rep(10,10),
  as.list(1:10),
  list(1:3, 4:7, 8:10),
  list(1:3)
)
S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
)

# ------------------------------------------------
# Tourism data -----------------------------------
# ------------------------------------------------

visnights_full <- readRDS("tourism/data/visnights_full.rds")

# Eigenvalues of raw data
y_mat <- visnights_full |>
  filter(is_aggregated(Purpose)) |> 
  as_tibble() |>
  # pivot_wider(id_cols = Month, names_from = c(State, Purpose, Zone, Region), values_from = Nights) |>
  pivot_wider(id_cols = Month, names_from = c(State, Zone, Region), values_from = Nights) |>
  select(-Month) |>
  as.matrix()
bottom_idx <- which(
  grepl("<aggregated>", colnames(y_mat), fixed = TRUE) == FALSE
)
bottom_eigen <- eigen(cor(y_mat[, bottom_idx]))$values
hts_eigen <- eigen(cor(y_mat))$values

ggplot() +
  geom_point(aes(x = 1:length(bottom_eigen), y = bottom_eigen, color = "Bottom")) +
  geom_point(aes(x = 1:length(hts_eigen), y = hts_eigen, color = "HTS"))

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
  filter(is_aggregated(Purpose)) |> # remove PUrpose 
  select(Month, State, Zone, Region, .resid) %>%
  pivot_wider(names_from = c(State, Zone, Region), values_from = .resid) %>%
  select(-Month) %>%
  as.matrix()

bottom_idx <- which(
  grepl("<aggregated>", colnames(resid), fixed = TRUE) == FALSE
)
bottom_eigen <- eigen(cor(resid[, bottom_idx]))$values
hts_eigen <- eigen(cor(resid))$values
ggplot() +
  geom_point(aes(x = 1:length(bottom_eigen), y = bottom_eigen, color = "Bottom Residuals")) +
  geom_point(aes(x = 1:length(hts_eigen), y = hts_eigen, color = "HTS Residuals"))

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



## 4. Simulate data from Identity and Explore effects of Shrinkage ---------------------

set.seed(1)
sigma <- diag(rep(1, 100))
n <- 1000

{
bottom <- MASS::mvrnorm(n = n, mu = rep(0, ncol(sigma)), Sigma = sigma)
hts_mat <- bottom %*% t(S)

# shrinkage
bottom_cor <- cor(bottom)
bottom_corS <- shrinkage_est(bottom, zero_mean = FALSE)$cov |> cov2cor()
hts_cor <- cor(hts_mat)
hts_corS <- shrinkage_est(hts_mat, zero_mean = FALSE)$cov |> cov2cor()

# Eigenvalues
bottom_eigen <- eigen(bottom_cor)$values
bottom_eigenS <- eigen(bottom_corS)$values
hts_eigen <- eigen(hts_cor)$values
hts_eigenS <- eigen(hts_corS)$values

# points plot
m <- 100
ylim <- c(-0.01, 8)
print(
  ggplot() +
    geom_point(aes(x = 1:m, y = head(bottom_eigen, m), color = "Bottom Sample")) +
    geom_point(aes(x = 1:m, y = head(bottom_eigenS, m), color = "Bottom Shrinkage")) +
    # geom_point(aes(x = 1:m, y = head(hts_eigen, m), color = "HTS")) +
    # geom_point(aes(x = 1:m, y = head(hts_eigenS, m), color = "HTS Shrinkage"))
    ylim(ylim)
)

true_hts_eigen <- eigen(cov2cor(S %*% sigma %*% t(S)))$values
ylim <- c(-0.01, 9)
ggplot() +
  geom_point(aes(x = 1:m, y = head(true_hts_eigen, m), color = "True HTS")) +
  geom_point(aes(x = 1:m, y = head(hts_eigen, m), color = "HTS Sample")) +
  geom_point(aes(x = 1:m, y = head(hts_eigenS, m), color = "HTS Shrinkage")) +
  ylim(ylim)
}

## 5. Mimic hase forecasts ---------------------

set.seed(1)

### Cov: Slight uniform eigenvalues ---------------------
p <- 100
lambda <- runif(p, min = 1, max = 3)
# random orthogonal Q via QR
A <- matrix(rnorm(p * p), nrow = p, ncol = p)
qrA <- qr(A)
Q <- qr.Q(qrA)
# Sigma = Q Λ Q^T
Sigma <- Q %*% diag(lambda) %*% t(Q)

eigen(cov2cor(Sigma))$values |> plot(ylim = c(0, 3))

### Cov: factor structure ---------------------
n_factors <- 5
loadings <- matrix(rnorm(p * n_factors, sd = 0.3), nrow = p, ncol = n_factors)
Sigma <- loadings %*% t(loadings) + diag(rep(1, p))
eigen(cov2cor(Sigma))$values |> plot()


### Simulate data ---------------------

W <- S %*% Sigma %*% t(S)
# plot_heatmap((abs(W))^(1/5))
# add some noise to the covariance
noise <- matrix(rnorm(nrow(W)*ncol(W), sd = 0.05) |> abs(), nrow = nrow(W))
noise[upper.tri(noise)] <- t(noise)[upper.tri(noise)]
W <- W + noise + diag(rep(0.3, nrow(W))) 
W_eigen <- eigen(cov2cor(W))$values
plot(W_eigen)

n <- 1000

{
y <- MASS::mvrnorm(n = n, mu = rep(0, ncol(W)), Sigma = W)
W_sample <- cor(y)
W_sample_eigen <- eigen(W_sample)$values

W_S <- shrinkage_est(y, zero_mean = FALSE)$cov |> cov2cor()
W_S_eigen <- eigen(W_S)$values

m <- nrow(W)
ggplot() +
  geom_point(aes(x = 1:m, y = head(W_eigen, m), color = "True Covariance")) +
  geom_point(aes(x = 1:m, y = head(W_sample_eigen, m), color = "Sample Covariance")) +
  geom_point(aes(x = 1:m, y = head(W_S_eigen, m), color = "Shrinkage Covariance")) +
  # add geom_line for each
  geom_line(aes(x = 1:m, y = head(W_eigen, m), color = "True Covariance")) +
  geom_line(aes(x = 1:m, y = head(W_sample_eigen, m), color = "Sample Covariance")) +
  geom_line(aes(x = 1:m, y = head(W_S_eigen, m), color = "Shrinkage Covariance")) +
  ylim(c(-0.1, 13))
}

