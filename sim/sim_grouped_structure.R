library(MASS)
# library(matrixcalc)
library(Matrix)
library(tidyr)
library(ggplot2)

library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(dplyr)

library(ReconCov)



# Parameters -----------------------------------

# groups <- c(2,2)
groups <- c(4,4,4,4)

T <- 304
h <- 4
Tsplit <- T - h

# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2))
# )
structure <- list(
  groups,
  as.list(seq(1,length(groups))),
  lapply(0:3, \(x) c(1,5,9,13) + x), # this line not work
  list(c(1,2,3,4))
)
S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
)
mat <- matrix(diag(rep(1,4)), 4, 16)
S[2:5, ] <- mat
S %>% plot_heatmap()
order_S <- rownames(S)

# ranges for coefs in VAR
diag_range <- c(-0.6, 0.6)
offdiag_range <- c(-0.4, 0.4)

## VAR(1) block -------------------------
A <- generate_block_diag(
  groups = groups,
  diag_range = diag_range,
  offdiag_range = offdiag_range,
  stationary = TRUE,
)$A
A %>% plot_heatmap()

diag(A[5:16, 1:12]) <- runif(12, offdiag_range[1], offdiag_range[2])
diag(A[1:12, 5:16]) <- runif(12, offdiag_range[1], offdiag_range[2])
diag(A[9:16, 1:8]) <- runif(8, offdiag_range[1], offdiag_range[2])
diag(A[1:8, 9:16]) <- runif(8, offdiag_range[1], offdiag_range[2])
diag(A[13:16, 1:4]) <- runif(4, offdiag_range[1], offdiag_range[2])
diag(A[1:4, 13:16]) <- runif(4, offdiag_range[1], offdiag_range[2])

rho_min <- 0.4
rho_max <- 0.7
rho <- runif(length(groups), rho_min, rho_max)
Sigma <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0.2,
  # delta = 0.15,
  epsilon = (1-max(rho))*0.6,
  # epsilon = 0.15,
  eidim = length(groups)
)
Sigma %>% plot_heatmap()

# Correlations for 2nd attribute
diag(Sigma[5:16, 1:12]) <- rep(runif(3, rho_min, rho_max), each=4)
diag(Sigma[1:12, 5:16]) <- rep(runif(3, rho_min, rho_max), each=4)
diag(Sigma[9:16, 1:8]) <- rep(runif(2, rho_min, rho_max), each=4)
diag(Sigma[1:8, 9:16]) <- rep(runif(2, rho_min, rho_max), each=4)
diag(Sigma[13:16, 1:4]) <- rep(runif(1, rho_min, rho_max), each=4)
diag(Sigma[1:4, 13:16]) <- rep(runif(1, rho_min, rho_max), each=4)

# convert to cov using random sd
Sigma <- convert_cor2cov(
  Sigma,
  stdevs = runif(nrow(Sigma), 1*sqrt(2), 2*sqrt(3))
)
# flip signs
V <- diag(x = sample(
  c(-1,1),
  size = sum(groups), replace = TRUE,
  prob= c(0.5, 0.5)
))
Sigma <- V %*% Sigma %*% V

#enforce symmetric
Sigma <- (Sigma + t(Sigma)) / 2
Sigma <- nearPD(Sigma)$mat %>% as.matrix()
plot_heatmap(Sigma %>% cov2cor(), TRUE)

# Function -----------------------------------
run <- function(A = NULL, Sigma = NULL, message = F) {

  # # # # # #
  # generate bottom-up series and transforming data
  bottom <- simulate_bottom_var(groups, T, intercept = 100, A=A, Sig=Sigma)
  hts_mat <- bottom$Y %*% t(S)

  hts <- hts_mat %>%
    as_tibble() %>%
    mutate(time = seq(1, nrow(hts_mat))) %>%
    # select(time, everything()) %>%
    as_tsibble(index = time) %>%
    pivot_longer(
      cols = -time,
      names_to = "series",
      values_to = "value"
    )

  # # # # # #
  # Fit and base forecasts
  fit <- hts %>%
    filter(time <= Tsplit) %>%
    model(
      arima = ARIMA(value)
    )
  fc <- fit |>
    forecast(h = h) |>
    as_fable(response = "value", distribution = value)
  fc <- fc %>%
    mutate(h = time - Tsplit)

  # Get data
  y <- hts_mat[1:Tsplit, order_S]
  actual <- hts_mat[(Tsplit + 1):T, order_S]

  y_hat <- fit %>%
    augment() %>%
    select(series, .fitted) %>%
    pivot_wider(names_from = series, values_from = .fitted, names_sort = FALSE) %>%
    as_tibble() %>%
    select(-time) %>%
    as.matrix()

  base_fc <- fc %>%
    as_tibble() %>%
    select(series, .mean, time) %>%
    pivot_wider(names_from = series, values_from = .mean) %>%
    select(-time) %>%
    as.matrix()

  y_hat <- y_hat[, order_S]
  base_fc <- base_fc[, order_S]

  # # # # # #
  # Get covariance estimates
  W_shr <- shrinkage_est(
    y - y_hat
  )

  window <- round(Tsplit * 0.7)
  W_n <- novelist_cv(
    y,
    y_hat,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    ensure_PD = TRUE,
    message = message
  )

  # C++ version
  # W_n <- novelist_cv_cpp(
  #   y,
  #   y_hat,
  #   S,
  #   window = window,
  #   deltas = seq(0, 1, by = 0.05),
  #   ensure_PD = TRUE
  # )

  W_n_hstep <- novelist_cv(
    y,
    y_hat,
    S,
    window = window,
    deltas = seq(0, 1, by = 0.05),
    h = h,
    ensure_PD = TRUE,
    message = message
  )

  # # # # # #
  # Reconcile
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)
  recon_mint_n <- reconcile_mint(base_fc, S, W_n$cov)

  recon_mint_n_hstep <- reconcile_mint(base_fc, S, W_n_hstep$cov)

  sample_cov <- compute_cov_matrix(y - y_hat, zero_mean = T)
  if (any(eigen(sample_cov)$values < 1e-10)) {
    cat("Sample covariance for mint_sample is singular, using nearPD\n")
    sample_cov <- nearPD(sample_cov)$mat
  }
  recon_mint_sample <- reconcile_mint(base_fc, S, sample_cov)

  # Sigma_true <- S %*% Sigma %*% t(S)
  # Sigma_true <- nearPD(Sigma_true)$mat # ensure positive-definite
  # recon_mint_true <- reconcile_mint(base_fc, S, Sigma_true)

  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix


  # # # # # #
  # Return
  SSE <- list(
    base = ((actual - base_fc)^2),
    ols = ((actual - recon_ols)^2),
    mint_shr = ((actual - recon_mint_shr)^2),
    mint_n = ((actual - recon_mint_n)^2),
    mint_n_hstep = ((actual - recon_mint_n_hstep)^2),
    mint_sample = ((actual - recon_mint_sample)^2)
    # mint_true = ((actual - recon_mint_true)^2)
  )

  list(
    SSE = SSE,
    W_shr = W_shr$lambda,
    W_n = c(W_n$lambda, W_n$delta),
    W_n_hstep = c(W_n_hstep$lambda, W_n_hstep$delta)
    # W1_hat = compute_cov_matrix(y - y_hat, zero_mean = TRUE)
  )
}




# ---- PARALELL ----------------------------------
library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 1)

M <- 100

# PARALLEL
# res_list <- future_lapply(seq_len(M), function(i) run(), future.seed=TRUE)

# parallel with progress bar
with_progress({
  p <- progressor(along = 1:M)  # auto sets steps = length

  set.seed(2)
  res_list <- future_lapply(
    X = 1:M,
    FUN = function(i) {
      # this guarantees the bar advances even on error
      on.exit(p(sprintf("Sim %d", i)), add = TRUE)

      ## run the simulation but swallow any error
      tryCatch(
        run(A, Sigma, message = FALSE),

        ## put the error into the return value instead of stopping
        error = function(e) {
          structure(list(message = e$message,
                         call    = e$call,
                         sim_id  = i),
                    class = "sim_error")
        }
      )
    },
    future.seed = TRUE
  )
})

# remove any error simulation
cat("Any error in sim:", any(sapply(res_list, inherits, "sim_error")))
res_list <- res_list[!sapply(res_list, inherits, "sim_error")]

model_names <- c("base", "ols", "mint_shr", "mint_n", "mint_n_hstep", "mint_sample")
SSE_cum <- setNames(
  lapply(model_names, function(name) {
    matrix(0, h, length(order_S), dimnames = list(1:h, order_S))
  }),
  model_names
)

# check orders of results
res_list[[1]]$SSE |> names() == SSE_cum |> names()

# wrangle
SSE_cum <- Reduce(function(acc, res) Map(`+`, acc, res$SSE),
                  res_list, init = SSE_cum)
W_shr_store <- sapply(res_list, `[[`, "W_shr")
W_n_store <- t(sapply(res_list, `[[`, "W_n")) ; colnames(W_n_store) <- c("lambda", "delta")
W_n_hstep_store <- t(sapply(res_list, `[[`, "W_n_hstep")) ; colnames(W_n_hstep_store) <- c("lambda", "delta")

MSE <- lapply(SSE_cum, function(mat) mat / M)

plan(sequential) # Reset to sequential





# Save --------------------------------


# Combine all your outputs into a named list
results <- list(
  groups = groups,
  S      = S,
  A      = A,
  Sigma  = Sigma,
  MSE    = MSE,    # or averaged MSE
  W_shr  = W_shr_store,       # can be a list of matrices or one matrix
  W_n    = W_n_store          # same
)
error_list <- purrr::map(res_list, "SSE")
W1_hat_list <- purrr::map(res_list, "W1_hat")

# Save to file
S_string <- paste0("S", sum(groups))
for (i in 2:length(structure)) {
  S_string <- paste0(S_string, "-", length(structure[[i]]))
}
file <- paste0(
  S_string,
  "_T", T-h,
  "_M", M,
  "_sparsebwgrp_hstep"
)
saveRDS(results, file = paste("sim/sim_results/", file, ".rds", sep = ""))

saveRDS(error_list, file = paste("sim/sim_results/errorlist/", file, "_errorlist.rds", sep = ""))

saveRDS(W1_hat_list, file = paste("sim/sim_results/W1hat/", file, "_W1hat.rds", sep = ""))



# Inspect -----------------------------------------------


library(purrr)

## Line plot -------------------

MSE

for (model in names(MSE)) {
  MSE[[model]] <- as.matrix(MSE[[model]])
}
MSE_ts <- transform_sim_MSE(MSE, F)

MSE_ts |> group_by(.model, h) |>
  summarise(mse = mean(MSE)) |>
  ggplot(aes(x = h, y = mse, color = .model)) +
  geom_line() +
  labs(x = "Horizon", y = "MSE") +
  theme_minimal()


MSE$mint_shr - MSE$mint_n

MSE_ts |>
  group_by(series, h) |>
  mutate(
    base_MSE = MSE[.model == "base"]
  ) |>
  ungroup() |>
  group_by(.model, h) |>
  summarise(
    MSE = mean(MSE),
    base_MSE = mean(base_MSE),
    pct_change = (MSE - base_MSE) / base_MSE * 100
  ) |>
  filter(h <=16) |>
  ggplot(aes(x = h, y = pct_change, color = .model)) +
  geom_line() +
  labs(x = "Horizon", y = "% improvements",
       title = "% relative improvements in MSE compared to Base") +
  theme_minimal()


## Box plot -------------------

# transform into df
error_df <- transform_error_list(error_list)

# box plot of 1-step-ahead error2
error_df %>%
  filter(h == 1) %>%  # filter for 1-step-ahead errors
  group_by(.model, id) %>%
  summarise(MSE = mean(e2)) %>%
  ggplot(aes(x = .model, y = MSE, color = .model)) +
  geom_boxplot() +
  theme_minimal()

# box plot of 1-step-ahead relative improvement
error_df %>%
  # filter(h == 1) %>%
  group_by(.model, id) %>%
  summarise(MSE = mean(e2)) %>%
  # calculate relative improvement compared to base model
  group_by(id) %>%
  mutate(base_MSE = MSE[.model == "base"]) %>%
  ungroup() %>%
  mutate(pct_change = (MSE - base_MSE) / base_MSE * 100) %>%
  # remove outliers from mint
  filter(pct_change < 200) %>%
  # plot
  ggplot(aes(x = .model, y = pct_change, color = .model)) +
  geom_boxplot() +
  labs(x = "Model", y = "% relative improvements in MSE",
       title = "% relative improvements in MSE compared to base, 1-to-4-step-ahead forecasts") +
  theme_minimal()






