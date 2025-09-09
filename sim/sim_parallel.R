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

groups <- params$groups
S <- params$S
A <- params$A
Sigma <- params$Sigma
colnames(Sigma) <- rownames(Sigma) <- colnames(A) <- rownames(A) <- colnames(S)
order_S <- rownames(S)

# groups <- c(2,2)
# groups <- c(4,4,4,4)
# groups <- c(6,6,6,6,6,6)
# groups <- c(50,50)

T <- 54
h <- 4
Tsplit <- T - h

structure <- list(
  groups,
  as.list(seq(1,length(groups))),
  list(c(1,2))
)
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   list(c(1,2,3,4))
# )
# structure <- list(
#   groups,
#   as.list(seq(1,length(groups))),
#   # list(c(1,2,3), c(4,5,6)),
#   # list(c(1,2))
#   list(1:6)
# )
structure <- list(
  rep(10,10),
  as.list(1:10),
  list(1:4, 5:6, 7:10),
  list(1:3)
)

S <- construct_S(
  structure = structure,
  sparse = FALSE,
  ascending = FALSE
)
S |> plot_heatmap()
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

plot_heatmap(A, TRUE)
A <- edit(A) ;colnames(A) <- NULL # edit manually

# check stationary
for (block in seq_along(groups)) {
  size <- groups[block]
  index <- seq(sum(groups[1:(block-1)]) + 1, sum(groups[1:block]))
  print(
    any(abs(A[block, block]) > 0.9)
  )
}


## Sigma ---------------------------------
rho <- runif(length(groups), 0.5, 0.7)
Sigma <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0.5,
  # delta = 0.15,
  epsilon = (1-max(rho))*0.9,
  # epsilon = 0.15,
  eidim = length(groups)
)

plot_heatmap(Sigma, TRUE)

Sigma <- edit(Sigma) ; colnames(Sigma) <- NULL # edit manually
Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
any((eigen(Sigma)$values) < 1e-8)

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

plot_heatmap(Sigma |> cov2cor(), TRUE)

### modify a range of values---------
R <- cov2cor(Sigma)
# scale down by x
lower <- 0.66
upper <- 0.69
R[lower < abs(R) & abs(R) < upper] <-
  R[lower < abs(R) & abs(R) < upper] * 0.05   # scale x

plot_heatmap(R)
any((eigen(R)$values) < 1e-8)
R <- nearPD(R)$mat |> as.matrix()
D <- diag(sqrt(diag(Sigma)))
Sigma <- D %*% R %*% D

### formulate as NOVELIST ------------
delta <- 0.5
lambda <- 0.2
R <- cov2cor(Sigma)
R_thresh <- sign(R) * pmax(abs(R) - delta, 0)
diag(R_thresh) <- 1
R_novelist <- lambda * R_thresh + (1 - lambda) * R
D <- diag(sqrt(diag(Sigma)))
Sigma <- D %*% R_novelist %*% D


# temporary save
params <- list(
  groups = groups,
  S      = S,
  A      = A,
  Sigma  = Sigma
)
saveRDS(params, "sim/temp_params.rds")

# Function -----------------------------------
run <- function(A = NULL, Sigma = NULL, message = F) {

  # # # # # #
  # generate bottom-up series and transforming data
  bottom <- simulate_bottom_var(groups, T, intercept = 100, A=A, Sig=Sigma)
  hts_mat <- bottom$Y %*% t(S)

  hts <- hts_mat |>
    as_tibble() |>
    mutate(time = seq(1, nrow(hts_mat))) |>
    # select(time, everything()) |>
    as_tsibble(index = time) |>
    pivot_longer(
      cols = -time,
      names_to = "series",
      values_to = "value"
    )

  # # # # # #
  # Fit and base forecasts
  fit <- hts |>
    filter(time <= Tsplit) |>
    model(
      arima = ARIMA(value)
    )
  fc <- fit |>
    forecast(h = h) |>
    as_fable(response = "value", distribution = value)
  fc <- fc |>
    mutate(h = time - Tsplit)

  # Get data
  y <- hts_mat[1:Tsplit, order_S]
  actual <- hts_mat[(Tsplit + 1):T, order_S]

  y_hat <- fit |>
    augment() |>
    select(series, .fitted) |>
    pivot_wider(names_from = series, values_from = .fitted, names_sort = FALSE) |>
    as_tibble() |>
    select(-time) |>
    as.matrix()

  base_fc <- fc |>
    as_tibble() |>
    select(series, .mean, time) |>
    pivot_wider(names_from = series, values_from = .mean) |>
    select(-time) |>
    as.matrix()

  y_hat <- y_hat[, order_S]
  base_fc <- base_fc[, order_S]

  # # # # # #
  # Get covariance estimates
  Ks <- c(0, 1, 2)

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

  sample_cov <- compute_cov_matrix(y - y_hat, zero_mean = T)
  if (any(eigen(sample_cov)$values < 1e-10)) {
    cat("Sample covariance for mint_sample is singular, using nearPD\n")
    sample_cov <- as.matrix(
      nearPD(sample_cov)$mat
    )
  }

  # # # # # #
  # Reconcile
  recon_ols <- reconcile_mint(base_fc, S, diag(rep(1, nrow(S)))) # identity matrix
  recon_mint_shr <- reconcile_mint(base_fc, S, W_shr$cov)
  recon_mint_n <- reconcile_mint(base_fc, S, W_n)
  # pc versions
  recon_mint_shr_pc <- lapply(W_shr_pc, function(W) {
    reconcile_mint(base_fc, S, W$cov)
  })
  recon_mint_n_pc <- lapply(W_n_list$cov[-1], function(W) {
    reconcile_mint(base_fc, S, W)
  })
  names(recon_mint_shr_pc) <- paste0("mint_shr_pc_", names(recon_mint_shr_pc))
  names(recon_mint_n_pc) <- paste0("mint_n_pc_", names(recon_mint_n_pc))
  # mint sample
  recon_mint_sample <- reconcile_mint(base_fc, S, sample_cov)

  # # # # # #
  # Return vanilla
  e <- list(
    base = (actual - base_fc),
    ols = (actual - recon_ols),
    mint_shr = (actual - recon_mint_shr),
    mint_n = (actual - recon_mint_n),
    mint_sample = (actual - recon_mint_sample),
    
    pc = c( # sublist of all pc versions
      lapply(recon_mint_shr_pc, function(x) actual - x),
      lapply(recon_mint_n_pc, function(x) actual - x)
    )
  )

  # # # 
  # Construct cov from 1:h-step-ahead resid (_sv and _hcov)
  # initialise the error matrices
  e_hresid <- setNames(
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
  e_hresid[["mint_shr_sv"]][1, ] <- e_hresid[["mint_shr_hcov"]][1, ] <- e$mint_shr[1, ]
  e_hresid[["mint_n_sv"]][1, ] <- e_hresid[["mint_n_hcov"]][1, ] <- e$mint_n[1, ]

  for (h_i in 2:h) {
    
    fit_augment_h <- augment(fit, h = h_i) |> 
      select(series, .fitted)

    y_hat_h <- fit_augment_h |>
      pivot_wider(names_from = series, values_from = .fitted, names_sort = FALSE) |>
      as_tibble() |>
      select(-time) |>
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
    W_shr_sv_h <- D_half_h %*% cov2cor(W_shr$cov) %*% D_half_h
    W_n_sv_h <- D_half_h %*% cov2cor(W_n) %*% D_half_h

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
    e_hresid[["mint_shr_sv"]][h_i, ] <- actual_h - recon_shr_sv_h
    e_hresid[["mint_n_sv"]][h_i, ] <- actual_h - recon_n_sv_h
    e_hresid[["mint_shr_hcov"]][h_i, ] <- actual_h - recon_shr_hcov_h
    e_hresid[["mint_n_hcov"]][h_i, ] <- actual_h - recon_n_hcov_h
  }

  # combine but keep e_hresid as a sublist name "hresid"
  e <- c(e, list(hresid = e_hresid))

  # lambdas for shrinkage
  W_shr_lambdas <- c(
    K0 = W_shr$lambda,
    sapply(W_shr_pc, function(x) x$lambda)
  )

  list(
    e = e,
    resid = y - y_hat, # in-sample residuals
    W_shr = W_shr_lambdas,
    W_n = rbind(
      delta = W_n_list$delta,
      lambda = W_n_list$lambda
    )
  )
}


# ---- PARALELL ----------------------------------
library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 1)

M <- 500

# PARALLEL
# res_list <- future_lapply(seq_len(M), function(i) run(), future.seed=TRUE)

# parallel with progress bar
with_progress({
  p <- progressor(along = 1:M)  # auto sets steps = length

  set.seed(1)
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

model_names <- c("base", "ols", "mint_shr", "mint_n", "mint_sample")
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
# W_n_hstep_store <- t(sapply(res_list, `[[`, "W_n_hstep")) ; colnames(W_n_hstep_store) <- c("lambda", "delta")

MSE <- lapply(SSE_cum, function(mat) mat / M)

plan(sequential) # Reset to sequential


# Warning message:
# In sqrt(diag(best$var.coef)) : NaNs produced
#
# This happens in fitting ARIMA, when the sum of AR coefs is very
# close to 1, causing difficulties in computing the standard errors.


## benchmark ---------------------
# W_shr_store <- numeric(M)
# W_n_store   <- matrix(0, M, 2,
#                       dimnames = list(NULL, c("lambda", "delta")))
# run_sequential <- function() {
#   for(i in seq_len(M)) {
#     cat("Iteration ", i, "\n")
#
#     res <- run(A, Sigma)
#     # add up the SSE matrices
#     SSE_cum <- Map(`+`, SSE_cum, res$SSE)
#     # store the parameters
#     W_shr_store[i]         <- res$W_shr
#     W_n_store[i, ]         <- res$W_n
#   }
# }
#
# bench::mark(
#   future_lapply(seq_len(M), function(i) run(A, Sigma), future.seed=TRUE),
#   run_sequential(),
#   iterations = 1,
#   check = FALSE
# )



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
  # W_n_hstep = W_n_hstep_store # same
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
  "_M", M
  # "_sparsebwgrp_hstep"
)
saveRDS(results, file = paste("sim/sim_results/", file, ".rds", sep = ""))

saveRDS(error_list, file = paste("sim/sim_results/errorlist/", file, "_errorlist.rds", sep = ""))

saveRDS(W1_hat_list, file = paste("sim/sim_results/W1hat/", file, "_W1hat.rds", sep = ""))



# Inspect -----------------------------------------------


library(purrr)

## Line plot -------------------

MSE

# for (model in names(MSE)) {
#   MSE[[model]] <- as.matrix(MSE[[model]])
# }
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
error_df |>
  filter(h == 1) |>  # filter for 1-step-ahead errors
  group_by(.model, id) |>
  summarise(MSE = mean(e2)) |>
  ggplot(aes(x = .model, y = MSE, color = .model)) +
    geom_boxplot() +
    theme_minimal()

# box plot of 1-step-ahead relative improvement
error_df |>
  # filter(h == 1) |>
  group_by(.model, id) |>
  summarise(MSE = mean(e2)) |>
  # calculate relative improvement compared to base model
  group_by(id) |>
  mutate(base_MSE = MSE[.model == "base"]) |>
  ungroup() |>
  mutate(pct_change = (MSE - base_MSE) / base_MSE * 100) |>
  # remove outliers from mint
  filter(pct_change < 200) |>
  # plot
  ggplot(aes(x = .model, y = pct_change, color = .model)) +
    geom_boxplot() +
    labs(x = "Model", y = "% relative improvements in MSE",
         title = "% relative improvements in MSE compared to base, 1-to-4-step-ahead forecasts") +
    theme_minimal()



