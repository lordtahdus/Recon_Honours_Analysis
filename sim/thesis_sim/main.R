# ---------------------------------------------
# This script is for re-running the simulations for the thesis.
# Adapted from sim/sim_parallel.R
# ---------------------------------------------

library(MASS)
# library(matrixcalc)
library(Matrix)
library(tidyr)
library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(dplyr)
library(ReconCov)

# Parameters -----------------------------------

# params <- readRDS("sim/thesis_sim/params_S4_2_1.rds")

groups <- params$groups
S <- params$S
A <- params$A
Sigma <- params$Sigma
colnames(Sigma) <- rownames(Sigma) <- colnames(A) <- rownames(A) <- colnames(S)
order_S <- rownames(S)

T <- 54
h <- 4
Tsplit <- T - h

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
    actual = actual, # actual values
    W_shr = W_shr_lambdas,
    W_n = rbind(
      delta = W_n_list$delta,
      lambda = W_n_list$lambda
    )
  )
}

# Run PARALLELL -----------------------------------

library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 1)

M <- 500

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

# Checks
cat("Any error in sim:", any(sapply(res_list, inherits, "sim_error")))
# any NA or 0s or Inf in errors


# # amend the actual values each iter to old result list
# old <- readRDS("sim/thesis_sim/S100_10_3_1_dense_T50_M500.rds")
# for (i in 1:M) {
#   print(all.equal(old[[i]]$resid, res_list[[i]]$resid))
# } # check

# for (i in 1:10) {
#   print(all.equal(temp[[i]]$actual, res_list[[i]]$actual))
# } # check


# for (i in 1:500) {
#   old[[i]]$actual <- res_list[[i]]$actual
# }
# # save
# saveRDS(old, "sim/thesis_sim/S100_10_3_1_dense_T50_M500.rds")


# Run SEQUENTIAL -----------------------------------
# set.seed(1)
# res_list <- lapply(1:M, function(i) {
#   cat("Running sim ", i, "\n")
#   run(A, Sigma, message = TRUE)
# })
