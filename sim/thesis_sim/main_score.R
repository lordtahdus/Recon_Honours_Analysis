# ---------------------------------------------
# This script is for getting the scores for the reran simulations.
# Adapted from tourism/job/results/read_results_cluster.Rmd
# ---------------------------------------------

{library(MASS)
# library(matrixcalc)
library(Matrix)
library(tidyr)
library(fabletools)
library(fable)
library(feasts)
library(tsibble)
library(dplyr)
library(ReconCov)
library(distributional)}

# Helpers -----------------------------------

# flatten out the sublist in `e`
{flatten_e_list <- function(e_list) {
  out <- list()

  for (name in names(e_list)) {
    obj <- e_list[[name]]

    if (is.list(obj) && all(vapply(obj, is.matrix, logical(1)))) {
      ## e.g. name == "mint_shr_pc", obj contains K1 … K20
      for (k in names(obj)) {
        out[[k]] <- obj[[k]]   # “mint_shr_pc_K1”, …
      }
    } else {
      out[[name]] <- obj                               # “base”, “ols”, …
    }
  }
  out                                             # all matrices at top level
}
vectorise_winkler <- function(mu, sd, actual, interval) {
  if (length(mu) == 0 || length(sd) == 0) {
    stop("mu or sd empty")
  }
  sapply(
    rownames(S),
    function(series) winkler_score(
      dist_normal(mu[series], sd[series]),
      actual[series],
      interval
    )
  )
}
vectorise_crps <- function(mu, sd, actual) {
  if (length(mu) == 0 || length(sd) == 0) {
    stop("mu or sd empty")
  }
  sapply(
    rownames(S),
    function(series) CRPS(
      dist_normal(mu[series], sd[series]),
      actual[series]
    )
  )
}

coherent_draws <- function(mu_base, Sigma_base, S, G, m) {
  # mu: vector of base forecast means (all series)
  # Sigma: covariance matrix of base forecast (all series)
  mu_bot <- G %*% mu_base
  Sigma_bot <- G %*% Sigma_base %*% t(G)
  
  b_samples <- MASS::mvrnorm(m, mu_bot, Sigma_bot) # m x p_bottpm
  y_samples <- t(S %*% t(b_samples))           # m x p
  return(y_samples)
}
energy_score_from_samples <- function(samples, actual) {
  m <- nrow(samples)
  term1 <- mean(sqrt(rowSums(
    (sweep(samples, 2, actual, "-"))^2
  )))
  term2 <- mean(sqrt(rowSums(
    (samples[sample(1:m, m, replace = TRUE), ] - samples[sample(1:m, m, replace = TRUE), ])^2
  )))
  return(term1 - 0.5 * term2)
}}



# Load and wrangle data -------------------------------
path <- "sim/thesis_sim/"
file <- "S4_2_1_T300_M500.rds"

result_list <- readRDS(paste0(path, file))
params <- readRDS("sim/thesis_sim/params_S4_2_1.rds")

S <- params$S

M <- length(result_list)
all_model_names <- names(flatten_e_list(result_list[[1]]$e))
model_names <- c(
  all_model_names[!grepl("sv|hcov", all_model_names)],
  "ols_shr"
)
#   row - sim iterations
#   col - series
result_winklers_80 <- result_winklers_95 <- result_crps <- setNames(
  lapply(model_names, function(name) {
    matrix(0, M, length(rownames(S)), dimnames = list(1:M, rownames(S)))
  }),
  model_names
)
model_names <- c(model_names, "base_shr") # add base_shr for energy only
result_energy <- matrix(0, M, length(model_names),
  dimnames = list(1:M, model_names)
)
rep <- 10000

for (i in 1:M) {

  result <- result_list[[i]]
  result$e <- flatten_e_list(result$e)
  resid <- result$resid
  actual <- result$actual[1,]

  # Recover a list of W1_hat for each model
  W1_hats <- list()

  sample_cov <- compute_cov_matrix(resid, T)
  if (any(eigen(sample_cov)$values < 1e-6)) # check PD for mint_sample
    sample_cov <- as.matrix(nearPD(sample_cov)$mat)

  W1_hats[["base"]] <- W1_hats[["ols"]] <- W1_hats[["mint_sample"]] <- sample_cov # sample variance of 1-step-ahead residuals as a proxy


  W1_hats[["mint_shr"]] <- shrinkage_est(resid, result$W_shr[1])$cov
  W1_hats[["mint_n"]] <- novelist_est(resid, result$W_n["delta", 1], result$W_n["lambda", 1])$cov
  rownames(W1_hats[["mint_n"]]) <- colnames(W1_hats[["mint_n"]]) <- rownames(S)
  # pc models
  W1_hats[["mint_shr_pc_K1"]] <- shrinkage_pc_est(resid, K = 1, result$W_shr[2])$cov
  W1_hats[["mint_shr_pc_K2"]] <- shrinkage_pc_est(resid, K = 2, result$W_shr[3])$cov
  W1_hats[["mint_n_pc_K1"]] <- novelist_pc_est(resid, K = 1, result$W_n["delta", 2], result$W_n["lambda", 2])$cov
  W1_hats[["mint_n_pc_K2"]] <- novelist_pc_est(resid, K = 2, result$W_n["delta", 3], result$W_n["lambda", 3])$cov

  # base_shr FOR MULTIVARIATE ONLY
  W1_hats[["base_shr"]] <- W1_hats[["ols_shr"]] <- shrinkage_est(resid)$cov

  base_fc <- actual - result$e$base[1, ] # recover 1-step ahead base forecast
  
  for (model in model_names) {

    if (model %in% c("base", "base_shr")) {
      fc <- base_fc
      sd_fc <- sqrt(diag(W1_hats[["base"]]))
      samples <- MASS::mvrnorm(rep, fc, W1_hats[[model]])

    } else {
      # obtain mapping G
      if (model %in% c("ols", "ols_shr")) {
        G <- solve(t(S) %*% S) %*% t(S)
      } else {
        P <- t(S) %*% solve(W1_hats[[model]])
        G <- solve(P %*% S) %*% P
      }
      # get reconciled forecast and covariance
      if (model == "ols_shr") {
        # ols_shr uses ols reconciled forecast
        fc <- actual - result$e[["ols"]][1, ]
      } else {
        fc <- actual - result$e[[model]][1, ]
      }
      W1_recon <- S %*% G %*% W1_hats[[model]] %*% t(G) %*% t(S)
      sd_fc <- sqrt(diag(W1_recon))

      samples <- coherent_draws(base_fc, W1_hats[[model]], S, G, rep)

    }
    # Get scores
    # Energy
    result_energy[i, model] <- energy_score_from_samples(samples, actual)

    if (model == "base_shr") next # skip winkler and crps for base_shr
    # Winkler
    result_winklers_80[[model]][i, ] <- vectorise_winkler(
      fc, sd_fc, actual, 80
    )
    result_winklers_95[[model]][i, ] <- vectorise_winkler(
      fc, sd_fc, actual, 95
    )
    # CRPS
    result_crps[[model]][i, ] <- vectorise_crps(fc, sd_fc, actual)
  }
  cat("Done", i, "\n")
}


save_obj <- list(
  w80 = result_winklers_80,
  w95 = result_winklers_95,
  crps = result_crps,
  energy = result_energy
)
saveRDS(save_obj, paste0(path, "scores_S4_2_1_T300_M500.rds"))
