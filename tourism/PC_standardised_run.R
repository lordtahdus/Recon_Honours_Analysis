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



visnights_full <- readRDS("tourism/data/visnights_full.rds")
S <- readRDS("tourism/data/S.rds")
row_names <- readRDS("tourism/data/row_names.rds")
seq_dates <- readRDS("tourism/job/seq_dates.rds")

# Get actual data
y_full <- visnights_full |> 
  as_tibble() |> 
  arrange(State, Region) |>
  left_join(row_names, by = c("State", "Zone", "Region", "Purpose")) |>
  select(name, Nights, Month) |>
  pivot_wider(names_from = name, values_from = Nights) |>
  select(-Month) |>
  as.matrix()
y_full <- y_full[, rownames(S)]

path <- "D:/Github/tourism_cluster/results/"
start_date <- seq_dates$start[1]
h <- 12

run <- function(i) {
  result <- readRDS(paste0(path, "result_", i, ".rds")) 
  resid <- result$resid

  actual <- y_full[(seq_dates$end[i] + 1 - start_date) + 1:h, ] # 1-step ahead actuals
  base_fc <- actual - result$e$base # reconstruct base forecasts

  remove(result) # free memory

  # y <- y_full[i:(seq_dates$end[i] - start_date + 1), ]
  # y_hat <- y - resid
  
  Ks <- c(1, 2, 5, 10)

  # shrinkage estimator with PC adjustment
  W_shr_pc <- lapply(Ks, function(K) {
    shrinkage_pc_est(resid, K = K)
  })
  names(W_shr_pc) <- paste0("K", Ks)
  
  recon_mint_shr_pc <- lapply(W_shr_pc, function(W) {
    reconcile_mint(base_fc, S, W$cov)
  })

  # window <- round(dim(y)[1] * 0.7)
  # # NOVELIST estimator with PC adjustment
  # W_n_list <- novelist_pc_cv(
  #   y,
  #   y_hat,
  #   Ks = Ks,
  #   S,
  #   window = window,
  #   deltas = seq(0, 1, by = 0.05),
  #   standardise = TRUE,
  #   ensure_PD = TRUE,
  #   message = FALSE
  # )

  # recon_mint_n_pc <- lapply(W_n_list$cov, function(W) {
  #   reconcile_mint(base_fc, S, W)
  # })

  e <- lapply(recon_mint_shr_pc, function(x) actual - x) # sublist
  # e <- lapply(recon_mint_n_pc, function(x) actual - x) # sublist

  # lambdas for shrinkage
  W_shr_lambdas <- c(
    sapply(W_shr_pc, function(x) x$lambda)
  )

  return(list(
    e = e,
    W_shr = W_shr_lambdas
    # W_n = rbind(
    #   delta = W_n_list$delta,
    #   lambda = W_n_list$lambda
    # )
  ))
}


library(future.apply)
library(progressr)

handlers(global = TRUE) # Setup progress bar handler
handlers("txtprogressbar")  # or "progress" for a fancier bar

plan(multisession, workers = parallel::detectCores() - 1)

M <- 97
# parallel with progress bar
with_progress({
  p <- progressor(along = 1:M)  # auto sets steps = length
  
  result_list <- future_lapply(
    1:M, 
    function(i) {
      p(message = sprintf("Sim %d", i))  # advances safely
      run(i)
    },
    future.seed=TRUE
  )
})

