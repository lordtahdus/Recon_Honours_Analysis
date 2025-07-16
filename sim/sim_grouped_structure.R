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

rho <- runif(length(groups), 0.5, 0.7)
Sigma <- generate_cor(
  groups = groups,
  rho = rho,
  delta = min(rho) * 0,
  # delta = 0.15,
  epsilon = (1-max(rho))*0.6,
  # epsilon = 0.15,
  eidim = length(groups)
)
Sigma %>% plot_heatmap()


diag(Sigma[5:16, 1:12]) <- runif(12, offdiag_range[1], offdiag_range[2])
diag(Sigma[1:12, 5:16]) <- runif(12, offdiag_range[1], offdiag_range[2])
diag(Sigma[9:16, 1:8]) <- runif(8, offdiag_range[1], offdiag_range[2])
diag(Sigma[1:8, 9:16]) <- runif(8, offdiag_range[1], offdiag_range[2])
diag(Sigma[13:16, 1:4]) <- runif(4, offdiag_range[1], offdiag_range[2])
diag(Sigma[1:4, 13:16]) <- runif(4, offdiag_range[1], offdiag_range[2])
