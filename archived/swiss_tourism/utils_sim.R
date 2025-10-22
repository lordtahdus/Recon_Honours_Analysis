dir_path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

source("utils.R")

initialize_results_sim = function(N, N_exp, methods, type, specific_method,
                              confidence_levels, quantiles) {
  performance_measures = c("Fc", "MSE", "Extremes_Conf_Int_L", 
                           "Extremes_Conf_Int_R", "MIS", 
                           "COVERAGE", "QUANTILES_SCORE")
  results = list()
  
  for (measure in performance_measures) {
    results[[measure]] = list()
    for (met in methods) {
      if (met == specific_method) {
        results[[measure]][[met]] = list()
          if (measure == "MIS" | measure == "COVERAGE"| measure == "Extremes_Conf_Int_L"|measure == "Extremes_Conf_Int_R"){
            results[[measure]][[met]][[type]] = array(0, dim = c(N, 
                                                                 length(confidence_levels), 
                                                                 N_exp))
          }else if (measure == "QUANTILES_SCORE"){
            results[[measure]][[met]][[type]] = array(0, dim = c(N, 
                                                                 length(quantiles), 
                                                                 N_exp))
          }else{
            results[[measure]][[met]][[type]] = array(0, dim = c(N, N_exp))
          }
        
      } else if (met == "MinT"){
        results[[measure]][[met]] = list()
        for (shrink in c("shrink")) {
          if (measure == "MIS" | measure == "COVERAGE"|measure == "Extremes_Conf_Int_L"|measure == "Extremes_Conf_Int_R"){
            results[[measure]][[met]][[shrink]] = array(0, dim = c(N, 
                                                                   length(confidence_levels), 
                                                         N_exp))
          }else if (measure == "QUANTILES_SCORE"){
            results[[measure]][[met]][[shrink]] = array(0, dim = c(N, 
                                                                   length(quantiles), 
                                                         N_exp))
          }else{
            results[[measure]][[met]][[shrink]] = array(0, dim = c(N, N_exp))
          }
      }
      }else {
        if (measure == "MIS" | measure == "COVERAGE"| measure == "Extremes_Conf_Int_L"|measure == "Extremes_Conf_Int_R"){
          results[[measure]][[met]] = array(0, dim = c(N, 
                                                       length(confidence_levels), 
                                                       N_exp))
        }else if (measure == "QUANTILES_SCORE"){
          results[[measure]][[met]] = array(0, dim = c(N, 
                                                       length(quantiles), 
                                                       N_exp))
        }else{
          results[[measure]][[met]] = array(0, dim = c(N, N_exp))
        }
        
      }
    }
  }
  
  return(results)
}

incoherenceVSreconciledVariance = function(posterior.psi, posterior.nu, 
                                           base_forecasts.mu, A){
  k = nrow(A)
  m = ncol(A)
  n = k + m
  Psi_u = posterior.psi[1:k, 1:k]
  Psi_b = posterior.psi[(k + 1):n, (k + 1):n]
  Psi_ub = posterior.psi[1:k, (k + 1):n, drop = FALSE]
  mu_u = base_forecasts.mu[1:k]
  mu_b = base_forecasts.mu[(k + 1):n]
  inco = ((A %*% mu_b) - mu_u)
  Q = Psi_u - (Psi_ub %*% t(A)) - (A %*% t(Psi_ub)) + (A %*% 
                                                         Psi_b %*% t(A))
  
  g_u = Psi_u - sum(Psi_ub)
  g_b = colSums(Psi_b) - Psi_ub
  
  threshold_u = (Q*(((posterior.nu - 3)/(posterior.nu - 4))*(Psi_u/(Psi_u - (g_u^2/Q))) - 1))
  threshold_b1 = (Q*(((posterior.nu - 3)/(posterior.nu - 4))*(Psi_b[1]/(Psi_b[1] - (g_b[1]^2/Q))) - 1))
  threshold_b2 = (Q*(((posterior.nu - 3)/(posterior.nu - 4))*(Psi_b[2]/(Psi_b[2] - (g_b[2]^2/Q))) - 1))
  
  thresholds = c(threshold_u, threshold_b1, threshold_b2)
  
  return(list(thresholds = thresholds, inco = inco ))
}

simulate_trend_seasonal_ts = function(L, N_b, A, gen_model_par, err_cov){
  s = gen_model_par$seasons
  errors = MASS::mvrnorm(n = L, mu = rep(0,N_b), Sigma = err_cov)
  
  mu = matrix(0, nrow = L + 1, ncol = N_b)
  nu = matrix(0, nrow = L + 1, ncol = N_b)
  gamma = matrix(0, nrow = L + s, ncol = N_b)  # s lags for seasonal
  B = matrix(0, nrow = L, ncol = N_b)
  eta_t = matrix(0, nrow = L, ncol = N_b)
  
  #Simulating arima noise
  for (i in 1:N_b){
    eta_t[, i] = arima.sim(list(order(1,0,1), ar = 0.3, ma = 0.5), n = L, innov = errors[,i])
  }
  
  #Other noises
  zeta = matrix(rnorm(L * N_b, mean = 0, sd = sqrt(gen_model_par$sigma_err_trend2)), nrow = L, ncol = N_b)
  omega = matrix(rnorm(L * N_b, mean = 0, sd = sqrt(gen_model_par$sigma_err_seas)), nrow = L, ncol = N_b)
  epsilon = matrix(rnorm(L * N_b, mean = 0, sd = sqrt(gen_model_par$sigma_err_trend1)), nrow = L, ncol = N_b)
  
  #Initial values
  mu[1, ] = rnorm(N_b)
  nu[1, ] = rnorm(N_b)
  for (i in 1:s) {
    gamma[i, ] = rnorm(N_b)
  }
  
  #TS components
  for (t in 1:L) {
    nu[t + 1, ] = nu[t, ] + zeta[t, ] 
    mu[t + 1, ] = mu[t, ] + nu[t, ] + epsilon[t, ]
    gamma[t + s, ] = -colSums(gamma[(t + s - (s - 1)):(t + s - 1), , drop = FALSE]) + omega[t, ]
    B[t, ] = mu[t + 1, ] + gamma[t + s, ] + eta_t[t, ]
  }
  
  B = t(B)
  U = A %*% B
  Y = rbind(U,B)
  
  return(Y)
}

initialize_if_null_sim = function(results, measure, met, type, shrink, dims) {
  if (is.null(results[[measure]][[met]][[type]])) {
    results[[measure]][[met]][[type]] = array(0, dim = dims)
  }
  return(results)
}

update_results_matrix_sim = function(results, metrics_list, dist_matrix, met, type = NULL, shrink = NULL, t) {
  
  N = length(dist_matrix$mean)
  num_conf_levels = ncol(metrics_list$is)
  num_quantiles = ncol(metrics_list$quantile_score)
  
  # Define metric keys and their dimensions, aligning with `initialize_results`
  metric_dims = list(
    "Fc" = c(N, t), 
    "MIS" = c(N, num_conf_levels, t),
    "COVERAGE" = c(N, num_conf_levels, t),
    "MSE" = c(N, t),
    "QUANTILES_SCORE" = c(N, num_quantiles, t),
    "Extremes_Conf_Int_L" = c(N, num_conf_levels, t),
    "Extremes_Conf_Int_R" = c(N, num_conf_levels, t)
  )
  
  # Determine whether to use type/shrink as a sub-key
  sub_key = if (!is.null(shrink)) shrink else if (!is.null(type)) type else NULL
  
  # Initialize missing elements in results using `initialize_if_null_sim`
  for (measure in names(metric_dims)) {
    if (!is.null(sub_key)) {
      results = initialize_if_null_sim(results, measure, met, sub_key, shrink, metric_dims[[measure]])
    } else {
      if (is.null(results[[measure]][[met]])) {
        results[[measure]][[met]] = array(0, dim = metric_dims[[measure]])
      }
    }
  }
  
  # Assign new computed metrics in a vectorized way, ensuring indexing is correct
  if (!is.null(sub_key)) {
    results[["Fc"]][[met]][[sub_key]][, t] = dist_matrix$mean
    results[["MIS"]][[met]][[sub_key]][,, t] = metrics_list$is
    results[["COVERAGE"]][[met]][[sub_key]][,, t] = metrics_list$coverage
    results[["MSE"]][[met]][[sub_key]][, t] = metrics_list$se
    results[["QUANTILES_SCORE"]][[met]][[sub_key]][,, t] = metrics_list$quantile_score
    results[["Extremes_Conf_Int_L"]][[met]][[sub_key]][,, t] = metrics_list$Extremes_Conf_Int[[1]]
    results[["Extremes_Conf_Int_R"]][[met]][[sub_key]][,, t] = metrics_list$Extremes_Conf_Int[[2]]
    
  } else {
    results[["Fc"]][[met]][, t] = dist_matrix$mean
    results[["MIS"]][[met]][,, t] = metrics_list$is
    results[["COVERAGE"]][[met]][,, t] = metrics_list$coverage
    results[["MSE"]][[met]][, t] = metrics_list$se
    results[["QUANTILES_SCORE"]][[met]][,, t] = metrics_list$quantile_score
    results[["Extremes_Conf_Int_L"]][[met]][,, t] = metrics_list$Extremes_Conf_Int[[1]]
    results[["Extremes_Conf_Int_R"]][[met]][,, t] = metrics_list$Extremes_Conf_Int[[2]]
  }
  
  return(results)
}

#Function to compute and store the metrics
# j is the index for the rolling window
# recon_results is the reconciled distributions on the hierarchy
# results is the initialized empty list to be filled
compute_and_store_metrics_matrix_sim = function(met, recon_results, actuals, .scales, 
                                            confidence_levels, quantiles, results, type, shrink, j) {
  # Create distribution object based on the method
  .dist = if (met == "t_Rec") {
    list(mean = recon_results$mean, 
         scale_par = diag(recon_results$scale_par), 
         df = recon_results$df)
  } else {
    list(mean = recon_results$mean, 
         cov = diag(recon_results$cov))
  }
  
  # Determine distribution type
  dist_type = if (met == "t_Rec") "t_stud" else "norm"
  
  # Compute all performance metrics in a vectorized manner
  metrics = list(
    is = is_score_vectorized(.dist, actuals, confidence_levels, .scales$scale_mis, dist_type),
    coverage = coverage_vectorized(.dist, actuals, confidence_levels, dist_type),
    se = (.dist$mean - actuals)**2,
    quantile_score = quant_score_vectorized(.dist, actuals, quantiles, dist_type, .scales$scale_qt),
    Extremes_Conf_Int = confidence_intervals_vectorized(.dist, confidence_levels, dist_type)
  )

  results = update_results_matrix_sim(results, metrics, .dist, met, type, shrink, j)
  
  return(results)
}


