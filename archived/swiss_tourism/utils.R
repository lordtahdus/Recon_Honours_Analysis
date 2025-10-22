dir_path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

#Function for initializing the results
{# Initialization of scores and other measures to compute and save
  
  #INPUT
  # N: number of variables;
  # tot_windows: total number of rolling windows;
  # methods: methods to be tested;
  # types: types of prior matrix for t-Rec
  # specific_method: method in methods that has the types (in our case is "t_Rec")
  # confidence_levels: confidence levels for the MIS, the PI_WIDTH and the COVERAGE
  # quantiles: quantiles to be computed for the QUANTILES_SCORE
  
  #OUTPUT
  #results: an empty list with measures, methods and types levels
  initialize_results = function(N, tot_windows, methods, types, specific_method,
                                confidence_levels, quantiles) {
    performance_measures = c("Fc", "MSE", "CRPS", "LOG", "MIS", 
                             "COVERAGE", "QUANTILES_SCORE", "PI_WIDTH")
    results = list()
    
    for (measure in performance_measures) {
      results[[measure]] = list()
      for (met in methods) {
        if (met == specific_method) {
          results[[measure]][[met]] = list()
          for (type in types) {
            if (measure == "MIS" | measure == "COVERAGE" | measure == "PI_WIDTH"){
              results[[measure]][[met]][[type]] = array(0, dim = c(N, 
                                                                   length(confidence_levels), 
                                                                   tot_windows))
            }else if (measure == "QUANTILES_SCORE"){
              results[[measure]][[met]][[type]] = array(0, dim = c(N, 
                                                                   length(quantiles), 
                                                                   tot_windows))
            }else{
              results[[measure]][[met]][[type]] = array(0, dim = c(N, tot_windows))
            }
          }
        } else if (met == "MinT"){
          results[[measure]][[met]] = list()
          for (shrink in c("shrink", "no_shrink")) {
            if (measure == "MIS" | measure == "COVERAGE"| measure == "PI_WIDTH"){
              results[[measure]][[met]][[shrink]] = array(0, dim = c(N, 
                                                                     length(confidence_levels), 
                                                                     tot_windows))
            }else if (measure == "QUANTILES_SCORE"){
              results[[measure]][[met]][[shrink]] = array(0, dim = c(N, 
                                                                     length(quantiles), 
                                                                     tot_windows))
            }else{
              results[[measure]][[met]][[shrink]] = array(0, dim = c(N, tot_windows))
            }
          }
        }else {
          if (measure == "MIS" | measure == "COVERAGE"| measure == "PI_WIDTH"){
            results[[measure]][[met]] = array(0, dim = c(N, 
                                                         length(confidence_levels), 
                                                         tot_windows))
          }else if (measure == "QUANTILES_SCORE"){
            results[[measure]][[met]] = array(0, dim = c(N, 
                                                         length(quantiles), 
                                                         tot_windows))
          }else{
            results[[measure]][[met]] = array(0, dim = c(N, tot_windows))
          }
          
        }
      }
    }
    
    return(results)
  }}

#Functions for preparing the data (Swiss_tourism, Swiss_tourism_Bigger,
#Australian_tourism_zone, Australian_tourism_no_zone)

#You need to have saved the file "SwissTourism.csv" and "Regions.csv" in a folder 
#called data in the same folder where you saved this code.

{
  #Function for formatting the Swiss Tourism dataset
  
  #INPUT
  #file_name: the "SwissTourism.csv" file
  
  #OUTPUT
  #a list of two component:
  # data_ts: the data as time series object
  # mdata: the data as matrix (length_training x N)
  
  FormatSwissTourismData = function(file_name){
    current_year = format(Sys.Date(), "%Y")
    current_month = format(Sys.Date(), "%m")
    
    data_ = read.csv(file_name, check.names = F, fileEncoding = "latin1")
    years = unique(data_[, "Anno"])
    month_names = c("Gennaio",   "Febbraio",  "Marzo",     "Aprile",    "Maggio",
                    "Giugno",    "Luglio",    "Agosto",   "Settembre", "Ottobre",
                    "Novembre",  "Dicembre" )
    months = match(data_[, "Mese"], month_names)
    start_year = as.numeric(years[1])
    start_month = as.numeric(months[1])
    
    bottom = unique(data_$Cantone)[2:length(unique(data_$Cantone))]
    upper = unique(data_$Cantone)[1]
    
    mdata = matrix(t(as.numeric(data_$Pernottamenti[!is.na(as.numeric(data_$Pernottamenti))])), ncol = 27, byrow = T)
    colnames(mdata) = c(upper, bottom)
    
    dfData = data.frame(mdata)
    
    num_months <- nrow(dfData)
    
    time_labels <- format(seq(from = as.Date(paste(start_year, start_month, "01", sep = "-")), 
                              by = "month", length.out = num_months), 
                          "%Y-%m")
    
    rownames(dfData) = time_labels
    data_ts = ts(dfData[,2:ncol(dfData)], start = c(as.numeric(start_year), as.numeric(start_month)), frequency = 12)
    return(list(data_ts = data_ts, mdata = mdata))
  }
  
  
  #Function for building the hierarchy for a selected datasets among 
  #Swiss_tourism, Swiss_tourism_Bigger, Australian_tourism_zone and 
  #Australian_tourism_no_zone
  
  #INPUT
  #name_dataset: a string with the name of the dataset among Swiss_tourism, 
  #Swiss_tourism_Bigger, Australian_tourism_zone and Australian_tourism_no_zone
  
  #OUTPUT
  #out: a list of eight component:
  #     dataset: the raw data read from .csv
  #     n_bottom: number of bottom variable in the hierarchy
  #     n_upper: number of upper variable in the hierarchy
  #     ts_length: number of the observations in the ts
  #     agg_mat: the aggregation matrix
  #     sum_mat: the summation matrix
  #     ts: a matrix of dimension (length_training x N)
  #     freq: the frequency of the ts (for monthly data is 12, for quartely is 4)
  
  get_dataset = function(name_dataset){
    if (name_dataset == "Swiss_tourism"){
      freq = 12
      file_name = paste0(dir_path, "/data/SwissTourism.csv")
      data_ = read.csv(file_name)
      
      Y = t(FormatSwissTourismData(file_name)$mdata)
      n_bottom = 26
      n_upper = 1
      ts_length = nrow(data_)
      
      A = matrix(1, nrow = n_upper, ncol = n_bottom, dimnames = list(rownames(Y)[1],
                                                                     rownames(Y)[-1]))
      S = rbind(A, diag(rep(1, n_bottom)))
      
    }else if (name_dataset == "Swiss_tourism_Bigger"){
      freq = 12
      file_name = paste0(dir_path, "/data/SwissTourism.csv")
      data_ = read.csv(file_name)
      
      n_bottom = 26
      
      ts_length = nrow(data_)
      
      B = t(FormatSwissTourismData(file_name)$mdata)[-1,]
      
      franc = c("Genève", "Vaud", "Neuchâtel", "Jura")
      ita = c("Graubünden / Grigioni / Grischun", "Ticino")
      ted = setdiff(rownames(B), c(franc, ita))
      
      upper_names = c("Schweiz", "FZone", "IZone", "DZone")
      
      n_upper = 4
      
      A = matrix(0, nrow = n_upper, ncol = n_bottom, dimnames = list(upper_names,
                                                                     rownames(B)))
      for (i in 1:ncol(A)){
        for (j in 1:nrow(A)){
          if (rownames(A)[j] == "Schweiz"){
            A[j,i] = 1
          }else if ((rownames(A)[j] == "FZone")&(colnames(A)[i] %in% franc)){
            A[j,i] = 1
          }else if ((rownames(A)[j] == "IZone")&(colnames(A)[i] %in% ita)){
            A[j,i] = 1
          }else if ((rownames(A)[j] == "DZone")&(colnames(A)[i] %in% ted)){
            A[j,i] = 1
          }
        }
      }
      
      S = rbind(A, diag(rep(1, n_bottom)))
      Y = S %*% B
      dimnames(Y) = list(c(upper_names,rownames(B) ), NULL)
      
    }else if (name_dataset == "Australian_tourism_zone"){
      freq = 12
      data_ = read.csv(file = paste0(dir_path, "/data/Regions.csv"))
      n_bottom = 76
      n_upper = 35
      ts_length = nrow(data_)
      
      state_regions_map = list()
      zones_regions_map = list()
      
      state_list = LETTERS[1:7]  
      for (state in state_list) {
        
        state_data = data_[,2:ncol(data_)] %>% dplyr::select(starts_with(state))
        
        
        regions_per_state = unique(unlist(colnames(state_data)))
        state_regions_map[[state]] = regions_per_state
        
        
        state_zones = colnames(state_data)
        for (zone in state_zones) {
          zone_short = substr(zone, start = 1, stop = 2)
          zones_regions_map[[zone_short]] = colnames(data_[,2:ncol(data_)] %>% dplyr::select(starts_with(zone_short)))
        }
      }
      
     
      all_regions = unique(unlist(c(state_regions_map, zones_regions_map)))
      
      
      num_rows = n_upper-1 
      num_cols = n_bottom
      
      
      A = matrix(0, nrow = num_rows, ncol = num_cols, dimnames = list(c(names(state_regions_map), names(zones_regions_map)), all_regions))
      new_row = rep(1, num_cols)
      A = rbind(Total = new_row, A)
      
      for (state in names(state_regions_map)) {
        for (region in state_regions_map[[state]]) {
          A[state, region] = 1
        }
      }
      
      
      for (zone in names(zones_regions_map)) {
        for (region in zones_regions_map[[zone]]) {
          A[zone, region] = 1
        }
      }
      
      S = rbind(A, diag(rep(1, n_bottom)))
      B = matrix(nrow = n_bottom, ncol = ts_length)
      for (i_B in 1:n_bottom) {
        B[i_B,] = t(data_[, i_B + 1])
      }
      
      
      U = A %*% B
      Y = (rbind(U,B))
      
    }else if (name_dataset == "Australian_tourism_no_zone"){
    
      freq = 4
      
      data_ = tourism %>%
        index_by(Quarter) %>%
        group_by(Region, State) %>%
        summarise(Total_trips = sum(Trips))
      
      n_bottom = 77
      n_upper = 8 
      ts_length = length(unique(data_$Quarter)) 
      
      unique_states = setdiff(unique(data_$State), "ACT")
      unique_regions = c(unique(data_$Region), "ACT")
      
      A = matrix(0, nrow = n_upper, ncol = n_bottom,
                 dimnames = list(c("Australian Tourism",unique_states), unique_regions))
      A[1,] = rep(1, ncol(A))
      for (i_A in 2:n_upper) {
        for (j_A in 1:n_bottom) {
          if (unique_regions[j_A] %in% data_$Region[data_$State == unique_states[i_A-1]]) {
            A[i_A, j_A] = 1
          }
        }
      }
      
      S = rbind(A, diag(rep(1, n_bottom)))
      B = matrix(nrow = n_bottom, ncol=ts_length)
      for (i_B in 1:(n_bottom-1)) {
        select_region = data_ %>%
          filter(Region == unique_regions[i_B])
        B[i_B,] = select_region$Total_trips
      }
      select_region = data_ %>%
        filter(State == "ACT")
      B[n_bottom,] = select_region$Total_trips
      
      U = A %*% B
      Y = (rbind(U,B))
      
      
    }
  
    
    out = list(dataset = data_, n_bottom = n_bottom, n_upper = n_upper, 
              ts_length = ts_length, agg_mat = A, sum_mat = S, ts = Y, freq = freq)
    
    
    return(out)
  }
  }

#Function for computing base forecasts
{#Function that allow to choose the model for the base forecasts
  
  #INPUT
  #.model: univariate forecasting model in the forecast package among ets, 
  #        auto.arima,naive, seasonal_naive, mixed
  #.data: a vector representing a single time series
  #freq: frequency of the time series considered in .data
  #mask: boolean value resulted from the seasonal test for only one time series 
  
  #OUTPUT
  #base_fc_model: the function correspondent to the .model selected
  
  base_fc_model = function(.model, .data, freq, mask){
    if (.model == "ets"){
      return(forecast::ets(ts(.data, frequency = freq), model = "AZZ"))
    }else if (.model == "auto.arima"){
      return(forecast::auto.arima(ts(.data, frequency = freq)))
    }else if (.model == "naive"){
      return(forecast::naive(ts(.data)))
    }else if (.model == "seasonal_naive"){
      return(forecast::snaive(ts(.data, frequency = freq)))
    }else if (.model == "mixed"){
      return(if (mask) {forecast::snaive(ts(.data, frequency = freq))
      }else {forecast::naive(ts(.data))})
    }
  }
  
  
  # Function to Fit Models and Compute Base Forecasts
  
  #INPUT
  #train: a (N x length_training) matrix where N is the number of time series
  #h: forecast horizon
  #.model: forecasting model in the forecast package among ets, auto.arima,
  #        naive, seasonal_naive, mixed 
  #freq: frequency of the ts
  #mask: a boolean vector resulted from the seasonal test for all the time series
  
  #OUTPUT
  #out: list of two component: 
  #     base_fc: a (N x h) matrix with the point forecasts correspondent to the .model
  #              chosen
  #     res: a (length_training x N) matrix with the residuals correspondent to 
  #          the .model chosen
  
  compute_base_forecasts = function(train, h, .model, freq, mask = NULL) {
    N = nrow(train) #number of variable
    base_fc = matrix(NA, nrow = N, ncol = h)
    res = matrix(NA, nrow = ncol(train), ncol = N)
    models = vector("list", N)
    
    for (i in 1:N){
      models[[i]] = base_fc_model(.model, train[i, ], freq, mask[i])
      f = forecast(models[[i]], h = h)
      base_fc[i, ] = f$mean
      res[, i] = models[[i]]$residuals
    }
    
    if (((.model == "mixed"))|(.model == "seasonal_naive")){
      if (freq == 0){
        res = res[(freq + 2):ncol(train), ] #since when we have seasonality we have 
        #shorter residuals the dimension is adjusted
      }else{
        res = res[(freq + 1):ncol(train), ] #since when we have seasonality we have 
        #shorter residuals the dimension is adjusted
      }
      
    }else if (.model == "naive"){
      res = res[2:ncol(train), ]
    }
    
    
    out = list(base_fc = base_fc, res = res)
    return(out)
  }
}

# Function for reconciling hierarchical probabilistic forecasts
{
  # Check if it is a covariance matrix (i.e. symmetric p.d.)
  .check_cov = function(cov_matrix, Sigma_str,pd_check=FALSE,symm_check=FALSE) {
    # Check if the matrix is square
    if (!is.matrix(cov_matrix) || nrow(cov_matrix) != ncol(cov_matrix)) {
      stop(paste0(Sigma_str, " is not square"))
    }
    
    # Check if the matrix is positive semi-definite
    if(pd_check){
      eigen_values = eigen(cov_matrix, symmetric = TRUE)$values
      if (any(eigen_values <= 0)) {
        stop(paste0(Sigma_str, " is not positive semi-definite"))
      }
    }
    if(symm_check){
      # Check if the matrix is symmetric
      if (!isSymmetric(cov_matrix)) {
        stop(paste0(Sigma_str, " is not symmetric"))
      }
    }
    # Check if the diagonal elements are non-negative
    if (any(diag(cov_matrix) < 0)) {
      stop(paste0(Sigma_str, ": some elements on the diagonal are negative"))
    }
    # If all checks pass, return TRUE
    return(TRUE)
  }
  
  # Function to check aggregation matrix A
  .check_A <- function(A) {
    if (!all(A %in% c(0,1))) {
      stop("Input error in A: A must be a matrix containing only 0s and 1s.")
    }
    
    if(any(colSums(A)==0)){
      stop("Input error in A: some columns do not have any 1.
          All bottom level forecasts must aggregate into an upper.")
    }
    
    if(nrow(unique(A))!=nrow(A)){
      warning("A has some repeated rows.")
    }
  }
  
  #t-Reconciliation given prior or posterior IW parameters
  
  #INPUT
  #A: aggregation matrix
  #base.forecasts.mu: the mean of the base forecasts
  #prior.psi: the scale matrix of the prior IW distribution
  #prior.nu: the degree of freedom of the prior IW distribution
  #sample_cov: the sample covariance matrix of the residuals of the method with 
  #            which the base forecasts are obtained
  #n_oss: length_training
  #post.psi: the scale matrix of the posterior IW distribution
  #post.nu: the degree of freedom of the posterior IW distribution
  
  #OUTPUT
  #out: a list of six components:
  #     bottom_reconciled_mean: mean of the reconciled multivariate t-distribution 
  #                             of the bottom variables
  #     bottom_reconciled_scale_parameter: scale matrix of the reconciled 
  #                                        multivariate t-distribution of the 
  #                                        bottom variables
  #     bottom_reconciled_dof_parameter: degree of freedom of the reconciled 
  #                                      multivariate t-distribution of the 
  #                                      bottom variables
  #     posterior.psi: the scale matrix of the posterior IW distribution
  #     posterior.nu: the degree of freedom of the posterior IW distribution
  #     const: the constant that depends on the incoherence (A\hat{b} - \hat{u})
  
  reconc_t_Rec = function (A, base_forecasts.mu, prior.psi, prior.nu, sample_cov, 
                               n_oss, post.psi, post.nu){
    .check_A(A)
    k = nrow(A)
    m = ncol(A)
    n = length(base_forecasts.mu)
    if (!(nrow(prior.psi) == n) && !(nrow(post.psi) == n)) {
      stop("Input error: nrow(prior.psi) != length(base_forecasts.mu)")
    }
    if (!(k + m == n)) {
      stop("Input error: the shape of A is not correct")
    }
    .check_cov(prior.psi, "Psi", pd_check = FALSE, symm_check = TRUE)
    if (is.null(post.psi) | is.null(post.nu)){
      posterior.psi = (n_oss * sample_cov) + prior.psi
      posterior.nu = prior.nu + n_oss
    }else{
      posterior.psi = post.psi
      posterior.nu = post.nu
    }
    
    .check_cov(posterior.psi, "posterior", pd_check = T, symm_check = T)
    Psi_u = posterior.psi[1:k, 1:k]
    Psi_b = posterior.psi[(k + 1):n, (k + 1):n]
    Psi_ub = posterior.psi[1:k, (k + 1):n, drop = FALSE]
    mu_u = base_forecasts.mu[1:k]
    mu_b = base_forecasts.mu[(k + 1):n]
    inco = ((A %*% mu_b) - mu_u)
    Q = Psi_u - (Psi_ub %*% t(A)) - (A %*% t(Psi_ub)) + (A %*% 
                                                           Psi_b %*% t(A))
    .check_cov(Q, "Q", pd_check = TRUE, symm_check = FALSE)
    invQ = solve(Q)
    C = 1 + ((t(inco) %*% invQ %*% inco))
    nu_b_tilde = posterior.nu - m + 1
    mu_b_tilde = mu_b + (t(Psi_ub) - Psi_b %*% t(A)) %*% 
      invQ %*% inco
    Sigma_b_tilde = as.numeric(C/nu_b_tilde) * (Psi_b - ((t(Psi_ub) - (Psi_b %*% t(A))) %*% 
                                                           invQ %*% t(t(Psi_ub) - (Psi_b %*% t(A)))))
    out = list(bottom_reconciled_mean = mu_b_tilde, 
               bottom_reconciled_scale_parameter = Sigma_b_tilde, 
               bottom_reconciled_dof_parameter = nu_b_tilde,
               posterior.psi = posterior.psi,
               posterior.nu = posterior.nu,
               const = C)
    return(out)
  }
  
  #t-Reconciliation through the different types
  #(it is meant to apply this function for every rolling windows)
  
  #INPUT
  #A: aggregation matrix
  #S: summation matrix
  #h: forecast horizon
  #base_fc: the mean of the base forecasts
  #type: the selected type from types
  #res: a (length_training x N) matrix of the residuals of the method with 
  #     which the base forecasts are obtained
  #train: (N x length_training) matrix of the ts data
  #freq: frequency of the ts
  #mask: boolean value resulted from the seasonal test for only one time series
  
  #OUTPUT
  #out: a list of seven components:
  #     mean: mean of the reconciled multivariate t-distribution 
  #           of all the variables
  #     scale_par: scale matrix of the reconciled 
  #                multivariate t-distribution of all the variables
  #     df: degree of freedom of the reconciled 
  #         multivariate t-distribution of all the variables
  #     prior.psi: the scale matrix of the prior IW distribution
  #     prior.nu: the degree of freedom of the prior IW distribution
  #     posterior.psi: the scale matrix of the posterior IW distribution
  #     posterior.nu: the degree of freedom of the posterior IW distribution
  
  reconcile_t_forecasts = function(A, S, h, base_fc, type, res, train,
                                   freq = NULL, mask = NULL) {
    
    Sigma = crossprod(res)/nrow(res) 
    Sigmahat = (1 - 10^-4)*Sigma + 10^-4*diag(diag(Sigma))
    N = nrow(Sigmahat)
    
    #we want that the mean of the IW prior is the covariance matrix given by the 
    #mixed naive residuals (seasonal and not seasonal)
    
    #the following is for extracting the residuals, it is not for any base forecast
    if (ncol(train) <= 2*freq){
      fc_and_res = compute_base_forecasts(train, h, "naive", freq, mask = mask)
    }else{
      fc_and_res = compute_base_forecasts(train, h, "mixed", freq, mask = mask)
    }
    
    res_naive = fc_and_res$res
    
    naive_mat = schaferStrimmer_cov(res_naive)$shrink_cov
    
    if (type == "Full"){
    
      bayesian_LOO = multi_log_score_optimization(res, nrow(res), N, 
                                                     naive_mat)
        
      prior.nu = bayesian_LOO[[1]]$solution[1]
      prior.psi = (naive_mat)*(prior.nu - N - 1)
      
    }else if (type == "Diag"){
      
      naive_mat_diag = diag(diag(naive_mat))
      
      bayesian_LOO = multi_log_score_optimization(res, nrow(res), N, 
                                                     naive_mat_diag)
        
      prior.nu = bayesian_LOO[[1]]$solution[1]
      prior.psi = (naive_mat_diag)*(prior.nu - N - 1)

    }
    
    # t-Student reconciliation
    recon_t_stud = reconc_t_Rec(A, base_fc, prior.nu = prior.nu,
                                prior.psi = prior.psi, sample_cov = Sigmahat,
                                n_oss = nrow(res), post.psi = NULL, 
                                post.nu = NULL)
    
    reconciled_mean = S %*% recon_t_stud$bottom_reconciled_mean
    scale_par = S %*% recon_t_stud$bottom_reconciled_scale_parameter %*% t(S)
    df = recon_t_stud$bottom_reconciled_dof_parameter
    
    posterior.psi = recon_t_stud$posterior.psi
    posterior.nu = recon_t_stud$posterior.nu
    
    
    out = list(mean = reconciled_mean, scale_par = scale_par, df = df, 
               prior.psi = prior.psi, prior.nu = prior.nu,
               posterior.psi = posterior.psi, posterior.nu = posterior.nu)
    return(out)
  }

  
  #Bottom-Up reconciliation 
  
  #INPUT
  #A: aggregation matrix
  #S: summation matrix
  #base_fc: the mean of the base forecasts
  #res: a (length_training x N) matrix of the residuals of the method with 
  #     which the base forecasts are obtained
  
  #OUTPUT
  #out: list of two components:
  #     mean: mean of the BU reconciled multivariate Gaussian distribution 
  #           of all the variables
  #     cov: covariance matrix of the BU reconciled 
  #          Gaussian distribution of all the variables
  
  reconcile_BU_forecasts = function(A, S, base_fc, res){
    Sigmahat = crossprod(res)/nrow(res)
    N = nrow(Sigmahat)
    N_u = nrow(A)
    
    #BU reconciliation
    reconciled_mean = S %*% base_fc[(N_u + 1):N]
    cov_mat = S %*% Sigmahat[(N_u + 1):N, (N_u + 1):N] %*% t(S)
    out = list(mean = reconciled_mean, cov = cov_mat)
    return(out)
    
  }
  
  #MinT Gaussian reconciliation
  
  #INPUT
  #A: aggregation matrix
  #S: summation matrix
  #base_fc: the mean of the base forecasts
  #res: a (length_training x N) matrix of the residuals of the method with 
  #     which the base forecasts are obtained
  #shrink: boolean value --> T: the shrunk sample covariance matrix will be used
  #                          F: the sample covariance matrix will be used
  
  #OUTPUT
  #out: list of two components:
  #     mean: mean of the BU reconciled multivariate Gaussian distribution 
  #           of all the variables
  #     cov: covariance matrix of the BU reconciled 
  #          Gaussian distribution of all the variables
  
  reconcile_gauss_forecasts = function(A, S, base_fc, res, shrink){
    Sigmahat = crossprod(res)/nrow(res)
    SigmahatS = schaferStrimmer_cov(res)$shrink_cov
    
    if (shrink == "shrink"){
      #Gaussian reconciliation with shrunk matrix
      recon_gauss = reconc_gaussian(A, base_fc, SigmahatS)
      reconciled_mean = S %*% recon_gauss$bottom_reconciled_mean
      cov_mat = S %*% recon_gauss$bottom_reconciled_covariance %*% t(S)
    }else{
      #Gaussian reconciliation NOT shrunk 
      Sigma = (1-10^-4)*Sigmahat + 10^-4*diag(diag(Sigmahat)) #for stability
      recon_gauss = reconc_gaussian(A, base_fc, Sigma)
      reconciled_mean = S %*% recon_gauss$bottom_reconciled_mean
      cov_mat = S %*% recon_gauss$bottom_reconciled_covariance %*% t(S)
    }
    
    out = list(mean = reconciled_mean, cov = cov_mat)
    return(out)
  }
  
  #Function that select the right reconciliation according to the method 
  
  #INPUT
  #method: a method selected among the vector ("t-Rec", "BU", "MinT")
  #A: aggregation matrix
  #S: summation matrix
  #h: forecast horizon
  #base_fc: the mean of the base forecasts
  #type: the selected type from types
  #res: a (length_training x N) matrix of the residuals of the method with 
  #     which the base forecasts are obtained
  #train: (N x length_training) matrix of the ts data
  #freq: frequency of the ts
  #mask: boolean value resulted from the seasonal test for only one time series
  #shrink: boolean value --> T: the shrunk sample covariance matrix will be used
  #                          F: the sample covariance matrix will be used
  
  #OUTPUT
  #the output of the called functions
  
  reconcile_forecasts = function(method, A, S, h, base_fc, type, res, train, 
                                 freq, mask, shrink){
    if (method == "t_Rec"){
      return(reconcile_t_forecasts(A, S, h, base_fc, type, res, train, freq, mask))
    }else if (method == "BU"){
      return(reconcile_BU_forecasts(A, S, base_fc, res))
    }else if (method == "MinT"){
      return(reconcile_gauss_forecasts(A, S, base_fc, res, shrink))
    }
  }
}

#Functions for computing performance metrics
{
  
  #Read the "no_scaled" files and compute the geometric mean over the variables
  geom_mean_computing = function(dataset, metric_name, length_training, n_digit, conf_levs, 
                                 quants, methods){
    
    substrings = c("Base", methods)
    pattern <- paste(substrings, collapse = "|")
    
    
    # Define the folders where the measure are stored
    measure_folder_path = paste0("results_folder/",dataset, "/", metric_name, "/length_train_", length_training, "/")
    
    
    # Check if both folders exist
    if (!dir.exists(measure_folder_path)) {
      dir.create(measure_folder_path, recursive = TRUE)
    }
    
    
    # Get a list of all CSV files related to the specified measure
    measure_pattern = paste0("no_scaled_",metric_name, "_.*_length_training_", length_training, "\\.csv")
    
    
    measure_files = list.files(measure_folder_path, pattern = measure_pattern, full.names = TRUE)
    
    
    # Check if there are matching files for both measure 
    if (length(measure_files) == 0) {
      stop("No CSV files found for the specified measure: ", metric_name)
    }
    
    # Initialize lists to store data for measures
    measure_data_list = list()
    
    measure_files_ = grep(pattern, measure_files, value = TRUE)
    for (file in measure_files_) {
      # Extract the method name from the file name
      method_name = sub(paste0("no_scaled_",metric_name, "_(.*)_length_training_", length_training, "\\.csv"), "\\1", basename(file))
      
      # Read the CSV file for the measure
      .df = read.csv2(file, row.names = 1)
      
        if (metric_name == "MIS"| metric_name == "PI_WIDTH"){
          
          measure_data_list[[method_name]] = setNames(as.data.frame(round(exp(colMeans(log(.df))), n_digit), 
                                                                    row.names = paste0(as.character(100 * conf_levs), "%")), 
                                                      method_name)
        }else if (metric_name == "COVERAGE"){
          measure_data_list[[method_name]] = setNames(as.data.frame(round(colMeans(.df), n_digit), 
                                                                    row.names = paste0(as.character(100 * conf_levs), "%")), 
                                                      method_name)
        }else if (metric_name == "MSE"){
          measure_data_list[[method_name]] = setNames(as.data.frame(round(colMeans(.df), n_digit)), 
                                                      method_name)
        }else if (metric_name == "QUANTILES_SCORE"){
          
          measure_data_list[[method_name]] = setNames(as.data.frame(round(exp(colMeans(log(.df))), n_digit), 
                                                                    row.names = paste0(as.character(100 * quants), "th")), 
                                                      method_name)
        }else{
          measure_data_list[[method_name]] = setNames(as.data.frame(round(exp(colMeans(log(.df))), n_digit)), 
                                                      method_name)
        }
        
        combined_df = do.call(cbind, lapply(names(measure_data_list), function(method) {
          if (!is.null(measure_data_list[[method]])) {
            measure_values = measure_data_list[[method]][, 1]
          } else {
            NA
          }
        }))
        
        
        # Set the column names to the method names
        colnames(combined_df) = names(measure_data_list)
        
        # Optionally, set row names based on confidence levels or quantiles
        if (!is.null(confidence_levels) && metric_name %in% c("MIS", "COVERAGE","PI_WIDTH")) {
          rownames(combined_df) = paste0(as.character(100 * confidence_levels), "%")
        } else if (!is.null(quantiles) && metric_name %in% c("QUANTILES_SCORE")) {
          rownames(combined_df) = paste0(as.character(100 * quantiles), "th")
        }
      
    }
    
    return(combined_df)
    
    
  }
  
  #Function for computing the extremes of the prediction intervals

  #INPUT
  #.dist: a list with the parameters of the distributions
  #conf_levs: vector of the confidence levels
  #name: name of the distribution between "norm" and "t_stud"
  #i.e., for the Gaussian distribution the format is: 
  #.dist = list(mean = mean, cov = cov) and name = "norm";
  #for the t distribution the format is: 
  #.dist = list(mean = mean, scale_par = scale_par, df = df) and name = "t_stud"

  # mean is the vector of the distributions mean of all the variables
  # cov/scale_par is the vector of the distributions variances/scale parameters of all the variables

  #OUTPUT
  #a list of two components:
  #lower_bound: a matrix of dimension (N x length(conf_levs)) with the lower 
  #             bounds for all the variables at the different confidence levels
  #upper_bound: a matrix of dimension (N x length(conf_levs)) with the upper 
  #             bounds for all the variables at the different confidence levels
  
  confidence_intervals_vectorized = function(.dist, conf_levs, name) {
    predicted_mean = .dist$mean
    N = length(predicted_mean)
    m = length(conf_levs)
    alpha = 1 - conf_levs  # Vector of alpha values
    
    if (name == "t_stud") {
      predicted_sd = sqrt(.dist$scale_par)
      lower = outer(predicted_sd, qt(alpha / 2, df = .dist$df), `*`) + 
        matrix(predicted_mean, nrow = N, ncol = m)
      upper = outer(predicted_sd, qt(1 - (alpha / 2), df = .dist$df), `*`) + 
        matrix(predicted_mean, nrow = N, ncol = m)
    } else if (name == "norm") {
      predicted_sd = sqrt(.dist$cov)
      lower = outer(predicted_sd, qnorm(alpha / 2), `*`) + 
        matrix(predicted_mean, nrow = N, ncol = m)
      upper = outer(predicted_sd, qnorm(1 - (alpha / 2)), `*`) + 
        matrix(predicted_mean, nrow = N, ncol = m)
    }
    
    return(list(lower_bound = (lower), upper_bound = (upper)))
  }
  
  
  #Function for computing the coverage for each confidence levels 
  
  #INPUT
  #.dist: a list with the parameters of the distributions
  #.actual: vector of dimension N with the actual values for every variable
  #conf_levs: vector of the confidence levels
  #name: name of the distribution between "norm" and "t_stud"
  
  #OUTPUT
  #coverage: a matrix of dimension (N x length(conf_levs)) of boolean values 
  #          --> T: the actual value for the variable is inside the prediction interval
  #              F: the actual value for the variable is not inside the prediction interval
  
  coverage_vectorized = function(.dist, .actual, conf_levs, name) {
    predict_intervals = confidence_intervals_vectorized(.dist, conf_levs, name)
    
    coverage = (predict_intervals$lower_bound < .actual) & (.actual < predict_intervals$upper_bound)
    return(coverage)  # Boolean matrix of results
  }
  
  #Function to compute the interval score given the predictive distributions
  
  #INPUT
  #.dist: a list with the parameters of the distributions
  #.actual: vector of dimension N with the actual values for every variable
  #conf_levs: vector of the confidence levels
  #.scale: a matrix of dimension (N x length(conf_levs)) with the scaling terms
  #        for each variable and confidence levels
  #name: name of the distribution between "norm" and "t_stud"
  
  #OUTPUT
  #a matrix of dimension (N x length(conf_levs)) of scaled interval scores 
  #for each variable and confidence levels
  
  is_score_vectorized = function(.dist, .actual, conf_levs, .scale, name) {
    N = length(.dist$mean)
    m = length(conf_levs)
    predict_intervals = confidence_intervals_vectorized(.dist, conf_levs, name)
    
    lower = predict_intervals$lower_bound
    upper = predict_intervals$upper_bound
    alpha = 1 - conf_levs  # Vector of alpha values
    
    interval_scores = (upper - lower) + 
      (matrix(rep(2 / alpha, N), ncol = m, byrow = T) * 
         matrix(pmax(0, lower - .actual), nrow = N, byrow = F)) + 
      (matrix(rep(2 / alpha, N), ncol = m, byrow = T) * 
         matrix(pmax(0, .actual - upper), nrow = N, byrow = F))
    
    return(interval_scores / .scale)
  }
  
  #Function to compute the quantile score given the predictive distribution
  
  #INPUT
  #.dist: a list with the parameters of the distributions
  #.actual: vector of dimension N with the actual values for every variable
  #quants: vector of the quantiles
  #.scale: a matrix of dimension (N x length(quants)) with the scaling terms
  #        for each variable and confidence levels
  #name: name of the distribution between "norm" and "t_stud"
  
  #OUTPUT
  #a matrix of dimension (N x length(quants)) of scaled quantile scores 
  #for each variable and quantiles
  
  quant_score_vectorized = function(.dist, .actual, quants, name, .scale) {
    N = length(.dist$mean)
    m = length(quants)
    
    if (name == "norm") {
      q = outer(sqrt(.dist$cov), qnorm(quants), `*`) + 
        matrix(.dist$mean, nrow = N, ncol = m)
    } else if (name == "t_stud") {
      q = outer(sqrt(.dist$scale_par), qt(quants, .dist$df), `*`) + 
        matrix(.dist$mean, nrow = N, ncol = m)
    }
    
    quantile_scores = ifelse(.actual > q, 
                             matrix(rep(quants, N), nrow = N, byrow = T) * (.actual - q), 
                             matrix(rep(1 - quants, N), nrow = N, byrow = T) * (q - .actual))
    return(quantile_scores / .scale)
  }
  
  #Function to compute the CRPS given the predictive distribution
  
  #INPUT
  #.dist: a list with the parameters of the distributions
  #.actual: vector of dimension N with the actual values for every variable
  #name: name of the distribution between "norm" and "t_stud"
  #.scale: a vector of dimension N with the scaling terms
  #        for each variable and confidence levels
  
  
  #OUTPUT
  #a vetor of dimension N of scaled CRPS for each variable
  
  crps_score_vectorized = function(.dist, .actual, name, .scale) {
    if (name == "norm") {
      score = mapply(crps_norm, .actual, .dist$mean, sqrt(.dist$cov)) / .scale #if .actual is a matrix better using apply
    } else if (name == "t_stud") {
      score = mapply(crps_t, .actual, .dist$df, .dist$mean, sqrt(.dist$scale_par)) / .scale
    }
    return(score)
  }
  
  #Function to compute the log score given the predictive distribution
  
  #INPUT
  #.dist: a list with the parameters of the distributions
  #.actual: vector of dimension N with the actual values for every variable
  #name: name of the distribution between "norm" and "t_stud"
  
  #OUTPUT
  #a vetor of dimension N of log scores for each variable
  
  log_score_vectorized = function(.dist, .actual, name) {
    if (name == "norm") {
      score = mapply(logs_norm, .actual, .dist$mean, sqrt(.dist$cov)) #if .actual is a matrix better using apply
    } else if (name == "t_stud") {
      score = mapply(logs_t, .actual, .dist$df, .dist$mean, sqrt(.dist$scale_par))
    }
    return(score)
  }

  
}

#Functions for storing the results (for every rolling and variable and 
#confidence levels and quantiles)
{
  
  initialize_if_null = function(results, measure, met, type, shrink, dims) {
    if (is.null(results[[measure]][[met]][[type]])) {
      results[[measure]][[met]][[type]] = array(0, dim = dims)
    }
    return(results)
  }
  
  update_results_matrix = function(results, metrics_list, dist_matrix, met, type = NULL, shrink = NULL, t) {
    
    N = length(dist_matrix$mean)
    num_conf_levels = ncol(metrics_list$is)
    num_quantiles = ncol(metrics_list$quantile_score)
    
    # Define metric keys and their dimensions, aligning with `initialize_results`
    metric_dims = list(
      "Fc" = c(N, t), 
      "MIS" = c(N, num_conf_levels, t),
      "COVERAGE" = c(N, num_conf_levels, t),
      "CRPS" = c(N, t),
      "LOG" = c(N, t),
      "MSE" = c(N, t),
      "QUANTILES_SCORE" = c(N, num_quantiles, t),
      "PI_WIDTH" = c(N, num_conf_levels, t)
    )
    
    # Determine whether to use type/shrink as a sub-key
    sub_key = if (!is.null(shrink)) shrink else if (!is.null(type)) type else NULL
    
    # Initialize missing elements in results using `initialize_if_null`
    for (measure in names(metric_dims)) {
      if (!is.null(sub_key)) {
        results = initialize_if_null(results, measure, met, sub_key, shrink, metric_dims[[measure]])
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
      results[["CRPS"]][[met]][[sub_key]][, t] = metrics_list$crps
      results[["LOG"]][[met]][[sub_key]][, t] = metrics_list$log_score
      results[["MSE"]][[met]][[sub_key]][, t] = metrics_list$se
      results[["QUANTILES_SCORE"]][[met]][[sub_key]][,, t] = metrics_list$quantile_score
      results[["PI_WIDTH"]][[met]][[sub_key]][,, t] = metrics_list$pi_width
    } else {
      results[["Fc"]][[met]][, t] = dist_matrix$mean
      results[["MIS"]][[met]][,, t] = metrics_list$is
      results[["COVERAGE"]][[met]][,, t] = metrics_list$coverage
      results[["CRPS"]][[met]][, t] = metrics_list$crps
      results[["LOG"]][[met]][, t] = metrics_list$log_score
      results[["MSE"]][[met]][, t] = metrics_list$se
      results[["QUANTILES_SCORE"]][[met]][,, t] = metrics_list$quantile_score
      results[["PI_WIDTH"]][[met]][,, t] = metrics_list$pi_width
    }
    
    return(results)
  }
  
  #Function to compute and store the metrics
  
  #INPUT
  #met: a string in the vector methods = c("t_Rec", "MinT", "Base")
  #recon_results: the list in OUTPUT of the function reconcile_forecasts 
  #               (list with the parameters of the reconciled distribution)
  #actuals: vector of dimension N with the actual values for every variable
  #.scales: a matrix of dimension (N x length(confidence_levels)) or 
  #        (N x length(quantiles)) with the scaling terms for each variable 
  #        and confidence levels/quantiles
  #confidence_levels: vector of the confidence levels
  #quantiles: vector of the quantiles
  #results: the initialized empty list (OUTPUT of initialize_results function) to be filled
  #type: a string in the vector types = c("Diag", "Full") for the t-Rec method
  #shrink: a string in the vector c("shrink", "no_shrink") for the MinT method
  #j: is the index of the current rolling origin
  
  #OUTPUT
  #results: a nested list with all the metrics as in initialize_results function
  
  compute_and_store_metrics_matrix = function(met, recon_results, actuals, .scales, 
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
    
    
    pi_extremes = confidence_intervals_vectorized(.dist, confidence_levels, dist_type)
    # Compute all performance metrics/scores in a vectorized manner
    metrics = list(
      is = is_score_vectorized(.dist, actuals, confidence_levels, .scales$scale_mis, dist_type),
      coverage = coverage_vectorized(.dist, actuals, confidence_levels, dist_type),
      crps = crps_score_vectorized(.dist, actuals, dist_type, .scales$scale_crps),
      log_score = log_score_vectorized(.dist, actuals, dist_type),
      se = (.dist$mean - actuals)**2,
      quantile_score = quant_score_vectorized(.dist, actuals, quantiles, dist_type, .scales$scale_qt),
      pi_width = pi_extremes$upper_bound - pi_extremes$lower_bound
    )
    
    results = update_results_matrix(results, metrics, .dist, met, type, shrink, j)
    
    return(results)
  }
  
}

#Functions for aggregating the results (over rolling windows)
{# Helper function to aggregate metrics over rolling windows
  aggregate_metric = function(results, met, type = NULL, metric_name, n_digit, is_t_Rec, scaled, shrink) {
    r = list()
    metric_data = if (is_t_Rec) {
      results[[metric_name]][[met]][[type]]
    } else if (met == "MinT"){
      results[[metric_name]][[met]][[shrink]]
    }else{
      results[[metric_name]][[met]]
    }
    
    if (scaled == "no_scaled"){
      if (metric_name %in% c("MIS", "COVERAGE", "QUANTILES_SCORE", "PI_WIDTH")) {
        r$result = round(apply(metric_data, c(1, 2), mean), n_digit)
      } else if (metric_name == "Fc") {
        r$result = round(apply(metric_data, 1, mean), n_digit)
      } else if (metric_name %in% c("MSE", "LOG", "CRPS")){
        r$result = round(apply(metric_data, 1, mean), n_digit)
      }
    }
    
    return(r)
  }
  
  final_results = function(results, methods, types, n_digit, scaled){
    final_res = list()
    
    for (met in methods) {
      is_t_Rec = (met == "t_Rec")
      metric_names = c("CRPS", "Fc", "MSE", "LOG", "MIS", "COVERAGE", "QUANTILES_SCORE", "PI_WIDTH")
      
      if (is_t_Rec) {
        for (type in types) {
          for (metric_name in metric_names) {
            a = aggregate_metric(results, met, type, metric_name, n_digit, is_t_Rec, scaled = scaled, shrink = NULL)
            final_res[[metric_name]][[met]][[type]] = a$result
          }
        }
      } else if (met == "MinT"){
        for (shrink in c("shrink", "no_shrink")){
          for (metric_name in metric_names) {
            a = aggregate_metric(results, met, metric_name = metric_name, n_digit = n_digit, is_t_Rec = is_t_Rec, scaled = scaled, shrink = shrink)
            final_res[[metric_name]][[met]][[shrink]] = a$result
          }
        }
      }else{
        for (metric_name in metric_names) {
          a = aggregate_metric(results, met, metric_name = metric_name, n_digit = n_digit, is_t_Rec = is_t_Rec, scaled = scaled, shrink = NULL)
          final_res[[metric_name]][[met]] = a$result
        }
      }
    }
    return(final_res)
  }
}

#Main function parallelized over cores
{
  computing_results_par = function(methods, types, confidence_levels, 
                                   quantiles, name_dataset, h, 
                                   length_training, tot_windows, 
                                   specific_method){
    L = length_training + tot_windows #ts total length needed 
    
    dataset_info = get_dataset(name_dataset) #extract the data and the info about the dataset
    
    N_b = dataset_info$n_bottom #number of bottom variables
    
    N_u = dataset_info$n_upper #number of upper variables
    
    N = N_b + N_u #total number of variables
    
    A = dataset_info$agg_mat #aggregation matrix
    
    S = dataset_info$sum_mat #summation matrix
    
    Y = dataset_info$ts #ts
    
    freq = dataset_info$freq #frequency for seasonality
    
    #Seasonality test for each variable
    mask = as.logical(apply(Y, 1, function(row) nsdiffs(ts(row, frequency = freq))))
    
    globals = list(N = N, methods = methods, types = types, 
                   specific_method = specific_method, confidence_levels = confidence_levels, 
                   quantiles = quantiles, A = A, S = S, freq = freq, 
                   mask = mask, length_training = length_training, h = h, Y = Y)
    
    # ---------------------- #
    # Set Parallel Processing #
    # ---------------------- #
    plan(multisession, workers = detectCores() - 2)  # Use all available cores except two
    
    # ---------------------- #
    # Run the main with Progress Tracking #
    # ---------------------- #
    start_time = Sys.time()
    
    with_progress({
      progress = progressor(along = 1:tot_windows)  # Create progress bar
      
      results_list = future_lapply(1:tot_windows, function(j) {
        
        # Explicitly reference global variables
        N = globals$N
        methods = globals$methods
        types = globals$types
        specific_method = globals$specific_method
        confidence_levels = globals$confidence_levels
        quantiles = globals$quantiles
        A = globals$A
        S = globals$S
        freq = globals$freq
        mask = globals$mask
        length_training = globals$length_training
        h = globals$h
        Y = globals$Y
        
        # ---------------------- #
        # Progress Update
        progress()
        
        # ---------------------- #
        # Initialize Local Storage (1 window at a time)
        results_no_scaled_par = initialize_results(N, tot_windows = 1, methods, types, 
                                                   specific_method, confidence_levels, quantiles)
        
        parms = list()
        
        # Extract Training and Actual Data
        train = Y[, j:(j + length_training - 1)]
        actuals = Y[, (j + length_training):(j + length_training + h - 1)]
        
        # Compute Base Forecasts
        base_results = compute_base_forecasts(train, h, "ets", freq)
        base_fc = base_results$base_fc
        res = base_results$res 
        
        # Compute Scaling Parameters (we use one values because we don't scale since we consider the 
        # geometric mean over the time series for the scores)
        ones_list = list(
          scale_crps = rep(1, N),
          scale_mis = matrix(1, nrow = N, ncol = length(confidence_levels)),
          scale_qt = matrix(1, nrow = N, ncol = length(quantiles))
        )
        
        # Reconcile Forecasts and Compute Metrics
        for (met in methods) {
          if (met == "t_Rec") {

            for (type in types) {
              
              recon_results = reconcile_forecasts(met, A, S, h, base_fc , type, res, 
                                                  train, freq, mask, shrink = NULL)
              
              results_no_scaled_par = compute_and_store_metrics_matrix(met, recon_results, actuals, ones_list,
                                                                       confidence_levels, quantiles, results_no_scaled_par, 
                                                                       type, shrink = NULL, j = 1)
              parms[[met]][[type]] = list(mean = recon_results$mean, 
                                          scale_par = recon_results$scale_par,
                                          dof = recon_results$df )
            }
          } else if (met == "MinT") {
            for (shrink in c("shrink", "no_shrink")) {
              
              recon_results = reconcile_forecasts(met, A, S, h, base_fc , type = NULL, res, 
                                                  train, freq, mask, shrink = shrink)
  
              results_no_scaled_par = compute_and_store_metrics_matrix(met, recon_results, actuals, ones_list,
                                                                       confidence_levels, quantiles, results_no_scaled_par, 
                                                                       type = NULL, shrink = shrink, j = 1)
              
              parms[[met]][[shrink]] = list(mean = recon_results$mean, 
                                            cov = recon_results$cov)
              
            }
          } else if (met == "Base"){
            Sigma = (crossprod(res)/nrow(res))
            Sigmahat = (1 - 10^-4)*Sigma + 10^-4* diag(diag(Sigma))
            pred_dist = list(mean = base_fc, cov = Sigmahat)
            
            parms[[met]] = pred_dist
            
            
            results_no_scaled_par = compute_and_store_metrics_matrix(met, pred_dist, actuals, ones_list,
                                                                     confidence_levels, quantiles, results_no_scaled_par, 
                                                                     type = NULL, shrink = NULL, j = 1)
          }else{
            recon_results = reconcile_forecasts(met, A, S, h, base_fc , type = NULL, res, 
                                                train, freq, mask, shrink = NULL)
            
            results_no_scaled_par = compute_and_store_metrics_matrix(met, recon_results, actuals, ones_list,
                                                                     confidence_levels, quantiles, results_no_scaled_par, 
                                                                     type = NULL, shrink = NULL, j = 1)
            
          }
        }
        
        return(list(results_no_scaled_par = results_no_scaled_par,
                    base_fc = base_fc, dist_parms = parms))
      }, future.seed = TRUE)  # Ensures reproducibility
    })
    
    # ---------------------- #
    # Merge Results Back #
    # ---------------------- #
    
    final_results_no_scaled_par = initialize_results(N, tot_windows, methods, types, 
                                                     specific_method, confidence_levels, quantiles)
    
    # Iterate over the list to merge the results
    agg_results = aggregate_parallel_results(results_list, tot_windows, final_results_no_scaled_par)
    
    parms_list = Reduce(function(x, y) {
      Map(function(a, b) c(a, b), x, y)
    }, lapply(results_list, `[[`, "dist_parms"))

    # ---------------------- #
    # Reset Processing #
    # ---------------------- #
    plan(sequential)  # Reset back to single-thread processing
    
    
    out = list(result_par = agg_results, parms_list = parms_list)
    return(out)
  }
  
  
  aggregate_parallel_results = function(results_list, tot_windows, 
                                        final_results_no_scaled_par){
    for (j in 1:tot_windows) {
      for (measure in names(final_results_no_scaled_par)) {
        for (met in names(final_results_no_scaled_par[[measure]])) {
          if (is.list(final_results_no_scaled_par[[measure]][[met]])) {
            for (sub in names(final_results_no_scaled_par[[measure]][[met]])) {
              if (measure == "MIS" | measure == "COVERAGE" | measure == "QUANTILES_SCORE"| measure == "PI_WIDTH"){
                final_results_no_scaled_par[[measure]][[met]][[sub]][, ,j] = results_list[[j]]$results_no_scaled_par[[measure]][[met]][[sub]]
                
              }else{
                final_results_no_scaled_par[[measure]][[met]][[sub]][, j] = results_list[[j]]$results_no_scaled_par[[measure]][[met]][[sub]]
                
              }
            }
          } else {
            if (measure == "MIS" | measure == "COVERAGE" | measure == "QUANTILES_SCORE"| measure == "PI_WIDTH"){
              final_results_no_scaled_par[[measure]][[met]][, ,j] = results_list[[j]]$results_no_scaled_par[[measure]][[met]]
              
            }else{
              final_results_no_scaled_par[[measure]][[met]][, j] = results_list[[j]]$results_no_scaled_par[[measure]][[met]]
              
            }
          }
        }
      }
    }
    
    out = list(final_results_no_scaled_par = final_results_no_scaled_par)
    return(out)
  }
}

#Optimization for computing nu prior according to the data
{
  library(nloptr)
  
  multi_log_score_optimization = function(res, n_obs, n_var, psi){
    n_obs_effective = round(0.1*n_obs)
    inv_sigma_fun = function(x){solve(((psi * (x[1] - n_var - 1)) + (x[2]*crossprod(res))))}
    log_det_sigma_fun = function(x){sum(log(eigen(((psi * (x[1] - n_var - 1)) + 
                                                     (x[2]*crossprod(res))))$values))}
    
    log_den_tz_inv_sigma_z_fun = function(l, inv_sigma){log((1 - (((res[l, ])%*%inv_sigma%*%(res[l,])))))}
    
    objective_function = function(x){
      inv_sigma = inv_sigma_fun(x)
      
      log_det_sigma = log_det_sigma_fun(x)
      
      - ((n_obs - n_obs_effective)*(lgamma((x[1] + (x[2]*n_obs))/2)- 
                                      lgamma((x[1] + (x[2]*n_obs) - n_var)/2) -
                                      (log_det_sigma/2)) + ((x[1] + (x[2]*n_obs) - 1)/2)*
           sum(sort(vapply(1:n_obs, function(l){log_den_tz_inv_sigma_z_fun(l, inv_sigma)}, 
                           numeric(1)))[(n_obs_effective + 1):n_obs]))
    }
    
    initial_guess = c(n_var + 2, 1)
    
    lower_bounds = c(n_var + 2, 1)
    upper_bounds = c(5*n_var, 1)
    
    opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 10^-5, 
                "maxeval" = 1000)
    
    start_time = Sys.time()
    results = nloptr(x0 = initial_guess, eval_f = objective_function, 
                     lb = lower_bounds, ub = upper_bounds, opts = opts)
    end_time = Sys.time()
    
    lag = end_time - start_time
    obj_fun_eval = NULL
    return(list(results, obj_fun_eval))
  }
  
}

#Computing the geometric mean of a vector
{
  geom_mean = function(x){
    return((prod(x)**(1/length(x))))
  }
}
