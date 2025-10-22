dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

library(GGMncv)
library(bayesRecon)
library(MASS)
library(forecast)
library(scoringRules)
library(ggplot2)

source("utils_sim.R")
source("utils.R")
source("results_utils.R")

##############################################################################################################

###########################################################
#                 Set parameters                          #
###########################################################

set.seed(27)

h = 1 #forecast horizon

Train = 12 #training length

N_b = 2 #bottom

N = 3 #total time series

L = Train + h #length of the generated time series

N_exp = 1000 #number of independent experiments


confidence_levels = c(0.8, 0.95)
quantiles = c(0.9, 0.975)

sigma_err_trend1 = 2
sigma_err_trend2 = 0.007
sigma_err_seas = 7

freq = 4 #quartely ts

gen_model_par = list(sigma_err_trend1 = 2,
                     sigma_err_trend2 = 0.007,
                     sigma_err_seas = 7, seasons = freq)

A = rbind(c(1,1))

err_cov = matrix(c(5, 3, 
                   3, 4), ncol = N_b)

S = rbind(A, diag(rep(1, N_b)))
  

methods = c("t_Rec", "MinT", "Base")
type = c("Full")


#Just formal to reuse already existing code
{ mask = rep(F, N)}


########################################################
#              Initialization                          #
########################################################
incoherence = list()

recon_results_mint = list()
recon_results_t_student = list()
pred_dist_base = list()

results_no_scaled = initialize_results_sim(N, N_exp, methods, type, "t_Rec",
                                           confidence_levels, quantiles)


######################################################
#                   Simulations                      #
######################################################
{pb = txtProgressBar(min = 0, max = N_exp, style = 3)
for (j in 1:N_exp){
  
  setTxtProgressBar(pb, j)
  
  
  #Simulating ts from an ARIMA models with gen_model_par parameters 
  #and err_cov the covariance matrix of errors
  Y = simulate_trend_seasonal_ts(L, N_b, A, gen_model_par, err_cov)
  
  
  train = Y[,1:(L-h)]
  actuals = Y[,(L-h+1):L]
  
  # Compute Base Forecasts
  base_results = compute_base_forecasts(train, h, "ets", freq)
  base_fc = base_results$base_fc
  res = base_results$res 
  
  Sigmahat = crossprod(res) / nrow(res)
  Sigmahat = (1 - 10^-4)*Sigmahat + 10^-4*diag(diag(Sigmahat))
  
  ones_list = list(
    scale_crps = rep(1, N),
    scale_mis = matrix(1, nrow = N, ncol = length(confidence_levels)),
    scale_qt = matrix(1, nrow = N, ncol = length(quantiles))
  )
  
  for (met in methods){
    if (met == "t_Rec"){
      recon_results_t_student[[j]] = reconcile_forecasts(met, A, S, h, base_fc , type, res, 
                                                         train, freq, mask, shrink = NULL)
      
      inco_var_results = incoherenceVSreconciledVariance(recon_results_t_student[[j]]$posterior.psi, 
                                                         recon_results_t_student[[j]]$posterior.nu, 
                                                         base_fc, A)
      
      incoherence[[j]] = inco_var_results$inco
      
      results_no_scaled = compute_and_store_metrics_matrix_sim(met, recon_results_t_student[[j]], actuals, ones_list,
                                                               confidence_levels, quantiles, results_no_scaled, 
                                                               type, shrink = NULL, j)
      
      
    }else if (met == "MinT"){
      recon_results_mint[[j]] = reconcile_forecasts(met, A, S, h, base_fc , type = NULL, res, 
                                                    train, freq, mask, shrink = "shrink")
      
      results_no_scaled = compute_and_store_metrics_matrix_sim(met, recon_results_mint[[j]], actuals, ones_list,
                                                               confidence_levels, quantiles, results_no_scaled, 
                                                               type = NULL, shrink = "shrink", j)
    }else if (met == "Base"){
      pred_dist_base[[j]] = list(mean = base_fc, cov = Sigmahat)
      
      results_no_scaled = compute_and_store_metrics_matrix_sim(met, pred_dist_base[[j]], actuals, ones_list,
                                                               confidence_levels, quantiles, results_no_scaled, 
                                                               type = NULL, shrink = NULL, j)
    }
    
  }
  
}

close(pb)}

#Intervals width
{sharpness = list()
  for (met in methods){
    if (met == "Base"){
      Left = results_no_scaled$Extremes_Conf_Int_L[[met]]
      Right = results_no_scaled$Extremes_Conf_Int_R[[met]]
    }else{
      Left = results_no_scaled$Extremes_Conf_Int_L[[met]][[1]]
      Right = results_no_scaled$Extremes_Conf_Int_R[[met]][[1]]
    }
    sharpness[[met]] = Right - Left
  }}


#########################################################
#                 Saving the workspace                  #
#########################################################

if (!dir.exists(paste0(dir_path, "/results_folder/Simulations/"))) {
  dir.create(paste0(dir_path, "/results_folder/Simulations/"), recursive = TRUE)
}

file_name = paste0(dir_path, "/results_folder/Simulations/length_", L, "_", N, 
                   "_variables_", N_b,"_bottom_", N_exp,"exp_ets_seas_trend_min_hierarchy.RData")
save.image(file_name)

  


##########################################################
#              Visualization of results                  #
##########################################################
#Run the following to get the results in the paper.
#The values in the Table 4 will be displayed in the Console, while 
#all the plots will be saved in the "Paper_plots" folder.


{
  #load the workspace
dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

name = paste0(dir_path, "/results_folder/Simulations/length_13_3_variables_2_bottom_1000exp_ets_seas_trend_min_hierarchy.RData")
load(name)


#Table 1: Relative Interval Width with respect to the Base method for MinT and t-Rec
{
  new_sharp = lapply(sharpness, function(x) round(apply(x, c(1,2), mean), 2))
  levs = c("80%", "95%")
  variables = c("U", "B1", "B2")
  
  result_sharp = list()
  result_sharp$t_Rec = round(new_sharp$t_Rec/new_sharp$Base, 2)
  result_sharp$MinT = round(new_sharp$MinT/new_sharp$Base, 2)
  
  rownames(result_sharp$t_Rec) = variables
  colnames(result_sharp$t_Rec) = levs
  rownames(result_sharp$MinT) = variables
  colnames(result_sharp$MinT) = levs
  
  print(result_sharp)

}

#Figure 4: Plot the distribution at two different incoherence level
{
  
#Low incoherence
index = 976
x_vals = seq(-60, 100, length.out = 10000)
normal_density_mints = dnorm(x_vals, mean = recon_results_mint[[index]]$mean[variable],
                       sd = sqrt(recon_results_mint[[index]]$cov[variable]))
normal_density_base = dnorm(x_vals, mean = pred_dist_base[[index]]$mean[variable], 
                       sd = sqrt(pred_dist_base[[index]]$cov[variable, variable]))


mean_BU = S%*%pred_dist_base[[index]]$mean[2:3]
cov_BU = S %*%  pred_dist_base[[index]]$cov[2:3, 2:3] %*% t(S)
normal_density_BU = dnorm(x_vals, mean = mean_BU[variable], 
                            sd = sqrt(cov_BU[variable, variable]))

t_student_density = dt((x_vals - recon_results_t_student[[index]]$mean[variable])/
                         sqrt(recon_results_t_student[[index]]$scale_par[variable, variable]),
                       df = recon_results_t_student[[index]]$df)/
  sqrt(recon_results_t_student[[index]]$scale_par[variable, variable])


df <- data.frame(
  x = rep(x_vals, 4),
  density = c(normal_density_base, t_student_density, normal_density_BU, normal_density_mints),
  Method = rep(c("Base", "t-Rec", "BU", "MinT"), each = length(x_vals))
)



# Plot using ggplot2
p = ggplot(df, aes(x = x, y = density)) +
  # Area fill only for MinT and t-Rec
  geom_area(
    data = subset(df, Method == "MinT"),
    aes(fill = Method),
    alpha = 0.3,
    position = "identity",
    show.legend = F
  ) +
  geom_area(
    data = subset(df, Method == "t-Rec"),
    aes(fill = Method),
    alpha = 0.3,
    position = "identity",
    show.legend = F
  ) +
  # Contour lines for all methods
  geom_line(aes(color = Method), size = 3) +
  
  # Manual color scales
  scale_fill_manual(values = c(
    "MinT"  = "#440154",
    "t-Rec" = "#29AF7F"
  )) +
  scale_color_manual(values = c(
    "MinT"  = "#440154",
    "t-Rec" = "#29AF7F",
    "BU"    = "#FDE725",
    "Base"  = "#3B528B"
  )) +
  
  coord_cartesian(xlim = c(25, 80)) +
  theme_minimal() +
  labs(title = "Low incoherence", x = "", y = "") +
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_blank()
  )

#High incoherence
index = 292
x_vals = seq(-60, 100, length.out = 10000)
normal_density_mints = dnorm(x_vals, mean = recon_results_mint[[index]]$mean[variable],
                             sd = sqrt(recon_results_mint[[index]]$cov[variable]))
normal_density_base = dnorm(x_vals, mean = pred_dist_base[[index]]$mean[variable], 
                            sd = sqrt(pred_dist_base[[index]]$cov[variable, variable]))


mean_BU = S%*%pred_dist_base[[index]]$mean[2:3]
cov_BU = S %*%  pred_dist_base[[index]]$cov[2:3, 2:3] %*% t(S)
normal_density_BU = dnorm(x_vals, mean = mean_BU[variable], 
                          sd = sqrt(cov_BU[variable, variable]))

t_student_density = dt((x_vals - recon_results_t_student[[index]]$mean[variable])/
                         sqrt(recon_results_t_student[[index]]$scale_par[variable, variable]),
                       df = recon_results_t_student[[index]]$df)/
  sqrt(recon_results_t_student[[index]]$scale_par[variable, variable])


df <- data.frame(
  x = rep(x_vals, 4),
  density = c(normal_density_base, t_student_density, normal_density_BU, normal_density_mints),
  Method = rep(c("Base", "t-Rec", "BU", "MinT"), each = length(x_vals))
)



# Plot using ggplot2
q = ggplot(df, aes(x = x, y = density)) +
  # Area fill only for MinT and t-Rec
  geom_area(
    data = subset(df, Method == "MinT"),
    aes(fill = Method),
    alpha = 0.3,
    position = "identity",
    show.legend = F
  ) +
  geom_area(
    data = subset(df, Method == "t-Rec"),
    aes(fill = Method),
    alpha = 0.3,
    position = "identity",
    show.legend = F
  ) +
  # Contour lines for all methods
  geom_line(aes(color = Method), size = 3) +
  
  # Manual color scales
  scale_fill_manual(values = c(
    "MinT"  = "#440154",
    "t-Rec" = "#29AF7F"
  )) +
  scale_color_manual(values = c(
    "MinT"  = "#440154",
    "t-Rec" = "#29AF7F",
    "BU"    = "#FDE725",
    "Base"  = "#3B528B"
  )) +
  
  coord_cartesian(xlim = c(-55, 25)) +
  theme_minimal() +
  labs(title = "High incoherence", x = "", y = "") +
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_blank()
  )


if (!dir.exists(paste0(dir_path, "/Paper_plots"))) {
  dir.create(paste0(dir_path, "/Paper_plots"), recursive = TRUE)
}

ggsave(paste0(dir_path, "/Paper_plots/Simulation_low_incoherence.pdf" ),
       plot = p, device = "pdf", width = 10, height = 10)
ggsave(paste0(dir_path, "/Paper_plots/Simulation_high_incoherence.pdf" ),
       plot = q, device = "pdf", width = 10, height = 10)
}

#Figure 3: Relative 95% prediction interval width against incoherence
{ 
  # Prepare data
  library(ggExtra)
  library(tidyr)
  library(patchwork)
  library(dplyr)
  
  inco = (abs(unlist(incoherence)/unlist(lapply(pred_dist_base, 
                                                function(x) sqrt(x$cov[variable, variable])))))
  mask = inco < 3
  
  df_ <- data.frame(
    inco = inco[mask],
    "t-Rec" = (sharpness$t_Rec[1, 2, ] / sharpness$Base[1, 2, ])[mask],
    check.names = FALSE
  )
  
  df_$MinT <- (sharpness$MinT[1, 2, ] / sharpness$Base[1, 2, ])[mask]
  
  df_long <- df_ %>%
    pivot_longer(cols = c(`t-Rec`, MinT), names_to = "Method", values_to = "Value")
  
  main <- ggplot(data = df_long, aes(x = inco, y = Value, color = Method)) +
    geom_point(aes(color = Method), size = 4, alpha = 0.3, show.legend = FALSE) +
    geom_smooth(
      data = df_long %>% filter(Method == "t-Rec"),
      aes(x = inco, y = Value), 
      method = "lm",
      se = FALSE,
      linewidth = 3,
      color = "#29AF7F",
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 1, color = "#E69F00", linewidth = 3, show.legend = F) +
    coord_cartesian(xlim = c(0, 3), ylim = c(min(min(df_$`t-Rec`), min(df_$MinT)), 
                                                     max(max(df_$`t-Rec`), max(df_$MinT)))) +
    scale_y_continuous(
      position = "right", 
      breaks = round(c(seq(min(min(df_$`t-Rec`), min(df_$MinT)), max(max(df_$`t-Rec`), 
                                                             max(df_$MinT)), 
                   by = 0.2), 1),1)  # Customize breaks to include 1
    ) +
    theme_minimal() +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154" )) +
    labs(x = "Relative incoherence", y = "Relative 95% prediction interval width") +
    theme(
      plot.title = element_text(size = 30, face = "bold"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey80"),
      axis.title = element_text(size = 40),
      axis.title.x = element_text(size = 40, margin = margin(t = 40)),
      axis.text = element_text(size = 30),
      axis.title.y.right = element_text(size = 40, vjust = 5),
      axis.text.y = element_text(size = 30),
      legend.text = element_text(size = 40),
      legend.title = element_text(size = 40, face = "bold"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  # Density for t-Rec (already done)
  d_tRec <- density(df_$`t-Rec`)
  density_tRec <- data.frame(
    y = d_tRec$x,
    density = d_tRec$y,
    Method = "t-Rec"
  )
  density_tRec$neg_density <- -density_tRec$density
  
  # Density for MinT
  d_minT <- density(df_$MinT)
  density_minT <- data.frame(
    y = d_minT$x,
    density = d_minT$y,
    Method = "MinT"
  )
  density_minT$neg_density <- -density_minT$density
  
  density_all <- rbind(density_tRec, density_minT)
  
  density_polygons <- density_all %>%
    group_by(Method) %>%
    group_split() %>%
    lapply(function(df) {
      rbind(
        data.frame(y = min(df$y), neg_density = 0, Method = df$Method[1]),
        df[, c("y", "neg_density", "Method")],
        data.frame(y = max(df$y), neg_density = 0, Method = df$Method[1])
      )
    }) %>%
    do.call(rbind, .)
  
  y_density <- ggplot(density_polygons, aes(x = neg_density, y = y, group = Method, fill = Method, color = Method)) +
    geom_polygon(alpha = 0.3, size = 2, show.legend = F) +
    geom_hline(yintercept = 1, color = "#E69F00", linewidth = 3, show.legend = F) +
    scale_y_continuous() +  # match main plot y axis
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    
    theme_void() +
    labs(y = "") +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    scale_fill_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154")) +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154"))
  
  y_limits <- range(c(df_$`t-Rec`, df_$MinT))
  
  main <- main + coord_cartesian(ylim = y_limits)
  y_density <- y_density + scale_y_continuous(limits = y_limits)
  
  p = (y_density | main) + plot_layout(widths = c(3, 2, 1))
  
  ggsave(paste0(dir_path, "/Paper_plots/Simulation_relative_PI_width_VS_incoherence.pdf" ),
         plot = p, device = "pdf", width = 20, height = 11)
}
}
