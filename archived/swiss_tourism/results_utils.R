dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

#Function to save the results in the format ".csv" or ".RDS".
{
  save_nested_list_to_csv <- function(nested_list, metrics_name, dataset, quantiles, 
                                    confidence_levels, length_training, scaled,
                                    folder_name, form) {
  
  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  for (met in names(nested_list)) {
    # Handle the presence of 'type' only for 't_stud'
    if ((met == "t_Rec") | (met == "MinT")) {
      for (type in names(nested_list[[met]])) {
        # Extract the data
        result <- nested_list[[met]][[type]]
        
        # Convert to data frame for easier saving
        if (is.matrix(result) || is.array(result)) {
          .df <- as.data.frame(result, row.names = if (metrics_name %in% c("MIS", "COVERAGE")){as.character(confidence_levels)}else if (metrics_name == "QUANTILES_SCORE"){as.character(quantiles)})
        } else if (is.vector(result)) {
          .df <- setNames(data.frame(
            result,
            row.names = if (metrics_name %in% c("MIS", "COVERAGE")) {
              paste0(as.character(100*confidence_levels),"%")
            } else if (metrics_name == "QUANTILES_SCORE") {
              paste0(as.character(100*quantiles),"th")
            } else {
              NULL  
            }
          ), paste0(met,"-", type))
        } else {
          next 
        }
        
        # Create a file name dynamically and include the folder
        file_name <- paste0(folder_name, scaled, "_",metrics_name,"_",met,"_", type, "_length_training_", length_training)
        
        if (!dir.exists(folder_name)) {
          dir.create(folder_name, recursive = TRUE)
        }
        # Save the data frame to a CSV file
        if (form == "csv"){
          write.csv2(.df, paste0(file_name, ".csv"), row.names = T)
        }else if (form == "RDS"){
          saveRDS(.df, paste0(file_name, ".RDS"))
        }
        
      }
    }else {
      # No 'type' present for other methods
      # Extract the data
      result <- nested_list[[met]]
      
      # Convert to data frame for easier saving
      if (is.matrix(result) || is.array(result)) {
        .df <- as.data.frame(result, row.names = if (metrics_name %in% c("MIS", "COVERAGE")){as.character(confidence_levels)}else if (metrics_name == "QUANTILES_SCORE"){as.character(quantiles)})
      } else if (is.vector(result)) {
        .df <- setNames(data.frame(
          result,
          row.names = if (metrics_name %in% c("MIS", "COVERAGE")) {
            paste0(as.character(100*confidence_levels),"%")
          } else if (metrics_name == "QUANTILES_SCORE") {
            paste0(as.character(100*quantiles),"th")
          } else {
            NULL  
          }
        ), met)
      } else {
        next 
      }
      
      # Create a file name dynamically and include the folder
      file_name <- paste0(folder_name, scaled, "_",metrics_name,"_",met, "_length_training_", length_training)
      
      
      # Save the data frame to a CSV file
      if (form == "csv"){
        write.csv2(.df, paste0(file_name, ".csv"), row.names = T)
      }else if (form == "RDS"){
        saveRDS(.df, paste0(file_name, ".RDS"))
      }
    }
  }
  
}
}

#Table 4
{
  read_measure_to_dataframe <- function(measure, dataset, length_training,
                                      confidence_levels = NULL, quantiles = NULL, n_digit = 4) {
  # Define the folders where the measure files are stored
  measure_folder_path <- paste0("results_folder/",dataset, "/", measure, "/length_train_", length_training, "/")
  
  
  # Check if both folders exist
  if (!dir.exists(measure_folder_path)) {
    stop("The measure folder does not exist: ", measure_folder_path)
  }
  
  
  # Get a list of all CSV files related to the specified measure
  measure_pattern <- paste0(measure, "_.*_length_training_", length_training, "\\.csv")
  
  
  measure_files <- list.files(measure_folder_path, pattern = measure_pattern, full.names = TRUE)
  
  
  # Check if there are matching files for both measure
  if (length(measure_files) == 0) {
    stop("No CSV files found for the specified measure: ", measure)
  }
  
  # Initialize lists to store data for measures
  measure_data_list_all = list()
  
  for (file in measure_files) {
      method_name <- sub(paste0(measure, "_(.*)_length_training_", length_training, "\\.csv"), "\\1", basename(file))
      
      # Read the CSV file for the measure
      .df <- read.csv2(file, row.names = 1)
      
      if (measure == "COVERAGE" | measure == "MSE"){
        measure_data_list_all[[method_name]] = as.data.frame(round(colMeans(.df), n_digit))
      }else{
        measure_data_list_all[[method_name]] = as.data.frame(apply(.df, 2, geom_mean))
      }
    
  }
  
  combined_df = do.call(cbind, lapply(names(measure_data_list_all), function(method) {
    if (!is.null(measure_data_list_all[[method]])) {
      measure_values <- measure_data_list_all[[method]]
    } else {
      NA
    }
  }))
  
  # Set the column names to the method names
  colnames(combined_df) <- names(measure_data_list_all)
  
  if (measure != "COVERAGE"){
    combined_df = round(combined_df/combined_df[["no_scaled_Base"]], n_digit)
  }
  
  # Optionally, set row names based on confidence levels or quantiles
  if (!is.null(confidence_levels) && measure %in% c("MIS", "COVERAGE", "PI_WIDTH")) {
    rownames(combined_df) <- paste0(as.character(100 * confidence_levels), "%")
  } else if (!is.null(quantiles) && measure %in% c("QUANTILES_SCORE")) {
    rownames(combined_df) <- paste0(as.character(100 * quantiles), "th")
  }
  
  return(combined_df)
}
}

#Figure 7, C.12, C.13
{
  library(tsutils)
library(latex2exp)
#MCB nemenyi test for each of the metrics of interest
#data = matrix N x n_methods to compare 
MCB_test = function(dataset, metrics_name, length_training, lev){
  folder_name = paste0("results_folder/",dataset,"/",metrics_name ,"/length_train_",length_training,"/")
  
  file_name_minTS <- paste0(folder_name, "no_scaled_",metrics_name,"_MinT_shrink_length_training_", length_training, ".csv")
  
  file_name_t_stud_cov_naive = paste0(folder_name, "no_scaled_",metrics_name,"_t_Rec_Full_length_training_", length_training, ".csv")
  
  file_name_base_fc = paste0(folder_name, "no_scaled_",metrics_name,"_Base_length_training_", length_training, ".csv")
  
  v1 = read.csv2(file_name_minTS)
  v2 = read.csv2(file_name_t_stud_cov_naive)
  v5 = read.csv2(file_name_base_fc)
  
  if (metrics_name == "MIS"){
    
    if (lev == "80"){
      data = setNames(as.data.frame(cbind(v5[,-c(1,3)],v1[,-c(1,3)], 
                                          v2[,-c(1,3)])), 
                      c("Base", "MinT", "t-Rec"))
    }else if (lev == "95"){
      data = setNames(as.data.frame(cbind(v5[,-c(1,2)],v1[,-c(1,2)], 
                                          v2[,-c(1,2)])), 
                      c("Base", "MinT", "t-Rec"))
    }
    
  }else{
    
    data = setNames(as.data.frame(cbind(v5[,-1],v1[,-1], v2[,-1])), 
                    c("Base", "MinT", "t-Rec"))
    
  }
  
  if (!dir.exists(paste0(dir_path, "/Paper_plots"))) {
    dir.create(paste0(dir_path, "/Paper_plots"), recursive = TRUE)
  }
  
  pdf(if(metrics_name == "MIS"){paste0("Paper_plots/",dataset, "_", length_training,"_nemenyi_plot_",metrics_name,lev,".pdf")}else{paste0("Paper_plots/",dataset, "_", length_training,"_nemenyi_plot_",metrics_name,".pdf")}, width = 20, height = 16) 
  
  par(pch = 16, cex = 4, lwd = 2)
  results = nemenyi(data, sort = T, plottype = "vmcb", conf.level = 0.95, xlab = if (metrics_name == "MIS"){paste("Mean Rank", metrics_name, lev)}else{paste("Mean Rank", metrics_name)})
  
  dev.off()
}
}

#Figure 5, 10, 11, D.14
{
  TrainPlotScore = function(trains, metrics, dataset, n_digit, conf_levs, 
                          quants, grid_search, methods, levs){
  df = list()
  df_ = list()
  
  for (metric in metrics){
    folder_name = paste0(dir_path, "/results_folder/",dataset,"/", metric)
    
    
    if (metric == "MIS"| metric == "COVERAGE"){
      for (lev in levs){
        df[[dataset]][[paste(metric, lev)]] = matrix(0, nrow = length(trains), ncol = 3)
        
        for (i in seq_along(trains)){
          
          geom_mean_result = geom_mean_computing(dataset, metric, trains[i], 
                                                 5, conf_levs, NULL, methods)
          
          if ((metric == "MIS") & (lev == "95")){
            row1 = round(geom_mean_result[2, methods]/geom_mean_result[2, "Base"], 5)
          }else if ((metric == "MIS") & (lev == "80")){
            row1 = round(geom_mean_result[1, methods]/geom_mean_result[1, "Base"], 5)
          }else if ((metric == "COVERAGE") & (lev == "95")){
            row1 = round(geom_mean_result[2, methods], 2)
          }else if ((metric == "COVERAGE") & (lev == "80")){
            row1 = round(geom_mean_result[1, methods], 2)
          }
          names(row1) = c("MinT", "t-Rec-Diag", "t-Rec")
          df[[dataset]][[paste(metric, lev)]][i,] = row1
          
        }
        row_names = vector(length = length(trains))
        for (i in seq_along(trains)) {
          row_names[i] = paste("Train", trains[i])
        }
        df_[[dataset]][[paste(metric, lev)]] = setNames(as.data.frame(df[[dataset]][[paste(metric, lev)]], 
                                                                      row.names = row_names), names(row1))
      }
    }else{
      df[[dataset]][[metric]] = matrix(0, nrow = length(trains), ncol = 3)
      for (i in seq_along(trains)){
        geom_mean_result = geom_mean_computing(dataset, metric, trains[i], 
                                               5, conf_levs, NULL, methods)
        
        row1 = round(geom_mean_result[, methods]/geom_mean_result[, "Base"], 5)
        names(row1) = c("MinT", "t-Rec-Diag", "t-Rec")
        df[[dataset]][[metric]][i,] = row1
      }
      
      row_names = vector(length = length(trains))
      for (i in seq_along(trains)) {
        row_names[i] = paste("Train", trains[i])
      }
      df_[[dataset]][[metric]] = setNames(as.data.frame(df[[dataset]][[metric]], 
                                                        row.names = row_names), names(row1))
    } 
    
  } 
  
  scores = c("MSE", "CRPS", "MIS 80", "MIS 95")
  for (score in scores){
    data = data.frame(
      Train = trains,
      MinT = df_[[dataset]][[score]][["MinT"]],
      "t-Rec-Diag"  = df_[[dataset]][[score]][["t-Rec-Diag"]],
      "t-Rec" = df_[[dataset]][[score]][["t-Rec"]],
      check.names = F)
    
    # Reshape the data to long format
    data_long = pivot_longer(data, 
                             cols = c("MinT", "t-Rec-Diag", "t-Rec"),
                             names_to = "Method",
                             values_to = "Value")
    
    ylabel = if(score == "MSE"){
      expression(RelMSE)
    }else if (score == "CRPS"){
      expression(RelCRPS)
    }else if (score == "MIS 80"){
      expression(RelMIS^{80*"%"})
    }else if (score == "MIS 95"){
      expression(RelMIS^{95*"%"})
    }else if (score == "COVERAGE 80"){
      expression(RelCoverage^{80*"%"})
    }else if (score == "COVERAGE 95"){
      expression(RelCoverage^{95*"%"})
    }
    
    # Create the ggplot
    p = ggplot(data_long, aes(x = Train, y = Value, color = Method, fill = Method)) +
      geom_line(size = 4) +
      geom_point(size = 6) +
      labs(title = "",
           x = "Training length",
           y = ylabel,
           color = "Method") +
      scale_color_manual(values = c("MinT" = "#440154FF",  
                                    "t-Rec-Diag" = "#FDE725FF",
                                    "t-Rec" = "#29AF7F"
      )) +
      scale_fill_manual(values = c("MinT" = "#440154",
                                   "t-Rec-Diag" = "#FDE725FF",
                                   "t-Rec" = "#29AF7F")) +
      guides(color = guide_legend(override.aes = list(size = 15))) +
      scale_x_continuous(breaks = data$Train)+
      theme_minimal()+
      theme(
        plot.title = element_text(size = 80, face = "bold"),
        axis.title.y = element_text(size = 80, margin = margin(r = 65)),
        axis.title.x = element_text(size = 80, margin = margin(t = 65)),
        axis.text = element_text(size = 50),
        legend.title = element_text(size = 50, face = "bold"),
        legend.text = element_text(size = 50)
      )
    
    if (score == "COVERAGE 80"){
      p = p + 
        geom_hline(yintercept = 0.8, color = "red", linewidth = 1)
      
    }else if (score == "COVERAGE 95"){
      p = p + 
        geom_hline(yintercept = 0.95, color = "red", linewidth = 1)
      
    }
    
    if (!dir.exists(paste0(dir_path, "/Paper_plots"))) {
      dir.create(paste0(dir_path, "/Paper_plots"), recursive = TRUE)
    }
    
    filename = paste0(dir_path, "/Paper_plots/", dataset, "_",score,"_diag.pdf")
    height = 15
    
    ggsave(filename, plot = p, device = "pdf", width = 25, height = height)
    
  }
  
}
}

#Figure 6
{
  HistogramsCoverage = function(measure, dataset, length_training,
                              confidence_levels = NULL, quantiles = NULL, n_digit = 4, lev){
  library(forcats)
  
    # Define the folders where the measure files are stored
    measure_folder_path <- paste0("results_folder/",dataset, "/", measure, "/length_train_", length_training, "/")
    
    
    # Check if both folders exist
    if (!dir.exists(measure_folder_path)) {
      stop("The measure folder does not exist: ", measure_folder_path)
    }
    
    
    # Get a list of all CSV files related to the specified measure
    measure_pattern <- paste0(measure, "_.*_length_training_", length_training, "\\.csv")
    
    
    measure_files <- list.files(measure_folder_path, pattern = measure_pattern, full.names = TRUE)
    
    
    # Check if there are matching files for both measure
    if (length(measure_files) == 0) {
      stop("No CSV files found for the specified measure: ", measure)
    }
    
    # Initialize lists to store data for measures
    measure_data_list_all = list()
    
    for (file in measure_files) {
      method_name <- sub(paste0(measure, "_(.*)_length_training_", length_training, "\\.csv"), "\\1", basename(file))
      
      # Read the CSV file for the measure
      .df <- read.csv2(file, row.names = 1)
      
      measure_data_list_all[[method_name]] = as.data.frame(round(.df, n_digit))
      
    }
  
    if (lev == "80"){
      booo = do.call(cbind, lapply(measure_data_list_all, function(x) x[, 1]))[, c("no_scaled_MinT_shrink", 
                                                                                   "no_scaled_t_Rec_Diag",
                                                                                   "no_scaled_t_Rec_Full")]
      
    }else if (lev == "95"){
      booo = do.call(cbind, lapply(measure_data_list_all, function(x) x[, 2]))[, c("no_scaled_MinT_shrink", 
                                                                                   "no_scaled_t_Rec_Diag",
                                                                                   "no_scaled_t_Rec_Full")]
      
    }
  colnames(booo) = c("MinT", "t-Rec-Diag", "t-Rec")
  booo = booo[, c("MinT", "t-Rec")]
  booo = round(booo, 2)
  
  #This is the order as in the original dataset
  series_names  <- c(
    "CH", "ZH", "BE", "LU", "UR", "SZ", "OW", "NW", "GL", "ZG",
    "FR", "SO", "BS", "BL", "SH", "AR", "AI", "SG", "GR", "AG",
    "TG", "TI", "VD", "VS", "NE", "GE", "JU"
  )
  
  booo = as.data.frame(booo)
  
  library(tibble)
  
  data <- tibble(
    Series = series_names,
    MinT = booo$MinT,
    tRec = booo$"t-Rec"
  )
  
  # Pivot to long format
  long_data <- data %>%
    dplyr::select(Series, MinT, tRec) %>%
    pivot_longer(-Series, names_to = "Method", values_to = "Value") %>%
    mutate(Method = recode(Method,
                           "MinT" = "MinT",
                           "tRec" = "t-Rec"))
  
  new_order = c("CH", sort(setdiff(long_data$Series, "CH")))
  
  long_data$Series <- factor(long_data$Series, levels = new_order)
  
  if(lev == "80"){ yintercept = 0.80}else if(lev == "95"){ yintercept = 0.95}
  
  if(lev == "80"){breaks = c(0.80)}else if(lev == "95"){breaks = c(0.95)}
  
  c = ggplot(long_data, aes(x = Series, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("MinT" = "#440154FF", "t-Rec" = "#29AF7F")) +
    geom_hline(yintercept = yintercept, 
               color = "red", linewidth = 3)+
    scale_y_continuous(
      breaks = breaks
    )+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 30, face = "bold"),
          axis.title.y = element_text(size = 50, margin = margin(r = 65)),
          axis.title.x = element_text(size = 50, margin = margin(t = 65)),
          axis.text = element_text(size = 35),
          legend.title = element_text(size = 50, face = "bold"),
          legend.text = element_text(size = 50)) +
    labs(title = "",
         y = paste("Coverage", lev, "%"),
         x = "",
         fill = "Method")
  
  if (!dir.exists(paste0(dir_path, "/Paper_plots"))) {
    dir.create(paste0(dir_path, "/Paper_plots"), recursive = TRUE)
  }
  
  filename = paste0(dir_path, "/Paper_plots/",dataset,"_detailed_coverage_",lev,".pdf")
  ggsave(filename, plot = c, device = "pdf", width = 25, height = 11)
  
}
}

#Figure 9
{
  ScoresAllVariables = function(dataset, metric_name, length_training, n_digit, conf_levs, 
                              quants, methods, lev, form){
  substrings = c("Base", methods)
  pattern <- paste(substrings, collapse = "|")
  
  
  # Define the folders where the measure files are stored
  measure_folder_path = paste0("results_folder/",dataset, "/", metric_name, "/length_train_", length_training, "/")
  
  
  if (form == "RDS"){
    # Check if both folders exist
    if (!dir.exists(measure_folder_path)) {
      stop("The measure folder does not exist: ", measure_folder_path)
    }
    
    
    # Get a list of all CSV files related to the specified measure
    measure_pattern = paste0("no_scaled_",metric_name, "_.*_length_training_", length_training, "\\.RDS")
    
    
    measure_files = list.files(measure_folder_path, pattern = measure_pattern, full.names = TRUE)
    
    
    # Check if there are matching files for both measure
    if (length(measure_files) == 0) {
      stop("No CSV files found for the specified measure: ", metric_name)
    }
    
    # Initialize lists to store data for measures
    measure_data_list = list()
    
    measure_files_ = grep(pattern, measure_files, value = TRUE)
    df_list = list()
    for (file in measure_files_) {
      # Extract the method name from the file name
      method_name = sub(paste0("no_scaled_",metric_name, "_(.*)_length_training_", length_training, "\\.RDS"), "\\1", basename(file))
      
      # Read the CSV file for the measure
      .df = readRDS(file)
      
      df_list[[method_name]] = .df
    }
    
    even = seq(2, length(df_list[[method_name]]), by = 2)
    
    if (lev == "80"){
      df_complete = as.data.frame(do.call(cbind, lapply(df_list, 
                                                        function(x) x[, -even])))
      
    }else if (lev == "95"){
      df_complete = as.data.frame(do.call(cbind, lapply(df_list, 
                                                        function(x) x[, even])))
    }
    
    if (dataset == "Swiss_tourism"){
      series_names  <- c(
        "CH", "ZH", "BE", "LU", "UR", "SZ", "OW", "NW", "GL", "ZG",
        "FR", "SO", "BS", "BL", "SH", "AR", "AI", "SG", "GR", "AG",
        "TG", "TI", "VD", "VS", "NE", "GE", "JU"
      )
      rownames(df_complete) = series_names
    }
    
    colnames(df_complete) = c(rep("Base", length(df_list[[method_name]])/length(conf_levs)), 
                              rep("MinT", length(df_list[[method_name]])/length(conf_levs)), 
                              rep("t-Rec", length(df_list[[method_name]])/length(conf_levs)))
    
    
    
  }else if (form == "CSV"){
    
    # Check if both folders exist
    if (!dir.exists(measure_folder_path)) {
      stop("The measure folder does not exist: ", measure_folder_path)
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
    df_list = list()
    for (file in measure_files_) {
      # Extract the method name from the file name
      method_name = sub(paste0("no_scaled_",metric_name, "_(.*)_length_training_", length_training, "\\.csv"), "\\1", basename(file))
      
      # Read the CSV file for the measure
      .df = read.csv2(file, row.names = 1)
      
      df_list[[method_name]] = .df
    }
    if (lev == "80"){
      df_complete = as.data.frame(do.call(cbind, lapply(df_list, 
                                                        function(x) x[, 1])))
      
    }else if (lev == "95"){
      df_complete = as.data.frame(do.call(cbind, lapply(df_list, 
                                                        function(x) x[, 2])))
    }
    if (dataset == "Swiss_tourism"){
      series_names  <- c(
        "CH", "ZH", "BE", "LU", "UR", "SZ", "OW", "NW", "GL", "ZG",
        "FR", "SO", "BS", "BL", "SH", "AR", "AI", "SG", "GR", "AG",
        "TG", "TI", "VD", "VS", "NE", "GE", "JU"
      )
      rownames(df_complete) = series_names
    }
    colnames(df_complete) = c("Base", "MinT", "t-Rec")
    
  }
  
  return(df_complete)
}


PredictionIntervalOverRolling = function(dataset, metric_name, 
                                         length_training, n_digit, 
                                         conf_levs, quants, 
                                         methods, lev, form, index_top, index_bot){
  library(tidyr)
  library(patchwork)
  
  
  df_complete = ScoresAllVariables(dataset, metric_name, 
                                   length_training, n_digit, 
                                   conf_levs, quants, 
                                   methods, lev, form) 
  
  even = seq(2, ncol(df_complete)*length(conf_levs)/3, by = 2)
  
  measure_folder_path = paste0("/results_folder/",dataset,"/Fc/length_train_",length_training,"/")
  file_base = paste0(dir_path,measure_folder_path,"no_scaled_Fc_Base_length_training_",length_training,".RDS")
  r = readRDS(file_base)
  incoherence = r[index_top,] - colSums(r[index_bot, ])
  
  measure_folder_path = paste0("/results_folder/",dataset,"/dist_parms/")
  file_base = paste0(dir_path,measure_folder_path,dataset,"_", length_training,"_dist_parms.RDS")
  r = readRDS(file_base)
  sd_base = (lapply(r$Base[even], function(x) x[index_top,index_top]))
  inco = abs((incoherence)/sqrt(unlist((sd_base))))
  
  df_upper <- data.frame(
    inco = t(inco),
    "t-Rec" = t(df_complete[index_top, which(names(df_complete) == "t-Rec")]/df_complete[index_top, which(names(df_complete) == "Base")]), 
    MinT = t(df_complete[index_top, which(names(df_complete) == "MinT")] / df_complete[index_top, which(names(df_complete) == "Base")]),
    check.names = F
  )
  
  colnames(df_upper) = c("inco", "t-Rec", "MinT")
  
  df_bottom <- data.frame(
    inco = t(inco),
    "t-Rec" = apply(t(df_complete[index_bot, which(names(df_complete) == "t-Rec")]/df_complete[index_bot, which(names(df_complete) == "Base")]), 1, geom_mean),
    MinT = apply(t(df_complete[index_bot, which(names(df_complete) == "MinT")] / df_complete[index_bot, which(names(df_complete) == "Base")]), 1, geom_mean),
    check.names = F
  )
  colnames(df_bottom) = c("inco", "t-Rec", "MinT")
  
  df_long <- df_upper %>%
    pivot_longer(cols = c(`t-Rec`, MinT), names_to = "Method", values_to = "Value")
  
  df_long_bot <- df_bottom %>%
    pivot_longer(cols = c(`t-Rec`, MinT), names_to = "Method", values_to = "Value")
  
  
  y_lim_range <- range(c(df_upper$`t-Rec`, df_upper$MinT, df_bottom$`t-Rec`, df_bottom$MinT, 1.2))
  
  main_bot = ggplot(df_long_bot, aes(x = inco, y = Value, color = Method)) +
    geom_point(size = 4, alpha = 0.3, show.legend = F) +
    geom_smooth(
      data = df_long_bot %>% filter(Method == "t-Rec"),
      aes(x = inco, y = Value), 
      method = "lm",
      se = FALSE,
      linewidth = 3,
      color = "#29AF7F",
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 1, color = "#E69F00", linewidth = 3, show.legend = F) +
    coord_cartesian(xlim = c(0, max(inco)), ylim = y_lim_range) +
    scale_y_continuous(
      position = "right",        # Move y-axis to the right
      breaks = round(c(seq(y_lim_range[1], y_lim_range[2], 
                           by = 0.2), 1),1)  # Customize breaks to include 1
    ) +
    theme_minimal() +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154" )) +
    labs(x = "Relative incoherence", y = "Relative 95% prediction interval width", title = "Bottom time series") +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey80"),
      axis.title = element_text(size = 50),
      axis.title.x = element_text(size = 50, margin = margin(t = 40)),
      axis.text = element_text(size = 40),
      axis.title.y.right =  element_text(size = 50, vjust = 5),
      axis.text.y = element_text(size = 40),
      legend.text = element_text(size = 50),
      legend.title = element_text(size = 50, face = "bold"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  # Density for t-Rec (already done)
  d_tRec <- density(df_bottom$`t-Rec`)
  density_tRec <- data.frame(
    y = d_tRec$x,
    density = d_tRec$y,
    Method = "t-Rec"
  )
  density_tRec$neg_density <- -density_tRec$density
  
  # Density for MinT
  d_minT <- density(df_bottom$MinT)
  density_minT <- data.frame(
    y = d_minT$x,
    density = d_minT$y,
    Method = "MinT"
  )
  density_minT$neg_density <- -density_minT$density
  
  density_all <- rbind(density_tRec, density_minT)
  library(dplyr)
  
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
  
  y_density_bot <- ggplot(density_polygons, aes(x = neg_density, y = y, group = Method, fill = Method, color = Method)) +
    geom_polygon(alpha = 0.3, size = 3, show.legend = F) +
    geom_hline(yintercept = 1, color = "#E69F00", linewidth = 3, show.legend = F) +
    scale_y_continuous(limits = y_lim_range) +  # match main plot y axis
    
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    scale_fill_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154")) +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154"))
  
  
  main_up = ggplot(df_long, aes(x = inco, y = Value, color = Method)) +
    geom_point(size = 4, alpha = 0.3, show.legend = F) +
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
    coord_cartesian(xlim = c(0, max(inco)), ylim = c(min(min(df_upper[which(colnames(df_upper) == "t-Rec")]), 
                                                         min(df_upper[which(colnames(df_upper) == "MinT")])), 
                                                     max(max(df_upper[which(colnames(df_upper) == "t-Rec")]), 
                                                         max(df_upper[which(colnames(df_upper) == "MinT")])))) +
    scale_y_continuous(
      position = "right",        # Move y-axis to the right
      breaks = round(c(seq(min(min(df_upper[which(colnames(df_upper) == "t-Rec")]), 
                               min(df_upper[which(colnames(df_upper) == "MinT")])), 
                           max(max(df_upper[which(colnames(df_upper) == "t-Rec")]), 
                               max(df_upper[which(colnames(df_upper) == "MinT")])), 
                           by = 0.2), 1),1)  # Customize breaks to include 1
    ) +
    theme_minimal() +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154" )) +
    labs(x = "Relative incoherence", y = "", title = "Upper time series") +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey80"),
      axis.title = element_text(size = 50),
      axis.title.x = element_text(size = 50, margin = margin(t = 40)),
      axis.text = element_text(size = 40),
      axis.title.y.right =  element_text(size = 50, vjust = 5),
      axis.text.y = element_text(size = 40),
      legend.text = element_text(size = 50),
      legend.title = element_text(size = 50, face = "bold"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  # Density for t-Rec (already done)
  d_tRec <- density(df_upper$`t-Rec`)
  density_tRec <- data.frame(
    y = d_tRec$x,
    density = d_tRec$y,
    Method = "t-Rec"
  )
  density_tRec$neg_density <- -density_tRec$density
  
  # Density for MinT
  d_minT <- density(df_upper$MinT)
  density_minT <- data.frame(
    y = d_minT$x,
    density = d_minT$y,
    Method = "MinT"
  )
  density_minT$neg_density <- -density_minT$density
  
  density_all <- rbind(density_tRec, density_minT)
  library(dplyr)
  
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
  
  y_density_up <- ggplot(density_polygons, aes(x = neg_density, y = y, group = Method, fill = Method, color = Method)) +
    geom_polygon(alpha = 0.3, size = 3, show.legend = F) +
    geom_hline(yintercept = 1, color = "#E69F00", linewidth = 3, show.legend = F) +
    scale_y_continuous() +  # match main plot y axis
    
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    scale_fill_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154")) +
    scale_color_manual(values = c("t-Rec" = "#29AF7F", "MinT" = "#440154"))
  
  #check they are on the same scale
  y_limits <- range(c(df_upper$`t-Rec`, df_upper$MinT, df_bottom$`t-Rec`, df_bottom$MinT, 1))
  
  main_up <- main_up + coord_cartesian(ylim = y_limits)
  y_density_up <- y_density_up + scale_y_continuous(limits = y_limits)
  
  main_bot <- main_bot + coord_cartesian(ylim = y_limits)
  y_density_bot <- y_density_bot + scale_y_continuous(limits = y_limits)
  
  p = (y_density_up | main_up | y_density_bot | main_bot) + plot_layout(widths = c(2, 3, 2, 3, 0.5))
  
  if (!dir.exists(paste0(dir_path, "/Paper_plots/"))) {
    dir.create(paste0(dir_path, "/Paper_plots/"), recursive = TRUE)
  }
  
  ggsave(paste0(dir_path, "/Paper_plots/",dataset,"_scatter_plot_",metric_name,"_inco_density_up_bot_length_training_",length_training,".pdf" ),
         plot = p, device = "pdf", width = 25, height = 14)
  
}
}
