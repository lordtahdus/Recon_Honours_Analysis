dir_path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

# Load required libraries
library(dplyr)
library(bayesRecon)
library(forecast)
library(scoringRules)
library(mvtnorm)
library(profvis)
library(tsibble)
library(ggplot2)
library(tseries)
library(fpp3)

# Parallelization libraries
library(forecast)
library(future.apply)
library(progressr)
library(parallel)

#Source utilities
source("utils.R")
source("results_utils.R")


########################## Initialize Parameters and Data ####################################

#"Australian_tourism_no_zone" corresponds to "Australian tourism - Q" in the paper;
#"Australian_tourism_zone" corresponds to "Australian tourism - M" in the paper;
#"Swiss_tourism" corresponds to "Swiss tourism" in the paper.

#The following vector contains the datasets used in the paper, if you want to compute results for only one of them
#you can just remove the others
datasets = c("Swiss_tourism", "Australian_tourism_zone","Australian_tourism_no_zone")


# Confidence levels and quantiles
confidence_levels = c(0.8, 0.95)
quantiles = c(0.9, 0.975)

#These are the methods tested and the types of prior for t-Rec computed
methods = c("t_Rec", "MinT", "Base")
types = c("Diag", "Full")

h = 1 #forecast horizon

for (name_dataset in datasets){
  #These are the choices we made for the training lengths in the paper, you can change them
  if ((name_dataset == "Swiss_tourism")|(name_dataset == "Swiss_tourism_Bigger")){
    train = c(seq(30, 40, 2), seq(45, 60, 5)) 
  }else if (name_dataset == "Australian_tourism_zone"){
    train = seq(30, 125, 5) 
  }else if (name_dataset == "Australian_tourism_no_zone"){
    train = seq(15, 50, 5)
  }
  
  #According to the number of training lengths tested we made a choice on the number of rolling windows
  if ((name_dataset == "Swiss_tourism")|(name_dataset == "Swiss_tourism_Bigger")){
    tot_windows = rep(100, length(train)) 
  }else if (name_dataset == "Australian_tourism_zone"){
    tot_windows = rep(100, length(train))
  }else if (name_dataset == "Australian_tourism_no_zone"){
    tot_windows = (80 - train) #Maximum number of observations for this dataset is 80
  }
  
  
  for (length_training in train){
    main_function_results = computing_results_par(methods, types, confidence_levels, 
                                                  quantiles, name_dataset, h, 
                                                  length_training , tot_windows[match(length_training, train)], 
                                                  specific_method = "t_Rec")
    
    #Compute the final result averaged over the rolling windows
    final_res_no_scaled = final_results(main_function_results$result_par$final_results_no_scaled_par, methods, types, 4,
                                        "no_scaled") 
    
    #####################################Saving the results##############################################
    
    for (name in names(final_res_no_scaled)){
      folder_name = paste0("results_folder/",name_dataset,"/",name ,"/length_train_",length_training,"/")
      save_nested_list_to_csv(final_res_no_scaled[[name]], name, name_dataset,
                              quantiles, confidence_levels, length_training, 
                              "no_scaled", folder_name, "csv")
    }
    
    for (name in names(final_res_no_scaled)){
      folder_name = paste0("results_folder/",name_dataset,"/",name ,"/length_train_",length_training,"/")
      save_nested_list_to_csv(main_function_results$result_par$final_results_no_scaled_par[[name]], 
                              name, name_dataset,
                              quantiles, confidence_levels, length_training, 
                              "no_scaled", folder_name, "RDS")
    }
    
    
    folder = paste0(dir_path, "/results_folder/", name_dataset, "/dist_parms/")
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
    saveRDS(main_function_results$parms_list, file = paste0(folder, name_dataset,"_", length_training,"_dist_parms.RDS" ))
    
  }
}


####################################VISUALIZATION OF RESULTS##########################################################################################
#ATTENTION: TO OBTAIN ALL THE RESULTS ON THE PAPER YOU NEED TO RUN THE CODE ABOVE FOR EVERY DATASETS
#           IF YOU WANT TO VISUALIZE THE RESULTS OF ONLY ONE DATASET, SELECT IT BY RUNNING THE FOLLOWING 
#           CHUNK

#Running the code below all the plots of the paper related to the real datasets will be 
#computed and saved in a folder called "Paper_plots".
#The values in the Table 4 will be displayed in the Console.

{
  # Display options
options <- c("Swiss_tourism", "Australian_tourism_zone", "Australian_tourism_no_zone")
cat("Choose one or more datasets by number, separated by commas:\n")

for (i in seq_along(options)) {
  cat(i, ":", options[i], "\n")
}

# Get input and parse it
input <- readline(prompt = "Enter numbers (e.g., 1,2): ")
selected_indices <- as.integer(unlist(strsplit(input, ",")))

# Validate input and get selected fruits
datasets <- options[selected_indices[!is.na(selected_indices) & selected_indices %in% seq_along(options)]]


for (name_dataset in datasets){
  if (name_dataset == "Swiss_tourism"){
    train = c(40) 
  }else if (name_dataset == "Australian_tourism_zone"){
    train = c(55, 110)
  }else if (name_dataset == "Australian_tourism_no_zone"){
    train = c(25, 40)
  }
  
  
  #Visualizing the values for a single training length and a single score choice for all the methods.
  
  #For "COVERAGE" and "MSE" will be computed the arithmetic mean over the variables, while for the others
  #scores the geometric mean over the variables.
  #Furthermore, for all scores except "COVERAGE" will be displayed the relative measure with respect to the Base results.
  
  #Choices for name_dataset: "Swiss_tourism", "Australian_tourism_no_zone",
  #"Australian_tourism_zone".
  
  #Results of Table 4 of the paper: insert as measure "COVERAGE" or "PI_WIDTH" and length_training, 
  #name_dataset as in the paper.
  
  {
    for (length_training in train){
      read_measure_to_dataframe(measure = "COVERAGE", name_dataset, length_training,
                                confidence_levels, quantiles, 2)
      
      read_measure_to_dataframe(measure = "PI_WIDTH", name_dataset, length_training,
                                confidence_levels, quantiles, 2)
    }
    
  }
  
  #Computing Friedman and post-hoc Nemenyi tests results as in the paper.
  #All the plots will be saved in a folder called "Paper_plots".
  {
    #Other parameters
    metrics = c("MSE", "CRPS", "MIS")
    levs = c("80", "95")
    
    for (length_training in train){
      for (name in metrics){
        if (name == "MIS"){
          for (lev in levs) {
            MCB_test(name_dataset, name, length_training, lev)
          }
        }else{
          MCB_test(name_dataset, name, length_training, NULL)
        }
      }
    } 
  }
  
  
  #To plot the scores over the training set lengths (Figure 5, 10, 11, D.14).
  
  #The plots produced are with both t-Rec and t-Rec-Diag (i.e. Figure 11 and D.14), 
  #but Figure 5 and 11 are exactly the same without t-Rec-Diag. 
  
  #All the plots will be saved in a folder called "Paper_plots".
  {
    if ((name_dataset == "Swiss_tourism")|(name_dataset == "Swiss_tourism_Bigger")){
      train = c(seq(30, 40, 2), seq(45, 60, 5)) 
    }else if (name_dataset == "Australian_tourism_zone"){
      train = seq(30, 120, 10)
    }else if (name_dataset == "Australian_tourism_no_zone"){
      train = seq(15, 50, 5)
    }
    methods = c("MinT_shrink", "t_Rec_Diag", 
                "t_Rec_Full")
    scores <- c("MSE", "CRPS", "MIS")
    levs = c("80", "95")
    TrainPlotScore(train, scores, name_dataset, 4, confidence_levels, quantiles,
                   F, methods, levs)
  }
  
  #Plot of Figure 9 in the paper: relative 95% prediction interval width over incoherence
  
  #The plot will be saved in a folder called "Paper_plots", with the name 
  #paste0(dataset, "_scatter_plot_PI_WIDTH_inco_density_up_bot_length_training_", length_training, ".pdf").
  {
    if (name_dataset == "Swiss_tourism"){
      methods = c("Base", "MinT_shrink",
                  "t_Rec_Full")
      PredictionIntervalOverRolling(dataset = name_dataset, metric_name = "PI_WIDTH", 
                                    length_training = 40, n_digit = 4, 
                                    conf_levs = confidence_levels, quants = quantiles, 
                                    methods = methods, lev = "95", form = "RDS", 
                                    index_top = 1, index_bot = 2:27)
    }
  }
  
  
  #Plot of coverage MinT vs t-Rec (Figure 6).
  #(The code is specific for the Swiss_tourism dataset)
  {
    if (name_dataset == "Swiss_tourism"){
      HistogramsCoverage("COVERAGE", "Swiss_tourism", length_training = 40,
                         confidence_levels = confidence_levels, quantiles = NULL, n_digit = 4,
                         lev = "95")
    }
    
  }
  
}
}
