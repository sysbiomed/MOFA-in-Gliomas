# This function takes input data (a list of matrices, with samples as rows and features as columns)
# and an optional annotation (a data frame with additional information about the samples).
# 
# Arguments:
# - X: A list of matrices. Rows represent samples, and columns are features.
# - n_runs: Number of MOFA initializations.
# - num_factors: Number of factors to be found in X. If NULL, the factors explaining less than 5% will be dropped.
# - likelihoods: Vector with the desired likelihoods to train each matrix given in X.
# - cinical_data: (Optional) A data frame containing annotations for the rows of X.
#               The row names of this annotation should match the row names of input_data.
#
# Returns:
# - The trained MOFA model.
# - The time taken to train the model.
# - Saves the factors and loadings matrices in a folder "results" in the directory with the respective key in X
#

run_mofa <- function(X, n_runs = 25, num_factors = NULL, likelihoods = NULL, clinical_data = NULL){
  
  time_taken <- system.time({
    
    seeds <- n_runs
    mofa_models <- list()
    
    for (seed in 1: seeds){
      X_t <- lapply(X, t)
      
      mofa <- create_mofa(X_t)
      
      data_options <- get_default_data_options(mofa)
      
      model_options <- get_default_model_options(mofa) 
      model_options$spikeslab_weights <- T
      model_options$ard_weights <- T 
      if (is.null(likelihoods) == F){
        model_options$likelihoods <- likelihoods
      }
      
      training_options <- get_default_training_options(mofa)
      training_options$seed <- seed
      training_options$convergence_mode <- "medium" 
      
      if (is.null(num_factors) == F){
        model_options$num_factors <- num_factors
      }
      else {
        training_options$drop_factor_threshold <- 0.05 
      }
      
      mofa <- prepare_mofa(
        mofa, 
        data_options = data_options, 
        model_options =  model_options, 
        training_options = training_options)
      
      if (!"mofa_env" %in% reticulate::conda_list()$name) {
        message("Creating 'mofa_env' with necessary Python packages...")
        reticulate::conda_create("mofa_env", packages = c("python=3.10"))
        reticulate::conda_install("mofa_env", packages = c("numpy==1.25.2", "scipy==1.9.3", "mofapy2==0.7.0"), pip = TRUE)
      }
      
      reticulate::use_condaenv("mofa_env", required = T)
      mofa <- MOFA2::run_mofa(mofa, use_basilisk = F)
      
      mofa_models <- append(mofa_models, mofa)
    }

    mofa <- select_model(mofa_models, plot = F)

    if (is.null(clinical_data) == F){
      clinical_data_2_df <- cbind(mofa@samples_metadata, clinical_data[mofa@samples_metadata$sample, ])
      samples_metadata(mofa) <- clinical_data_2_df 
    }

    loadings <- get_weights(mofa)
    loadings <- lapply(loadings, function(mat) {
      mat <- as.matrix(mat)  
      colnames(mat) <- stringr::str_replace_all(colnames(mat), "Factor", "Factor ")
      return(mat)  
    })
   
    factors <- get_factors(mofa)[[1]]
    colnames(factors) <- gsub("Factor", "Factor ", colnames(factors))
    result <- list(model_output = mofa)
  })

  write.csv(factors, "mofa-output-nomiRNA/factors.csv", row.names = T)
  for (i in 1:length(loadings)){
    write.csv(loadings[[i]], paste0("mofa-output-nomiRNA/loadings_", tolower(names(X)[i]), ".csv"), row.names = T)
  }

  result$time_taken <- time_taken["elapsed"]
  return (result)
}


