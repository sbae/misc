# Load necessary libraries
require(xgboost)
require(rBayesianOptimization)

#' @title Tune XGBoost Hyperparameters with Bayesian Optimization
#'
#' @description
#' Performs hyperparameter tuning for an XGBoost model using Bayesian Optimization
#' and k-fold cross-validation. Handles both maximization and minimization of metrics.
#'
#' @param data A numeric matrix or data frame of features.
#' @param label A numeric vector of target labels.
#' @param param_bounds A named list defining the search space for hyperparameters.
#'        Each element should be a numeric vector of length 2: c(min, max).
#' @param nfold Integer. The number of folds for cross-validation.
#' @param init_points Integer. The number of random exploration steps for Bayesian Optimization.
#' @param n_iter Integer. The number of Bayesian optimization iterations.
#' @param objective Character. The learning objective for XGBoost (e.g., "binary:logistic").
#' @param eval_metric Character. The evaluation metric to optimize (e.g., "auc", "logloss", "rmse").
#' @param maximize Logical. If TRUE, the metric is maximized. If FALSE, it is minimized.
#'        If NULL (default), the function will infer the direction. Metrics in
#'        `c("auc", "aucpr", "ndcg", "map")` will be maximized; all others will be minimized.
#'
#' @return A list containing the best score and the best set of hyperparameters found.
#'
tune_xgb_bayes <- function(data, label, param_bounds, nfold = 5,
                               init_points = 5, n_iter = 15,
                               objective = "binary:logistic",
                               eval_metric = "auc",
                               subsample = 0.1,
                               seed = 256,
                               maximize = NULL) {
  set.seed(seed)  # Set seed for reproducibility

  # --- 1. Determine Optimization Direction (NEW LOGIC) ---
  if (is.null(maximize)) {
    # Define common metrics that are typically maximized
    max_metrics <- c("auc", "aucpr", "ndcg", "map")
    
    if (eval_metric %in% max_metrics) {
      maximize <- TRUE
      cat(sprintf("Inferred optimization direction: Maximizing '%s'.\n", eval_metric))
    } else {
      maximize <- FALSE
      cat(sprintf("Inferred optimization direction: Minimizing '%s'.\n", eval_metric))
    }
  }

  # --- 2. Prepare data for xgboost ---
  if (is.null(subsample)){
    data_sub <- data
    label_sub <- label
  } else {
    subsample_n <- floor(subsample*nrow(data))

    # If samsumple will give a sample size below 100, display a warning.
    if (is.numeric(subsample) && subsample_n < 100) {
      warning("Subsample size is less than 100. This may lead to unreliable results.")
    }

    # Randomly sample a subset of the data if subsample is specified
    if (is.numeric(subsample) && subsample < 1) {
      sample_indices <- sample(seq_len(nrow(data)), size = subsample_n, replace = FALSE)
      data_sub <- data[sample_indices, , drop = FALSE]
      label_sub <- label[sample_indices]
    } else {
      stop("subsample must be a numeric value less 1.")
    }
  }

  dtrain <- xgb.DMatrix(data = data_sub, label = label_sub)

  # --- 3. Define the Objective Function for Bayesian Optimization ---
  xgb_cv_bayes <- function(max_depth, eta, subsample, colsample_bytree, nrounds, gamma) {

    params <- list(
      objective = objective,
      eval_metric = eval_metric,
      booster = "gbtree",
      max_depth = floor(max_depth),
      eta = eta,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      gamma = gamma
    )

    cv_result <- xgb.cv(
      params = params,
      data = dtrain,
      nfold = nfold,
      nrounds = floor(nrounds),
      verbose = FALSE,
      early_stopping_rounds = 10
    )

    # Extract the best score from the evaluation log
    if (maximize) {
        best_score <- max(cv_result$evaluation_log[[paste0("test_", eval_metric, "_mean")]])
    } else {
        best_score <- min(cv_result$evaluation_log[[paste0("test_", eval_metric, "_mean")]])
    }
    
    # The optimizer always maximizes. If we want to minimize a metric, we return its negative value.
    score_for_optimizer <- if (maximize) best_score else -best_score

    list(
      Score = score_for_optimizer,
      Pred = NA # Not used, but required by the package
    )
  }

  # --- 4. Run Bayesian Optimization ---
  cat("--- Starting Bayesian Optimization ---\n")
  optimizer <- BayesianOptimization(
    FUN = xgb_cv_bayes,
    bounds = param_bounds,
    init_points = init_points,
    n_iter = n_iter,
    acq = "ucb",
    verbose = TRUE
  )

  # --- 5. Format and return the results ---
  cat("\n--- Optimization Finished ---\n")
  
  best_params <- as.list(optimizer$Best_Par)
  best_params$max_depth <- floor(best_params$max_depth)
  best_params$nrounds <- floor(best_params$nrounds)

  # Return the true best score (not the potentially negative one fed to the optimizer)
  actual_best_score <- if (maximize) optimizer$Best_Value else -optimizer$Best_Value

  result <- list(
    best_score = actual_best_score,
    best_params = best_params,
    optimization_details = optimizer
  )

  return(result)
}


# # --- USAGE EXAMPLE ---
# # 1. Load data from the xgboost package
# # The 'agaricus' dataset is a classic binary classification problem.
# data(agaricus.train, package = "xgboost")
# data(agaricus.test, package = "xgboost") # We'll just use the training set for tuning

# # Separate features and labels from the loaded list
# features <- agaricus.train$data # This is a sparse matrix (dgCMatrix)
# labels <- agaricus.train$label

# # 2. Define search space
# param_bounds_r <- list(
#   max_depth = c(2L, 8L), # Smaller max_depth for this simpler dataset
#   eta = c(0.01, 0.3),
#   subsample = c(0.6, 1.0),
#   colsample_bytree = c(0.6, 1.0),
#   nrounds = c(50L, 500L),
#   gamma = c(0, 1)
# )

# # --- EXAMPLE 1: MAXIMIZING AUC ---
# cat("\n\n--- EXAMPLE 1: Maximizing AUC ---\n")
# tuning_results_auc <- tune_xgb_bayes(
#   data = features,
#   label = labels,
#   param_bounds = param_bounds_r,
#   nfold = 5,
#   init_points = 5,
#   n_iter = 10, # Reduced iterations for a quicker demo
#   eval_metric = "auc"
# )
# cat(sprintf("\nBest AUC found: %.4f\n", tuning_results_auc$best_score))
# cat("Best parameters for AUC:\n")
# print(tuning_results_auc$best_params)

# # --- EXAMPLE 2: MINIMIZING LOGLOSS ---
# cat("\n\n--- EXAMPLE 2: Minimizing LogLoss ---\n")
# tuning_results_logloss <- tune_xgb_bayes(
#   data = features,
#   label = labels,
#   param_bounds = param_bounds_r,
#   nfold = 5,
#   init_points = 5,
#   n_iter = 10,
#   eval_metric = "logloss"
# )
# cat(sprintf("\nBest (minimum) LogLoss found: %.4f\n", tuning_results_logloss$best_score))
# cat("Best parameters for LogLoss:\n")
# print(tuning_results_logloss$best_params)

# # --- EXAMPLE 3: MINIMIZING RMSE (Regression) ---
# # This example remains self-contained as xgboost doesn't have a standard built-in regression dataset.
# cat("\n\n--- EXAMPLE 3: Minimizing RMSE (Regression) ---\n")
# set.seed(123)
# reg_features <- matrix(rnorm(500 * 10), ncol = 10)
# reg_labels <- reg_features[, 1] * 2 + reg_features[, 2] * 0.5 + rnorm(500)

# tuning_results_rmse <- tune_xgb_bayes(
#   data = reg_features,
#   label = reg_labels,
#   param_bounds = param_bounds_r,
#   nfold = 5,
#   init_points = 5,
#   n_iter = 10,
#   objective = "reg:squarederror", # Change objective for regression
#   eval_metric = "rmse"            # Change metric for regression
# )
# cat(sprintf("\nBest (minimum) RMSE found: %.4f\n", tuning_results_rmse$best_score))
# cat("Best parameters for RMSE:\n")
# print(tuning_results_rmse$best_params)
