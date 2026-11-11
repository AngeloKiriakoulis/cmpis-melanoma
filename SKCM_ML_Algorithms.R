# Complete SKCM Analysis Pipeline
# Data Preparation, Univariate Cox Analysis, and 286 ML Models

# Load required libraries
library(survival)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(plsRcox)
library(survivalsvm)
library(superpc)
library(gbm)
library(pec)
library(Hmisc)
library(dplyr)
library(foreach)
library(tibble)
library(readr)
library(limma)
library(pheatmap)
library(survminer)
library(parallel)
library(doParallel)

library(ggplot2)
library(gridExtra)
library(viridis)
library(patchwork)
library(tidyr)
Sys.setenv(RENV_WATCHDOG_ENABLED = FALSE)

# cores <- detectCores()
# cl <- makeCluster(8L)
# registerDoParallel(cl)

# ============================================================================
# PART 1: DATA PREPARATION FUNCTIONS
# ============================================================================

#' Prepare TCGA expression data
#' @param expression_file Path to TCGA expression file (genes x samples)
#' @param clinical_file Path to TCGA clinical data file
#' @param pcdrg_genes Vector of PCDRG gene names (1,254 genes)
#' @param cancer_type Cancer type code (default: "SKCM")
prepare_tcga_data <- function(expression_file, clinical_file, pcdrg_genes, cancer_type = "SKCM") {
  
  # Load expression data
  cat("Loading TCGA expression data...\n")
  if(grepl("\\.csv$", expression_file)) {
    expr_data <- read.csv(expression_file, row.names = 1, check.names = FALSE)
  } else if(grepl("\\.txt$|tsv$", expression_file)) {
    expr_data <- read.table(expression_file, header = TRUE, row.names = 1, 
                            sep = "\t", check.names = FALSE)
  } else {
    stop("Expression file must be .csv, .txt, or .tsv format")
  }
  
  # Load clinical data
  cat("Loading TCGA clinical data...\n")
  if(grepl("\\.csv$", clinical_file)) {
    clinical_data <- read.csv(clinical_file, check.names = FALSE)
  } else {
    clinical_data <- read.table(clinical_file, header = TRUE, sep = "\t", check.names = FALSE)
  }

  # Standardize column names for clinical data
  colnames(clinical_data) <- tolower(colnames(clinical_data))
  
  # Extract sample IDs (first 12 characters for TCGA barcodes)
  sample_ids <- substr(colnames(expr_data), 1, 12)
  colnames(expr_data) <- sample_ids
  
  # Match clinical and expression data
  clinical_ids <- clinical_data$submitter_id %in% colnames(expr_data)
  clinical_matched <- clinical_data[clinical_ids, ]
  
  expr_matched <- expr_data[, colnames(expr_data) %in% clinical_matched$submitter_id]
  
  # Ensure same order
  common_samples <- intersect(colnames(expr_matched), clinical_matched$submitter_id)
  expr_matched <- expr_matched[, common_samples]
  clinical_matched <- clinical_matched[match(common_samples, clinical_matched$submitter_id), ]
  
  # Filter for PCDRG genes
  pcdrg_in_data <- intersect(rownames(expr_matched), pcdrg_genes)
  cat("Found", length(pcdrg_in_data), "PCDRG genes in expression data\n")
  
  expr_matched <- expr_matched[pcdrg_in_data, ]
  
  
  # Prepare survival data
  survival_data <- prepare_survival_data(clinical_matched)
  
  # Combine expression and survival data
  final_data <- merge(survival_data, data.frame(submitter_id = colnames(expr_matched), t(expr_matched)), by = "submitter_id")
  
  # Remove samples with missing survival data
  complete_cases <- complete.cases(final_data[, c("time", "status")])
  final_data <- final_data[complete_cases, ]
  
  cat("Final TCGA dataset: ", nrow(final_data), " samples, ", 
      length(pcdrg_in_data), " genes\n")
  
  return(list(
    data = final_data,
    genes = pcdrg_in_data,
    sample_info = list(total_samples = nrow(final_data))
  ))
}

#' Prepare META cohort data
#' @param meta_files List of paths to META cohort files
#' @param pcdrg_genes Vector of PCDRG gene names
prepare_meta_data <- function(meta_files, pcdrg_genes) {
  
  meta_datasets <- list()
  
  for(i in 1:length(meta_files)) {
    cat("Processing META cohort", i, "...\n")
    
    # Load data (assuming combined expression and clinical data)
    if(grepl("\\.csv$", meta_files[[i]])) {
      data <- read.csv(meta_files[[i]], check.names = FALSE)
    } else {
      data <- read.table(meta_files[[i]], header = TRUE, sep = "\t", check.names = FALSE)
    }
    
    # Identify survival columns (common names)
    time_cols <- grep("time|survival|os|dfs|rfs", colnames(data), ignore.case = TRUE)
    status_cols <- grep("status|event|death|recur", colnames(data), ignore.case = TRUE)
    
    if(length(time_cols) == 0 || length(status_cols) == 0) {
      warning(paste("Could not identify survival columns in META cohort", i))
      next
    }
    
    # Use first identified survival columns
    survival_data <- data.frame(
      time = data[, time_cols[1]],
      status = data[, status_cols[1]]
    )
    
    # Filter for PCDRG genes
    gene_cols <- intersect(colnames(data), pcdrg_genes)
    expr_data <- data[, gene_cols]
    
    # Combine
    final_data <- cbind(survival_data, expr_data)
    
    # Remove missing cases
    complete_cases <- complete.cases(final_data[, c("time", "status")])
    final_data <- final_data[complete_cases, ]
    
    meta_datasets[[i]] <- list(
      data = final_data,
      genes = gene_cols
    )
    
    cat("META cohort", i, ": ", nrow(final_data), " samples, ", 
        length(gene_cols), " genes\n")
  }
  
  return(meta_datasets)
}

#' Prepare survival data from clinical information
prepare_survival_data <- function(clinical_data) {
  
  # Handle different survival time formats
  if("days_to_death" %in% colnames(clinical_data)) {
    time_to_death <- as.numeric(clinical_data$days_to_death)
  } else if("overall_survival_months" %in% colnames(clinical_data)) {
    time_to_death <- as.numeric(clinical_data$overall_survival_months) * 30.44
  } else {
    # Look for any time column
    time_cols <- grep("time|survival|days|months", colnames(clinical_data), ignore.case = TRUE)
    if(length(time_cols) > 0) {
      time_to_death <- as.numeric(clinical_data[, time_cols[1]])
    } else {
      stop("Could not identify survival time column")
    }
  }
  
  if("days_to_last_follow_up" %in% colnames(clinical_data)) {
    time_to_last_followup <- as.numeric(clinical_data$days_to_last_follow_up)
  } else if("last_follow_up_months" %in% colnames(clinical_data)) {
    time_to_last_followup <- as.numeric(clinical_data$last_follow_up_months) * 30.44
  } else {
    time_to_last_followup <- rep(NA, nrow(clinical_data))
  }
  
  # Calculate overall survival time and status
  time <- ifelse(!is.na(time_to_death), time_to_death, time_to_last_followup)
  status <- ifelse(!is.na(time_to_death), 1, 0)
  
  # Convert months to days if necessary
  if(max(time, na.rm = TRUE) < 100) {  # Likely in months
    time <- time * 30.44
  }
  
  # Remove invalid times
  valid_time <- !is.na(time) & time > 0
  
  survival_data <- data.frame(
    submitter_id = clinical_data$submitter_id,
    time = time,
    status = status
  )
  
  return(survival_data[valid_time, ])
}

# ============================================================================
# PART 2: MACHINE LEARNING ALGORITHMS (286 COMBINATIONS)
# ============================================================================

# [Previous ML algorithm functions remain the same - stepwise_cox, gbm_cox, etc.]

# Function to perform Stepwise Cox regression with direction parameter
stepwise_cox <- function(train_data, test_data, candidate_genes, direction = "both") {
  formula_str <- paste("Surv(time, status) ~", paste(candidate_genes, collapse = " + "))
  full_model <- coxph(as.formula(formula_str), data = train_data)
  step_model <- step(full_model, direction = direction, trace = FALSE)
  
  # Predict risk scores
  risk_scores <- predict(step_model, newdata = test_data, type = "risk")
  return(list(model = step_model, risk_scores = risk_scores))
}

# Function for Generalized Boosted Regression Modeling (GBM)
gbm_cox <- function(train_data, test_data, candidate_genes) {
  # Suppress all output using capture.output
  invisible(capture.output({
    gbm_model <- gbm(Surv(time, status) ~ ., 
                     data = train_data[, c("time", "status", candidate_genes)],
                     distribution = "coxph",
                     n.trees = 500,
                     n.cores = 1,
                     shrinkage = 0.05,
                     interaction.depth = 3,
                     bag.fraction = 0.5,
                     train.fraction = 0.8,
                     cv.folds = 5,
                     verbose = FALSE)
    
    # This is what was printing the CV fold numbers
    best_iter <- gbm.perf(gbm_model, method = "cv", plot.it = FALSE)
  }))
  
  risk_scores <- predict(gbm_model, newdata = test_data, n.trees = best_iter, type = "link")
  return(list(model = gbm_model, risk_scores = risk_scores))
}

# Function for Random Survival Forest (RSF)
rsf_model <- function(train_data, test_data, candidate_genes) {
  formula_str <- paste("Surv(time, status) ~", paste(candidate_genes, collapse = " + "))
  
  rsf_fit <- rfsrc(as.formula(formula_str), 
                   data = train_data,
                   ntree = 1000,
                   nodesize = 10,
                   mtry = max(1, floor(sqrt(length(candidate_genes)))),
                   importance = TRUE,
                   seed = 123)
  
  pred <- predict(rsf_fit, newdata = test_data)
  risk_scores <- pred$predicted
  return(list(model = rsf_fit, risk_scores = risk_scores))
}

# Function for CoxBoost
coxboost_model <- function(train_data, test_data, candidate_genes) {
  tryCatch({
    # Basic data preparation with complete cases only
    complete_rows <- complete.cases(train_data[, c("time", "status", candidate_genes)])
    if (sum(complete_rows) < 10) {
      return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
    }
    
    train_clean <- train_data[complete_rows, ]
    
    train_x <- as.matrix(train_clean[, candidate_genes, drop = FALSE])
    train_time <- train_clean$time
    train_status <- train_clean$status
    
    # Use simple CoxBoost without cross-validation
    cb_fit <- CoxBoost(
      time = train_time,
      status = train_status,
      x = train_x,
      stepno=100  # Fixed number of steps
    )
    
    test_x <- as.matrix(test_data[, candidate_genes, drop = FALSE])
    
    # Handle missing values in test data
    test_complete <- complete.cases(test_x)
    risk_scores <- rep(0, nrow(test_data))
    
    if (sum(test_complete) > 0) {
      pred_scores <- predict(cb_fit, newx = test_x[test_complete, , drop = FALSE], type = "lp")
      risk_scores[test_complete] <- pred_scores
    }
    
    return(list(model = cb_fit, risk_scores = risk_scores))
    
  }, error = function(e) {
    warning(paste("Simple CoxBoost error:", e$message))
    return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
  })
}

# Function for Elastic Net (Enet)
enet_model <- function(train_data, test_data, candidate_genes, alpha = 0.5) {
  train_x <- as.matrix(train_data[, candidate_genes])
  train_y <- Surv(train_data$time, train_data$status)
  
  cv_fit <- cv.glmnet(x = train_x, y = train_y, 
                      family = "cox", alpha = alpha, nfolds = 3, nlambda=100)
  
  test_x <- as.matrix(test_data[, candidate_genes])
  risk_scores <- predict(cv_fit, newx = test_x, s = "lambda.min", type = "link")
  return(list(model = cv_fit, risk_scores = as.vector(risk_scores)))
}

# Function for PLS Cox regression
plsrcox_model <- function(train_data, test_data, candidate_genes) {
  train_x <- as.matrix(train_data[, candidate_genes])
  train_time <- train_data$time
  train_status <- train_data$status
  
  pls_fit <- plsRcox(Xplan = train_x, 
                     time = train_time, 
                     event = train_status,
                     nt = min(6, ncol(train_x)),
                     verbose = FALSE)
  
  test_x <- as.matrix(test_data[, candidate_genes])
  risk_scores <- predict(pls_fit, newdata = test_x, type = "lp")
  return(list(model = pls_fit, risk_scores = risk_scores))
}

# Function for Survival SVM
survival_svm_model <- function(train_data, test_data, candidate_genes) {
  # Combine survival object with candidate gene data
  train_df <- data.frame(time = train_data$time,
                         status = train_data$status,
                         train_data[, candidate_genes, drop = FALSE])
  train_df$SurvObj <- with(train_df, Surv(time, status))
  
  # survivalsvm expects formula like Surv(time, status) ~ predictors, so use SurvObj here
  svm_fit <- survivalsvm(SurvObj ~ ., data = train_df[, -c(1,2)],  # Exclude original time/status columns
                         gamma.mu = 0.1, kernel = "rbf_kernel")
  
  # For prediction, create test_df similarly
  test_df <- test_data[, candidate_genes, drop = FALSE]
  
  risk_scores <- predict(svm_fit, newdata = test_df)
  
  return(list(model = svm_fit, risk_scores = risk_scores$predicted))
}

superpc_model_improved <- function(train_data, test_data, candidate_genes) {
  tryCatch({
    # Ensure data is properly formatted
    train_x <- t(as.matrix(train_data[, candidate_genes]))
    train_y <- train_data$time
    train_censoring.status <- train_data$status
    
    # Check for valid data
    if (length(train_y) == 0 || length(train_censoring.status) == 0) {
      return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
    }
    
    # Remove any NA values
    valid_idx <- complete.cases(t(train_x)) & !is.na(train_y) & !is.na(train_censoring.status)
    if(sum(valid_idx) < 10) {  # Need minimum samples
      return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
    }
    
    train_x <- train_x[, valid_idx, drop = FALSE]
    train_y <- train_y[valid_idx]
    train_censoring.status <- train_censoring.status[valid_idx]
    
    data_obj <- list(x = train_x, 
                     y = train_y, 
                     censoring.status = train_censoring.status,
                     featurenames = candidate_genes)
    
    # Train model
    train_obj <- superpc.train(data_obj, type = "survival")
    
    # Simple cross-validation (reduce complexity)
    cv_obj <- superpc.cv(train_obj, data_obj, n.threshold = 20, n.fold = 3)
    
    # Use first component and find best threshold
    if(is.null(cv_obj$scor) || nrow(cv_obj$scor) == 0) {
      return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
    }
    
    component <- 1
    scores <- cv_obj$scor[component, ]
    
    # Handle case where all scores are NA
    if(all(is.na(scores))) {
      threshold <- cv_obj$thresholds[1]
    } else {
      best_idx <- which.max(scores)
      threshold <- cv_obj$thresholds[best_idx]
    }
    
    # Prepare test data
    test_x <- t(as.matrix(test_data[, candidate_genes]))
    test_y <- test_data$time
    test_censoring.status <- test_data$status
    
    test_data_obj <- list(
      x = test_x,
      y = test_y,
      censoring.status = test_censoring.status,
      featurenames = candidate_genes
    )
    
    # Predict with error handling
    pred_result <- superpc.predict(train_obj, data_obj, test_data_obj, 
                                   threshold = threshold, 
                                   n.components = 1)
    
    if(is.null(pred_result$v.pred) || all(is.na(pred_result$v.pred))) {
      return(list(model = train_obj, risk_scores = rep(0, nrow(test_data))))
    }
    
    return(list(model = train_obj, risk_scores = pred_result$v.pred))
    
  }, error = function(e) {
    cat("SuperPC error:", e$message, "\n")
    return(list(model = NULL, risk_scores = rep(NA, nrow(test_data))))
  })
}


# Function for Lasso regression
lasso_model <- function(train_data, test_data, candidate_genes) {
  return(enet_model(train_data, test_data, candidate_genes, alpha = 1))
}

# Function for Ridge regression
ridge_model <- function(train_data, test_data, candidate_genes) {
  return(enet_model(train_data, test_data, candidate_genes, alpha = 0))
}

# Function to calculate C-index
calculate_cindex <- function(risk_scores, time, status, higher_risk_higher_value = FALSE) {
  if(length(risk_scores) != length(time) || length(time) != length(status)) {
    return(NA)
  }
  
  valid_idx <- !is.na(risk_scores) & !is.na(time) & !is.na(status)
  if(sum(valid_idx) < 10) return(NA)
  
  risk_scores <- risk_scores[valid_idx]
  time <- time[valid_idx]
  status <- status[valid_idx]
  
  if (!higher_risk_higher_value) {
    risk_scores <- -risk_scores
  }
  
  tryCatch({
    concordance(risk_scores ~ Surv(time, status))$concordance
  }, error = function(e) {
    NA
  })
}


# Main function to run all 286 algorithm combinations
run_286_ml_combinations <- function(train_data, validation_cohorts, candidate_genes) {

  # Define base algorithms (8 with default parameters)
  base_algorithms <- list(
    "SuperPC" = superpc_model_improved,
    "GBM" = gbm_cox,
    "RSF" = rsf_model,
    "CoxBoost" = coxboost_model,
    "plsRcox" = plsrcox_model,
    "Survival-SVM" = survival_svm_model,
    "Lasso" = lasso_model,
    "Ridge" = ridge_model
  )
  
  # Create all algorithm combinations with parameters
  algorithm_combinations <- list()
  combination_names <- character()
  
  # 1. StepCox with 3 directions (3 combinations)
  stepcox_directions <- c("forward", "backward", "both")
  for(direction in stepcox_directions) {
    algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
      algorithms = list(stepwise_cox),
      params = list(list(direction = direction)),
      names = paste0("StepCox_", direction)
    )
    combination_names <- c(combination_names, paste0("StepCox_", direction))
  }
  
  # 2. Enet with 9 alpha values (9 combinations)
  enet_alphas <- seq(0.1, 0.9, by = 0.1)
  for(alpha in enet_alphas) {
    algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
      algorithms = list(enet_model),
      params = list(list(alpha = alpha)),
      names = paste0("Enet_alpha", alpha)
    )
    combination_names <- c(combination_names, paste0("Enet_alpha", alpha))
  }
  
  # 3. Base algorithms with default parameters (8 combinations)
  for(i in 1:length(base_algorithms)) {
    algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
      algorithms = list(base_algorithms[[i]]),
      params = list(list()),
      names = names(base_algorithms)[i]
    )
    combination_names <- c(combination_names, names(base_algorithms)[i])
  }
  
  # 4. Pairwise combinations
  # StepCox (3) x Base algorithms (8) = 24 combinations
  for(direction in stepcox_directions) {
    for(i in 1:length(base_algorithms)) {
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(stepwise_cox, base_algorithms[[i]]),
        params = list(list(direction = direction), list()),
        names = paste0("StepCox_", direction, "+", names(base_algorithms)[i])
      )
      combination_names <- c(combination_names, paste0("StepCox_", direction, "+", names(base_algorithms)[i]))
    }
  }
  
  # StepCox (3) x Enet (9) = 27 combinations
  for(direction in stepcox_directions) {
    for(alpha in enet_alphas) {
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(stepwise_cox, enet_model),
        params = list(list(direction = direction), list(alpha = alpha)),
        names = paste0("StepCox_", direction, "+Enet_alpha", alpha)
      )
      combination_names <- c(combination_names, paste0("StepCox_", direction, "+Enet_alpha", alpha))
    }
  }
  
  # Enet (9) x Base algorithms (8) = 72 combinations
  for(alpha in enet_alphas) {
    for(i in 1:length(base_algorithms)) {
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(enet_model, base_algorithms[[i]]),
        params = list(list(alpha = alpha), list()),
        names = paste0("Enet_alpha", alpha, "+", names(base_algorithms)[i])
      )
      combination_names <- c(combination_names, paste0("Enet_alpha", alpha, "+", names(base_algorithms)[i]))
    }
  }
  
  # Base algorithms pairwise combinations: C(8,2) = 28 combinations
  for(i in 1:(length(base_algorithms)-1)) {
    for(j in (i+1):length(base_algorithms)) {
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(base_algorithms[[i]], base_algorithms[[j]]),
        params = list(list(), list()),
        names = paste0(names(base_algorithms)[i], "+", names(base_algorithms)[j])
      )
      combination_names <- c(combination_names, paste0(names(base_algorithms)[i], "+", names(base_algorithms)[j]))
    }
  }
  
  # Additional combinations to reach exactly 286
  # StepCox x StepCox combinations
  for(i in 1:length(stepcox_directions)) {
    for(j in 1:length(stepcox_directions)) {
      # Avoid duplicates if you want, or keep all combinations (including same-direction pairs)
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(stepwise_cox, stepwise_cox),
        params = list(list(direction = stepcox_directions[i]), list(direction = stepcox_directions[j])),
        names = paste0("StepCox_", stepcox_directions[i], "+StepCox_", stepcox_directions[j])
      )
      combination_names <- c(combination_names, paste0("StepCox_", stepcox_directions[i], "+StepCox_", stepcox_directions[j]))
    }
  }
  
  # Add Enet_alpha + StepCox_direction pairs
  for(i in 1:length(enet_alphas)) {
    for(j in 1:length(stepcox_directions)) {
      algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
        algorithms = list(enet_model, stepwise_cox),
        params = list(list(alpha = enet_alphas[i]), list(direction = stepcox_directions[j])),
        names = paste0("Enet_alpha", enet_alphas[i], "+StepCox_", stepcox_directions[j])
      )
      combination_names <- c(combination_names, paste0("Enet_alpha", enet_alphas[i], "+StepCox_", stepcox_directions[j]))
    }
  }
  
  # Enet x Enet combinations with different alphas to reach 286
  enet_pairs_needed <- 286 - length(algorithm_combinations)
  pair_count <- 0
  
  # Add Enet pairs including same-alpha pairs (i == j allowed)
  for(i in 1:length(enet_alphas)) {
    for(j in 1:length(enet_alphas)) {
      # Only add if we still need combos
      if(pair_count < enet_pairs_needed) {
        algorithm_combinations[[length(algorithm_combinations) + 1]] <- list(
          algorithms = list(enet_model, enet_model),
          params = list(list(alpha = enet_alphas[i]), list(alpha = enet_alphas[j])),
          names = paste0("Enet_alpha", enet_alphas[i], "+Enet_alpha", enet_alphas[j])
        )
        combination_names <- c(combination_names, paste0("Enet_alpha", enet_alphas[i], "+Enet_alpha", enet_alphas[j]))
        pair_count <- pair_count + 1
      } else {
        break  # Stop early if we hit needed combos
      }
    }
  }
  cat("Total combinations generated:", length(algorithm_combinations), "\n")
  
  # Enhanced version with detailed C-index storage and logging
  
  # Add these columns to your results dataframe initialization
  results <- data.frame(
    combination = combination_names,
    mean_cindex = numeric(length(combination_names)),
    train_cindex = numeric(length(combination_names)),
    stringsAsFactors = FALSE
  )
  
  # Add columns for each validation cohort
  validation_names <- paste0("validation_", 1:length(validation_cohorts))
  for(val_name in validation_names) {
    results[[val_name]] <- numeric(length(combination_names))
  }
  
  # Create a detailed results list for further analysis
  detailed_results <- list()
  
  # Run evaluation for each combination
  for(combo_idx in 1:length(algorithm_combinations)) {
    cat("Processing combination", combo_idx, "of", length(algorithm_combinations), 
        ":", combination_names[combo_idx], "\n")
    
    current_combo <- algorithm_combinations[[combo_idx]]
    
    # Evaluate on all datasets (training + validation)
    all_datasets <- c(list(train_data), validation_cohorts)
    all_cindices <- numeric(length(all_datasets))
    dataset_names <- c("training", validation_names)
    
    for(dataset_idx in 1:length(all_datasets)) {
      dataset <- all_datasets[[dataset_idx]]
      ensemble_scores <- numeric(nrow(dataset))
      n_models <- 0
      
      for(algo_idx in 1:length(current_combo$algorithms)) {
        tryCatch({
          algorithm_func <- current_combo$algorithms[[algo_idx]]
          params <- current_combo$params[[algo_idx]]
          
          if(length(params) == 0) {
            model_result <- algorithm_func(train_data, dataset, candidate_genes)
          } else {
            model_result <- do.call(algorithm_func, c(list(train_data, dataset, candidate_genes), params))
          }
          
          if(!is.null(model_result$risk_scores) && length(model_result$risk_scores) == nrow(dataset)) {
            if(n_models == 0) {
              ensemble_scores <- model_result$risk_scores
            } else {
              ensemble_scores <- ensemble_scores + model_result$risk_scores
            }
            n_models <- n_models + 1
          }
        }, error = function(e) {
          cat("Error in algorithm", algo_idx, ":", e$message, "\n")
        })
      }
      
      if(n_models > 0) {
        ensemble_scores <- ensemble_scores / n_models
        all_cindices[dataset_idx] <- calculate_cindex(ensemble_scores, 
                                                      dataset$time, 
                                                      dataset$status)
      } else {
        all_cindices[dataset_idx] <- NA
      }
      
      # Enhanced logging with dataset names
      cat("  ", dataset_names[dataset_idx], "C-index:", 
          round(all_cindices[dataset_idx], 4), "\n")
    }
    
    # Store detailed results for this combination
    detailed_results[[combination_names[combo_idx]]] <- list(
      combination_name = combination_names[combo_idx],
      train_cindex = all_cindices[1],
      validation_cindices = all_cindices[-1],
      validation_names = validation_names,
      all_cindices = all_cindices,
      dataset_names = dataset_names
    )
    
    # Store in results dataframe
    results$train_cindex[combo_idx] <- all_cindices[1]
    for(i in 1:length(validation_cohorts)) {
      results[[validation_names[i]]][combo_idx] <- all_cindices[i + 1]
    }
    
    # Calculate mean C-index across all datasets (training + validation)
    results$mean_cindex[combo_idx] <- mean(all_cindices, na.rm = TRUE) 
    
    # Print summary for this combination
    cat("  Mean C-index across all datasets:", 
        round(results$mean_cindex[combo_idx], 4), "\n") 
    cat("  Validation C-indices:", 
        paste(round(all_cindices[-1], 4), collapse = ", "), "\n\n")
  }
  
  # Return both summary and detailed results
  return_object <- list(
    summary_results = results,
    detailed_results = detailed_results
  )
  
  return(return_object)
}
# Additional analysis functions you can use after running:

# Function to extract validation results for analysis
extract_validation_results <- function(results_object) {
  detailed <- results_object$detailed_results
  
  # Create matrix of validation C-indices
  validation_matrix <- matrix(NA, 
                              nrow = length(detailed), 
                              ncol = length(detailed[[1]]$validation_cindices))
  
  rownames(validation_matrix) <- names(detailed)
  colnames(validation_matrix) <- detailed[[1]]$validation_names
  
  for(i in 1:length(detailed)) {
    validation_matrix[i, ] <- detailed[[i]]$validation_cindices
  }
  
  return(validation_matrix)
}
  
  # Function to print detailed summary
print_detailed_summary <- function(results_object) {
  cat("=== DETAILED C-INDEX SUMMARY ===\n\n")
  
  for(combo_name in names(results_object$detailed_results)) {
    combo_data <- results_object$detailed_results[[combo_name]]
    cat("Combination:", combo_name, "\n")
    cat("Training C-index:", round(combo_data$train_cindex, 4), "\n")
    
    for(i in 1:length(combo_data$validation_cindices)) {
      cat("Validation", i, "C-index:", 
          round(combo_data$validation_cindices[i], 4), "\n")
    }
    
    cat("Mean C-index:", 
        round(mean(combo_data$all_cindices, na.rm = TRUE), 4), "\n")
    cat("---\n")
  }
}

# ============================================================================
# PART 3: MODEL CONSTRUCTION AND EVALUATION
# ============================================================================
build_final_model <- function(train_data, best_combination, candidate_genes) {
  
  cat("Building final model with combination:", best_combination$algorithms_used, "\n")
  cat("Training StepCox forward selection + Random Survival Forest ensemble\n")
  
  # Train final models
  final_model_stepcox <- stepwise_cox(train_data, train_data, candidate_genes, direction = "forward")
  final_model_rsf <- rsf_model(train_data, train_data, candidate_genes)
  
  # StepCox coefficients
  stepcox_coef <- coef(final_model_stepcox$model)
  selected_genes <- names(stepcox_coef)
  
  # RSF importance (subset to selected genes)
  rsf_importance <- final_model_rsf$model$importance
  rsf_importance_selected <- rsf_importance[selected_genes]
  
  # ---------- Simple Average Ensemble ----------
  # Just average the standardized coefficients and importance
  
  # Standardize coefficients and importance to same scale
  stepcox_coef_std <- stepcox_coef / max(abs(stepcox_coef))
  rsf_importance_std <- rsf_importance_selected / max(abs(rsf_importance_selected))
  
  # Simple average
  ## CORRECTED (Phase 2): averaging a signed Cox coefficient with an unsigned
  ## RSF importance corrupted gene directions. Sign comes from Cox; RSF
  ## importance (clipped at 0, scaled to [0,1]) modulates magnitude.
  ## Downstream scores must NOT be inverted with these weights.
  rsf_imp01 <- pmax(rsf_importance_selected, 0) / max(pmax(rsf_importance_selected, 0))
  ensemble_coefficients <- stepcox_coef * rsf_imp01
  ## alternative ensemble formula, kept for reference:
  ## ensemble_coefficients <- (stepcox_coef_std + rsf_importance_std) / 2
  
  cat("Final StepCox model includes", length(selected_genes), "genes\n")
  cat("RSF model trained on", length(candidate_genes), "genes\n")
  cat("Ensemble average coefficients created for", length(ensemble_coefficients), "genes\n")
  
  return(list(
    stepcox_model = final_model_stepcox,
    rsf_model = final_model_rsf,
    selected_genes = selected_genes,
    stepcox_coefficients = stepcox_coef,
    rsf_importance = rsf_importance,
    rsf_importance_selected = rsf_importance_selected,
    ensemble_coefficients = ensemble_coefficients   # use this for MPCDS
  ))
}

#' Calculate risk scores for all samples using StepCox + RSF ensemble
calculate_risk_scores <- function(final_model, all_datasets) {
  
  risk_scores_list <- list()
  
  for(i in 1:length(all_datasets)) {
    dataset <- all_datasets[[i]]
    
    # Calculate risk scores from both models
    
    # StepCox scores (linear predictor)
    stepcox_scores <- predict(final_model$stepcox_model$model, 
                              newdata = dataset, 
                              type = "lp")
    
    # RSF scores (mortality risk)
    # Note: RSF predict returns survival probabilities, we need mortality risk
    rsf_pred <- predict(final_model$rsf_model$model, 
                        newdata = dataset[, final_model$selected_genes],type = "lp") #results$ or nothing
    
    # Convert RSF survival to mortality risk (1 - survival probability at median time)
    # Assuming we want risk at median follow-up time
    median_time <- median(dataset$time[dataset$status == 1], na.rm = TRUE)
    
    # Extract mortality risk from RSF (this may need adjustment based on randomForestSRC version)
    if("mortality" %in% names(rsf_pred)) {
      rsf_scores <- rsf_pred$mortality
    } else {
      # Alternative: use ensemble mortality or calculate from survival
      rsf_scores <- 1 - rsf_pred$survival[, which.min(abs(rsf_pred$time.interest - median_time))]
    }
    
    # Standardize scores to have similar scales
    stepcox_scores_std <- scale(stepcox_scores)[,1]
    rsf_scores_std <- scale(rsf_scores)[,1]
    
    # Calculate ensemble risk scores (average of standardized scores)
    ensemble_scores <- (stepcox_scores_std + rsf_scores_std) / 2
    
    # Classify as high/low risk based on median
    risk_group <- ifelse(ensemble_scores > median(ensemble_scores), "High", "Low")
    
    risk_scores_list[[i]] <- data.frame(
      sample_id = rownames(dataset),
      stepcox_score = stepcox_scores,
      rsf_score = rsf_scores,
      ensemble_score = ensemble_scores,
      risk_group = risk_group,
      time = dataset$time,
      status = dataset$status
    )
    
    cat("Cohort", i, "- High risk:", sum(risk_group == "High"), 
        "Low risk:", sum(risk_group == "Low"), "\n")
  }
  
  names(risk_scores_list) <- paste0("Cohort_", 1:length(all_datasets))
  return(risk_scores_list)
}

#' Perform survival analysis and generate Kaplan-Meier plots
perform_survival_analysis <- function(risk_scores, output_dir = "results") {
  
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  survival_results <- list()
  
  for(i in 1:length(risk_scores)) { #results$ or nothing
    cohort_name <- names(risk_scores)[i] #results$ or nothing
    data <- risk_scores[[i]] #results$ or nothing
    
    # Fit survival model
    surv_fit <- survfit(Surv(time, status) ~ risk_group, data = data)
    
    # Log-rank test
    surv_diff <- survdiff(Surv(time, status) ~ risk_group, data = data)
    p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    
    # Cox regression
    cox_model <- coxph(Surv(time, status) ~ risk_group, data = data)
    cox_summary <- summary(cox_model)
    
    # Store results
    survival_results[[cohort_name]] <- list(
      surv_fit = surv_fit,
      p_value = p_value,
      cox_hr = cox_summary$coefficients[1, "exp(coef)"],
      cox_pvalue = cox_summary$coefficients[1, "Pr(>|z|)"],
      cox_ci_lower = cox_summary$conf.int[1, "lower .95"],
      cox_ci_upper = cox_summary$conf.int[1, "upper .95"]
    )
    
    # Generate Kaplan-Meier plot
    p <- ggsurvplot(surv_fit,
                    data = data,
                    pval = TRUE,
                    conf.int = FALSE,
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    linetype = "strata",
                    surv.median.line = "hv",
                    ggtheme = theme_bw(),
                    palette = c("#E7B800", "#2E9FDF"),
                    title = paste("Kaplan-Meier Survival Curve -", cohort_name),
                    xlab = "Time (days)",
                    ylab = "Overall Survival Probability")
    
    # Save plot
    ggsave(file.path(output_dir, paste0("KM_plot_", cohort_name, ".png")),
           plot = p$plot, width = 10, height = 8, dpi = 300)
    
    cat("Cohort", cohort_name, ": HR =", round(cox_summary$coefficients[1, "exp(coef)"], 3),
        ", p =", format(cox_summary$coefficients[1, "Pr(>|z|)"], scientific = TRUE), "\n")
  }
  
  return(survival_results)
}

# ============================================================================
# PART 4: MAIN ANALYSIS PIPELINE
# ============================================================================

#' Main analysis pipeline for SKCM
run_skcm_analysis <- function(tcga_expression_file, 
                              tcga_clinical_file,
                              meta_files_list,
                              pcdrg_genes_file,
                              output_dir = "skcm_results") {
  
  cat("=== Starting SKCM Analysis Pipeline ===\n")
  
  # Create output directory
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  # Load PCDRG genes
  cat("Loading PCDRG genes...\n")
  if(grepl("\\.csv$", pcdrg_genes_file)) {
    pcdrg_genes <- read.csv(pcdrg_genes_file)$gene_name  # Adjust column name as needed
  } else {
    pcdrg_genes <- read.table(pcdrg_genes_file, header = TRUE)$gene_name
  }
  cat("Loaded", length(pcdrg_genes), "PCDRG genes\n")
  
  # Step 1: Prepare datasets
  cat("\n=== Step 1: Data Preparation ===\n")
  tcga_data <- prepare_tcga_data(tcga_expression_file, tcga_clinical_file, pcdrg_genes)
  meta_data <- prepare_meta_data(meta_files_list, pcdrg_genes)

  candidate_genes <- pcdrg_genes
  
  if(length(candidate_genes) < 5) {
    stop("Too few candidate genes identified. Consider relaxing p-value threshold.")
  }
  
  # Step 3: Prepare data for ML algorithms
  cat("\n=== Step 2: Preparing ML Data ===\n")
  
  # Filter datasets to include only significant genes
  tcga_ml_data <- tcga_data$data[, c("time", "status", candidate_genes)]

  validation_cohorts <- list()
  for(i in 1:length(meta_data)) {
    common_genes <- intersect(candidate_genes, colnames(meta_data[[i]]$data))
    if(length(common_genes) > 0) {
      validation_cohorts[[i]] <- meta_data[[i]]$data[, c("time", "status", common_genes)]
    }
  }
  
  # Update candidate genes to only include those available in all cohorts
  final_genes <- candidate_genes
  for(cohort in validation_cohorts) {
    final_genes <- intersect(final_genes, colnames(cohort))
  }
  
  cat("Final candidate genes for ML:", length(final_genes), "\n")
  
  # Step 4: Run 286 ML algorithm combinations
  cat("\n=== Step 3: Running 286 ML Combinations ===\n")
  ml_results <- run_286_ml_combinations(tcga_ml_data, validation_cohorts, final_genes)
  
  # Now this will work correctly
  best_idx <- which.max(as.numeric(ml_results$summary_results$mean_cindex)) #results$ or ml_results$
  best_combination <- ml_results$summary_results[best_idx, ] #results$ or ml_results$
  
  cat("Best combination:", best_combination$combination, "\n") #results$ or nothing
  cat("Best C-index:", best_combination$mean_cindex, "\n") #results$ or nothing
  
  # Step 5: Build final model
  cat("\n=== Step 4: Building Final Model ===\n")
  final_model <- build_final_model(tcga_ml_data, best_combination, final_genes)
  
  # Save model information
  model_info <- list(
    best_combination = best_combination, # results$ or nothing
    selected_genes = final_model$selected_genes, # results$ or nothing
    coefficients = final_model$ensemble_coefficients # results$ or nothing
  )
  cat("Coefficients:", final_model$coefficients)
  saveRDS(model_info, file.path(output_dir, "final_model_info.rds"))
  
  # Step 6: Calculate risk scores for all cohorts
  cat("\n=== Step 5: Calculating Risk Scores ===\n")
  all_datasets <- c(list(tcga_ml_data), validation_cohorts)
  risk_scores <- calculate_risk_scores(final_model, all_datasets)
  
  # Save risk scores
  for(i in 1:length(risk_scores)) {
    write.csv(risk_scores[[i]], 
              file.path(output_dir, paste0("risk_scores_", names(risk_scores)[i], ".csv")), 
              row.names = FALSE)
  }
  
  # Step 7: Survival analysis
  cat("\n=== Step 6: Survival Analysis ===\n")
  survival_results <- perform_survival_analysis(risk_scores, output_dir)
  
  # Save survival results summary
  survival_summary <- data.frame(
    cohort = names(survival_results),
    hr = sapply(survival_results, function(x) x$cox_hr),
    p_value = sapply(survival_results, function(x) x$cox_pvalue),
    ci_lower = sapply(survival_results, function(x) x$cox_ci_lower),
    ci_upper = sapply(survival_results, function(x) x$cox_ci_upper),
    km_p_value = sapply(survival_results, function(x) x$p_value)
  )
  
  write.csv(survival_summary, file.path(output_dir, "survival_analysis_summary.csv"), row.names = FALSE)
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results saved in:", output_dir, "\n")
  
  return(list(
    ml_results = ml_results, #results$ or ml_results$ (2nd one)
    best_combination = best_combination,
    final_model = final_model,
    risk_scores = risk_scores,
    survival_results = survival_results
  ))
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

#Example of how to run the complete analysis:

# Define file paths
tcga_expr_file <- "tcga_skcm_expression.txt"
tcga_clinical_file <- "tcga_skcm_clinical.txt"
meta_files <- c("gse19234_cohort.txt", "gse98394_cohort.txt")
pcdrg_genes_file <- "pcdrg_genes.txt"

# Run complete analysis
results <- run_skcm_analysis(
  tcga_expression_file = tcga_expr_file,
  tcga_clinical_file = tcga_clinical_file,
  meta_files_list = meta_files,
  pcdrg_genes_file = pcdrg_genes_file,
  output_dir = "skcm_analysis_results"
)

# View results
print(results$best_combination) 
print(results$final_model$selected_genes)

#stopCluster(cl)


scores_df_wide <- results$ml_results$summary_results

# 2. Separate out mean_cindex for sorting and for rejoining later
mean_cindex_df <- scores_df_wide %>%
  dplyr::select(combination, mean_cindex)

# 3. Pivot to long format (exclude mean_cindex)
scores_df_long <- scores_df_wide %>%
  dplyr::select(-mean_cindex) %>%
  pivot_longer(
    cols = -combination,
    names_to = "Cohort",
    values_to = "Score"
  )

# 4. Join mean_cindex as a separate column (for sorting purposes)
scores_df_long <- scores_df_long %>%
  left_join(mean_cindex_df, by = "combination")

# 5. Append mean_cindex as its own "Cohort" column (to show as last column in plot)
scores_combined <- bind_rows(
  scores_df_long,
  mean_cindex_df %>%
    dplyr::rename(Score = mean_cindex) %>%
    mutate(Cohort = "mean_cindex")
)

# 6. Sort combinations by mean_cindex (descending so best at top)
sorted_combinations <- mean_cindex_df %>%
  arrange(desc(mean_cindex)) %>%
  pull(combination)

# 7. Define cohort column order (mean_cindex last)
cohort_levels <- scores_combined %>%
  filter(Cohort != "mean_cindex") %>%
  distinct(Cohort) %>%
  pull(Cohort) %>%
  sort() %>%
  c("mean_cindex")  # mean_cindex goes last

# 8. Apply correct factor levels for sorting
scores_combined <- scores_combined %>%
  mutate(
    combination = factor(combination, levels = rev(sorted_combinations)),  # so best models are at the top
    Cohort = factor(Cohort, levels = cohort_levels)
  )

# Keep only the top half of the models based on mean_cindex
top_n <- floor(length(sorted_combinations) / 2)

top_combinations <- head(sorted_combinations, top_n)

# Filter the combined dataframe to include only those top combinations
scores_combined_top <- scores_combined %>%
  filter(combination %in% top_combinations) %>%
  mutate(combination = factor(combination, levels = rev(top_combinations)))

scores_combined_top <- scores_combined_top %>%
  mutate(
    fill_type = ifelse(Cohort == "mean_cindex", "mean", "score")
  )


ggplot(scores_combined_top, aes(x = Cohort, y = combination)) +
  geom_tile(aes(
    fill = ifelse(fill_type == "mean", NA, Score)
  ), color = "white", linewidth = 0.1) +
  geom_tile(
    data = subset(scores_combined_top, fill_type == "mean"),
    fill = "lightblue", color = "white", linewidth = 0.1
  ) +
  geom_text(aes(label = round(Score, 3)), size = 1.5) +
  scale_fill_gradient2(
    low = "turquoise",
    mid = "white",
    high = "orange",
    midpoint = 0.5,
    na.value = "lightblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(size = 5)
  ) + 
  labs(
    title = "Top-performing models (with mean C-index in blue)",
    x = "Cohort",
    y = "Model",
    fill = "Score"
  )

coef_df <- data.frame(
  selected_genes = results$final_model$selected_genes,
  coefficients = results$final_model$ensemble_coefficients)

ggplot(coef_df, aes(x = coefficients, y = reorder(selected_genes, coefficients))) +
  geom_point(aes(color = coefficients), size = 5) +
  scale_color_gradient2(low = "orange", mid = "white", high = "purple", midpoint = -1) +
  theme_minimal() +
  labs(x = "Coefficient", y = "Gene", title = "Top gene coefficients") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")

current_names <- names(results[["survival_results"]][["Cohort_1"]][["surv_fit"]][["strata"]])

# Replace prefix
new_names <- sub("^risk_group=", "CMPIS=", current_names)

# Assign back
names(results[["survival_results"]][["Cohort_1"]][["surv_fit"]][["strata"]]) <- new_names

sf <- results$survival_results$Cohort_3$surv_fit

# Convert to data frame manually
df_plot <- data.frame(
  time  = sf$time,
  surv  = sf$surv,
  strata = rep(names(sf$strata), sf$strata),
  n_risk = sf$n.risk,
  n_event = sf$n.event,
  n_censor = sf$n.censor
)

medians <- summary(sf)$table[, "median"]

median_df <- data.frame(
  medians = medians,
  colors = c("red", "#2E9FDF")
)


p_text <- paste0("Log-rank\np ", format.pval(results$survival_results$Cohort_1$cox_pvalue))

# Plot with ggplot
ggplot(df_plot, aes(x = time, y = surv, color = strata)) +
  geom_step(size = 1) +
  # Horizontal line at y=0.5
  geom_segment(aes(x = 0, xend = 2184, y = 0.5, yend = 0.5), 
               color = "black", linetype = "dashed") + 
  geom_segment(data = median_df, aes(x = medians, xend = medians, y = 0, yend = 0.5), 
               color = c("red", "blue"), linetype = "dashed") +
  labs(
    title = "Kaplan-Meier Survival Curve - TCGA Data",
    x = "Time (days)",
    y = "Overall Survival Probability",
    color = "Risk Group"
  ) +
  # Annotate p-value
  annotate("text", x = 0.1, y = 0.1, label = p_text, hjust = 0) +
  scale_color_manual(values = c("red", "blue")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = c(0.85, 0.85),  # move legend inside top right
    legend.background = element_rect(fill = "white", color = "black")
  )


