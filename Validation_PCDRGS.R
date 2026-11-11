# Load required libraries
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(ggplot2)
library(pheatmap)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridExtra)
library(forestplot)
library(reshape2)
library(MOVICS)
library(limma)

options(warn=-1)

tcga_skcm <- load("skcm.tcga.RData")  # TCGA-SKCM cohort
meta_cohort <- read.csv("scaled_expression_data.csv")    # META validation cohort
# gse_cohorts <- list(
#   GSE54467 = read.csv("GSE54467.csv"),    # Example GEO datasets for SKCM
#   GSE65904 = read.csv("GSE65904.csv"),
#   GSE59455 = read.csv("GSE59455.csv")
# )


pcdrg_genes <- c("CFLAR", "CYBB", "CYLD", "GYS1", "IRF1", "NFKBIA", "PIM2", "VCAM1")

# Function to calculate PCDRG signature score (PCDRS)
# Improved Function to calculate PCDRG signature score (PCDRS)
calculate_pcdrs <- function(expr_data, signature_genes, coefficients = NULL, method = "weighted", 
                            survival_data = NULL, risk_direction = "auto") {
  # Ensure genes are present in expression data
  common_genes <- intersect(signature_genes, rownames(expr_data))
  
  if(length(common_genes) < length(signature_genes)) {
    warning(paste("Missing genes:", setdiff(signature_genes, common_genes)))
    cat("Available genes:", paste(common_genes, collapse = ", "), "\n")
    cat("Missing genes:", paste(setdiff(signature_genes, common_genes), collapse = ", "), "\n")
  }
  
  # Extract expression data for signature genes
  sig_expr <- expr_data[common_genes, , drop = FALSE]
  
  # Check if expression data needs log transformation
  max_expr <- max(sig_expr, na.rm = TRUE)
  if(max_expr > 50) {
    cat("High expression values detected (max =", max_expr, "), applying log2 transformation\n")
    sig_expr <- log2(sig_expr + 1)
  }  
  if(method == "weighted" && !is.null(coefficients)) {
    # Use provided coefficients from ML model
    if(length(coefficients) != length(common_genes)) {
      stop("Number of coefficients must match number of available genes")
    }
    
    # Calculate weighted risk score
    pcdrs <- colSums(sig_expr * coefficients)
    cat("Using weighted method with provided coefficients\n")
    
    # Check coefficient directions
    positive_coefs <- sum(coefficients > 0)
    negative_coefs <- sum(coefficients < 0)
    cat("Coefficient summary: Positive =", positive_coefs, ", Negative =", negative_coefs, "\n")
    
  } else if(method == "cox_based" && !is.null(survival_data)) {
    # Calculate coefficients using Cox regression
    library(survival)
    
    # Match samples between expression and survival data
    common_samples <- intersect(colnames(sig_expr), rownames(survival_data))
    if(length(common_samples) < 10) {
      stop("Not enough overlapping samples for Cox regression")
    }
    
    sig_expr_matched <- sig_expr[, common_samples, drop = FALSE]
    surv_matched <- survival_data[common_samples, ]
    
    # Calculate Cox coefficients for each gene
    cox_coefs <- numeric(length(common_genes))
    names(cox_coefs) <- common_genes
    
    for(i in 1:length(common_genes)) {
      gene <- common_genes[i]
      tryCatch({
        cox_model <- coxph(Surv(futime, fustat) ~ gene_expr, 
                           data = data.frame(
                             futime = surv_matched$futime,
                             fustat = surv_matched$fustat,
                             gene_expr = as.numeric(sig_expr_matched[gene, ])
                           ))
        cox_coefs[gene] <- coef(cox_model)
      }, error = function(e) {
        warning(paste("Cox regression failed for gene", gene))
        cox_coefs[gene] <- 0
      })
    }
    
    # Calculate risk score using Cox coefficients
    pcdrs <- colSums(sig_expr * cox_coefs[common_genes])
    coefficients <- cox_coefs  # Store for direction checking
    cat("Using Cox-based method\n")
    
  } else {
    # Simple average method - but this might be problematic!
    cat("WARNING: Using simple mean expression method.\n")
    cat("This assumes all genes are risk-promoting with equal weights.\n")
    cat("Consider using weighted method with proper coefficients.\n")
    pcdrs <- colMeans(sig_expr, na.rm = TRUE)
    coefficients <- rep(1, length(common_genes))  # All positive weights
  }
  
  # Direction checking and correction
  if(risk_direction == "auto" && !is.null(survival_data) && !is.null(coefficients)) {
    # Match samples for direction checking
    common_samples <- intersect(names(pcdrs), rownames(survival_data))
    if(length(common_samples) > 20) {  # Need sufficient samples
      
      pcdrs_matched <- pcdrs[common_samples]
      surv_matched <- survival_data[common_samples, ]
      
      # Check correlation between PCDRS and survival time
      valid_idx <- !is.na(pcdrs_matched) & !is.na(surv_matched$futime) & !is.na(surv_matched$fustat)
      
      if(sum(valid_idx) > 20) {
        # Correlation test
        cor_test <- cor.test(pcdrs_matched[valid_idx], surv_matched$futime[valid_idx], 
                             method = "spearman")
        
        # Risk-benefit analysis
        high_risk <- pcdrs_matched >= median(pcdrs_matched, na.rm = TRUE)
        low_risk <- pcdrs_matched < median(pcdrs_matched, na.rm = TRUE)
        
        if(sum(high_risk & valid_idx) > 5 && sum(low_risk & valid_idx) > 5) {
          high_risk_events <- mean(surv_matched$fustat[high_risk & valid_idx])
          low_risk_events <- mean(surv_matched$fustat[low_risk & valid_idx])
          
          cat("\nDirection Analysis:\n")
          cat("Correlation with survival time:", round(cor_test$estimate, 3), 
              "(p =", round(cor_test$p.value, 4), ")\n")
          cat("High PCDRS group event rate:", round(high_risk_events, 3), "\n")
          cat("Low PCDRS group event rate:", round(low_risk_events, 3), "\n")
          
          # If high PCDRS group has fewer events or positive correlation with survival time
          # Then the score might need inversion
          if(high_risk_events < low_risk_events || cor_test$estimate > 0.1) {
            cat("WARNING: PCDRS appears to be inversely related to risk!\n")
            cat("Consider setting risk_direction = 'invert' or checking your coefficients.\n")
            
            if(risk_direction == "auto") {
              cat("Auto-inverting PCDRS...\n")
              pcdrs <- -pcdrs
              cat("PCDRS has been inverted. Higher values now indicate higher risk.\n")
            }
          } else {
            cat("PCDRS direction appears correct (higher values = higher risk).\n")
          }
        }
      }
    }
  } else if(risk_direction == "invert") {
    cat("Manually inverting PCDRS...\n")
    pcdrs <- -pcdrs
  }
  
  # Final summary
  cat("\nFinal PCDRS summary:\n")
  print(summary(pcdrs))
  cat("PCDRS standard deviation:", sd(pcdrs, na.rm = TRUE), "\n")
  
  # Return additional information
  result <- list(
    pcdrs = pcdrs,
    coefficients_used = if(exists("coefficients")) coefficients else NULL,
    genes_used = common_genes,
    method = method
  )
  
  return(result)
}

# Function for time-dependent ROC analysis
perform_timeROC_analysis <- function(risk_score, clinical_data, time_points = c(1, 3, 5, 10)) {
  common_ids <- intersect(names(risk_score), rownames(clinical_data))
  if (length(common_ids) == 0) {
    stop("No overlapping patient IDs between risk_score and clinical_data.")
  }
  
  risk_score <- risk_score[common_ids]
  clinical_data <- clinical_data[common_ids, , drop = FALSE]
  
  # Convert time from days to years
  time_years <- clinical_data$futime / 365.25
  
  # Check PCDRS direction before analysis
  # timeROC expects higher values = higher risk
  # If PCDRS is inverted, we need to detect and correct it
  
  # Quick direction check for PCDRS
  valid_for_check <- !is.na(risk_score) & !is.na(time_years) & !is.na(clinical_data$fustat)
  if(sum(valid_for_check) > 20) {
    # Check correlation with survival time
    cor_test <- cor.test(risk_score[valid_for_check], time_years[valid_for_check], method = "spearman")
    
    # Check event rates in high vs low risk groups
    high_risk <- risk_score >= median(risk_score, na.rm = TRUE)
    low_risk <- risk_score < median(risk_score, na.rm = TRUE)
    
    high_risk_events <- mean(clinical_data$fustat[high_risk & valid_for_check])
    low_risk_events <- mean(clinical_data$fustat[low_risk & valid_for_check])
    
    # If PCDRS shows wrong direction, invert it
    if(cor_test$estimate > 0.1 || high_risk_events < low_risk_events) {
      ## Phase 2: report only - no outcome-driven sign flips (they bias AUC/C upward)
      message("NOTE: PCDRS direction looks inverted in this cohort (no flip applied)")
      risk_score_corrected <- risk_score
    } else {
      risk_score_corrected <- risk_score
    }
  } else {
    risk_score_corrected <- risk_score
  }
  
  # Prepare predictors list
  predictors <- list(
    CMPIS = risk_score_corrected,
    Age   = if("Age" %in% colnames(clinical_data)) clinical_data$Age else NULL,
    
    # Stage predictor (convert to numeric)
    Stage = if("Stage" %in% colnames(clinical_data)) {
      stage_map <- c("Stage 0"=0, "Stage I"=1, "Stage II"=2, "Stage III"=3, "Stage IV"=4)
      stage_num <- stage_map[as.character(clinical_data$Stage)]
      stage_num
    } else NULL,
    
    # TNM staging components
    T_stage = if("T" %in% colnames(clinical_data)) {
      t_map <- c("T1"=1, "T2"=2, "T3"=3, "T4"=4)
      t_num <- t_map[as.character(clinical_data$T)]
      t_num
    } else NULL,
    
    N_stage = if("N" %in% colnames(clinical_data)) {
      n_map <- c("NX"=0, "N0"=0, "N1"=1, "N2"=2, "N3"=3)
      n_num <- n_map[as.character(clinical_data$N)]
      n_num
    } else NULL,
    
    M_stage = if("M" %in% colnames(clinical_data)) {
      m_map <- c("M0"=0, "M1"=1)
      m_num <- m_map[as.character(clinical_data$M)]
      m_num
    } else NULL,
    
    # Gender (convert to numeric: 0 = female, 1 = male)
    Gender = if("Gender" %in% colnames(clinical_data)) {
      gender_num <- ifelse(clinical_data$Gender == "male", 1, 0)
      gender_num
    } else NULL
  )
  
  # Remove NULL predictors
  predictors <- predictors[!sapply(predictors, is.null)]
  
  # Initialize results
  roc_results <- list()
  auc_results <- data.frame(Time = time_points)
  
  # Initialize AUC columns for all predictors
  for(pred_name in names(predictors)) {
    auc_results[[paste0(pred_name, "_AUC")]] <- numeric(length(time_points))
  }
  
  # Loop over time points
  for(i in seq_along(time_points)) {
    t <- time_points[i]
    roc_results[[paste0("Year_", t)]] <- list()
    
    for(pred_name in names(predictors)) {
      marker <- predictors[[pred_name]]
      
      # Skip if marker is NULL
      if(is.null(marker)) next
      
      # Handle NA values - only use complete cases for this predictor
      valid_idx <- !is.na(marker) & !is.na(time_years) & !is.na(clinical_data$fustat)
      
      if(sum(valid_idx) < 10) {  # Need minimum samples for ROC analysis
        warning(paste("Not enough valid samples for", pred_name, "at time", t))
        auc_results[i, paste0(pred_name, "_AUC")] <- NA
        next
      }
      
      tryCatch({
        roc_obj <- timeROC(
          T = time_years[valid_idx],
          delta = clinical_data$fustat[valid_idx],
          marker = marker[valid_idx],
          cause = 1,
          times = t,
          iid = TRUE
        )
        
        auc_results[i, paste0(pred_name, "_AUC")] <- roc_obj$AUC[2]
        roc_results[[paste0("Year_", t)]][[pred_name]] <- roc_obj
        
      }, error = function(e) {
        warning(paste("Error calculating ROC for", pred_name, "at time", t, ":", e$message))
        auc_results[i, paste0(pred_name, "_AUC")] <- NA
      })
    }
  }
  
  return(list(auc_results = auc_results, roc_results = roc_results))
}

# Plot time-dependent ROC curves
# Plot time-dependent ROC curves based on your perform_timeROC_analysis function
plot_timeROC <- function(roc_results, time_points = c(1, 3, 5, 10)) {
  
  # Extract AUC results
  auc_data <- roc_results$auc_results
  
  # Check auc_data structure
  cat("AUC data structure:\n")
  str(auc_data)
  
  if(!is.data.frame(auc_data)) {
    stop("auc_results should be a data.frame")
  }
  
  # Create AUC comparison plot
  auc_melted <- melt(auc_data, id.vars = "Time", variable.name = "Feature", value.name = "AUC")
  
  p1 <- ggplot(auc_melted, aes(x = factor(Time), y = AUC, fill = Feature)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(title = "Time-dependent AUC Comparison",
         x = "Years", y = "AUC Value",
         fill = "Features") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom") +
    ylim(0, 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7)
  
  # Create individual ROC curve plots for each time point
  roc_plots <- list()
  
  for(i in 1:length(time_points)) {
    t <- time_points[i]
    roc_key <- paste0("Year_", t)
    
    roc_data <- roc_results$roc_results[[roc_key]]
    
    if(!is.null(roc_data)) {
      
      # Extract ROC curve data from timeROC object
      if(!is.null(roc_data)) {
        roc_df <- data.frame()
        
        for(pred_name in names(roc_data)) {
          roc_obj <- roc_data[[pred_name]]
          temp_df <- data.frame(FPR = roc_obj$FP[,2], 
                                TPR = roc_obj$TP[,2], 
                                Predictor = paste0(pred_name, " (AUC=", round(roc_data[[pred_name]]$AUC[2], 3), ")"))
          roc_df <- rbind(roc_df, temp_df)
        }
        
        # Create ROC plot
        roc_plots[[paste0(t, "_year")]] <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Predictor)) +
          geom_line(size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", alpha = 0.7) +
          labs(title = paste0(t, "-Year ROC Curves"), x = "False Positive Rate (1-Specificity)", 
               y = "True Positive Rate (Sensitivity)", color = "Predictors (AUC)") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                legend.position = "inside",
                legend.justification = c("right", "bottom")) +
          xlim(0, 1) + ylim(0, 1) +
          coord_fixed()
      } else {
        cat("Missing required components in ROC data. Available:", names(roc_data), "\n")
      }
    } else {
      cat("No ROC data found for", roc_key, "\n")
    }
  }
  
  # Remove NULL elements from roc_plots
  roc_plots <- roc_plots[!sapply(roc_plots, is.null)]
  
  # Combine all plots
  if(length(roc_plots) > 0) {
    return(list(
      auc_comparison = p1,
      individual_plots = roc_plots,
      n_plots = length(roc_plots)
    ))
  } else {
    warning("No ROC curves could be created")
    return(list(auc_comparison = p1))
  }
}

# Univariate and Multivariate Cox regression analysis
perform_cox_analysis <- function(clinical_data, risk_score) {
  common_ids <- intersect(names(risk_score), rownames(clinical_data))
  if (length(common_ids) == 0) {
    stop("No overlapping patient IDs between risk_score and clinical_data.")
  }
  
  risk_score <- risk_score[common_ids]
  clinical_data <- clinical_data[common_ids, , drop = FALSE]
  
  # Add risk score to clinical data
  clinical_data$PCDRS <- risk_score
  clinical_data$CMPIS <- risk_score
  clinical_data$PCDRS_group <- ifelse(risk_score > median(risk_score), "High", "Low")
  
  # Convert time from days to years for Cox analysis (optional - can keep in days)
  # clinical_data$futime_years <- clinical_data$futime / 365.25
  
  # Define all possible variables (adjust based on your column names)
  variables <- c("CMPIS", "Age", "Gender", "Stage", "T", "N")
  
  # Univariate Cox regression
  univariate_results <- data.frame(
    Variable = character(),
    HR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(var in variables) {
    if(var %in% colnames(clinical_data)) {
      tryCatch({
        # Handle different variable types
        if(var %in% c("Stage", "T", "N", "M", "Gender", "CS")) {
          # For categorical variables, convert to factor
          clinical_data[[paste0(var, "_factor")]] <- as.factor(clinical_data[[var]])
          formula <- as.formula(paste("Surv(futime, fustat) ~", paste0(var, "_factor")))
        } else {
          # For continuous variables (PCDRS, Age)
          formula <- as.formula(paste("Surv(futime, fustat) ~", var))
        }
        
        # Remove rows with NA for this specific variable
        complete_data <- clinical_data[!is.na(clinical_data[[var]]), ]
        
        if(nrow(complete_data) < 10) {
          warning(paste("Not enough complete cases for", var))
          next
        }
        
        cox_model <- coxph(formula, data = complete_data)
        summary_cox <- summary(cox_model)
        
        # For categorical variables with multiple levels, get overall p-value
        if(var %in% c("Stage", "T", "N", "Gender", "CS") && nrow(summary_cox$coefficients) > 1) {
          # Use likelihood ratio test for overall significance
          null_model <- coxph(Surv(futime, fustat) ~ 1, data = complete_data)
          lr_test <- anova(null_model, cox_model, test = "Chisq")
          p_value <- lr_test$`Pr(>Chi)`[2]
          
          # Use first level's HR as representative (or could use median)
          hr_value <- summary_cox$coefficients[1, "exp(coef)"]
          ci_lower <- summary_cox$conf.int[1, "lower .95"]
          ci_upper <- summary_cox$conf.int[1, "upper .95"]
        } else {
          # Single coefficient
          hr_value <- summary_cox$coefficients[1, "exp(coef)"]
          ci_lower <- summary_cox$conf.int[1, "lower .95"]
          ci_upper <- summary_cox$conf.int[1, "upper .95"]
          p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
        }
        
        univariate_results <- rbind(univariate_results, data.frame(
          Variable = var,
          HR = hr_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper,
          P_value = p_value,
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        warning(paste("Error in univariate Cox analysis for", var, ":", e$message))
      })
    }
  }
  
  # Multivariate Cox regression
  available_vars <- intersect(variables, colnames(clinical_data))
  
  # Remove variables with too many missing values
  vars_to_include <- c()
  for(var in available_vars) {
    if(sum(!is.na(clinical_data[[var]])) >= 10) {
      vars_to_include <- c(vars_to_include, var)
    }
  }
  
  if(length(vars_to_include) > 1) {
    # Create complete case dataset for multivariate analysis
    vars_used <- c("futime", "fustat", vars_to_include)
    clinical_complete <- clinical_data[complete.cases(clinical_data[, vars_used]), ]
    
    if(nrow(clinical_complete) >= 20) {  # Minimum cases for multivariate
      tryCatch({
        # Convert categorical variables to factors for multivariate
        for(var in vars_to_include) {
          if(var %in% c("Stage", "T", "N", "Gender")) {
            clinical_complete[[var]] <- as.factor(clinical_complete[[var]])
          }
        }
        
        multivariate_formula <- as.formula(paste("Surv(futime, fustat) ~", 
                                                 paste(vars_to_include, collapse = " + ")))
        
        multivariate_cox <- coxph(multivariate_formula, data = clinical_complete)
        
        # Global likelihood ratio tests
        global_tests <- drop1(multivariate_cox, test = "Chisq")
        
        # Extract summary
        cox_sum <- summary(multivariate_cox)
        
        # Build results table
        coef_names <- rownames(cox_sum$coefficients)
        multivariate_summary <- data.frame(
          Variable = coef_names,
          HR = cox_sum$coefficients[, "exp(coef)"],
          CI_lower = cox_sum$conf.int[, "lower .95"],
          CI_upper = cox_sum$conf.int[, "upper .95"],
          P_value = cox_sum$coefficients[, "Pr(>|z|)"],
          stringsAsFactors = FALSE
        )
        
        # For categorical variables, replace individual level p-values with global p-values
        for(var in c("Stage", "T", "N", "M", "Gender", "CS")) {
          if(var %in% rownames(global_tests)) {
            # Find rows that belong to this variable
            var_rows <- grep(paste0("^", var), multivariate_summary$Variable)
            if(length(var_rows) > 1) {
              # Replace with single row showing global test
              other_rows <- multivariate_summary[-var_rows, ]
              global_row <- data.frame(
                Variable = var,
                HR = NA,
                CI_lower = NA,
                CI_upper = NA,
                P_value = global_tests[var, "Pr(>Chi)"],
                stringsAsFactors = FALSE
              )
              multivariate_summary <- rbind(other_rows, global_row)
            } else if(length(var_rows) == 1) {
              # Update p-value for single level variables
              multivariate_summary$P_value[var_rows] <- global_tests[var, "Pr(>Chi)"]
            }
          }
        }
        
      }, error = function(e) {
        warning(paste("Error in multivariate Cox analysis:", e$message))
        multivariate_summary <- data.frame(
          Variable = character(),
          HR = numeric(),
          CI_lower = numeric(),
          CI_upper = numeric(),
          P_value = numeric(),
          stringsAsFactors = FALSE
        )
      })
      
    } else {
      warning("Not enough complete cases for multivariate analysis")
      multivariate_summary <- data.frame(
        Variable = character(),
        HR = numeric(),
        CI_lower = numeric(),
        CI_upper = numeric(),
        P_value = numeric(),
        stringsAsFactors = FALSE
      )
    }
  } else {
    warning("Not enough variables available for multivariate analysis")
    multivariate_summary <- data.frame(
      Variable = character(),
      HR = numeric(),
      CI_lower = numeric(),
      CI_upper = numeric(),
      P_value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    univariate = univariate_results,
    multivariate = if(exists("multivariate_summary")) multivariate_summary else NULL,
    multivariate_model = if(exists("multivariate_cox")) multivariate_cox else NULL,
    clinical_data = clinical_data,
    n_complete_cases = if(exists("clinical_complete")) nrow(clinical_complete) else 0
  ))
}

c_index_calculation <- function(clinical_data, predictors = NULL, risk_score, time_points = c(1, 2, 3, 4, 5, 10, 15)) {
  library(survcomp)
  
  # Get common patient IDs
  common_ids <- intersect(names(risk_score), rownames(clinical_data))
  if (length(common_ids) == 0) {
    stop("No overlapping patient IDs between risk_score and clinical_data.")
  }
  
  # Subset to common IDs
  risk_score <- risk_score[common_ids]
  clinical_data <- clinical_data[common_ids, , drop = FALSE]
  
  # Add risk score to clinical data
  clinical_data$PCDRS <- risk_score
  
  # Define default predictors if not provided
  if(is.null(predictors)) {
    available_cols <- colnames(clinical_data)
    predictors <- c("PCDRS")
    
    # Add available clinical predictors
    if("Age" %in% available_cols) predictors <- c(predictors, "Age")
    if("Gender" %in% available_cols) predictors <- c(predictors, "Gender")
    if("Stage" %in% available_cols) predictors <- c(predictors, "Stage")
    if("T" %in% available_cols) predictors <- c(predictors, "T")
    if("N" %in% available_cols) predictors <- c(predictors, "N")
    if("M" %in% available_cols) predictors <- c(predictors, "M")
    if("CS" %in% available_cols) predictors <- c(predictors, "CS")
  }
  
  # Convert survival time to years if in days
  time_years <- clinical_data$futime / 365.25
  fustat <- clinical_data$fustat
  
  # Check the direction of PCDRS association with survival
  # Higher risk scores should be associated with shorter survival times
  if("PCDRS" %in% predictors) {
    valid_pcdrs <- !is.na(clinical_data$PCDRS) & !is.na(time_years) & !is.na(fustat)
    if(sum(valid_pcdrs) > 10) {
      cor_test <- cor.test(clinical_data$PCDRS[valid_pcdrs], time_years[valid_pcdrs], method = "spearman")
      if(cor_test$estimate > 0) {
        message("PCDRS shows positive correlation with survival time - this may indicate the risk score direction needs adjustment")
      }
    }
  }
  
  # Prepare output data frame
  cindex_df <- data.frame()
  
  # Create numeric versions of categorical predictors
  clinical_data_numeric <- clinical_data
  
  # Convert categorical variables to numeric for C-index calculation
  if("Gender" %in% predictors && "Gender" %in% colnames(clinical_data)) {
    clinical_data_numeric$Gender_numeric <- ifelse(clinical_data$Gender == "male", 1, 0)
  }
  
  if("Stage" %in% predictors && "Stage" %in% colnames(clinical_data)) {
    stage_map <- c("Stage 0"=0, "Stage I"=1, "Stage II"=2, "Stage III"=3, "Stage IV"=4)
    clinical_data_numeric$Stage_numeric <- stage_map[as.character(clinical_data$Stage)]
  }
  
  if("T" %in% predictors && "T" %in% colnames(clinical_data)) {
    t_map <- c("T1"=1, "T2"=2, "T3"=3, "T4"=4)
    clinical_data_numeric$T_numeric <- t_map[as.character(clinical_data$T)]
  }
  
  if("N" %in% predictors && "N" %in% colnames(clinical_data)) {
    n_map <- c("NX"=0, "N0"=0, "N1"=1, "N2"=2, "N3"=3)
    clinical_data_numeric$N_numeric <- n_map[as.character(clinical_data$N)]
  }
  
  if("M" %in% predictors && "M" %in% colnames(clinical_data)) {
    m_map <- c("M0"=0, "M1"=1)
    clinical_data_numeric$M_numeric <- m_map[as.character(clinical_data$M)]
  }
  
  if("CS" %in% predictors && "CS" %in% colnames(clinical_data)) {
    cs_map <- c("CS1"=1, "CS2"=2)
    clinical_data_numeric$CS_numeric <- cs_map[as.character(clinical_data$CS)]
  }
  
  # Calculate C-index for each predictor at each time point
  for(pred in predictors) {
    # Determine which column to use
    if(pred %in% c("PCDRS", "Age")) {
      x <- clinical_data_numeric[[pred]]
      pred_name <- pred
    } else if(pred %in% c("Gender", "Stage", "T", "N", "M", "CS")) {
      numeric_col <- paste0(pred, "_numeric")
      if(numeric_col %in% colnames(clinical_data_numeric)) {
        x <- clinical_data_numeric[[numeric_col]]
        pred_name <- pred
      } else {
        warning(paste("Could not find numeric version of", pred))
        next
      }
    } else {
      # Try original column
      if(pred %in% colnames(clinical_data_numeric)) {
        x <- clinical_data_numeric[[pred]]
        if(is.factor(x) || is.character(x)) {
          x <- as.numeric(as.factor(x))
        }
        pred_name <- pred
      } else {
        warning(paste("Predictor", pred, "not found in clinical data"))
        next
      }
    }
    
    for(t in time_points) {
      tryCatch({
        # Create indices for complete cases
        idx <- !is.na(x) & !is.na(time_years) & !is.na(fustat)
        
        if(sum(idx) < 15) {  # Increased minimum for more stable estimates
          warning(paste("Not enough complete cases for", pred_name, "at time", t, "years"))
          next
        }
        
        # Improved time-specific C-index calculation
        # Only include patients who either:
        # 1. Had an event by time t, OR
        # 2. Were followed for at least time t without an event
        
        time_subset <- time_years[idx]
        event_subset <- fustat[idx]
        x_subset <- x[idx]
        
        # Create time-restricted dataset
        # Patients who die by time t: keep their actual time and event
        # Patients who survive past t: censor them at time t
        time_restricted <- pmin(time_subset, t)
        event_restricted <- ifelse(time_subset <= t, event_subset, 0)
        
        # For more stable C-index, require minimum follow-up or event
        min_followup <- min(0.25, t/4)  # At least 3 months or 1/4 of target time
        include_patient <- (time_subset >= min_followup) | (event_subset == 1)
        
        final_time <- time_restricted[include_patient]
        final_event <- event_restricted[include_patient]
        final_x <- x_subset[include_patient]
        
        # Check we have enough events and variability
        if(length(final_x) < 15 || sum(final_event) < 3 || var(final_x) == 0) {
          warning(paste("Insufficient data for", pred_name, "at time", t, "years:",
                        "N =", length(final_x), "Events =", sum(final_event)))
          next
        }
        
        # For PCDRS, we might need to invert if it's a "protective" score
        # Most risk scores: higher values = higher risk = lower survival
        # But some scores might be inverted
        x_for_cindex <- final_x
        
        # Optional: Auto-detect and correct direction for PCDRS
        if(pred_name == "PCDRS") {
          # Check if higher scores are associated with better outcomes (wrong direction)
          high_risk_group <- final_x >= median(final_x)
          low_risk_group <- final_x < median(final_x)
          
          if(sum(high_risk_group) > 2 && sum(low_risk_group) > 2) {
            high_risk_events <- mean(final_event[high_risk_group])
            low_risk_events <- mean(final_event[low_risk_group])
            
            # If high risk group has fewer events, the score might be inverted
            if(high_risk_events < low_risk_events) {
              message(paste("PCDRS appears inverted at time", t, "- consider using -PCDRS"))
              # Uncomment next line to auto-correct:
              # x_for_cindex <- -final_x
            }
          }
        }
        
        ci <- concordance.index(
          x = x_for_cindex,
          surv.time = final_time,
          surv.event = final_event,
          method = "noether"
        )
        
        cindex_df <- rbind(cindex_df, data.frame(
          Predictor = pred_name,
          Time = t,
          C_index = ci$c.index,
          C_index_lower = ci$lower,
          C_index_upper = ci$upper,
          P_value = ci$p.value,
          N_cases = length(final_x),
          N_events = sum(final_event),
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        warning(paste("Error calculating C-index for", pred_name, "at time", t, ":", e$message))
      })
    }
  }
  
  return(cindex_df)
}

plot_cindex_over_time <- function(cindex_df, title = "C-index Performance Over Time") {
  
  # Ensure we have the required columns
  if(!all(c("Predictor", "Time", "C_index") %in% colnames(cindex_df))) {
    stop("DataFrame must contain columns: Predictor, Time, C_index")
  }
  
  # Create the main plot
  p1 <- ggplot(cindex_df, aes(x = Time, y = C_index, color = Predictor)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    
    # Add confidence intervals if available
    {if("C_index_lower" %in% colnames(cindex_df) & "C_index_upper" %in% colnames(cindex_df)) {
      geom_ribbon(aes(ymin = C_index_lower, ymax = C_index_upper, fill = Predictor), 
                  alpha = 0, color = NA)
    }} +
    
    # Styling
    scale_color_brewer(type = "qual", palette = "Set2") +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    
    # Add reference line at C-index = 0.5 (random performance)
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
    
    # Customize axes
    scale_x_continuous(breaks = unique(cindex_df$Time), 
                       name = "Time (Years)") +
    scale_y_continuous(limits = c(0.4, 1.0), 
                       breaks = seq(0.4, 1.0, 0.1),
                       name = "C-index") +
    
    # Theme and labels
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.5),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    
    ggtitle(title)
  
  # Create summary table
  summary_table <- cindex_df %>%
    group_by(Predictor) %>%
    summarise(
      Mean_Cindex = round(mean(C_index, na.rm = TRUE), 3),
      Max_Cindex = round(max(C_index, na.rm = TRUE), 3),
      Min_Cindex = round(min(C_index, na.rm = TRUE), 3),
      Stability = round(max(C_index, na.rm = TRUE) - min(C_index, na.rm = TRUE), 3),
      .groups = 'drop'
    ) %>%
    arrange(desc(Mean_Cindex))
  
  print("C-index Summary Table:")
  print(summary_table)
  
  return(p1)
}

create_nomogram_analysis <- function(clinical_data, risk_score, 
                                     significant_predictors = c("PCDRS", "Age", "Stage"),
                                     time_points = c(1, 3, 5),
                                     title_prefix = "MPCDS") {
  
  # Get common patient IDs
  common_ids <- intersect(names(risk_score), rownames(clinical_data))
  if (length(common_ids) == 0) {
    stop("No overlapping patient IDs between risk_score and clinical_data.")
  }
  
  # Subset to common IDs
  risk_score <- risk_score[common_ids]
  clinical_data <- clinical_data[common_ids, , drop = FALSE]
  
  # Add risk score to clinical data
  clinical_data$PCDRS <- risk_score
  clinical_data$CMPIS <- risk_score
  
  # Convert survival time to years
  clinical_data$futime_years <- clinical_data$futime / 365.25
  
  # Prepare data for modeling
  model_vars <- c("futime_years", "fustat", significant_predictors)
  complete_data <- clinical_data[complete.cases(clinical_data[, model_vars]), ]
  complete_data$Age <- max(complete_data$Age, na.rm=TRUE) - complete_data$Age
  
  if(nrow(complete_data) < 50) {
    stop("Not enough complete cases for nomogram development")
  }
  
  cat("Using", nrow(complete_data), "complete cases for nomogram development\n")
  
  # DIAGNOSTIC: Check predictor distributions
  cat("\n=== DIAGNOSTIC: Predictor Variable Summaries ===\n")
  for(var in significant_predictors) {
    cat(paste("\n", var, ":\n"))
    if(is.numeric(complete_data[[var]])) {
      summary_stats <- summary(complete_data[[var]])
      print(summary_stats)
      cat("Range:", diff(range(complete_data[[var]], na.rm = TRUE)), "\n")
      cat("SD:", sd(complete_data[[var]], na.rm = TRUE), "\n")
      
      # Check for potential issues
      if(diff(range(complete_data[[var]], na.rm = TRUE)) < 0.001) {
        warning(paste("Variable", var, "has very small range - this may cause nomogram issues"))
      }
    } else {
      freq_table <- table(complete_data[[var]])
      print(freq_table)
      if(length(freq_table) < 2) {
        warning(paste("Variable", var, "has only one level - remove from predictors"))
      }
    }
  }
  
  # Convert categorical variables to factors with proper handling
  if("Stage" %in% significant_predictors) {
    complete_data$Stage <- as.factor(complete_data$Stage)
    # Ensure at least 2 levels
    if(length(levels(complete_data$Stage)) < 2) {
      warning("Stage has < 2 levels, removing from predictors")
      significant_predictors <- significant_predictors[significant_predictors != "Stage"]
    }
  }
  
  # Handle other categorical variables similarly
  cat_vars <- c("Gender", "T", "N", "M", "CS")
  for(var in cat_vars) {
    if(var %in% significant_predictors) {
      complete_data[[var]] <- as.factor(complete_data[[var]])
      if(length(levels(complete_data[[var]])) < 2) {
        warning(paste(var, "has < 2 levels, removing from predictors"))
        significant_predictors <- significant_predictors[significant_predictors != var]
      }
    }
  }
  complete_data$Age <- as.numeric(complete_data$Age)
  dd <- datadist(complete_data[, significant_predictors, drop = FALSE])
  dd$Age$Low <- min(complete_data$Age, na.rm = TRUE)
  dd$Age$High <- max(complete_data$Age, na.rm = TRUE)
  options(datadist = dd)
  
  # Create survival object
  complete_data$surv_obj <- with(complete_data, Surv(futime_years, fustat))
  
  # Build the Cox model with updated predictors
  if(length(significant_predictors) == 0) {
    stop("No valid predictors remaining after filtering")
  }
  
  formula_str <- paste("surv_obj ~", paste(significant_predictors, collapse = " + "))
  cox_formula <- as.formula(formula_str)
  
  cat("\nBuilding Cox model with formula:", formula_str, "\n")
  
  cox_model <- cph(cox_formula, data = complete_data, 
                   x = TRUE, y = TRUE, surv = TRUE, 
                   time.inc = max(time_points))
  
  # Print model summary
  cat("\n=== Cox Model Summary ===\n")
  print(summary(cox_model))
  
  # DIAGNOSTIC: Check model coefficients
  cat("\n=== DIAGNOSTIC: Model Coefficients ===\n")
  print(cox_model$coefficients)
  
  # Create nomogram
  cat("\n=== Creating Nomogram ===\n")
  
  surv_fun <- Survival(cox_model)
  
  # Create survival probability functions
  surv_probs <- lapply(time_points, function(t) {
    function(x) surv_fun(t, x)
  })
  names(surv_probs) <- paste0("S", time_points, "yr")
  
  # Create the nomogram with better parameters
  nom <- nomogram(cox_model, 
                  fun = surv_probs,
                  fun.at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                  funlabel = paste("Survival Prob.", time_points, "Years"),
                  lp = TRUE,  # Include linear predictor
                  vnames = "labels")
  
  # DIAGNOSTIC: Print nomogram object
  cat("\n=== DIAGNOSTIC: Nomogram Object ===\n")
  str(nom)
  
  # Create improved plot
  png(paste0(title_prefix, "_nomogram_improved.png"), width = 1400, height = 1000, res = 150)
  par(mar = c(5, 8, 4, 2) + 0.1)  # Increase left margin
  
  plot(nom, 
       xfrac = 0.5,           # More space for variable names
       cex.axis = 0.85,       # Readable axis labels
       cex.var = 0.95,        # Readable variable names  
       col.grid = "gray85",   # Subtle grid
       lmgp = 0.2,           # Line margin spacing
       points.at = TRUE)      # Show points on scales
  
  title(paste(title_prefix, "Nomogram for Survival Prediction"), 
        cex.main = 1.3, font.main = 2, line = 2)
  dev.off()
  
  # Display nomogram in current device
  par(mar = c(5, 8, 4, 2) + 0.1)
  plot(nom, 
       xfrac = 0.5,           
       cex.axis = 0.85,       
       cex.var = 0.95,        
       col.grid = "gray85",   
       lmgp = 0.2,           
       points.at = TRUE)
  title(paste(title_prefix, "Nomogram for Survival Prediction"), 
        cex.main = 1.3, font.main = 2, line = 2)
  
  # Calibration analysis
  cat("\n=== Performing Calibration Analysis ===\n")
  
  calibration_results <- list()
  calibration_plots <- list()
  
  for(t in time_points) {
    cat("Calibration for", t, "year survival...\n")
    
    tryCatch({
      # Calculate predicted survival probabilities
      pred_surv <- survest(cox_model, newdata = complete_data, times = t)$surv
      
      # Create calibration plot
      cal_data <- data.frame(
        predicted = pred_surv,
        observed_time = complete_data$futime_years,
        observed_event = complete_data$fustat
      )
      
      # Group patients by predicted risk (deciles)
      cal_data$risk_group <- cut(1 - pred_surv, quantile(1 - pred_surv, probs = seq(0, 1, 0.1)), 
                                 include.lowest = TRUE, labels = 1:10)
      
      # Calculate observed survival in each group
      cal_summary <- cal_data %>%
        group_by(risk_group) %>%
        summarise(
          n = n(),
          predicted_surv = mean(predicted),
          observed_surv = sum(observed_time >= t | (observed_time < t & observed_event == 0)) / n(),
          mean_predicted_risk = mean(1 - predicted),
          .groups = 'drop'
        ) %>%
        filter(n >= 5)  # Only groups with at least 5 patients
      
      # Create calibration plot
      cal_plot <- ggplot(cal_summary, aes(x = predicted_surv, y = observed_surv)) +
        geom_point(aes(size = n), alpha = 0.7, color = "blue") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
        geom_smooth(method = "loess", se = TRUE, color = "darkblue", alpha = 0.3) +
        scale_x_continuous(limits = c(0, 1), name = "Nomogram-Predicted Survival Probability") +
        scale_y_continuous(limits = c(0, 1), name = "Observed Survival Probability") +
        scale_size_continuous(name = "N patients", range = c(2, 8)) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom"
        ) +
        ggtitle(paste("Calibration Plot -", t, "Year Survival")) +
        annotate("text", x = 0.2, y = 0.8, 
                 label = paste("Perfect calibration"), 
                 color = "red", size = 3.5, hjust = 0)
      
      calibration_plots[[paste0("year_", t)]] <- cal_plot
      calibration_results[[paste0("year_", t)]] <- cal_summary
      
      # Print calibration statistics
      cat("Calibration summary for", t, "years:\n")
      cat("Mean absolute error:", round(mean(abs(cal_summary$predicted_surv - cal_summary$observed_surv)), 3), "\n")
      cat("Calibration slope (should be ~1):", 
          round(coef(lm(observed_surv ~ predicted_surv, data = cal_summary))[2], 3), "\n\n")
      
    }, error = function(e) {
      warning(paste("Error in calibration for", t, "years:", e$message))
    })
  }
  
  # Combine calibration plots
  if(length(calibration_plots) > 1) {
    combined_cal_plot <- do.call(grid.arrange, c(calibration_plots, ncol = length(calibration_plots)))
    
    # Save combined plot
    ggsave(paste0(title_prefix, "_calibration_plots.png"), combined_cal_plot, 
           width = 4 * length(time_points), height = 4, dpi = 300)
  } else if(length(calibration_plots) == 1) {
    ggsave(paste0(title_prefix, "_calibration_plot.png"), calibration_plots[[1]], 
           width = 6, height = 5, dpi = 300)
  }
  
  # Display calibration plots
  if(length(calibration_plots) > 0) {
    for(i in 1:length(calibration_plots)) {
      print(calibration_plots[[i]])
    }
  }
  
  # Model performance metrics
  cat("\n=== Model Performance Metrics ===\n")
  
  # C-index
  c_index <- cox_model$stats["C Index"]
  cat("C-index:", round(c_index, 3), "\n")
  
  # AIC and likelihood ratio test
  cat("AIC:", round(AIC(cox_model), 2), "\n")
  cat("Likelihood ratio chi-square:", round(cox_model$stats["Model L.R."], 2), "\n")
  cat("p-value:", format.pval(cox_model$stats["P"], digits = 3), "\n")
  
  # Return results
  return(list(
    cox_model = cox_model,
    nomogram = nom,
    calibration_results = calibration_results,
    calibration_plots = calibration_plots,
    model_data = complete_data,
    c_index = c_index
  ))
}


## Phase 2: weights + sign come from the single CMPIS helper (corrected weights need no inversion)
source("R/cmpis_weights.R"); .cmpis <- load_cmpis()
pcdrs_scores <- calculate_pcdrs(skcm.tcga$mRNA.expr, pcdrg_genes, .cmpis$coef,
                                risk_direction = if (.cmpis$sign < 0) "invert" else "none")
roc_analysis <- perform_timeROC_analysis(pcdrs_scores$pcdrs, skcm.tcga$clin.info)
plots <- plot_timeROC(roc_analysis)
plots$auc_comparison
grid.arrange(
  grobs = plots$individual_plots,
  nrow = 2,
  ncol = 2
)

cox_results <- perform_cox_analysis(skcm.tcga$clin.info, pcdrs_scores$pcdrs)

df <- cox_results$univariate
tabletext <- cbind(
  c("Variable", df$Variable),
  c("p-value", round(df$P_value,6)),
  c("Hazard ratio (95% CI)", paste0(round(df$HR,3), " (", round(df$CI_lower,3), "-", round(df$CI_upper,3), ")"))
)

forestplot(tabletext,
           mean  = c(NA, df$HR),
           lower = c(NA, df$CI_lower),
           upper = c(NA, df$CI_upper),
           title = "Univariate Cox Regression Analysis",
           xlab = "Hazard Ratio",
           zero = 1,
           cex  = 2,
           lineheight = "auto",
           boxsize = 0.1,
           vertices = TRUE,
           col = fpColors(box="royalblue", line="darkblue", summary="blue"))

df <- cox_results$multivariate

# Check if multivariate results exist and have data
if(is.null(df) || nrow(df) == 0) {
  stop("No multivariate Cox results available")
}

# Handle rows with NA HR values (categorical variables with global p-values)
# Option 1: Remove rows with NA HR values
df_plot <- df[!is.na(df$HR), ]

# Or Option 2: Keep all rows and handle NAs in plotting
# df_plot <- df
# df_plot$HR[is.na(df_plot$HR)] <- 1  # Set to 1 for display
# df_plot$CI_lower[is.na(df_plot$CI_lower)] <- 1
# df_plot$CI_upper[is.na(df_plot$CI_upper)] <- 1

if(nrow(df_plot) == 0) {
  stop("No plottable results (all HR values are NA)")
}

# Create table text with proper column names from your updated function
tabletext <- cbind(
  c("Variable", df_plot$Variable),
  c("p-value", formatC(df_plot$P_value, format = "e", digits = 2)),  # Use scientific notation for very small p-values
  c("Hazard ratio (95% CI)", 
    ifelse(is.na(df_plot$HR), 
           "Global test", 
           paste0(round(df_plot$HR, 3), " (", 
                  round(df_plot$CI_lower, 3), "-", 
                  round(df_plot$CI_upper, 3), ")")))
)

# Load required library
library(forestplot)

# Create forest plot
forestplot(tabletext,
           mean  = c(NA, df_plot$HR),
           lower = c(NA, df_plot$CI_lower),
           upper = c(NA, df_plot$CI_upper),
           title = "Multivariate Cox Regression Analysis",
           xlab = "Hazard Ratio",
           zero = 1,
           cex  = 1.2,  # Reduced from 2 for better readability
           lineheight = "auto",
           boxsize = 0.15,  # Slightly larger boxes
           vertices = TRUE,
           col = fpColors(box="royalblue", 
                          line="darkblue", 
                          summary="blue"),
           # Add confidence interval lines
           ci.vertices = TRUE,
           ci.vertices.height = 0.1,
           # Format x-axis
           xticks = c(0.5, 1, 2),  # Adjust based on your HR range
           grid = TRUE)

cindex_results <- c_index_calculation(clinical_data = skcm.tcga$clin.info, 
                                      risk_score = pcdrs_scores$pcdrs)  # Phase 2: higher = higher risk, no negation

filtered_results <- cindex_results %>%
  dplyr::filter(Predictor %in% c("PCDRS", "Age", "Gender", "Stage"))

cindex_plot <- plot_cindex_over_time(filtered_results, "SKCM C-index Performance Over Time")
cindex_plot

nomogram_results <- create_nomogram_analysis(
  clinical_data = skcm.tcga$clin.info,
  risk_score = pcdrs_scores$pcdrs,  # Use your PCDRS scores
  significant_predictors = c("CMPIS", "Age", "Stage", "Gender", "CS"),
  time_points = c(1, 3, 5),  # 1, 3, and 5-year survival
  title_prefix = "MPCDS"
)

expr <- skcm.tcga$mRNA.expr
clin <- cox_results$clinical_data

# 1) Align samples between expression and clinical data
common <- intersect(colnames(expr), rownames(clin))
expr   <- expr[, common, drop = FALSE]
clin   <- clin[common, , drop = FALSE]

# 2) Build design (two-group example; change to your column name)
clin$group <- factor(clin$PCDRS_group)  # e.g., levels: "MPCDS_low","MPCDS_high"
design <- model.matrix(~ 0 + clin$group)
colnames(design) <- levels(clin$group)

# 3) Optional: add covariates if you have them (age, sex, batch)
# design <- model.matrix(~ 0 + clin$group + clin$age + clin$sex + clin$batch)
# colnames(design)[1:2] <- levels(clin$group)

# 4) Log-transform if not already log-scale
#logexpr <- log2(expr + 1)

# 5) Filter lowly expressed genes (keeps genes with reasonable signal)
keep <- rowMeans(expr) > 1  # tweak threshold as needed
expr <- expr[keep, ]

# 6) Fit limma with empirical Bayes *trend* (key for non-count log-intensity data)
fit <- lmFit(expr, design)

# 7) Set up contrasts (example: High vs Low)
contrast.matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# 8) eBayes with trend=TRUE for TPM/FPKM
fit2 <- eBayes(fit2, trend = TRUE)

# 9) Results
res <- topTable(fit2, coef = "High_vs_Low", number = Inf, adjust.method = "BH")
# Filter DEGs
deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)

# 10) Vectors of up/down genes if needed
up_genes   <- rownames(subset(deg, logFC > 0))
down_genes <- rownames(subset(deg, logFC < 0))

# 11) Save for enrichment
deg_genes <- rownames(deg)

library(clusterProfiler)
library(org.Hs.eg.db) # Human gene annotation
library(enrichplot)

# Assuming DEGs are in a vector called 'deg_genes' with Entrez IDs
ego_bp <- enrichGO(gene = deg_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

df_bp <- data.frame(
  category = 'BP',
  ID = ego_bp@result$ID,
  term = ego_bp@result$Description,
  count = ego_bp@result$Count,
  gene = gsub("/", ", ", ego_bp@result$geneID),
  logFC = log2(ego_bp@result$FoldEnrichment),
  adj_pval = ego_bp@result$p.adjust,
  zscore = ego_bp@result$zScore,
  stringsAsFactors = FALSE  # optional but avoids factor conversion
)

ego_cc <- enrichGO(gene = deg_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

df_cc <- data.frame(
  category = 'CC',
  ID = ego_cc@result$ID,
  term = ego_cc@result$Description,
  count = ego_cc@result$Count,
  gene = gsub("/", ", ", ego_cc@result$geneID),
  logFC = log2(ego_cc@result$FoldEnrichment),
  adj_pval = ego_cc@result$p.adjust,
  zscore = ego_cc@result$zScore,
  stringsAsFactors = FALSE  # optional but avoids factor conversion
)

ego_mf <- enrichGO(gene = deg_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

df_mf <- data.frame(
  category = 'MF',
  ID = ego_mf@result$ID,
  term = ego_mf@result$Description,
  count = ego_mf@result$Count,
  gene = gsub("/", ", ", ego_mf@result$geneID),
  logFC = log2(ego_mf@result$FoldEnrichment),
  adj_pval = ego_mf@result$p.adjust,
  zscore = ego_mf@result$zScore,
  stringsAsFactors = FALSE  # optional but avoids factor conversion
)

df_all <- rbind(df_bp, df_cc, df_mf)
df_all <- df_all[df_all$logFC > 1, ]

entrez_ids <- bitr(deg_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID), ]

ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                    organism = 'hsa', 
                    keyType = "kegg",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

df_ekegg <- data.frame(
  category = 'KEGG',
  ID = ekegg@result$ID,
  term = ekegg@result$Description,
  count = ekegg@result$Count,
  gene = gsub("/", ", ", ekegg@result$geneID),
  logFC = log2(ekegg@result$FoldEnrichment),
  adj_pval = ekegg@result$p.adjust,
  zscore = ekegg@result$zScore,
  stringsAsFactors = FALSE  # optional but avoids factor conversion
)


library(GOplot)
library(gridExtra)
library(grid)

# Create all titled plots
bp_titled <- arrangeGrob(GOCircle(df_bp,label.size = 3 ), 
                         top = textGrob("Biological Process", gp = gpar(fontsize = 18, fontface = "bold")))

cc_titled <- arrangeGrob(GOCircle(df_cc,label.size = 3 ), 
                         top = textGrob("Cellular Component", gp = gpar(fontsize = 18, fontface = "bold")))

mf_titled <- arrangeGrob(GOCircle(df_mf,label.size = 3 ), 
                         top = textGrob("Molecular Function", gp = gpar(fontsize = 18, fontface = "bold")))

kegg_titled <- arrangeGrob(GOCircle(df_ekegg,label.size = 3 ), 
                           top = textGrob("KEGG Pathways", gp = gpar(fontsize = 18, fontface = "bold")))

# Arrange all four in a 2x2 grid
grid.arrange(bp_titled, cc_titled, ncol = 2)
grid.arrange(mf_titled, kegg_titled, ncol = 2)

library(maftools)

# -----------------------------
# 1. Extract clinical sample IDs
# -----------------------------
# Assume rownames are case-level barcodes (first 12 chars)
high_samples <- rownames(cox_results$clinical_data)[
  cox_results$clinical_data$PCDRS_group == "High"
]
low_samples <- rownames(cox_results$clinical_data)[
  cox_results$clinical_data$PCDRS_group == "Low"
]

# Trim to first 12 characters just in case
high_samples <- substr(high_samples, 1, 12)
low_samples  <- substr(low_samples, 1, 12)

# -----------------------------
# 2. Read MAF file
# -----------------------------
maf <- read.maf(maf = skcm.tcga$maf)

# -----------------------------
# 3. Harmonize MAF barcodes to 12 chars
# -----------------------------
maf@data$Tumor_Sample_Barcode <- substr(maf@data$Tumor_Sample_Barcode, 1, 12)
maf@clinical.data$Tumor_Sample_Barcode <- substr(maf@clinical.data$Tumor_Sample_Barcode, 1, 12)

# -----------------------------
# 4. Check overlap before subsetting
# -----------------------------
cat("High group overlap:", length(intersect(high_samples, maf@data$Tumor_Sample_Barcode)), "\n")
cat("Low group overlap: ", length(intersect(low_samples, maf@data$Tumor_Sample_Barcode)), "\n")

# -----------------------------
# 5. Subset MAF by sample groups
# -----------------------------
high_maf <- subsetMaf(maf = maf, tsb = high_samples)
low_maf  <- subsetMaf(maf = maf, tsb = low_samples)

# -----------------------------
# 6. Summary plots (optional)
# -----------------------------
plotmafSummary(high_maf, dashboard = TRUE, addStat = "median")
plotmafSummary(low_maf, dashboard = TRUE, addStat = "median")


#-------------------------------------------------------------------------------


library(IOBR)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(tidyr)
library(ggpubr)
library(tibble)
mpcds_groups <- cox_results$clinical_data$PCDRS_group
mpcds_df <- data.frame(
  Sample = rownames(cox_results$clinical_data),
  MPCDS_Group  = mpcds_groups
)


expr_data <- skcm.tcga$mRNA.expr

# Clean gene names and remove duplicates FIRST
expr_data <- expr_data[!is.na(rownames(expr_data)), ]
expr_data <- expr_data[rownames(expr_data) != "", ]

# Aggregate duplicates by taking mean
if(any(duplicated(rownames(expr_data)))) {
  print("Found duplicate genes, aggregating...")
  expr_data <- aggregate(. ~ rownames(expr_data), 
                             data = as.data.frame(expr_data), 
                             FUN = mean)
  rownames(expr_data) <- expr_data[,1]
  expr_data <- expr_data[,-1]
}

# Filter low-expressed genes (keep genes expressed in at least 10% samples)
keep_genes <- rowSums(expr_data > 0) >= ncol(expr_data) * 0.1
expr_data_filtered <- expr_data[keep_genes, ]
print(paste("Kept", sum(keep_genes), "genes out of", nrow(expr_data)))

if(max(expr_data_filtered) > 1000) {
  print("Data appears to be counts, converting to TPM...")
  expr_data_tpm <- count2tpm(
    expr_data_filtered,
    idType = "symbol",
    org = "hsa",
    source = "local"
  )
} else {
  print("Data appears already normalized, using as-is...")
  expr_data_tpm <- expr_data_filtered
}

expr_data_clean <- expr_data_tpm

common_samples <- intersect(colnames(expr_data_clean), mpcds_df$ID)
expr_data_final <- expr_data_clean[, common_samples]
mpcds_df_final <- mpcds_df[mpcds_df$ID %in% common_samples, ]


print("Running deconvolution methods...")

# CIBERSORT - most widely validated
cibersort <- deconvo_tme(
  eset = expr_data_final,  # Use non-log transformed data
  method = "cibersort",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE,  # Let it handle normalization
  perm = 1000,  # Increased for better p-values,
)

# EPIC - good for tumor samples
epic <- deconvo_tme(
  eset = expr_data_final,
  method = "epic",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE
)

# MCPcounter - enrichment scores
mcp <- deconvo_tme(
  eset = expr_data_final,
  method = "mcpcounter",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE
)

# xCell - comprehensive cell types
xcell <- deconvo_tme(
  eset = expr_data_final,  # FIXED: was 'expr'
  method = "xcell",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE
)

# ESTIMATE - overall immune/stromal scores
estimate <- deconvo_tme(
  eset = expr_data_final,
  method = "estimate",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE
)

# TIMER - needs correct tissue type
timer <- deconvo_tme(
  eset = expr_data_final,
  method = "timer",
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE,
  group_list = rep("skcm", ncol(expr_data_final))  # FIXED: was 'stad'
)

# quanTIseq - absolute cell fractions
quantiseq <- deconvo_tme(
  eset = expr_data_final,
  tumor = TRUE,
  arrays = FALSE,
  scale_mrna = TRUE,
  method = "quantiseq"
)


tme_combine<-cibersort %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,epic,by      = "ID") %>% 
  inner_join(.,estimate,by  = "ID") %>% 
  inner_join(.,timer,by     = "ID") %>% 
  inner_join(.,quantiseq,by = "ID")


mpcds_df <- mpcds_df %>%
  dplyr::rename(ID = Sample)

tme_combine_filtered <- tme_combine[tme_combine$ID %in% mpcds_df$ID, ]
immune_clinical <- merge(mpcds_df, tme_combine_filtered, by = "ID")

print("Performing immune cell enrichment analysis between MPCDS groups...")

# Use iobr_cor_plot function for comprehensive analysis


res <- iobr_cor_plot(
  pdata_group           = mpcds_df,
  id1                   = "ID",
  feature_data          = tme_combine_filtered,
  id2                   = "ID",
  target                = NULL,
  group                 = "MPCDS_Group",
  is_target_continuous  = FALSE,
  padj_cutoff           = 1,
  index                 = 1,
  category              = "signature",
  signature_group       = sig_group,
  ProjectID             = "SKCM_MPCDS",
  palette_box           = "paired1",
  palette_corplot       = "pheatmap",
  palette_heatmap       = 2,
  feature_limit         = 50,
  character_limit       = 30,
  show_heatmap_col_name = TRUE,
  show_col              = TRUE,
  show_plot             = TRUE
)

res <- iobr_cor_plot(
  pdata_group           = mpcds_df,
  id1                   = "ID",
  feature_data          = tme_combine_filtered,
  id2                   = "ID",
  target                = NULL,
  group                 = "MPCDS_Group",
  is_target_continuous  = FALSE,
  padj_cutoff           = 1,
  index                 = 1,
  category              = "signature",
  signature_group       = sig_group[c(2,5,6,7,29)],
  ProjectID             = "SKCM_MPCDS",
  palette_box           = "paired1",
  palette_corplot       = "pheatmap",
  palette_heatmap       = 2,
  feature_limit         = 50,
  character_limit       = 30,
  show_heatmap_col_name = TRUE,
  show_col              = TRUE,
  show_plot             = TRUE
)

print("Calculating TMB from mutation data...")

if(exists("skcm.tcga") && "maf" %in% names(skcm.tcga)) {
  
  # Calculate TMB (mutations per sample)
  mutation_counts <- skcm.tcga$maf %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mutation_count = n(), .groups = "drop")
  
  # Convert to mutations per megabase (assuming ~30 Mb exome)
  mutation_counts$TMB_score <- mutation_counts$mutation_count / 30
  
  # Match sample IDs to expression data format
  mutation_counts$ID <- substr(mutation_counts$Tumor_Sample_Barcode, 1, 12)
  #mutation_counts$ID <- gsub("-", ".", mutation_counts$ID)
  
  # Create TMB dataframe
  tmb_data <- mutation_counts %>%
    dplyr::select(ID, TMB_score)
  
  names(mpcds_df)[1] <- "ID"
  
  mpcds_df <- mpcds_df %>%
    left_join(tmb_data, by = "ID")
  
  print(paste("TMB calculated for", nrow(tmb_data), "samples"))
  
} else {
  print("Mutation data (maf) not found in skcm.tcga")
}


# 2. TUMOR NEOANTIGEN BURDEN (TNB) CALCULATION
print("Calculating TNB proxy from HLA and antigen presentation signatures...")

# Calculate HLA expression and antigen presentation scores
hla_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1")
antigen_genes <- c("TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "CIITA")

# Calculate HLA and antigen processing scores
tnb_related_genes <- c(hla_genes, antigen_genes)

tnb_signature <- IOBR::calculate_sig_score(
  pdata = NULL,
  eset = expr_data,
  signature = list(TNB_signature = tnb_related_genes),
  method = "ssgsea",
  mini_gene_count = 3
)

mpcds_df <- mpcds_df %>%
  left_join(tnb_signature, by = "ID")

# 3. MICROSATELLITE INSTABILITY (MSI) SCORE
if(exists("signature_collection")) {
  # First, let's see what signatures are available
  
  # Let's also check their lengths to see which ones have genes
  sig_lengths <- sapply(signature_collection, length)
  non_empty_sigs <- sig_lengths[sig_lengths > 0]
  
  # Look for any signatures that might relate to MSI/DNA repair
  potential_msi <- grep("repair|damage|mismatch|homolog|DNA", names(signature_collection), ignore.case = TRUE, value = TRUE)
  
  if(length(potential_msi) > 0) {
    print(paste("Found potential MSI-related signatures:", paste(potential_msi, collapse = ", ")))
    
    # Check if they have genes
    for(sig_name in potential_msi) {
      sig_genes <- signature_collection[[sig_name]]
      overlap <- intersect(sig_genes, rownames(expr_data))
      cat(sig_name, " | Genes in sig:", length(sig_genes),
          " | Overlap:", length(overlap), "\n")
      if(length(signature_collection[[sig_name]]) > 0) {
        
        msi_signature <- IOBR::calculate_sig_score(
          pdata = mpcds_df,
          eset = expr_data,
          signature = signature_collection[sig_name],
          method = "ssgsea",
          mini_gene_count = 2
        )
        
        msi_data <- data.frame(
          ID = colnames(msi_signature),
          MSI_score = as.numeric(msi_signature[1, ]),
          stringsAsFactors = FALSE
        )
        
        mpcds_df <- mpcds_df %>%
          left_join(msi_data, by = "ID")
        
        print("MSI score calculated")
    
      }
    }
  } else {
    print("No MSI-related signatures found - skipping MSI calculation")
  }
} else {
  print("signature_collection not found - skipping MSI calculation")
}

data("sig_group", package = "IOBR")

if (exists("sig_group")) {
  
  # 1. Calculate Exhaustion Score
  exhaustion_res <- IOBR::calculate_sig_score(
    pdata = NULL,
    eset = expr_data,
    signature = IOBR::signature_collection[sig_group$immu_exhaustion],
    method = "ssgsea"
  )
  # Note: calculate_sig_score returns a dataframe with ID and the scores
  
  # 2. Calculate Exclusion Score
  exclusion_res <- IOBR::calculate_sig_score(
    pdata = NULL,
    eset = expr_data,
    signature = IOBR::signature_collection[sig_group$immu_exclusion],
    method = "ssgsea"
  )
  
  # 3. Merge and Normalize
  proxy_df <- exhaustion_res %>%
    inner_join(exclusion_res, by = "ID") %>%
    tibble::column_to_rownames("ID") %>% # Move ID to rownames for scaling
    scale() %>% # Z-score standardization (Crucial if combining different signatures)
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>% 
    rowwise() %>%
    # Taking the MEAN is usually standard for composite scores, but Median is okay too
    mutate(Combined_Evasion_Score = mean(c_across(!ID), na.rm = TRUE)) %>%
    ungroup()
  
  # 4. Extract Final Data
  proxy_data <- proxy_df %>%
    select(ID, TIDE_Proxy = Combined_Evasion_Score)
  
  # 5. Merge with Metadata
  mpcds_df <- mpcds_df %>%
    left_join(proxy_data, by = "ID")
  
  print("Proxy Evasion score calculated. High Score = High Evasion/Dysfunction.")
}

print("Comparing immunotherapy biomarkers between MPCDS groups...")

# Create comparison dataframe
biomarker_cols <- c("TMB_score", "TNB_signature", "MSI_score", "TIDE_score")
available_cols <- intersect(biomarker_cols, colnames(mpcds_df))

# Statistical comparisons
biomarker_comparisons <- list()

# Get available biomarker columns
actual_biomarker_cols <- intersect(c("TMB_score", "TNB_signature", "MSI_score", "TIDE_score"), 
                                   colnames(mpcds_df))

if(length(actual_biomarker_cols) == 0) {
  print("No biomarker scores found in data")
} else {
  print(paste("Found biomarkers:", paste(actual_biomarker_cols, collapse = ", ")))
  
  for(biomarker in actual_biomarker_cols) {
    # Remove NA values for comparison
    low_values <- mpcds_df[[biomarker]][mpcds_df$MPCDS_Group == "Low"]
    high_values <- mpcds_df[[biomarker]][mpcds_df$MPCDS_Group == "High"]
    
    low_values <- low_values[!is.na(low_values)]
    high_values <- high_values[!is.na(high_values)]
    
    if(length(low_values) > 0 && length(high_values) > 0) {
      # Wilcoxon test
      test_result <- wilcox.test(low_values, high_values, alternative = "two.sided")
      
      biomarker_comparisons[[biomarker]] <- list(
        biomarker = biomarker,
        p_value = test_result$p.value,
        low_median = median(low_values, na.rm = TRUE),
        high_median = median(high_values, na.rm = TRUE),
        low_n = length(low_values),
        high_n = length(high_values)
      )
    }
  }
}

# Print results
if(length(biomarker_comparisons) > 0) {
  cat("\nImmunotherapy Biomarker Comparison Results:\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  for(comp in biomarker_comparisons) {
    cat(sprintf("%s:\n", comp$biomarker))
    cat(sprintf("  Low MPCDS median: %.3f (n=%d)\n", comp$low_median, comp$low_n))
    cat(sprintf("  High MPCDS median: %.3f (n=%d)\n", comp$high_median, comp$high_n))
    cat(sprintf("  P-value: %.3e\n", comp$p_value))
    cat(sprintf("  Significant: %s\n\n", ifelse(comp$p_value < 0.05, "Yes", "No")))
  }
} else {
  cat("No biomarker comparisons could be performed\n")
}

# ===============================================================================
# VISUALIZATION - FIGURE 11A EQUIVALENT
# ===============================================================================

print("Creating comprehensive immunotherapy biomarker visualization...")

# Prepare data for plotting - only use available biomarkers
if(length(actual_biomarker_cols) > 0) {
  plot_data <- mpcds_df %>%
    select(ID, MPCDS_Group, all_of(actual_biomarker_cols)) %>%
    pivot_longer(cols = all_of(actual_biomarker_cols), 
                 names_to = "Biomarker", 
                 values_to = "Score") %>%
    filter(!is.na(Score))
  
  # Create the multi-panel plot
  p1 <- ggplot(plot_data, aes(x = MPCDS_Group, y = Score, fill = MPCDS_Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
    facet_wrap(~ Biomarker, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Low" = "#E31A1C", "High" = "#1F78B4")) +
    labs(
      title = "Immunotherapy Response Biomarkers by MPCDS Group",
      subtitle = "Analysis using actual TCGA mutation data and IOBR signatures",
      x = "MPCDS Group",
      y = "Biomarker Score",
      fill = "MPCDS Group"
    ) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Add statistical significance annotations
  p1 <- p1 + stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = list(c("Low", "High"))
  )
  
  print(p1)
} else {
  print("No biomarker data available for plotting")
}


library(oncoPredict)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tibble)


# Load the necessary libraries
library(depmap)
library(ExperimentHub) # This is used by depmap

# --- 1. Load the Expression Data ---
# Replace with your actual file path
expr_file_path <- "Expression_Public_25Q3_subsetted.csv"

expression_data_df <- read.csv(
  expr_file_path, 
  header = TRUE,         # The first row is a header
  row.names = 1,       # The first column is row names (genes)
  check.names = FALSE  # IMPORTANT: Keeps cell line names as-is
)

# --- 2. Load the Drug Sensitivity Data ---
# Replace with your actual file path
drug_file_path <- "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv"

drug_data_df <- read.csv(
  drug_file_path, 
  header = TRUE,         # The first row is a header
  row.names = 1,       # The first column is row names (drugs)
  check.names = FALSE  # IMPORTANT: Keeps cell line names as-is
)

# --- 3. Convert to Matrix ---
# The calcPhenotype function needs matrices
trainExpr <- as.matrix(expression_data_df)
drugData <- as.matrix(drug_data_df)

# --- 4. Check the results ---
cat("Expression Matrix (trainExpr):\n")
print(dim(trainExpr))
print(class(trainExpr))

cat("\nDrug Data Matrix (drugData):\n")
print(dim(drugData))
print(class(drugData))

# 1. Find the cell lines (colnames) that are in BOTH datasets
common_cell_lines <- intersect(rownames(trainExpr), rownames(drugData))

# 2. Filter both matrices to keep *only* these common cell lines
trainExpr_final <- t(trainExpr[common_cell_lines,7:ncol(trainExpr)])
drugData_final <- t(drugData[common_cell_lines,7:ncol(drugData)])

calcPhenotype_FIXED <- function (trainingExprData, trainingPtype, testExprData, batchCorrect = "eb", 
                                 powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2, 
                                 minNumSamples = 10, selection = -1, printOutput = TRUE, 
                                 removeLowVaringGenesFrom = "homogenizeData") 
{
  # --- FIX 1: Replaced buggy 'class()' check with 'is.matrix()' ---
  if (!is.matrix(testExprData)) 
    stop("ERROR: \"testExprData\" must be a matrix.")
  
  # --- FIX 2: Replaced buggy 'class()' check with 'is.matrix()' ---
  if (!is.matrix(trainingExprData)) 
    stop("ERROR: \"trainingExprData\" must be a matrix.")
  
  # --- FIX 3: Replaced buggy 'class()' check with 'is.numeric()' ---
  if (!is.numeric(trainingPtype)) 
    stop("ERROR: \"trainingPtype\" must be a numeric vector.")
  
  # --- ALL OTHER CODE IS THE SAME ---
  
  if (ncol(trainingExprData) != length(trainingPtype)) 
    stop("The training phenotype must be of the same length as the number of columns of the training expressin matrix.")
  if ((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < 
                                                   minNumSamples)) {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  homData <- homogenizeData(testExprData, trainingExprData, 
                            batchCorrect = batchCorrect, selection = selection, 
                            printOutput = printOutput)
  if (!(removeLowVaringGenesFrom %in% c("homogenizeData", 
                                        "rawData"))) {
    stop("\"removeLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }
  keepRows <- seq(1:nrow(homData$train))
  if (removeLowVaryingGenes > 0 && removeLowVaryingGenes < 
      1) {
    if (removeLowVaringGenesFrom == "homogenizeData") {
      keepRows <- doVariableSelection(cbind(homData$test, 
                                            homData$train), removeLowVaryingGenes = removeLowVaryingGenes)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) 
        cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."))
    }
    else if (removeLowVaringGenesFrom == "rawData") {
      evaluabeGenes <- rownames(homData$test)
      keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes, 
      ], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes, 
      ], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) 
        cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."))
    }
  }
  offset = 0
  if (powerTransformPhenotype) {
    if (min(trainingPtype) < 0) {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  if (printOutput) 
    cat("\nFitting Ridge Regression model... ")
  trainFrame <- data.frame(Resp = trainingPtype, t(homData$train[keepRows, 
  ]))
  rrModel <- ridge::linearRidge(Resp ~ ., data = trainFrame)
  if (printOutput) 
    cat("Done\n\nCalculating predicted phenotype...")
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # --- FIX 4: Proactively fixed another buggy 'class()' check ---
  if (is.numeric(homData$test) && !is.matrix(homData$test)) {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol = 1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata = rbind(testFrame, 
                                              testFrame))[1]
  }
  else {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata = testFrame)
  }
  if (powerTransformPhenotype) {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if (printOutput) 
    cat("Done\n\n")
  return(preds)
}

# --- 1. Create the full phenotype vector (as you did) ---
phenotype_vec_full <- as.numeric(drugData_final[1,])

# --- 2. Find which samples are *not* NA ---
# This creates a TRUE/FALSE vector
is_valid <- !is.na(phenotype_vec_full)

# --- 3. Create the FINAL, CLEAN training data ---
# Filter the phenotype vector to remove NAs
trainPtype_clean <- phenotype_vec_full[is_valid]

# Filter the expression matrix to keep *only* the matching columns
trainExpr_clean <- trainExpr_final[, is_valid]

# --- 4. Ensure your test data is numeric (as before) ---
storage.mode(gset_filtered) <- "numeric"

# --- 5. Run the function with the CLEAN data ---
predicted_phenotypes <- calcPhenotype_FIXED(
  trainingExprData = trainExpr_clean,   # <-- Use the clean matrix
  trainingPtype = trainPtype_clean,      # <-- Use the clean vector
  testExprData = gset_filtered,
  batchCorrect = "eb",
  powerTransformPhenotype = FALSE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  selection = 1,
  printOutput = TRUE,
  removeLowVaringGenesFrom = "homogenizeData"
)

# ===============================================================================
# MELANOMA-SPECIFIC DRUG SELECTION
# ===============================================================================

# Get all available drugs from your data
available_drugs <- rownames(drugData_final)
print(paste("Total drugs available:", length(available_drugs)))
print("First 20 drugs:")
print(head(available_drugs, 20))

# Define melanoma-relevant drugs by category
# These are FDA-approved or commonly used drugs for melanoma

melanoma_drugs <- list(
  
  # BRAF inhibitors (for BRAF-mutant melanoma)
  BRAF_inhibitors = c(
    "Dabrafenib", "Vemurafenib", "Encorafenib"
  ),
  
  # MEK inhibitors (often combined with BRAF inhibitors)
  MEK_inhibitors = c(
    "Trametinib", "Cobimetinib", "Binimetinib"
  ),
  
  # Immunotherapy checkpoint inhibitors (if available as proxy)
  Checkpoint_related = c(
    "Pictilisib"  # PI3K inhibitor that may relate to immune response
  ),
  
  # Chemotherapy agents
  Chemotherapy = c(
    "Cisplatin", "Carboplatin", "Paclitaxel", "Docetaxel", 
    "Dacarbazine", "Temozolomide", "5-Fluorouracil"
  ),
  
  # Targeted therapy - KIT inhibitors (for KIT-mutant melanoma)
  KIT_inhibitors = c(
    "Imatinib", "Nilotinib", "Dasatinib"
  ),
  
  # mTOR/PI3K/AKT pathway inhibitors
  PI3K_AKT_mTOR = c(
    "Everolimus", "Temsirolimus", "Alpelisib", "AZD8055",
    "BEZ235", "GDC-0941", "MK-2206"
  ),
  
  # CDK inhibitors (cell cycle)
  CDK_inhibitors = c(
    "Palbociclib", "Ribociclib", "Abemaciclib"
  ),
  
  # Other targeted therapies
  Other_targeted = c(
    "Sorafenib", "Sunitinib", "Pazopanib",  # Multi-kinase inhibitors
    "Selumetinib"  # MEK inhibitor
  )
)

# Find which drugs are actually available in your dataset
available_melanoma_drugs <- list()
for(category in names(melanoma_drugs)) {
  found_drugs <- intersect(melanoma_drugs[[category]], available_drugs)
  if(length(found_drugs) > 0) {
    available_melanoma_drugs[[category]] <- found_drugs
    cat(sprintf("\n%s (%d drugs found):\n", category, length(found_drugs)))
    cat(paste("  -", found_drugs, collapse = "\n"))
  }
}


#--------------LMS RMS RESPONDERS-------------#
# Create survival object
surv_obj <- Surv(time = clin$futime/30.44, event = clin$fustat)

fit_km <- survfit(surv_obj ~ group, data = clin)

# Plot with time points at 12 and 24 months
plot_B <- ggsurvplot(
  fit_km,
  data = clin,
  pval = TRUE,                    # Show p-value
  pval.method = TRUE,             # Show test method
  conf.int = FALSE,               # Don't show confidence intervals
  risk.table = TRUE,              # Show risk table
  risk.table.height = 0.25,       # Risk table height
  xlab = "Time(months)",
  ylab = "Overall survival",
  xlim = c(0, 30),                # Show up to 24 months
  break.time.by = 6,              # Tick marks every 6 months
  palette = c("red", "blue"),     # High = red, Low = blue
  legend.title = "",
  legend.labs = c("High MPCDS", "Low MPCDS"),
  title = "Comparison of RMS",
  ggtheme = theme_classic()
)

# Add vertical lines at 12 and 24 months
plot_B$plot <- plot_B$plot + 
  geom_vline(xintercept = 12, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 24, linetype = "dashed", color = "gray40") +
  annotate("text", x = 12, y = 0.05, label = "Comparison of RMS at 12-month", 
           size = 3, hjust = 0.5) +
  annotate("text", x = 24, y = 0.05, label = "Comparison of RMS at 24-month", 
           size = 3, hjust = 0.5)

# Print and save
plot_B
ggsave("Plot_B_RMS_12_24_months.pdf", plot = plot_B$plot, width = 8, height = 6)

plot_C <- ggsurvplot(
  fit_km,
  data = clin,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  xlab = "Time(months)",
  ylab = "Overall survival",
  xlim = c(0, 24),
  break.time.by = 3,
  palette = c("red", "blue"),
  legend.title = "",
  legend.labs = c("High MPCDS", "Low MPCDS"),
  title = "Comparison of LTS",
  ggtheme = theme_classic()
)

# Add vertical lines at 6 and 18 months
plot_C$plot <- plot_C$plot + 
  geom_vline(xintercept = 6, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 18, linetype = "dashed", color = "gray40") +
  annotate("text", x = 6, y = 0.05, label = "Comparison of LTS after 6-month", 
           size = 3, hjust = 0.5) +
  annotate("text", x = 18, y = 0.05, label = "Comparison of LTS after 18-month", 
           size = 3, hjust = 0.5)

print(plot_C)
ggsave("Plot_C_LTS_6_18_months.pdf", plot = plot_C$plot, width = 8, height = 6)

calculate_response_from_survival <- function(data, time_col, event_col, 
                                             cutoff_time = 12) {
  # Initialize Response column
  data$Response <- NA_character_
  
  for (i in seq_len(nrow(data))) {
    time <- data[[time_col]][i]
    event <- data[[event_col]][i]
    
    # Skip if missing values
    if (is.na(time) || is.na(event)) {
      data$Response[i] <- NA
      next
    }
    
    # Convert time to months (if in days)
    time_months <- time / 30.44
    
    if (time_months >= cutoff_time && event == 0) {
      # Alive at cutoff
      if (time_months >= 24) {
        data$Response[i] <- "CR"  # Complete Response
      } else {
        data$Response[i] <- "PR"  # Partial Response
      }
    } else if (time_months >= cutoff_time && event == 1) {
      # Event after cutoff
      data$Response[i] <- "SD"    # Stable Disease
    } else if (time_months < cutoff_time) {
      # Early event or short survival
      data$Response[i] <- "PD"    # Progressive Disease
    }
  }
  
  return(data)
}

# Apply the function
# Adjust 'OS.time' and 'OS' to match your actual column names
resp_data <- calculate_response_from_survival(
  clin, 
  time_col = "futime",  # Change this to your survival time column name
  event_col = "fustat",      # Change this to your event status column name
  cutoff_time = 12       # 12 months is common cutoff for immunotherapy
)

# ===============================================================================
# PLOT D: MPCDS scores by response level
# ===============================================================================

# You need a response column: CR, PD, PR, SD
# If you don't have this, you'll need to add it to resp_data
# resp_data$Response <- your_response_column

# Check if Response column exists
if("Response" %in% colnames(resp_data)) {
  
  # Make sure Response is a factor with correct order
  resp_data$Response <- factor(resp_data$Response, 
                              levels = c("CR", "PD", "PR", "SD"))
  
  # Create violin plot
  plot_D <- ggplot(resp_data, aes(x = Response, y = PCDRS, fill = Response)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    scale_fill_manual(values = c("CR" = "#D8494C", "PD" = "#3F76A9", 
                                 "PR" = "#E8B959", "SD" = "#5FAF80")) +
    labs(x = "Response", y = "MPCDS", title = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    stat_compare_means(comparisons = list(c("CR", "PD"), c("CR", "PR"), 
                                          c("CR", "SD"), c("PD", "PR")),
                       method = "wilcox.test", label = "p.format")
  
  print(plot_D)
  ggsave("Plot_D_MPCDS_by_response.pdf", plot = plot_D, width = 6, height = 6)
  
} else {
  cat("Response column not found in mpcds_df. Please add it to create Plot D.\n")
}

if("Response" %in% colnames(resp_data)) {
  
  resp_data$Binary_Response <- ifelse(
    resp_data$Response %in% c("CR", "PR"), "R", "NR"
  )
  
  # Make it a factor
  resp_data$Binary_Response <- factor(resp_data$Binary_Response, 
                                     levels = c("NR", "R"))
  
  # Create violin plot
  plot_G <- ggplot(resp_data, aes(x = Binary_Response, y = PCDRS, 
                                 fill = Binary_Response)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    scale_fill_manual(values = c("NR" = "#D8494C", "R" = "#5FA4D1")) +
    labs(x = "response", y = "MPCDS", title = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    stat_compare_means(method = "wilcox.test", label = "p.format")
  
  print(plot_G)
  ggsave("Plot_G_MPCDS_responder_vs_nonresponder.pdf", plot = plot_G, 
         width = 5, height = 6)
  
} else {
  cat("Response column not found. Please add it to create Plot G.\n")
}

library("IMvigor210CoreBiologies")
data(cds)

# Explore the dataset
print(dim(cds))  # Check dimensions
print(head(pData(cds)))  # Clinical data
print(head(fData(cds)))  # Gene features

# Key clinical variables in IMvigor210:
# - binaryResponse: CR/PR vs SD/PD
# - censOS: Overall survival censoring
# - os: Overall survival time
# - FMOne mutation burden per MB: TMB score
# - Immune phenotype: Desert/Excluded/Inflamed

# Extract relevant data
imvigor_clinical <- pData(cds)
imvigor_expr <- counts(cds)  # Expression data

# Key response variable
table(imvigor_clinical$binaryResponse)

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

# Assuming you successfully loaded the data (using Option 1 or 3)

# Step 1: Extract clinical and expression data
if(exists("cds")) {
  # From IMvigor210CoreBiologies package
  imvigor_clinical <- pData(cds)
  imvigor_expr <- counts(cds)
  
  # Key columns to extract
  imvigor_analysis <- imvigor_clinical %>%
    select(
      Response = binaryResponse,  # CR/PR vs SD/PD
      OS_time = os,               # Overall survival
      OS_event = censOS,          # Death event
      TMB = 'FMOne mutation burden per MB',  # Tumor mutation burden
      Immune_phenotype = `Immune phenotype`,  # Desert/Excluded/Inflamed
      IC_Level = 'IC Level'         # PD-L1 expression level
    ) %>%
    mutate(
      # Binary response for easy analysis
      Responder = ifelse(Response == "CR/PR", "Responder", "Non-Responder"),
      # Numeric response
      Response_Binary = ifelse(Response == "CR/PR", 1, 0)
    )
  
  print("IMvigor210 clinical data prepared")
  print(table(imvigor_analysis$Responder))
}

# Convert Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(imvigor_expr),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Replace rownames
rownames(imvigor_expr) <- gene_symbols

# Remove NAs (genes that couldn't be mapped)
imvigor_expr <- imvigor_expr[!is.na(rownames(imvigor_expr)), ]

# Remove duplicates (take mean if duplicated)
imvigor_expr <- imvigor_expr[!duplicated(rownames(imvigor_expr)), ]


available_genes <- intersect(pcdrg_genes, rownames(imvigor_expr))
print(paste("Available MPCDS genes:", length(available_genes), "out of", length(pcdrg_genes)))

# Extract expression for MPCDS genes
mpcds_expr_imvigor <- imvigor_expr[available_genes, ]

# METHOD 1: Using IOBR (RECOMMENDED - handles everything)
mpcds_scores_imvigor <- IOBR::calculate_sig_score(
  pdata = NULL,
  eset = mpcds_expr_imvigor,
  signature = list(MPCDS_signature = available_genes),
  method = "ssgsea",
  mini_gene_count = 3
)
