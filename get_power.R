########### Get power ##################

get_power <- function(n1, n2 , icc, eff.size, b, BF_thres, ndatasets = 5000, pair, batch_size) {
    # Libraries
    library(lme4)
    library(dplyr)
    
    # Source
    source("data_generation.R")
    source("small_functions.R")
    source("aafbf.R")
    
    total_var <- 1
    var_u <- icc * total_var
    var_e <- total_var - var_u
    eff_size0 <- 0
    
    if (pair == 1) {
        df_H0 <- matrix(NA, nrow = ndatasets, ncol = 4)
        df_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
        
        # Generate data
        data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u, var_e,
                                              mean_interv = eff_size0, batch_size = batch_size))
        data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u, var_e,
                                              mean_interv = eff.size, batch_size = batch_size))
        
        # Compute Bayes factor
        eff_n_H0 <- ((n1 * n2) / (1 + (n1 - 1) *  data_H0$rho_data)) / 2
        aafbf_H0 <- Map(calc_aafbf, type = "Equality", data_H0$estimates, data_H0$cov_list, list(b), eff_n_H0)
        eff_n_H1 <- ((n1 * n2) / (1 + (n1 - 1) *  data_H1$rho_data)) / 2
        aafbf_H1 <- Map(calc_aafbf, type = "Equality", data_H1$estimates, data_H1$cov_list, list(b), eff_n_H1)
        
        # Extract results
        df_H0[, 1] <- unlist(lapply(aafbf_H0, extract_res, 1)) # Bayes factor H1vsH0
        df_H0[, 2] <- unlist(lapply(aafbf_H0, extract_res, 4)) # Posterior model probabilities H1
        df_H0[, 3] <- unlist(lapply(aafbf_H0, extract_res, 2)) # Bayes factor H0vsH1
        df_H0[, 4] <- unlist(lapply(aafbf_H0, extract_res, 3)) # POsterior model probabilities H0
        colnames(df_H0) <- c("BF_10", "PMP_1", "BF_01", "PMP_0")
        
        df_H1[, 1] <- unlist(lapply(aafbf_H1, extract_res, 1)) # Bayes factor H1vsH2
        df_H1[, 2] <- unlist(lapply(aafbf_H1, extract_res, 4)) # Posterior model probabilities H1
        df_H1[, 3] <- unlist(lapply(aafbf_H1, extract_res, 2)) # Bayes factor H2vsH1
        df_H1[, 4] <- unlist(lapply(aafbf_H1, extract_res, 3)) # Posterior model probabilities H0
        colnames(df_H1) <- c("BF_10", "PMP_1", "BF_01", "PMP_0")
        
        # Bayesian power
        prop_BF01 <- length(which(df_H0[, "BF_01"] > BF_thres)) / ndatasets
        prop_BF10 <- length(which(df_H1[, "BF_10"] > BF_thres)) / ndatasets
        
        # Output
        results_matrix <- matrix(NA, nrow = b, ncol = 5)
        results_matrix[, 1] <- seq(b)
        results_matrix[, 2] <- n1 #n1
        results_matrix[, 3] <- n2 #n2
        results_matrix[, 4] <- round(prop_BF01, 3) #eta01
        results_matrix[, 5] <- round(prop_BF10, 3) #eta10
        colnames(results_matrix) <- c("b", "n1", "n2", paste("P(BF.01 >", BF_thres, "| H0)", sep = " "), 
                                      paste("P(BF.10 >", BF_thres, "| H1)", sep = " "))
        cat("Hypotheses:", "\n")
        cat("    H0:", "Intervention = Control", "\n")
        cat("    H1:", "Intervention > Control", "\n")
        cat("\n")
        
        cat("***********************************************************************", "\n")
        print(format(results_matrix, justify = "centre"))
        cat("***********************************************************************", "\n")
        cat("n1: Cluster sizes", "\n")
        cat("n2: Number of clusters", "\n")
        results <- list(data_H0 = df_H0,
                        data_H1 = df_H1)
        
    } else if (pair == 2) {
        df_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
        
        # Generate data
        data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u, var_e,
                                              mean_interv = eff.size, batch_size = 1000))
        
        # Compute Bayes factor
        eff_n_H1 <- ((n1 * n2) / (1 + (n1 - 1) *  data_H1$rho_data)) / 2
        aafbf_H1 <- Map(calc_aafbf, type = "Equality", data_H1$estimates, data_H1$cov_list, list(b), eff_n_H1)
        
        # Extract results
        df_H1[, 1] <- unlist(lapply(aafbf_H1, extract_res, 1)) # Bayes factor H1vsH1
        df_H1[, 2] <- unlist(lapply(aafbf_H1, extract_res, 4)) # Posterior model probabilities H1
        df_H1[, 3] <- unlist(lapply(aafbf_H1, extract_res, 2)) # Bayes factor H1vsH1
        df_H1[, 4] <- unlist(lapply(aafbf_H1, extract_res, 3)) # Posterior model probabilities H0
        colnames(df_H1) <- c("BF_12", "PMP_1", "BF_21", "PMP_2")
        
        # Bayesian power
        prop_BF12 <- length(which(df_H1[, "BF_12"] > BF_thres)) / ndatasets
        
        #Output
        results_matrix <- matrix(NA, nrow = b, ncol = 4)
        results_matrix[, 2] <- n1 #n1
        results_matrix[, 3] <- n2 #n2
        results_matrix[, 4] <- round(prop_BF12, 3) #eta12
        results_matrix[, 5] <- round((1 - prop_BF12), 3) #eta21
        colnames(results_matrix) <- c("n1", "n2", paste("P(BF.12 >", BF_thres, "| H1)", sep = " "), 
                                      paste("P(BF.21 >", BF_thres, "| H1)", sep = " "))
        cat("Hypotheses:", "\n")
        cat("    H1:", "Intervention > Control", "\n")
        cat("    H2:", "Intervention < Control", "\n")
        cat("\n")
        
        cat("***********************************************************************", "\n")
        print(format(results_matrix, justify = "centre"))
        cat("***********************************************************************", "\n")
        cat("n1: Cluster sizes", "\n")
        cat("n2: Number of clusters", "\n")
        
        results <- list(data_H1 = df_H1)
    }
    return(results)
}
