# Imports
library("R6")
library(splines)

get_clinical_cv <- function(snfl_value) {
  cv <- case_when(
    snfl_value < 7.2   ~ 0.052,
    snfl_value >= 7.2   & snfl_value < 13.25  ~ 0.086,
    snfl_value >= 13.25 & snfl_value < 20.65  ~ 0.051,
    snfl_value >= 20.65 & snfl_value < 54.7   ~ 0.023,
    snfl_value >= 54.7  & snfl_value < 153.6  ~ 0.015,
    snfl_value >= 153.6                       ~ 0.016
  )
  return(cv)
}


# Base Class
sNfL_Model <- R6Class("sNfL_Model",
                      public = list(
                        model_name = NULL,
                        
                        initialize = function(name = "Generic Model") {
                          self$model_name <- name
                        },

                        calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                          stop("Implement this method!")
                        },
                        
                        predict_relapse = function(new_value, patient_data = NULL, new_age = NULL) {
                          threshold <- self$calculate_threshold(patient_data, new_age)
                          is_relapse <- new_value > threshold 
                          
                          return(list(
                            is_relapse = is_relapse,
                            calculated_threshold = threshold
                          ))
                        }
                      )
)

# Child Class: Fixed Cut-off
FixedCutoff_Model <- R6Class("FixedCutoff_Model",
                             inherit = sNfL_Model, 
                             
                             public = list(
                               cutoff_value = NULL,
                               
                               initialize = function(cutoff = 12.9) {
                                 super$initialize(name = paste("Absolute Cut-off (", cutoff, " pg/mL)", sep=""))
                                 self$cutoff_value <- cutoff
                               },
                               
                               calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                                 return(self$cutoff_value)
                               }
                             )
)


# Personalized Models
# Child Class: Within-Subject Model
WithinSubject_Model <- R6Class("WithinSubject_Model",
                               inherit = sNfL_Model,
                               
                               public = list(
                                 cv_i = NULL,    
                                 cv_a = NULL,    
                                 z_value = NULL, 
                                 
                                 initialize = function(cv_i = 0.086, alpha = 0.05) {
                                   super$initialize(name = "Within-Subject")
                                   self$cv_i <- cv_i
                                   self$z_value <- qnorm(1 - alpha)
                                 },
                                 
                                 calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                                   snfl_hist <- if (!is.null(patient_data) && is.list(patient_data)) patient_data$snfl_history else patient_data
                                   if (is.null(snfl_hist) || length(snfl_hist) < 2) return(12.9)
                                   
                                   n <- length(snfl_hist)
                                   hsp <- mean(snfl_hist)
                                   
                                   current_cv_a <- get_clinical_cv(hsp) # cv_a dynamically for HSP calculated
                                   
                                   sd_i <- self$cv_i * hsp
                                   sd_a <- current_cv_a * hsp
                                   
                                   margin <- self$z_value * sqrt(((n + 1) / n) * (sd_i^2 + sd_a^2))
                                   threshold <- hsp + margin
                                   
                                   return(threshold)
                                 }
                               )
)

# Child Class : Within-Person Model
WithinPerson_Model <- R6Class("WithinPerson_Model",
                              inherit = sNfL_Model,
                              
                              public = list(
                                cv_a = NULL,    
                                alpha = NULL,   
                                
                                initialize = function(alpha = 0.05) {
                                  super$initialize(name = "Within-Person")
                                  self$alpha <- alpha
                                },
                                
                                calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                                  snfl_hist <- if (!is.null(patient_data) && is.list(patient_data)) patient_data$snfl_history else patient_data
                                  if (is.null(snfl_hist) || length(snfl_hist) < 2) return(12.9)
                                  
                                  n <- length(snfl_hist)
                                  hsp <- mean(snfl_hist) 
                                  
                                  current_cv_a <- get_clinical_cv(hsp)
                                  
                                  sd_tv <- sd(snfl_hist)
                                  sd_a <- current_cv_a * hsp
                                  sd_p <- sqrt(max(0, sd_tv^2 - sd_a^2))
                                  t_value <- qt(1 - self$alpha, df = n - 1)
                                  
                                  # prRI
                                  margin <- t_value * sqrt(((n + 1) / n) * (sd_p^2 + sd_a^2))
                                  threshold <- hsp + margin
                                  
                                  return(threshold)
                                }
                              )
)

# Child Class: Penalized Joint Quantile Model (PJQM2)
# helper function: check function
check_function <- function(w, tau) {
  # I(w <= 0) is 1, if w <= 0, otherwise 0
  indicator <- ifelse(w <= 0, 1, 0)
  return(w * (tau - indicator))
}

# helper function: penalty loss
pjqm2_loss <- function(params, residuals, beta1, beta2, tau1, tau2, lambda_u, lambda_z) {
  u_i <- params[1]
  z_i <- params[2]
  
  if (z_i <= 0) return(Inf) 
  
  # first error term
  w1 <- residuals - u_i - (z_i * beta1)
  loss1 <- sum(check_function(w1, tau1))
  
  #  second error term
  w2 <- residuals - u_i - (z_i * beta2)
  loss2 <- sum(check_function(w2, tau2))
  
  # penalty term
  penalty_u <- lambda_u * (u_i^2)
  penalty_z <- lambda_z * ((z_i - 1)^2)
  
  # final loss
  return(loss1 + loss2 + penalty_u + penalty_z)
}

# Final PJQM2 Child Class
PJQM2_Model <- R6Class("PJQM2_Model",
                       inherit = sNfL_Model,
                       
                       public = list(
                         beta_bmi = NULL,
                         beta_hba1c = NULL,
                         beta_crea = NULL,
                         tau_target = NULL,     # quantil for relapse identification
                         lambda_u = NULL,       # penalty for intercept
                         lambda_z = NULL,       # Penalty for variance
                         
                         initialize = function(beta_bmi = -0.020, 
                                               beta_hba1c = 0.090, 
                                               beta_crea = 0.467, 
                                               alpha = 0.05,
                                               lambda_u = 0.001,   
                                               lambda_z = 0.1) 
                           
                           {
                             super$initialize(name = "PJQM2")
                             self$beta_bmi <- beta_bmi
                             self$beta_hba1c <- beta_hba1c
                             self$beta_crea <- beta_crea
                             self$tau_target <- 1 - alpha
                             self$lambda_u <- lambda_u
                             self$lambda_z <- lambda_z
                           },
                         
                           calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                             
                             snfl_hist <- if (!is.null(patient_data) && is.list(patient_data)) patient_data$snfl_history else patient_data
                             if (is.null(snfl_hist) || length(snfl_hist) < 2) return(12.9)
                             
                             n <- length(snfl_hist)
                             
                             # current covariates
                             bmi_curr   <- if (!is.null(patient_data$current_bmi)) patient_data$current_bmi else 25.0
                             hba1c_curr <- if (!is.null(patient_data$current_hba1c)) patient_data$current_hba1c else 5.0
                             crea_curr  <- if (!is.null(patient_data$current_crea)) patient_data$current_crea else 0.8
                             age_curr   <- if (!is.null(new_age)) new_age else 40.0
                             
                             # historical covariates (constant)
                             age_hist   <- if (!is.null(patient_data$age_history)) patient_data$age_history else rep(age_curr, n)
                             bmi_hist   <- if (!is.null(patient_data$bmi_history)) patient_data$bmi_history else rep(bmi_curr, n)
                             hba1c_hist <- if (!is.null(patient_data$hba1c_history)) patient_data$hba1c_history else rep(hba1c_curr, n)
                             crea_hist  <- if (!is.null(patient_data$crea_history)) patient_data$crea_history else rep(crea_curr, n)
                             
                             # Covariate-adjusted baseline 
                             y_pop_hist <- private$calc_age_effect(age_hist) + 
                               self$beta_bmi * bmi_hist + 
                               self$beta_hba1c * hba1c_hist + 
                               self$beta_crea * crea_hist
                             
                             # residuals
                             log_y_hist <- log(snfl_hist)
                             residuals_hist <- log_y_hist - y_pop_hist
                             
                             # population noise
                             median_snfl <- median(snfl_hist)
                             current_cv_a <- get_clinical_cv(median_snfl)
                             sigma_eps <- sqrt(log(1 + current_cv_a^2))
                             
                             tau1 <- 0.05
                             tau2 <- 0.95
                             beta1 <- qnorm(tau1) * sigma_eps
                             beta2 <- qnorm(tau2) * sigma_eps
                             
                             # PJQM2 fitting (optimizer, start with 0 and 1)
                             initial_params <- c(0, 1) 
                             
                             # help function for optimizer
                             optim_loss <- function(params) {
                               u_i <- params[1]
                               z_i <- params[2]
                               
                               eff_lambda_u <- self$lambda_u / sqrt(n)
                               eff_lambda_z <- self$lambda_z / sqrt(n)
                               
                               if (z_i <= 0.01) return(1e9) # prevents for crashs at z_i <= 0
                               
                               w1 <- residuals_hist - u_i - (z_i * beta1)
                               w2 <- residuals_hist - u_i - (z_i * beta2)
                               
                               loss1 <- sum(w1 * (tau1 - ifelse(w1 <= 0, 1, 0)))
                               loss2 <- sum(w2 * (tau2 - ifelse(w2 <= 0, 1, 0)))
                               
                               return(loss1 + loss2 + (eff_lambda_u * u_i^2) + (eff_lambda_z * (z_i - 1)^2))
                             }
                             
                             # optimizer (with L-BFGS-B, which allows bounds)
                             fit <- optim(par = initial_params, 
                                          fn = optim_loss, 
                                          method = "L-BFGS-B", 
                                          lower = c(-Inf, 0.01), 
                                          upper = c(Inf, Inf),
                                          control = list(maxit = 500))
                             
                             u_i_opt <- fit$par[1]
                             z_i_opt <- fit$par[2]
                             
                             # prediction
                             y_pop_curr <- private$calc_age_effect(age_curr) + 
                               self$beta_bmi * bmi_curr + 
                               self$beta_hba1c * hba1c_curr + 
                               self$beta_crea * crea_curr
                             
                             # target quantile
                             beta_tau_target <- qnorm(self$tau_target) * sigma_eps
                             
                             # final PJQM2 formula 
                             log_threshold <- y_pop_curr + u_i_opt + (z_i_opt * beta_tau_target)
                             
                             return(exp(log_threshold))
                           }
                       ),
                       
                       private = list(
                         calc_age_effect = function(age) {
                           effect <- numeric(length(age))
                           for (i in 1:length(age)) {
                             a <- age[i]
                             if (a < 30) {
                               effect[i] <- 0
                             } else if (a >= 30 & a < 50) {
                               effect[i] <- 0.021 * (a - 30)
                             } else if (a >= 50 & a <= 60) {
                               effect[i] <- (0.021 * 20) + 0.031 * (a - 50)
                             } else {
                               effect[i] <- (0.021 * 20) + (0.031 * 10) + 0.033 * (a - 60)
                             }
                           }
                           return(effect)
                         }
                       )
)
# Child Class: GAMLSS Model, retrieved values from Excel table
GAMLSS_Model <- R6Class("GAMLSS_Model",
                        inherit = sNfL_Model,
                        
                        public = list(
                          alpha = NULL,
                          z_cutoff = NULL,
                          
                          initialize = function(alpha = 0.05) {
                            super$initialize(name = "GAMLSS (Benkert et al.)")
                            self$alpha <- alpha
                            self$z_cutoff <- qnorm(1 - alpha)
                          },
                          
                          calculate_threshold = function(patient_data = NULL, new_age = NULL) {
                            return(NA) 
                          },
                          
                          predict_relapse = function(new_value, patient_data = NULL, new_age = NULL) {
                            z_score <- patient_data$current_zscore
                            
                            if (is.null(z_score) || is.na(z_score)) {
                              return(list(is_relapse = FALSE, calculated_threshold = NA))
                            }
                            
                            is_relapse <- z_score > self$z_cutoff
                            
                            return(list(
                              is_relapse = is_relapse,
                              calculated_threshold = NA
                            ))
                          }
                        )
)
