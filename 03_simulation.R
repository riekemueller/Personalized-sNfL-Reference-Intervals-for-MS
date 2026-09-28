# Imports
library(dplyr)
library(tidyr)
library(ggplot2)

# Calculate linear age effect
#' @param ages Vector with age in years
#' @return cumulative age effect
calc_age_effect <- function(age) {
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

# calculates relapse effect
#' @param age current age
#' @param relapse_ages
#' @param t_b time for build-up (5 months)
#' @param K_0 calculated in thesis
#' @param peak_amp relapse amplitude in log scale (log(1.43))
#' @param K calculated in thesis
calc_relapse_effect <- function(age, relapse_ages, t_b = 5/12, K_0 = 0.86, peak_amp = log(1.43), K = 8.32) {
  effect <- numeric(length(age))
  for (i in 1:length(age)) {
    current_age <- age[i]
    total_relapse_activity <- 0
    for (tau_j in relapse_ages) {
      delta_t <- current_age - tau_j
      # Phase 1: Zero-order Build-up (-t_b <= delta_t <= 0)
      if (delta_t >= -t_b && delta_t <= 0) {
        total_relapse_activity <- total_relapse_activity + K_0 * (delta_t + t_b)
      } 
      # Phase 2: First-order Wash-out (delta_t > 0)
      else if (delta_t > 0) {
        total_relapse_activity <- total_relapse_activity + peak_amp * exp(-K * delta_t)
      }
    }
    effect[i] <- total_relapse_activity
  }
  return(effect)
}

#' calculates analytical noise, based on concentration
#' @param true_snfl_values true snfl concentration
add_clinical_assay_noise <- function(true_snfl_values) {

  cv_steps <- case_when(
    true_snfl_values < 7.2    ~ 0.052,  # mean: 4.4
    true_snfl_values >= 7.2   & true_snfl_values < 13.25  ~ 0.086,  # mean: 10.0
    true_snfl_values >= 13.25 & true_snfl_values < 20.65  ~ 0.051,  # mean: 16.5
    true_snfl_values >= 20.65 & true_snfl_values < 54.7   ~ 0.023,  # mean: 24.8
    true_snfl_values >= 54.7  & true_snfl_values < 153.6  ~ 0.015,  # mean: 84.6
    true_snfl_values >= 153.6                             ~ 0.016   # mean: 222.6
  )
  
  analytical_sd <- cv_steps * true_snfl_values
  
  noisy_snfl <- rnorm(n = length(true_snfl_values), 
                      mean = true_snfl_values, 
                      sd = analytical_sd)
  
  noisy_snfl <- pmax(noisy_snfl, 0.1) 
  
  return(noisy_snfl)
}


# Phase 1: Parameter estimation, build the regression model
#' @param cohort_data cleaned cohort NHANES dataset
#' @param follow_up_years How many years are simulated?
#' @param freq_years Intervals of follow-ups in years
#' @param start_age optional start data, if NULL, the actual NHANES age is used
#' @param var_beta0i variance of random intercept
#' @param force_beta0i optional: fixed beta0i value (useful for MS patient 0)
generate_longitudinal_baseline <- function(cohort_data, follow_up_years=35, freq_years = 0.25, start_age = NULL, var_beta0i = 0.20,
                                           force_beta0i = NULL){
  
  # beta coefficients
  beta_bmi   <- -0.020  # BMI
  beta_hba1c <- 0.090   # HbA1c
  beta_crea  <- 0.467   # Creatinine
  
  # get std from variance for beta_0i
  sd_beta0i <- sqrt(var_beta0i)
  
  # individual random intercept (beta_0i)
  baseline_calc <- cohort_data %>%
    mutate(
      beta_0i = if (!is.null(force_beta0i)) {
        force_beta0i
      } else {
        rnorm(n(), mean = 0, sd = sd_beta0i)
      },
      
      confounder_effects_nhanes = calc_age_effect(RIDAGEYR) + 
        (beta_bmi * BMXBMI) + 
        (beta_hba1c * LBXGH) + 
        (beta_crea * LBXSCR),
      
      beta_0 = log(SSSNFL) - beta_0i - confounder_effects_nhanes,
     
      Sim_Base_Age = if (!is.null(start_age)) start_age else RIDAGEYR
  )
  
  # create latent grid
  time_points_latent <- seq(0, follow_up_years, by = 0.05)
  latent_data <- baseline_calc %>%
    crossing(Time_Years = time_points_latent) %>%
    mutate(
      Current_Age = Sim_Base_Age + Time_Years,
      True_Log_sNfL = beta_0 + beta_0i + calc_age_effect(Current_Age) + 
        (beta_bmi * BMXBMI) + (beta_hba1c * LBXGH) + (beta_crea * LBXSCR),
      sNfL_Latent = exp(True_Log_sNfL) # Reiner Basiswert
    ) %>% select(SEQN, Status, Time_Years, Current_Age, True_Log_sNfL, sNfL_Latent)
  
  # create observed grid (with correct follow-up intervals)
  time_points_obs <- seq(0, follow_up_years, by = freq_years)
  
  obs_data <- baseline_calc %>%
    crossing(Time_Years = time_points_obs) %>%
    mutate(
      Current_Age = Sim_Base_Age + Time_Years,
      True_Log_sNfL = beta_0 + 
                      beta_0i + 
                      calc_age_effect(Current_Age) + 
                      (beta_bmi * BMXBMI) + 
                      (beta_hba1c * LBXGH) + 
                      (beta_crea * LBXSCR),
      
      True_sNfL = exp(True_Log_sNfL)
    ) %>%
    select(SEQN, Status, Gender, MS_Medication, Time_Years, Current_Age, 
           BMXBMI, LBXGH, LBXSCR, beta_0, beta_0i, True_Log_sNfL, True_sNfL) %>%
    arrange(SEQN, Time_Years)
  
  return(list(latent = latent_data, observed = obs_data))
}

# Phase 2: Longitudinal Data Generation using Relapse Engine
#' @param sim_data_list List from phase 1, having latent data and the observed data as tables
#' @param arr Annualized Relapse Rate (set to 0.3)
#' @param kappa_disp Dispersions parameter for Negative Binomial (set to 0.5)
inject_relapses_and_noise <- function(sim_data_list, 
                                      arr = 0.3,              
                                      kappa_disp = 0.5) {   
  #set.seed(42) # for reproducability, uncommend for plot patient 0 prediction interval as no monte carlo simulation is possible
  
  latent_list <- list()
  obs_list <- list()
  patient_ids <- unique(sim_data_list$observed$SEQN)
  
  for(id in patient_ids) {
    p_lat <- sim_data_list$latent %>% filter(SEQN == id)
    p_obs <- sim_data_list$observed %>% filter(SEQN == id)
    
    is_ms <- unique(p_obs$Status) == "MS Patient"
    max_time <- max(p_obs$Time_Years)
    start_age <- min(p_obs$Current_Age)
    
    relapse_times <- c()
    relapse_ages <- c()
    
    if(is_ms && max_time > 0) {
      n_relapses <- rnbinom(1, size = kappa_disp, mu = arr * max_time)
      
      if(n_relapses > 0) {
        relapse_times <- sort(runif(n_relapses, min = 0, max = max_time))
        relapse_ages <- start_age + relapse_times
      }
    }
    
    # relapse on latent data
    p_lat <- p_lat %>%
      mutate(
        Relapse_Effect = calc_relapse_effect(Current_Age, relapse_ages = relapse_ages),
        Log_sNfL_Latent = True_Log_sNfL + Relapse_Effect,
        sNfL_Latent = exp(Log_sNfL_Latent)
      )
    
    # relapse on observed data
    p_obs <- p_obs %>%
      mutate(
          Relapse_Effect = calc_relapse_effect(Current_Age, relapse_ages = relapse_ages),
          Log_sNfL_True_Total = True_Log_sNfL + Relapse_Effect,
          sNfL_True_Raw = exp(Log_sNfL_True_Total),
          sNfL_Obs = add_clinical_assay_noise(sNfL_True_Raw)
      )
    
    latent_list[[as.character(id)]] <- p_lat
    obs_list[[as.character(id)]] <- p_obs
  }
  
  return(list(
    latent = bind_rows(latent_list),
    observed = bind_rows(obs_list)
  ))
}


plot_data <- function(sim_data_list, p_id = NULL, cutoff_value = 12.9, filename = NULL) {
  
  sim_latent <- sim_data_list$latent
  sim_obs    <- sim_data_list$observed
  
  if(is.null(p_id)) p_id <- unique(sim_obs$SEQN)[1]
  
  sim_latent <- sim_latent %>% filter(SEQN == p_id)
  sim_obs    <- sim_obs %>% filter(SEQN == p_id)
  
  plot_relapse <- ggplot() +
    geom_line(data = sim_latent, aes(x = Current_Age, y = sNfL_Latent, color = "Latent disease trajectory"), linewidth = 0.8) +
    geom_point(data = sim_obs, aes(x = Current_Age, y = sNfL_Obs, fill = "Observed clinical measurements (3-month interval)"), 
               shape = 21, color = "black", size = 2, stroke = 0.3, alpha = 0.8) +
    geom_hline(aes(yintercept = cutoff_value, color = "Cutoff value"), linewidth = 0.8, linetype = "solid") +
    
    scale_color_manual(name = NULL, values = c("Latent disease trajectory" = "turquoise4", "Cutoff value" = "grey50")) +
    scale_fill_manual(name = NULL, values = c("Observed clinical measurements (3-month interval)" = "black")) +
    
    labs(x = "Age (Years)", y = "sNfL concentration (pg/mL)") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom", 
      legend.box = "vertical",
      legend.margin = margin(t = -0.5, unit = "cm"),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(filename)) {
    require(tikzDevice)
    cat(paste("\nSpeichere Plot als LaTeX TikZ Datei:", filename, "\n"))
    tikz(file = filename, width = 6, height = 5)
    print(plot_relapse)
    dev.off()
  }
  
  print(plot_relapse)
  return(plot_relapse)
}