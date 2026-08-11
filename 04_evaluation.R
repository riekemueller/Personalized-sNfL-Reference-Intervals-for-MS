evaluate_patient_trajectory <- function(pat_data, models, baseline_n) {
  
  if (is.null(pat_data) || !is.data.frame(pat_data)) return(NULL)
  
  if (nrow(pat_data) <= baseline_n) {
    return(NULL)
  }
  
  results <- list()
  
  for (t in (baseline_n + 1):nrow(pat_data)) {
    
    history_data <- pat_data[1:(t-1), ]
    current_visit <- pat_data[t, ]
    
    patient_history_list <- list(
      snfl_history  = history_data$sNfL_Obs,  
      age_history   = history_data$Current_Age,
      bmi_history   = history_data$BMXBMI,
      hba1c_history = history_data$LBXGH,
      crea_history  = history_data$LBXSCR,
      current_bmi   = current_visit$BMXBMI,
      current_hba1c = current_visit$LBXGH,
      current_crea  = current_visit$LBXSCR
    )
    
    actual_relapse <- current_visit$Relapse_Effect > 0.01 
    
    for (mod_name in names(models)) {
      model <- models[[mod_name]]
      
      pred <- model$predict_relapse(
        new_value = current_visit$sNfL_Obs,
        patient_data = patient_history_list, 
        new_age = current_visit$Current_Age
      )
      
      alarm <- pred$is_relapse
      
      tp <- ifelse(actual_relapse == TRUE  && alarm == TRUE,  1, 0)
      fp <- ifelse(actual_relapse == FALSE && alarm == TRUE,  1, 0)
      tn <- ifelse(actual_relapse == FALSE && alarm == FALSE, 1, 0)
      fn <- ifelse(actual_relapse == TRUE  && alarm == FALSE, 1, 0)
      
      results[[length(results) + 1]] <- data.frame(
        patient_id = current_visit$SEQN,
        model = mod_name,
        time_point = t,
        true_relapse = actual_relapse,
        predicted_alarm = alarm,
        threshold = pred$calculated_threshold,
        TP = tp, FP = fp, TN = tn, FN = fn
      )
    }
  }
  
  return(do.call(rbind, results))
}


# grid search
run_full_evaluation <- function(cohort_data, config, output_file) {
  
  grid <- config
  all_metrics <- list()
  
  cat("\n=======================================================\n")
  cat("Start Simulation & Evaluation. In total", nrow(grid), "scenarios\n")
  cat("=======================================================\n")
  
  eval_cohort <- cohort_data
  
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  
  for (i in 1:nrow(grid)) {
    params <- grid[i, ]
    cat(sprintf("[%d/%d] Evaluation -> Baseline: %d | Freq: %.2f | Alpha: %.3f\n", 
                i, nrow(grid), params$baseline_length, params$sampling_freq, params$alpha))
    
    scenario_results <- list()
    
    models <- list(
      FixedCutoff   = FixedCutoff_Model$new(cutoff = 12.9),
      WithinSubject = WithinSubject_Model$new(cv_i = 0.086, alpha = params$alpha),
      WithinPerson  = WithinPerson_Model$new(alpha = params$alpha),
      PJQM2         = PJQM2_Model$new(alpha = params$alpha)
    )
    
    for (p in 1:nrow(eval_cohort)) {
      pat <- eval_cohort[p, ]
      
      sim_data <- generate_longitudinal_baseline(pat, follow_up_years = 35, freq_years = params$sampling_freq) %>%
        inject_relapses_and_noise(arr = 0.3, kappa_disp = 0.5)
      
      pat_res <- evaluate_patient_trajectory(sim_data$observed, models, baseline_n = params$baseline_length)
      
      if (!is.null(pat_res)) {
        scenario_results[[length(scenario_results) + 1]] <- pat_res
      }
    }
    
    if (length(scenario_results) > 0) {
      scenario_df <- do.call(rbind, scenario_results)
      
      metrics_df <- scenario_df %>%
        group_by(model) %>%
        summarise(
          Total_Visits = n(),
          Total_TP = sum(TP),
          Total_FP = sum(FP),
          Total_TN = sum(TN),
          Total_FN = sum(FN),
          TPR_Sensitivity = Total_TP / max(1, (Total_TP + Total_FN)),
          FPR_False_Alarm = Total_FP / max(1, (Total_FP + Total_TN)),
          TNR_Specificity = Total_TN / max(1, (Total_TN + Total_FP)),
          NPV = Total_TN / max(1, (Total_TN + Total_FN)),
          .groups = "drop"
        ) %>%
        mutate(
          baseline_length = params$baseline_length,
          sampling_freq = params$sampling_freq,
          alpha = params$alpha
        )
      
      all_metrics[[i]] <- metrics_df
      
      if (i == 1) {
        write.table(metrics_df, file = output_file, sep = ",", 
                    row.names = FALSE, col.names = TRUE)
      } else {
        write.table(metrics_df, file = output_file, sep = ",", 
                    row.names = FALSE, col.names = FALSE, append = TRUE)
      }
      # ==========================================
    }
  }
  
  cat("\n>>> All experiments finished!\n")
  return(do.call(rbind, all_metrics))
}