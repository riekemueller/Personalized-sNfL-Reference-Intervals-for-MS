# hyperparam optimization, to get FPR for n=3 to 5%
library(here)
library(dplyr)
library(tidyr)

source(here("thesis_experiments/00_config.R"))
source(here("thesis_experiments/01_data.R"))
source(here("thesis_experiments/02_models.R"))
source(here("thesis_experiments/03_simulation.R"))
source(here("thesis_experiments/04_evaluation.R"))

set.seed(42)
cohort_data <- get_nhanes_cohort()

# take 100 healthy patients
tuning_cohort <- cohort_data %>%
  filter(Status != "MS Patient") %>% 
  sample_n(100)

cat(sprintf("Start Tuning with %d healthy patients \n", nrow(tuning_cohort)))

simulated_patients <- list()
for (i in 1:nrow(tuning_cohort)) {
  pat <- tuning_cohort[i, ]
  
  # 10 years, quarterly frequency
  sim_data <- generate_longitudinal_baseline(pat, follow_up_years = 10, freq_years = 0.25) %>%
    inject_relapses_and_noise() 
  
  simulated_patients[[i]] <- sim_data$observed
}

# hyperparam grid
test_baseline_length <- 3
test_alpha <- 0.05

hyper_grid <- expand.grid(
  lambda_u = c(0.001, 0.01, 0.03, 0.05, 0.1),
  lambda_z = c(0.1, 0.5, 1.0, 5.0)
)

tuning_results <- list()

cat("Start Grid-Search over", nrow(hyper_grid), "combinations (n =", test_baseline_length, ")...\n")
cat("----------------------------------------------------------\n")

# grid-search
for (g in 1:nrow(hyper_grid)) {
  l_u <- hyper_grid$lambda_u[g]
  l_z <- hyper_grid$lambda_z[g]
  
  cat(sprintf("[%2d/%d] Test lambda_u = %4.1f | lambda_z = %4.1f ... ", g, nrow(hyper_grid), l_u, l_z))
  
  pjqm2_model <- PJQM2_Model$new(
    alpha = test_alpha,
    lambda_u = l_u,
    lambda_z = l_z
  )
  
  models_list <- list(PJQM2_Tuning = pjqm2_model)
  
  grid_eval_results <- list()
  for (i in seq_along(simulated_patients)) {
    pat_data <- simulated_patients[[i]]
    
    res <- evaluate_patient_trajectory(
      pat_data = pat_data, 
      models = models_list, 
      baseline_n = test_baseline_length
    )
    
    if (!is.null(res)) {
      grid_eval_results[[length(grid_eval_results) + 1]] <- res
    }
  }
  
  # FPR calculation
  all_res_df <- do.call(rbind, grid_eval_results)
  
  metrics <- all_res_df %>%
    summarise(
      Total_TN = sum(TN),
      Total_FP = sum(FP),
      FPR = Total_FP / max(1, (Total_FP + Total_TN)) 
    ) %>%
    mutate(
      lambda_u = l_u,
      lambda_z = l_z
    )
  
  tuning_results[[g]] <- metrics
  cat(sprintf("FPR: %.4f\n", metrics$FPR))
}

final_tuning_df <- do.call(rbind, tuning_results)

final_tuning_df <- final_tuning_df %>%
  mutate(Abweichung_vom_Ziel = abs(FPR - test_alpha)) %>%
  arrange(Abweichung_vom_Ziel)

cat("\Results:\n")
print(head(final_tuning_df, 5))

# best combination
best_combo <- final_tuning_df %>% slice(1)

cat("\n>>> Best results:\n")
cat(sprintf("lambda_u = %.1f\n", best_combo$lambda_u))
cat(sprintf("lambda_z = %.1f\n", best_combo$lambda_z))
cat(sprintf("To get for n=%d a FPR of %.2f %%\n", test_baseline_length, best_combo$FPR * 100))