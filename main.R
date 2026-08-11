# imports
library(here)

# load modules
source(here("thesis_experiments/00_config.R"))
source(here("thesis_experiments/01_data.R"))
source(here("thesis_experiments/02_models.R"))
source(here("thesis_experiments/03_simulation.R"))
source(here("thesis_experiments/04_evaluation.R"))

# load cohorts
cohort_data <- get_nhanes_cohort()
#cohort_data %>% View()

# estimate variance of healthy groups
#empiric_var <- estimate_control_variance(cohort_data)
#cat("Calculated variance of beta_0i:", round(empiric_var, 4), "\n\n")

# one healthy patient
#healthy_patient <- cohort_data %>% filter(Status == "Control") %>% slice(1)
#sim_healthy <- generate_longitudinal_baseline(healthy_patient, follow_up_years = 5)
#sim_healthy <- inject_relapses_and_noise(sim_healthy)

# for patient 0
#pat_0 <- cohort_data %>% filter(SEQN == 74929)
#sim_pat0 <- generate_longitudinal_baseline(pat_0, follow_up_years = 50, start_age = 30, force_beta0i = 0.15) %>%
 # inject_relapses_and_noise()

# Plot erstellen
#plot_data(sim_pat0)

# Start experiments
csv_pfad <- here("thesis_experiments/results_metrics_20260810.csv")
final_experiment_results <- run_full_evaluation(cohort_data, experiment_config, output_file = csv_pfad)
print(head(final_experiment_results))
