# imports
library(here)
library(readxl)

# load modules
source(here("thesis_experiments/00_config.R"))
source(here("thesis_experiments/01_data.R"))
source(here("thesis_experiments/02_models.R"))
source(here("thesis_experiments/03_simulation.R"))
source(here("thesis_experiments/04_evaluation.R"))

# load cohorts
cohort_data <- get_nhanes_cohort()
#export_simulated_values_kuhle(cohort_data)

#cohort_data %>% View()

# estimate variance of healthy groups
#empiric_var <- estimate_control_variance(cohort_data)
#cat("Calculated variance of beta_0i:", round(empiric_var, 4), "\n\n")

# one healthy patient
#healthy_patient <- cohort_data %>% filter(Status == "Control") %>% slice(1)
#sim_healthy <- generate_longitudinal_baseline(healthy_patient, follow_up_years = 5)
#sim_healthy <- inject_relapses_and_noise(sim_healthy)

# for patient 0
#pat_ms <- cohort_data %>% filter(Status == "MS Patient")
#sim_pat0 <- generate_longitudinal_baseline(pat_0, follow_up_years = 50, start_age = 30, force_beta0i = 0.15) %>%
# inject_relapses_and_noise()

# Plot erstellen
#plot_data(sim_pat0)

# Start experiments
gamlss_results <- read_excel(here("thesis_experiments", "data", "Kuhle_results.xlsx"))
csv_pfad <- here("thesis_experiments/results_metrics_with_gamlss.csv")
final_experiment_results <- run_full_evaluation(pat_ms, experiment_config, output_file = csv_pfad, gamlss_data = gamlss_results)
print(head(final_experiment_results))
