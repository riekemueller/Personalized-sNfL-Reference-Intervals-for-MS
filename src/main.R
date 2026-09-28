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

# Start experiments, if GAMLSS_data is available, other wise comment first row and do not include it in final_experiments_results
gamlss_results <- read_excel(here("thesis_experiments", "data", "GAMLSS_data.xlsx"))
csv_pfad <- here("thesis_experiments/results/results_metrics_final.csv")
final_experiment_results <- run_full_evaluation(cohort_data, experiment_config, output_file = csv_pfad, gamlss_data = gamlss_results)
print(head(final_experiment_results))