# Configurable parameters for experiments
set.seed(42) #  for reproducability

sim_baseline_lengths <- c(2, 3, 5, 7, 10, 20, 30)
sim_sampling_frequencies <- c(0.25, 0.5, 1.0)
sim_alpha_levels <- c(0.10, 0.05, 0.025, 0.01)

experiment_config <- expand.grid(
  baseline_length = sim_baseline_lengths,
  sampling_freq = sim_sampling_frequencies,
  alpha = sim_alpha_levels
)