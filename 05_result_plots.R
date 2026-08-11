# ==============================================================================
# 05_results_plots.R
# ==============================================================================

# Imports
library(dplyr)
library(ggplot2)
library(tidyr) 

# Daten laden
df <- read.csv("C:/Users/Rieke/Documents/thesis_experiments/results_metrics_20260810_server.csv", stringsAsFactors = FALSE)

model_colors <- c("FixedCutoff" = "#E41A1C", 
                  "PJQM2" = "#4DAF4A", 
                  "WithinPerson" = "#377EB8", 
                  "WithinSubject" = "#984EA3")

# baseline length
p1_data <- df %>%
  group_by(model, baseline_length) %>%
  summarise(Mean_TPR = mean(TPR_Sensitivity), .groups = "drop")

p1 <- ggplot(p1_data, aes(x = baseline_length, y = Mean_TPR, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = model_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Impact of Baseline Length on Sensitivity",
       x = "Number of Baseline Measurements (n)",
       y = "Mean Sensitivity (TPR)",
       color = "Model") +
  theme_minimal(base_size = 14)

print(p1)
ggsave("C:/Users/Rieke/Documents/thesis_experiments/plots/plot_1_baseline_length.png", plot = p1, width = 8, height = 5)

# sampling frequency
p2_data <- df %>%
  group_by(model, sampling_freq) %>%
  summarise(Mean_TPR = mean(TPR_Sensitivity), .groups = "drop")

p2 <- ggplot(p2_data, aes(x = as.factor(sampling_freq), y = Mean_TPR, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = model_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Impact of Sampling Frequency on Relapse Detection",
       x = "Sampling Frequency (Years)",
       y = "Mean Sensitivity (TPR)",
       fill = "Model") +
  theme_minimal(base_size = 14)

print(p2)
ggsave("C:/Users/Rieke/Documents/thesis_experiments/plots/plot_2_sampling_freq.png", plot = p2, width = 8, height = 5)

# alpha
p3_data <- df %>%
  group_by(model, alpha) %>%
  summarise(Mean_FPR = mean(FPR_False_Alarm), 
            Mean_TPR = mean(TPR_Sensitivity), 
            .groups = "drop")

p3_data_long <- p3_data %>%
  pivot_longer(cols = c(Mean_TPR, Mean_FPR), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "Mean_TPR", "1. Sensitivity (TPR)", "2. False Positive Rate (FPR)"))

p3 <- ggplot(p3_data_long, aes(x = as.factor(alpha), y = Value, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Metric, scales = "free_y") + 
  scale_color_manual(values = model_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Clinical Trade-off: Sensitivity vs. False Alarms across \u03B1 levels",
       x = "Significance Level (\u03B1)",
       y = "Metric Value",
       color = "Model") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold", size = 12))

print(p3)
ggsave("C:/Users/Rieke/Documents/thesis_experiments/plots/plot_3_alpha_tradeoff.png", plot = p3, width = 10, height = 5)

# roc
p4 <- ggplot(df, aes(x = FPR_False_Alarm, y = TPR_Sensitivity, color = model)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = model_colors) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Clinical Utility: Overall Sensitivity vs. False Positive Rate",
       x = "False Positive Rate (1 - Specificity)",
       y = "Sensitivity (TPR)",
       color = "Model") +
  theme_minimal(base_size = 14) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")

print(p4)
ggsave("C:/Users/Rieke/Documents/thesis_experiments/plots/plot_4_overall_roc.png", plot = p4, width = 8, height = 6)