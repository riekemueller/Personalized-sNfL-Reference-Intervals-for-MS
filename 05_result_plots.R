# ==============================================================================
# Master Thesis: Chapter 5 Plot & Metrics Generation Script
# Author: Rieke Müller
# Project: Personalized sNfL Reference Intervals for Multiple Sclerosis
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Imports & Global Setup
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr) 
library(tikzDevice)

# Eigene Prozent-Formatierung für LaTeX (maskiertes %)
tex_percent <- scales::label_percent(suffix = "\\%")

# Konsistente Farbpalette für die Modelle
model_colors <- c(
  "FixedCutoff"   = "#E41A1C", # Rot (Klinischer Standard)
  "PJQM2"         = "#4DAF4A", # Grün (Penalized Joint Quantile)
  "WithinPerson"  = "#377EB8", # Blau (Personalisierte Biologische Varianz)
  "WithinSubject" = "#984EA3", # Violett (Personalisierte Population-Varianz)
  "GAMLSS"        = "#FF7F00"  # Orange
)

# Saubere Metric-Beschriftungen für Facets (erweitert um Balanced Accuracy)
metric_labels <- c(
  "Sensitivity"       = "Sensitivity (TPR)",
  "Specificity"       = "Specificity (TNR)",
  "Precision"         = "Precision (PPV)",
  "NPV"               = "Negative Pred. Value",
  "F1_Score"          = "F1-Score",
  "F2_Score"          = "F2-Score",
  "F0_5_Score"        = "F0.5-Score",
  "Balanced_Accuracy" = "Balanced Accuracy",
  "Error_Rate"        = "Overall Error Rate"
)

# ------------------------------------------------------------------------------
# 1. Dateneingabe & Metrikberechnung
# ------------------------------------------------------------------------------
data_path  <- "C:/Users/Rieke/Documents/thesis_experiments/results/results_metrics_final.csv"
output_dir <- "C:/Users/Rieke/Documents/thesis_experiments/plots/"

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Berechnung aller klinischen & diagnostischen Metriken
df <- df %>%
  mutate(
    Specificity       = ifelse((Total_TN + Total_FP) == 0, NA, Total_TN / (Total_TN + Total_FP)),
    Sensitivity       = ifelse((Total_TP + Total_FN) == 0, NA, Total_TP / (Total_TP + Total_FN)),
    Precision         = ifelse((Total_TP + Total_FP) == 0, NA, Total_TP / (Total_TP + Total_FP)),
    NPV               = ifelse((Total_TN + Total_FN) == 0, NA, Total_TN / (Total_TN + Total_FN)),
    Error_Rate        = (Total_FP + Total_FN) / (Total_TP + Total_FP + Total_FN + Total_TN),
    F1_Score          = ifelse(is.na(Precision) | is.na(Sensitivity) | (Precision + Sensitivity) == 0, NA, 
                               2 * (Precision * Sensitivity) / (Precision + Sensitivity)),
    F2_Score          = ifelse(is.na(Precision) | is.na(Sensitivity) | (4 * Precision + Sensitivity) == 0, NA, 
                               5 * (Precision * Sensitivity) / (4 * Precision + Sensitivity)),
    F0_5_Score        = ifelse(is.na(Precision) | is.na(Sensitivity) | (0.25 * Precision + Sensitivity) == 0, NA, 
                               1.25 * (Precision * Sensitivity) / (0.25 * Precision + Sensitivity)),
    Balanced_Accuracy = (Sensitivity + Specificity) / 2,
    Youden_Index      = Sensitivity + Specificity - 1
  )

# ==============================================================================
# SECTION 5.1: Initialization Stability across Baseline Lengths (n)
# ==============================================================================
p1_data <- df %>%
  group_by(model, baseline_length) %>%
  summarise(across(c(Sensitivity, Specificity, Precision, NPV, F1_Score, F2_Score, F0_5_Score, Balanced_Accuracy, Error_Rate), 
                   ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = Sensitivity:Error_Rate, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Sensitivity", "Specificity", "Precision", "NPV", "F1_Score", "F2_Score", "F0_5_Score", "Balanced_Accuracy", "Error_Rate")))

# Visuelle Markierung der harten klinischen Grenze (95% Spezifität)
hline_p1 <- data.frame(Metric = factor("Specificity", levels = levels(p1_data$Metric)), Value = 0.95)

# ------------------------------------------------------------------------------
# 5.1a Plot: jetzt im 3x3-Layout (9 Metriken -> ncol = 3)
# ------------------------------------------------------------------------------
# Y-Achsen-Grenzen pro Metrik berechnen: 5% Puffer unten, 10% Puffer oben,
# aber niemals unter 0% bzw. über 100%
y_limits_data <- p1_data %>%
  group_by(Metric) %>%
  summarise(
    data_min = min(Value, na.rm = TRUE),
    data_max = max(Value, na.rm = TRUE),
    range    = data_max - data_min,
    y_min    = pmax(0, data_min - 0.25 * range),
    y_max    = pmin(1, data_max + 0.25 * range),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(y_min, y_max), names_to = "bound", values_to = "Value") %>%
  select(Metric, Value)

p1 <- ggplot(p1_data, aes(x = baseline_length, y = Value, color = model, group = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_blank(data = y_limits_data, 
             aes(x = min(p1_data$baseline_length), y = Value),
             inherit.aes = FALSE) +
  geom_hline(data = hline_p1, aes(yintercept = Value), color = "red", linetype = "dotted", linewidth = 1) +
  facet_wrap(~ Metric, scales = "free_y", labeller = as_labeller(metric_labels), ncol = 3) + 
  scale_color_manual(values = model_colors) +
  scale_y_continuous(labels = tex_percent, expand = expansion(mult = 0)) +
  labs(x = "Number of Historical Baseline Measurements ($n$)",
       y = "Metric Value",
       color = "Model") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

# Breite/Höhe an quadratischeres 3x3-Layout angepasst (vorher 14x7 für 3x4)
ggsave(paste0(output_dir, "plot_5_1_metrics_vs_baseline.png"), plot = p1, width = 10.5, height = 10)

tikz(paste0(output_dir, "plot_5_1_metrics_vs_baseline.tex"), width = 10.5, height = 10)
print(p1)
dev.off()

# ------------------------------------------------------------------------------
# 5.1b Stabilitätsanalyse: exakte Berechnung basierend auf der Corridor-of-
# Stability-Definition (Schönbrodt & Perugini, 2013), erweitert um ein
# monotones Kriterium: Eine Metrik gilt bei n als stabil, wenn für ALLE
# nachfolgend getesteten n' > n gilt, dass sich der Wert monoton (gleiches
# Vorzeichen der Änderung, oder keine Änderung) und um weniger als 2
# Prozentpunkte relativ zum jeweils vorherigen n verändert.
# ------------------------------------------------------------------------------

#' Findet das kleinste n, ab dem eine Metrikkurve als "stabil" gilt.
#'
#' @param baseline_lengths Numerischer Vektor der getesteten n-Werte, aufsteigend sortiert.
#' @param values Numerischer Vektor der Metrikwerte (gleiche Reihenfolge wie baseline_lengths), als Proportion (0-1).
#' @param threshold Schwelle in Prozentpunkten (Default: 0.02 = 2 pp).
#' @return Liste mit stable (logical), stable_at_n (numeric oder NA)
find_stability_point <- function(baseline_lengths, values, 
                                 abs_threshold = 0.02, 
                                 min_remaining_steps = 2) {
  ord <- order(baseline_lengths)
  n_sorted <- baseline_lengths[ord]
  v_sorted <- values[ord]
  
  k <- length(v_sorted)
  if (k < 2) return(list(stable = FALSE, stable_at_n = NA_real_))
  
  diffs <- diff(v_sorted)
  
  # Nur Startpunkte prüfen, bei denen noch mindestens 
  # min_remaining_steps Schritte folgen (verhindert Trivialität am Rand)
  max_i <- k - min_remaining_steps
  if (max_i < 1) return(list(stable = FALSE, stable_at_n = NA_real_))
  
  for (i in seq_len(max_i)) {
    remaining_diffs <- diffs[i:(k - 1)]
    below_threshold <- all(abs(remaining_diffs) < abs_threshold)
    signs <- sign(remaining_diffs)
    nonzero_signs <- signs[signs != 0]
    is_monotonic <- length(unique(nonzero_signs)) <= 1
    
    if (below_threshold && is_monotonic) {
      return(list(stable = TRUE, stable_at_n = n_sorted[i]))
    }
  }
  
  return(list(stable = FALSE, stable_at_n = NA_real_))
}

# Stabilität nur für die vier aussagekräftigen Metriken berechnen
metrics_for_stability <- c("Sensitivity", "Specificity", "Balanced_Accuracy", "Error_Rate")

stability_results <- p1_data %>%
  filter(Metric %in% metrics_for_stability) %>%
  group_by(model, Metric) %>%
  arrange(baseline_length, .by_group = TRUE) %>%
  summarise(
    result = list(find_stability_point(baseline_length, Value, 
                                       abs_threshold = 0.02, 
                                       min_remaining_steps = 2)),
    .groups = "drop"
  ) %>%
  mutate(
    stable      = purrr::map_lgl(result, "stable"),
    stable_at_n = purrr::map_dbl(result, "stable_at_n")
  ) %>%
  select(-result)

# In breites Format für die LaTeX-Tabelle bringen (Modelle als Zeilen, Metriken als Spalten)
stability_table <- stability_results %>%
  mutate(Label = ifelse(stable, paste0("Yes, at $n=", stable_at_n, "$"), "No")) %>%
  select(model, Metric, Label) %>%
  pivot_wider(names_from = Metric, values_from = Label)

write.csv(stability_table, paste0(output_dir, "stability_table_5_1.csv"), row.names = FALSE)

# Zusätzlich: Rohwerte (exakte Prozentwerte je n) für Prüfzwecke/Anhang exportieren
write.csv(p1_data, paste0(output_dir, "5_1_metrics_vs_baseline_raw_values_.csv"), row.names = FALSE)

cat("\nStability table:\n")
print(stability_table)

# ==============================================================================
# SECTION 5.2: Impact of Monitoring Frequency on Detection Performance
# ==============================================================================
p3_data <- df %>%
  group_by(model, sampling_freq) %>%
  summarise(Mean_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(p3_data, aes(x = as.factor(sampling_freq), y = Mean_Balanced_Accuracy, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = model_colors) +
  scale_x_discrete(labels = c("0.0833333333333333" = "Monthly", 
                              "0.166666666666667" = "Every 2 months", 
                              "0.25" = "Quarterly", 
                              "0.5" = "Biannual", 
                              "1" = "Annual")) +
  scale_y_continuous(labels = tex_percent, limits = c(0, 0.8)) +
  labs(x = "Sampling Frequency",
       y = "Mean Balanced Accuracy",
       fill = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(paste0(output_dir, "plot_5_2_sampling_frequency.png"), plot = p3, width = 7.5, height = 4.5)

tikz(paste0(output_dir, "plot_5_2_sampling_frequency.tex"), width = 7.5, height = 4.5)
print(p3)
dev.off()

# Exakte Prozentwerte für den Fließtext exportieren
p3_export <- p3_data %>%
  mutate(Balanced_Accuracy_Pct = sprintf("%.2f%%", Mean_Balanced_Accuracy * 100)) %>%
  arrange(model, sampling_freq)

write.csv(p3_export, paste0(output_dir, "5_2_sampling_frequency_values.csv"), row.names = FALSE)

cat("\nExact Balanced Accuracy values by model and sampling frequency:\n")
print(p3_export)

# ==============================================================================
# SECTION 5.3: Influence of Significance Levels (alpha) on Diagnostic Trade-offs
# ==============================================================================
p2_data <- df %>%
  group_by(model, alpha) %>%
  summarise(
    Mean_Sens = mean(Sensitivity, na.rm = TRUE),
    Mean_Spec = mean(Specificity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(FPR = 1 - Mean_Spec) %>%
  select(model, alpha, Sensitivity = Mean_Sens, FPR) %>%
  pivot_longer(cols = c(Sensitivity, FPR), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Sensitivity", "FPR")))

# Visuelle Markierung der harten klinischen Grenze (max 5% FPR)
hline_p2 <- data.frame(Metric = factor("FPR", levels = levels(p2_data$Metric)), Value = 0.05)

p2 <- ggplot(p2_data, aes(x = as.factor(alpha), y = Value, color = model, 
                          linetype = Metric, group = interaction(model, Metric))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2, aes(shape = Metric)) +
  facet_wrap(~ model) + 
  scale_color_manual(values = model_colors) +
  scale_linetype_manual(values = c("Sensitivity" = "solid", "FPR" = "dashed"),
                        labels = c("Sensitivity (TPR)", "False Positive Rate (FPR)")) +
  scale_shape_manual(values = c("Sensitivity" = 16, "FPR" = 17), 
                     labels = c("Sensitivity (TPR)", "False Positive Rate (FPR)")) +
  scale_y_continuous(labels = tex_percent) +
  labs(x = "Significance Level ($\\alpha$)",
       y = "Percentage",
       linetype = "Metric",
       shape = "Metric") +
  guides(color = "none") + 
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(paste0(output_dir, "plot_5_3_metrics_vs_alpha.png"), plot = p2, width = 9, height = 6.5)

tikz(paste0(output_dir, "plot_5_3_metrics_vs_alpha.tex"), width = 9, height = 6.5)
print(p2)
dev.off()

# Exakte Prozentwerte für den Fließtext exportieren
p2_export <- p2_data %>%
  mutate(Value_Pct = sprintf("%.2f%%", Value * 100)) %>%
  arrange(model, Metric, alpha)

write.csv(p2_export, paste0(output_dir, "5_3_metrics_vs_alpha_values.csv"), row.names = FALSE)

cat("\nExact Sensitivity/FPR values by model and alpha:\n")
print(p2_export)

# ==============================================================================
# SECTION 5.4: Global Diagnostic Accuracy — exact values export
# ==============================================================================

# 5.4a Exportiere alle einzelnen Konfigurationen (roh) für Nachschlagezwecke
roc_raw_export <- df %>%
  filter(!is.na(Sensitivity), !is.na(Specificity)) %>%
  mutate(FPR = 1 - Specificity) %>%
  select(model, baseline_length, sampling_freq, alpha, Sensitivity, Specificity, FPR, Balanced_Accuracy) %>%
  arrange(model, FPR, Sensitivity)

write.csv(roc_raw_export, paste0(output_dir, "5_4_roc_raw_configurations.csv"), row.names = FALSE)

# 5.4b Zusammenfassende Kennzahlen pro Modell: Bereich von FPR & Sensitivity,
# sowie tatsächliche Werte in der Low-FPR-Region (<= 25%) — relevant für
# den Vergleich mit PJQM2 und die spätere 5%-Constraint-Section
model_ranges <- df %>%
  filter(!is.na(Sensitivity), !is.na(Specificity)) %>%
  mutate(FPR = 1 - Specificity) %>%
  group_by(model) %>%
  summarise(
    FPR_min           = min(FPR, na.rm = TRUE),
    FPR_max           = max(FPR, na.rm = TRUE),
    Sensitivity_min   = min(Sensitivity, na.rm = TRUE),
    Sensitivity_max   = max(Sensitivity, na.rm = TRUE),
    n_configs         = n(),
    n_configs_FPR_le5 = sum(FPR <= 0.05, na.rm = TRUE),
    n_configs_FPR_le25 = sum(FPR <= 0.25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(FPR_min, FPR_max, Sensitivity_min, Sensitivity_max), ~sprintf("%.2f%%", .x * 100))) %>%
  arrange(FPR_min)

write.csv(model_ranges, paste0(output_dir, "5_4_model_ranges.csv"), row.names = FALSE)

# 5.4c Innerhalb der klinisch relevanten Low-FPR-Region (<=25%): 
# Sensitivity-Spanne jedes Modells, das dort überhaupt Konfigurationen hat
low_fpr_detail <- df %>%
  filter(!is.na(Sensitivity), !is.na(Specificity)) %>%
  mutate(FPR = 1 - Specificity) %>%
  filter(FPR <= 0.25) %>%
  group_by(model) %>%
  summarise(
    n_configs_in_region = n(),
    Sensitivity_min = min(Sensitivity, na.rm = TRUE),
    Sensitivity_max = max(Sensitivity, na.rm = TRUE),
    Sensitivity_median = median(Sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(Sensitivity_min, Sensitivity_max, Sensitivity_median), ~sprintf("%.2f%%", .x * 100))) %>%
  arrange(desc(n_configs_in_region))

write.csv(low_fpr_detail, paste0(output_dir, "5_4_low_fpr_region.csv"), row.names = FALSE)

cat("\nModel ranges (full FPR/Sensitivity spans):\n")
print(model_ranges)
cat("\nLow-FPR region (<=25%) detail:\n")
print(low_fpr_detail)

# Exakte Median-/Mittelwerte pro Modell für den Fließtext
ba_summary <- df %>%
  filter(!is.na(Balanced_Accuracy)) %>%
  group_by(model) %>%
  summarise(
    Median_BA = median(Balanced_Accuracy, na.rm = TRUE),
    Mean_BA   = mean(Balanced_Accuracy, na.rm = TRUE),
    SD_BA     = sd(Balanced_Accuracy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Median_BA)) %>%
  mutate(across(c(Median_BA, Mean_BA, SD_BA), ~sprintf("%.2f%%", .x * 100)))

write.csv(ba_summary, paste0(output_dir, "5_4_balanced_accuracy_summary.csv"), row.names = FALSE)
print(ba_summary)

# ==============================================================================
# SECTION 5.5: Clinical Utility & Optimization Tables Export
# ==============================================================================
# balanced accuracy optimum
optimal_balanced <- df %>%
  filter(!is.na(Balanced_Accuracy)) %>%
  group_by(model) %>%
  slice_max(order_by = Balanced_Accuracy, n = 1, with_ties = FALSE) %>%
  select(model, baseline_length, sampling_freq, alpha, Sensitivity, Specificity, Balanced_Accuracy, F2_Score) %>%
  arrange(desc(Balanced_Accuracy))

write.csv(optimal_balanced, paste0(output_dir, "5_5_optimal_models_balanced_accuracy.csv"), row.names = FALSE)

# F2 score optimum
optimal_f2 <- df %>%
  filter(!is.na(F2_Score)) %>%
  group_by(model) %>%
  slice_max(order_by = F2_Score, n = 1, with_ties = FALSE) %>%
  select(model, baseline_length, sampling_freq, alpha, Sensitivity, Specificity, Precision, F2_Score, Balanced_Accuracy) %>%
  arrange(desc(F2_Score))

write.csv(optimal_f2, paste0(output_dir, "5_5_optimal_models_f2.csv"), row.names = FALSE)

all_clinical_95 <- df %>%
  filter(!is.na(Sensitivity), !is.na(Specificity), Specificity >= 0.95) %>%
  select(model, baseline_length, sampling_freq, alpha, Sensitivity, Specificity, Precision, Balanced_Accuracy, F2_Score) %>%
  arrange(model, desc(Sensitivity), desc(Specificity))

write.csv(all_clinical_95, paste0(output_dir, "5_5_all_models_clinical_95.csv"), row.names = FALSE)

# Best-of-Zusammenfassung (eine Zeile pro Modell) bleibt zusätzlich erhalten,
# für den schnellen Überblick / die Haupttabelle im Text
optimal_clinical_95 <- all_clinical_95 %>%
  group_by(model) %>%
  slice_max(order_by = Sensitivity, n = 1, with_ties = TRUE) %>%
  slice_max(order_by = Specificity, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(optimal_clinical_95, paste0(output_dir, "5_5_optimal_models_clinical_95.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Clinical constraint: Specificity >= 90%
# ALLE Konfigurationen, die den Constraint erfüllen, nicht nur die beste pro Modell
# ------------------------------------------------------------------------------
all_clinical_90 <- df %>%
  filter(!is.na(Sensitivity), !is.na(Specificity), Specificity >= 0.90) %>%
  select(model, baseline_length, sampling_freq, alpha, Sensitivity, Specificity, Precision, Balanced_Accuracy, F2_Score) %>%
  arrange(model, desc(Sensitivity), desc(Specificity))

write.csv(all_clinical_90, paste0(output_dir, "5_5_all_models_clinical_90.csv"), row.names = FALSE)

optimal_clinical_90 <- all_clinical_90 %>%
  group_by(model) %>%
  slice_max(order_by = Sensitivity, n = 1, with_ties = TRUE) %>%
  slice_max(order_by = Specificity, n = 1, with_ties = FALSE) %>%  ungroup() %>%
  ungroup()

write.csv(optimal_clinical_90, paste0(output_dir, "5_5_optimal_models_clinical_90.csv"), row.names = FALSE)

cat("\nAll configurations meeting Specificity >= 95%:\n")
print(all_clinical_95)
cat("\nAll configurations meeting Specificity >= 90%:\n")
print(all_clinical_90)

# ==============================================================================
# SECTION 5.5: Pareto Frontier of PJQM2 under the 95% Specificity Constraint
# ==============================================================================

pjqm2_feasible_95 <- df %>%
  filter(model == "PJQM2", !is.na(Sensitivity), !is.na(Specificity), Specificity >= 0.95) %>%
  mutate(
    sampling_freq_label = factor(sampling_freq,
                                 levels = c(0.0833333333333333, 0.166666666666667, 0.25, 0.5, 1),
                                 labels = c("Monthly", "Every 2 months", "Quarterly", "Biannual", "Annual")),
    alpha_label = factor(alpha, levels = c(0.01, 0.025, 0.05, 0.1))
  )

p7 <- ggplot(pjqm2_feasible_95, aes(x = Specificity, y = Sensitivity, 
                                    color = alpha_label, shape = sampling_freq_label)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_vline(xintercept = 0.95, color = "red", linetype = "dotted", linewidth = 1) +
  scale_color_brewer(palette = "Dark2", name = expression(alpha ~ "Level")) +
  scale_shape_manual(values = c(15, 16, 17, 18, 8), name = "Sampling Frequency") +
  scale_x_continuous(labels = tex_percent) +
  scale_y_continuous(labels = tex_percent) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(paste0(output_dir, "plot_5_5_pjqm2_pareto_95.png"), plot = p7, width = 8.5, height = 6)

tikz(paste0(output_dir, "plot_5_5_pjqm2_pareto_95.tex"), width = 8.5, height = 6)
print(p7)
dev.off()

cat("\nNumber of PJQM2 configurations meeting Specificity >= 95%:", nrow(pjqm2_feasible_95), "\n")

cat("\n========================================================================\n")
cat("SUCCESS: All plots (.png and .tex) and summary tables generated cleanly!\n")
cat("========================================================================\n")