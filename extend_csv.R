library(dplyr)

metrics_df <- read.csv("C:/Users/Rieke/Documents/thesis_experiments/results_metrics_alldata.csv", stringsAsFactors = FALSE)

metrics_extended <- metrics_df %>%
  mutate(
    TNR_Specificity = ifelse((Total_TN + Total_FP) > 0, 
                             Total_TN / (Total_TN + Total_FP), 
                             NA),
    NPV = ifelse((Total_TN + Total_FN) > 0, 
                 Total_TN / (Total_TN + Total_FN), 
                 NA)
  )

write.csv(metrics_extended, "C:/Users/Rieke/Documents/thesis_experiments/results_metrics_extended.csv", row.names = FALSE)

print("Successful. New data:")
head(metrics_extended %>% select(model, Total_TN, Total_FP, Total_FN, TNR_Specificity, NPV))