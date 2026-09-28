# Imports
library(dplyr)
library(stringr)
library(haven)

#' Retrieves full NHANES data, with all necessary data
get_nhanes_cohort <- function() {
  path = "C:/Users/Rieke/Documents/thesis_experiments/data"
  
  rx_data     <- read_xpt(file.path(path, "RXQ_RX_H.xpt")) 
  snfl_data   <- read_xpt(file.path(path, "SSSNFL_H.xpt")) 
  demo_data   <- read_xpt(file.path(path, "DEMO_H.xpt")) 
  mcq_data    <- read_xpt(file.path(path, "MCQ_H.xpt")) 
  pfq_data    <- read_xpt(file.path(path, "PFQ_H.xpt")) 
  ghb_data    <- read_xpt(file.path(path, "GHB_H.xpt")) 
  biopro_data <- read_xpt(file.path(path, "BIOPRO_H.xpt")) 
  bmx_data    <- read_xpt(file.path(path, "BMX_H.xpt")) 
  
  # define MS patients
  ms_drugs <- c("GLATIRAMER", "INTERFERON BETA", "FINGOLIMOD", 
                "DIMETHYL FUMARATE", "NATALIZUMAB", "OCRELIZUMAB", 
                "TERIFLUNOMIDE", "ALEMTUZUMAB")

  ms_medication_df <- rx_data %>%
    filter(str_detect(toupper(RXDDRUG), paste(ms_drugs, collapse = "|"))) %>%
    select(SEQN, RXDDRUG) %>%
    group_by(SEQN) %>%
    summarise(MS_Medication = paste(unique(toupper(RXDDRUG)), collapse = " + ")) %>%
    ungroup()
  
  ms_seqns <- ms_medication_df %>% pull(SEQN) %>% unique() 
  
  # define control group using exclusion criteria
  icd_exclude <- "I63|G45\\.9|G40|G30\\.9|G35|G31\\.84|G20|G31\\.9|G81\\.1|G50\\.0|G47\\.41"
  
  rx_exclude_seqns <- rx_data %>%
    filter(str_detect(RXDRSD1, icd_exclude) | str_detect(RXDRSD2, icd_exclude) | str_detect(RXDRSD3, icd_exclude)) %>%
    pull(SEQN) %>% unique() 
  
  mcq_exclude_seqns <- mcq_data %>% filter(MCQ160F == 1) %>% pull(SEQN) %>% unique()
  
  pfq_exclude_seqns <- pfq_data %>%
    filter(PFQ061B %in% c(3, 4) | (PFQ054 == 1 & PFQ057 %in% c(1, 2, 3, 4))) %>%
    pull(SEQN) %>% unique()
  
  preg_exclude_seqns <- demo_data %>% filter(RIDEXPRG %in% 1) %>% pull(SEQN) %>% unique()
  
  pool_0 <- snfl_data %>% filter(!is.na(SSSNFL) & !(SEQN %in% ms_seqns)) %>% pull(SEQN)
  pool_1 <- setdiff(pool_0, preg_exclude_seqns)
  pool_2 <- setdiff(pool_1, unique(c(rx_exclude_seqns, mcq_exclude_seqns))) 
  pool_3 <- setdiff(pool_2, pfq_exclude_seqns) 
  outliers_dropped <- snfl_data %>% filter(SEQN %in% pool_3 & SSSNFL > 66) %>% pull(SEQN)
  pool_4 <- setdiff(pool_3, outliers_dropped)
  
  # Final cohort
  final_cohort <- snfl_data %>% 
    select(SEQN, SSSNFL) %>%      
    filter(!is.na(SSSNFL)) %>%    
    filter(SEQN %in% ms_seqns | SEQN %in% pool_4) %>% 
    
    left_join(demo_data %>% select(SEQN, RIDAGEYR, RIAGENDR), by = "SEQN") %>%  
    left_join(ghb_data %>% select(SEQN, LBXGH), by = "SEQN") %>%   
    left_join(biopro_data %>% select(SEQN, LBXSCR), by = "SEQN") %>%
    left_join(bmx_data %>% select(SEQN, BMXBMI), by = "SEQN") %>% 
    left_join(ms_medication_df, by = "SEQN") %>%
    
    mutate(
      Status = ifelse(SEQN %in% ms_seqns, "MS Patient", "Control"),
      Gender = case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      MS_Medication = ifelse(is.na(MS_Medication), "None", MS_Medication)
    ) %>%
    
    select(-RIAGENDR) %>% 
    
    filter(!is.na(RIDAGEYR) & !is.na(BMXBMI) & !is.na(LBXGH) & !is.na(LBXSCR))
  
  cat("Successfully loaded cohort. Total MS:", sum(final_cohort$Status == "MS Patient"), 
      "| Total Control:", sum(final_cohort$Status == "Control"), "\n")
  
  return(final_cohort)
}

#' estimate variance of control group
#' @param cohort_data cleaned cohort NHANES dataset
estimate_control_variance <- function(cohort_data) {
  beta_bmi   <- -0.020 
  beta_hba1c <- 0.090  
  beta_crea  <- 0.467  
  
  healthy_data <- cohort_data %>%
    filter(Status == "Control") %>%
    mutate(
      log_sNfL_measured = log(SSSNFL),
      
      base_value = log_sNfL_measured - (
        calc_age_effect(RIDAGEYR) + 
          (beta_bmi * BMXBMI) +
          (beta_hba1c * LBXGH) +
          (beta_crea * LBXSCR)
      )
    )
  
  true_variance <- var(healthy_data$base_value, na.rm = TRUE)
  
  return(true_variance)
}