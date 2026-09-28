# Personalized sNfL Reference Intervals for Multiple Sclerosis: A Covariate-Adjusted Longitudinal Approach
This project is under the MIT license.
Author: Rieke Müller, Goethe University, Germany, 2026

## 0. Prerequisites
Install the following packages:
- here
- dplyr
- tidyr
- readxl
- stringr
- haven
- R6
- splines
- ggplot2
- tikzDevice

## 1. NHANES Dataset (Folder /data)
Data is downloaded from NHANES 2013-2024 cohort: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2013

The following datasets are relevant and should be stored in the data folder:
- Demographics Data: 
     - DEMO_H.xpt contains age and gender
- Examination Data:
     - BMX_H.xpt contains Body Mass Index in kg/m^2 in variable BMXBMI
- Laboratory Data: 
     - SSSNFL_H.xpt contains Serum Neurofilament Light Chain in pg/mL in variable SSSNFL
     - BIOPRO_H.xpt contains Serum Creatinine in mg/dL in variable LBXSCR
     - GHB_H.xpt contains Glycohemoglobin (HbA1c) in % in variable LBXGH
- Questionnaire Data:
     - RXQ_RX_H.xpt contains prescription medications in variables RXDDRUG (drug name) and corresponding ICD-code in variable RXDRSC1, RXDRSC2 and RXDRSC3 to detect MS patients
     - MCQ_H.xpt contains medical conditions
     - PFQ_H.xpt contains physical functioning data (e.g. walking problems)


## 2. Source Code (Folder /src)
1. in 00_config.R: All experimental parameter configurations are implemented. The parameters are: 
- sim_baseline_lengths: The number of historical baseline measurements n, that are necessary to initialize the prRI boundary are varied. To determine the minimum number of measurements and to observe asymptotic convergence, n is tested with n ∈ {2, 3, 5, 7, 10, 20, 30}.
- sim_sampling_frequencies: The follow-up periods between two measurements are varied to assess diagnostic delay in relapse detection. The simulated follow-up intervals are set to tfollow-up ∈ { 1/12 , 1/6 , 0.25, 0.5, 1} years.
- sim_alpha_levels: The significance level α describes the probability of rejecting the null hypothesis when it is actually true [14]. It is set to α ∈{0.10, 0.05, 0.025, 0.01}.

2. in 01_data.R: The data from NHANES dataset is processed

3. in 02_models: Five models are implemented:
- Fixed Cut-Off model with a static 12.9 threshold
- Within-Subject Model of Coskun et al.
- Within-Person Model of Coskun et al.
- PJQM2 Model of Pusparum et al.
- GAMLSS Model of Benkert et al.: here the model is not implemented itself but its results are processed.

4. in 03_simulation.R: The longitudinal simulation based on the Data Generation Process is implemented.

5. in 04_evaluation.R: The evaluation with all parameter configurations and all models is implemented.

6. in 05_results_plots.R: Generation of result plots in terms of initialization stability, impact of monitoring frequency, influence of significance level, diagnostic accuracy and clinical utility are implemented

## 3. Results (Folder /results)
All results (in csv, png and .tex) generated of the full experiment are stored here.
The detailed results of each configuration run is in the file: results_metrics_final.csv

## 4. PJQM2 Hyperparameter (Folder /PJQM2)
PJQM2 was optimized for n=3, FPR <= 5% with a grid search for hyperparameter setting. This is implemented in 06_hyperparam_PJQM2.R with the results in the .txt file in Folder PJQM2.

## Citation:
If you use this code or build on this work, please cite:
Müller, R. (2026). *Personalized sNfL Reference Intervals for Multiple Sclerosis: A Covariate-Adjusted Longitudinal Approach* [Goethe University Frankfurt]. GitHub. https://github.com/riekemueller/Personalized-sNfL-Reference-Intervals-for-MS