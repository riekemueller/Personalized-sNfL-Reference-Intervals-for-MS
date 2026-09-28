# Personalized sNfL Reference Intervals for Multiple Sclerosis: A Covariate-Adjusted Longitudinal Approach
This project is under the MIT license.
Author: Rieke Müller, Goethe University, Germany, 2026

## 0. Prerequisites
Install the following packages: here, dplyr, tidyr, readxl, stringr, haven, R6, splines, ggplot2, tikzDevice

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
     - sim_sampling_frequencies: The follow-up periods between two measurements are varied to assess diagnostic delay in relapse detection. The simulated follow-up intervals are set to t_follow-up ∈ {1/12, 1/6, 0.25, 0.5, 1} years.
     - sim_alpha_levels: The significance level α describes the probability of rejecting the null hypothesis when it is actually true. It is set to α ∈{0.10, 0.05, 0.025, 0.01}.
2. in 01_data.R: The data from NHANES dataset is processed
3. in 02_models: Five models are adapted to sNfL use case:
     - Fixed Cut-Off model with a static 12.9 threshold 
     - Within-Subject Model of Coskun et al. [1]
     - Within-Person Model of Coskun et al. [2, 3]
     - PJQM2 Model of Pusparum et al. [4]
     - GAMLSS Model of Benkert et al. [5]: here the model is not implemented itself but its results are processed.
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

## References
- [1] Coskun, A., Sandberg, S., Unsal, I., Cavusoglu, C., Serteser, M., Kilercik, M., & Aarsand, A. K. (2021). Personalized reference intervals in laboratory medicine: A new model based on within-subject biological variation. *Clinical Chemistry*, 67, 374–384. https://doi.org/10.1093/clinchem/hvaa233
- [2] Coskun, A., Sandberg, S., Unsal, I., Cavusoglu, C., Serteser, M., Kilercik, M., & Aarsand, A. K. (2022). Personalized reference intervals: Using estimates of within-subject or within-person biological variation requires different statistical approaches. *Clinica Chimica Acta*, 524, 201–202. https://doi.org/10.1016/j.cca.2021.10.034
- [3] Coskun, A., Sandberg, S., Unsal, I., Yavuz, F. G., Cavusoglu, C., Serteser, M., Kilercik, M., & Aarsand, A. K. (2022). Personalized reference intervals - statistical approaches and considerations. *Clinical Chemistry and Laboratory Medicine*, 60, 629–635. https://doi.org/10.1515/cclm-2021-1066
- [4] Pusparum, M., Ertaylan, G., & Thas, O. (2022). Individual reference intervals for personalised interpretation of clinical and metabolomics measurements. *Journal of Biomedical Informatics*, 131, 104111.
- [5] Benkert, P., Meier, S., Schaedelin, S., Manouchehrinia, A., Yaldizli, Ö., Maceski, A., et al. (2022). Serum neurofilament light chain for individual prognostication of disease activity in people with multiple sclerosis: a retrospective modelling and validation study. *The Lancet Neurology*, 21, 246–257.
