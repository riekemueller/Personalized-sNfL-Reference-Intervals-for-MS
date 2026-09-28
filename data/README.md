## NHANES Dataset

Data is downloaded from NHANES 2013-2024 cohort: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2013

The following datasets are relevant:
- Demographics Data: DEMO_H.xpt contains age and gender
- Laboratory Data: 
     - SSSNFL_H.xpt contains Serum Neurofilmanet Lightchain in pg/mL in variable SSSNFL
     - BIOPRO_H.xpt contains Serum Creatinine in mg/dL in variable LBXSCR
     - GHB_H.xpt contains Glycohemoglobin (HbA1c) in % in variable LBXGH
- Examination Data:
     - BMX_H.xpt containt Body Mass Index in kg/m^2 in variable BMXBMI
- Questionnaire Data:
     - RXG_PX contains prescription medications in variables RXDDRUG (drug name) and corresponding ICD-code in variable RXDRSC1, RXDRSC2 and RXDRSC3 to detect MS patients
     - MCQ_H.xpt contains medical conditions
