## NHANES Dataset

Data has to be downloaded from NHANES 2013-2024 cohort: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2013

The following datasets are relevant:
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
