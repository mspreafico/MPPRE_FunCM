# MPPRE_FunCM
Code for (i) modelling of Marked-Point Process compensators related to Recurrent Events (MPPRE) of interest with cross-validation; (ii) fitting a Functional linear Cox Model (FunCM) with the obtained compensators as multivariate functional covariates (model selection by cross-validation).
[This code is also freely available as Supporting Information in the main reference paper.]

### Reference
Spreafico, M. & Ieva, F. (2021). Functional modeling of recurrent events on time-to-event processes. *Biometrical Journal*, **63**(5):948-967. https://onlinelibrary.wiley.com/doi/10.1002/bimj.202000374


### Data Availability
We cannot provide the original administrative data due to confidentiality.
We provide some fake datasets in order to allow researchers who want to replicate the same analysis to properly get how the code has to be run and how results are displayed and should be read.

## Description

- Files:
  - **1_clinical_history.R**: Code for Table 3 (Appendix A.2 - Results of "Step 1: Data preprocessing and clinical history")
  - **2_modelling_compensators.R**: Code for "Step 2: Modelling compensators of marked point processes" (Table 1 included). It generates results 'best_AGmodels.Rdata', 'baseline_haz.R' and 'compensators.R' that will be saved in folder '/results'.
  - **3_fpca_compensators**: Code for "Step 3: Summarize compensators through Functional Principal Component Analysis". It generates results 'fpca_compensators.R' that will be saved in folder /results.
  - **4_functional_Cox.R**: Code for "Step 4: Predictive functional Cox model for overall survival".
  - **figures.R**: Code for Figures 3-4-5-6-7-8.
      
- Sub-folder **./utils_functions/** contains some auxiliary functions to run the code.
  
- Sub-folder **./data/** contains fake datasets along with their legends: 
  - **fake_recurrent.Rdata**: It contains four datasets (once for each process h = {ACE, BB, AA, HF hosp}), namely 'ACE_df', 'BB_df', 'AA_df' and 'HF_df', related to the clinical history of the events of interest for 500 fake patients. Datasets have been reformatted as explained in Appendix A.2. 
	- **recurrent_legend.txt**: Variables legend of dataset 'fake_recurrent.Rdata'.
	- **fake_survival.Rdata**: Survival data 'survival_df' related to the 500 fake patients.
	- **survival_legend.txt**: Variables legend of datasets 'fake_survival.Rdata'
 
 - Sub-folder **./results** is an empty folder where intermediate results will be saved.

## Software
- R software.
- File **sessionInfo.txt** contains a list of code configurations [software version (incl. package versions), platform].

(Last update: July 20th, 2023)
