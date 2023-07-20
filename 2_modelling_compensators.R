#############################################################
# Step 2: Modelling compensators of marked point processes #
############################################################
rm(list = ls())

# Load packages
library(data.table)
library(survival)

# Set directory
# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/Code" ## change to your working directory
setwd(directory_path)

# Load datasets
load('data/fake_recurrent.RData')

# In order to use some function in the following cells, we build features 
# corresponding to the transformed variables enum:marks, log(enum+1), log(marks+1), 
# log(enum+1):marks, enum:log(marks+1) and log(enum+1):log(marks+1).
source('utils_functions/transformed_variables.R')

ACE_df <- transformed_variables(ACE_df)
BB_df <- transformed_variables(BB_df)
AA_df <- transformed_variables(AA_df)
HF_df <- transformed_variables(HF_df)


#----------------------------------------------------#
# 2.1 Features selection and coefficients estimation #
#----------- ----------------------------------------#
source('utils_functions/AG_crossvalidation.R')

# AG-models selection trough 10-fold cross-validation minimizing MAMR
best_ACE <- Kfold_CV_MAMR(ACE_df)
best_BB <- Kfold_CV_MAMR(BB_df)
best_AA <- Kfold_CV_MAMR(AA_df)
best_HF <- Kfold_CV_MAMR(HF_df)

# Selected AG-models for recurrent events
model_ACE = coxph(best_ACE, cluster=id, data = ACE_df)
model_BB = coxph(best_BB, cluster=id, data = BB_df)
model_AA = coxph(best_AA, cluster=id, data = AA_df)
model_HF = coxph(best_HF, cluster=id, data = HF_df)

# Table 1
tableHR <- function(model, digits=3){
  table = cbind('HR'=round(exp(coef(model)), digits),
                '2.5%CI'=round(exp(confint(model)[,1]), digits),
                '97.5%CI'=round(exp(confint(model)[,2]), digits) )
  return(table)
}
tableHR(model_ACE)
tableHR(model_BB)
tableHR(model_AA)
tableHR(model_HF)

save(list=c('model_ACE', 'model_BB', 'model_AA', 'model_HF'), file='results/best_AGmodels.Rdata')


#-----------------------------------------------#
# 2.2 Fit and smooth cumulative baseline hazard #
#-----------------------------------------------#
#load('results/best_AGmodels.Rdata')
source('utils_functions/fit_smooth_cumbas_hazard.R')

ACE_bashaz <- fit_smooth_Lambda0(model_ACE)
BB_bashaz <- fit_smooth_Lambda0(model_BB)
AA_bashaz <- fit_smooth_Lambda0(model_AA)
HF_bashaz <- fit_smooth_Lambda0(model_HF)

save(list=c('ACE_bashaz', 'BB_bashaz', 'AA_bashaz', 'HF_bashaz'), file='results/baseline_haz.Rdata')


#------------------------------#
# 2.3 Reconstruct compensators #
#------------------------------#
#load('results/baseline_haz.Rdata')
source('utils_functions/compute_cumulative_hazards.R')

t_max <- 365
times <- seq(0,t_max,by=1)
ACE_cumhaz = compute_cumulative_hazard(model = model_ACE,
                                       sel_df = ACE_df,
                                       smoothed_baseline = ACE_bashaz$Lambda0S,
                                       times, verbose = TRUE)
BB_cumhaz = compute_cumulative_hazard(model = model_BB,
                                      sel_df = BB_df,
                                      smoothed_baseline = BB_bashaz$Lambda0S,
                                      times, verbose = TRUE)
AA_cumhaz = compute_cumulative_hazard(model = model_AA,
                                      sel_df = AA_df,
                                      smoothed_baseline = AA_bashaz$Lambda0S,
                                      times, verbose = TRUE)
HF_cumhaz = compute_cumulative_hazard(model = model_HF,
                                      sel_df = HF_df,
                                      smoothed_baseline = HF_bashaz$Lambda0S,
                                      times, verbose = TRUE)

save(list=c('ACE_cumhaz', 'BB_cumhaz', 'AA_cumhaz', 'HF_cumhaz'), file='results/compensators.Rdata')

