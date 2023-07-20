#########################################################################################
# Step 3: Summarize compensators through Functional Principal Component Analysis (FPCA) #
#########################################################################################
rm(list = ls())

# Load packages
library(data.table)
library(survival)

# Set directory
# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/Code" ## change to your working directory
setwd(directory_path)

# Load compensators
load('results/compensators.Rdata')

# Load util functions for FPCA
source('utils_functions/fpca_prcomp.R')

fpca_ACE <- fpca_prcomp(ACE_cumhaz, name='ACE')
fpca_BB <- fpca_prcomp(BB_cumhaz, name='BB')
fpca_AA <- fpca_prcomp(AA_cumhaz, name='AA')
fpca_HF <- fpca_prcomp(HF_cumhaz, name='HF')

save(list=c('fpca_ACE', 'fpca_BB', 'fpca_AA', 'fpca_HF'), file='results/fpca_compensators.Rdata')

