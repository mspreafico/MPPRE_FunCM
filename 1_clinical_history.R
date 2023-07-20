#################################################
# Appendix A.2 - Clinical history data (Step 1) #
#################################################
rm(list = ls())
library(data.table)

# Set directory
# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/Code" ## change to your working directory
setwd(directory_path)

# Load datasets
load('data/fake_recurrent.RData')

# Warning: if you use your own data, at the end of "Data preprocessing and 
# clinical history data" selection (Step 1) data must be reformatted in the 
# form explained in Appendix A.2 and shown in Table 3.

#---------#
# Table 3 #
#---------#
ACE_df[id==500]




