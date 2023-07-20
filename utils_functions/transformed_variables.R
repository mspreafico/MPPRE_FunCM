
# Auxiliary functions for file 02_modelling_compensators.R

library(data.table)

# Function to build features corresponding to the transformed variables 
# enum:marks, log(enum+1), log(marks+1), log(enum+1):marks, 
# enum:log(marks+1) and log(enum+1):log(marks+1).
transformed_variables <- function(dataset){
  sel_df <- data.table(dataset)
  
  sel_df[, enum_marks := enum*marks]
  sel_df[, logp1_enum := log(enum+1)]
  sel_df[, logp1_marks := log(marks+1)]
  sel_df[, logp1_enum_marks := log(enum+1)*marks]
  sel_df[, enum_logp1_marks := enum*log(marks+1)]
  sel_df[, logp1_enum_logp1_marks := log(enum+1)*log(marks+1)]
  sel_df[, gender_bin:=ifelse(gender=='M',1,0)]
  
  return(sel_df)
}

