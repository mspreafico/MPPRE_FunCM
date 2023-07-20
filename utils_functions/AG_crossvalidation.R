
# Auxiliary functions for file 02_modelling_compensators.R
# Step 2.1: Features selection and coefficients estimation

library(data.table)
library(survival)

# Kfold_CV_MAMR function 
#   - fits different AG models on trainingData;
#   - computes the Mean Absolute Martingale (MAMR) residual (see Appendix A.3)
#     for the validation data (k-th fold);
#   - returns the formula related to the AG model with minimum MAMR.

Kfold_CV_MAMR <- function(sel_df, K=10){
  
  set.seed(1234)
  # Randomly shuffle the patients
  patients = sample(unique(sel_df$id))
  # Create K equally size folds
  folds <- cut(seq(1,length(patients)),breaks=K,labels=FALSE)
  
  # Perform K-fold cross validation
  scores = NULL
  pb <- txtProgressBar(min = 0, max = K, style = 3)
  for(i in 1:K){
    setTxtProgressBar(pb, i)
    # Segement patients by fold 
    validIndexes <- which(folds==i,arr.ind=TRUE)
    valid_patients <- patients[validIndexes]
    train_patients <- patients[-validIndexes]
    # Split train, valid
    train = sel_df[id %in% train_patients]
    valid = sel_df[id %in% valid_patients]
    
    # Fit AG models on train data
    models = list(
      # AGE
      coxph(Surv(start,stop,status) ~ age + enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + marks + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + marks + enum_marks, cluster=id, data = train),
      
      # log(enum+1)
      coxph(Surv(start,stop,status) ~ age + logp1_enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + marks + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + marks + logp1_enum_marks, cluster=id, data = train),
      
      # log(marks+1)
      coxph(Surv(start,stop,status) ~ age + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_marks + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + enum + logp1_marks + enum_logp1_marks, cluster=id, data = train),
      
      # log(enum+1), log(marks+1)
      coxph(Surv(start,stop,status) ~ age + logp1_enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_marks + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ age + logp1_enum + logp1_marks + logp1_enum_logp1_marks, cluster=id, data = train),
      
      # GENDER
      coxph(Surv(start,stop,status) ~ gender_bin + enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + marks + enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + marks + enum_marks, cluster=id, data = train),
      
      # log(enum+1)
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + marks + logp1_enum_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + marks + logp1_enum_marks, cluster=id, data = train),
      
      # log(marks+1)
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_marks + enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + enum + logp1_marks + enum_logp1_marks, cluster=id, data = train),
      
      # log(enum+1), log(marks+1)
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_marks + logp1_enum_logp1_marks, cluster=id, data = train),
      coxph(Surv(start,stop,status) ~ gender_bin + logp1_enum + logp1_marks + logp1_enum_logp1_marks, cluster=id, data = train),
      
      # AGE + GENDER
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + marks + enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + marks + enum_marks, cluster= id, data = train),
      
      # log(enum+1)
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + logp1_enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + marks + logp1_enum_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + marks + logp1_enum_marks, cluster= id, data = train),
      
      # log(marks+1)
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_marks + enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + enum + logp1_marks + enum_logp1_marks, cluster= id, data = train),
      
      # log(enum+1), log(marks+1)
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + logp1_enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_marks + logp1_enum_logp1_marks, cluster= id, data = train),
      coxph(Surv(start,stop,status) ~ age + gender_bin + logp1_enum + logp1_marks + logp1_enum_logp1_marks, cluster= id, data = train)
      
    )
    # Evaluate martingale residuals on validation events
    fold_scores = c()
    for(model in models){
      expected_N <- predict(model, valid, type='expected')
      score <- mean(abs(valid$status - expected_N)) # mean absolute residual
      fold_scores = c(fold_scores,score)
    }
    scores = rbind(scores,fold_scores)
  }
  
  # summarize CV scores
  mean_scores = colMeans(scores)
  best_score = min(mean_scores)
  worst_score = max(mean_scores)
  
  # print cv scores
  print('****************** Cross validation Mean Absolute Martingale residual ******************')
  cat('\n')
  print('|------------|-----------|--------------------------------------|')
  print('| K-fold CV  |   MAMR    |               FORMULA                |')
  print('|------------|-----------|--------------------------------------|')
  for(i in 1:length(models)){
    score = mean_scores[i]
    if(score != best_score & score != worst_score){
      print(paste('|            |  ',round(mean_scores[i],3),'  | ',models[[i]]$formula[3]))
    }else if (score == best_score){
      print(paste('|   BEST --> |  ',round(mean_scores[i],3),'  | ',models[[i]]$formula[3]))
    }else{
      print(paste('|  WORST --> |  ',round(mean_scores[i],3),'  | ',models[[i]]$formula[3]))
    }
  }
  
  # Return formula related to the best AG model (mean MAMR)
  best_model <- models[[which(mean_scores==min(mean_scores))]]$formula
  return(best_model)
}