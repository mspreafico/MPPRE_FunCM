################################################################
# Step 4: Predictive functional Cox model for overall survival #
################################################################
rm( list = ls() )
library(data.table)
library(survival)

# Set directory
# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/Code" ## change to your working directory
setwd(directory_path)

# Load data
load('results/fpca_compensators.Rdata')
load('data/fake_survival.Rdata')

# K-FOLD CROSS VALIDATION #
K = 10 
# Create K equally size folds
set.seed(1234)
patients = sample(unique(survival_df$id))
folds <- cut(seq(1,length(patients)), breaks=K,labels=FALSE)
# Baseline covariates list
baseline_cov <- list(c('id', 'death_time', 'death_status'),
                     c('id', 'death_time', 'death_status', 'gender', 'age'),
                     c('id', 'death_time', 'death_status','age'),
                     c('id', 'death_time', 'death_status', 'gender', 'age'))
# Functional predictors: Consider scores such that cumPVE>=99%
K99_ACE <- which(fpca_ACE$cumPVE<0.99)+1
K99_BB <- which(fpca_BB$cumPVE<0.99)+1
K99_AA <- which(fpca_AA$cumPVE<0.99)+1
K99_HF <- which(fpca_HF$cumPVE<0.99)+1

# set progress bar
pb <- txtProgressBar(min = 0, max = K99_ACE*K99_BB*K99_AA*K99_HF*length(baseline_cov), style = 3)
pb_index = 0

best_CI = 0
best_components <- list('baseline' = NA,
                        'K_ACE' = NA,
                        'K_BB' = NA,
                        'K_AA' = NA,
                        'K_HF' = NA)

for(i in 1:length(baseline_cov)){
  cols = as.vector(baseline_cov[[i]])
  baseline = survival_df[, ..cols]
  for(Kace in 1:K99_ACE){
    data_cox_1 = merge(baseline, fpca_ACE$PCscores[,1:(Kace+1)],by='id')
    for(Kbb in 1:K99_BB){
      data_cox_12 = merge(data_cox_1, fpca_BB$PCscores[,1:(Kbb+1)],by='id')
      for(Kaa in 1:K99_AA){
        data_cox_123 = merge(data_cox_12, fpca_AA$PCscores[,1:(Kaa+1)],by='id')
        for(Khf in 1:K99_HF){
          pb_index = pb_index+1
          setTxtProgressBar(pb, pb_index)
          data_cox_1234 = merge(data_cox_123, fpca_HF$PCscores[,1:(Khf+1)],by='id')
          
          # Perform 10 fold cross validation
          fold_score = NULL
          for(i in 1:K){
            validIndexes <- which(folds==i,arr.ind=TRUE)
            valid_patients <- patients[validIndexes]
            train_patients <- patients[-validIndexes]
            
            train = data_cox_1234[id %in% train_patients]
            valid = data_cox_1234[id %in% valid_patients]
            model = coxph(Surv(death_time,death_status) ~ ., data=train[,-1])
            prediction = predict(model, valid)
            score = survConcordance(Surv(death_time, death_status) ~ prediction, data = valid)$concordance[[1]]
            fold_score = c(fold_score, score)
          }
          if(mean(fold_score)>best_CI){
            best_CI <- mean(fold_score)
            best_components <- list('baseline' = cols,
                                    'K_ACE' = Kace,
                                    'K_BB' = Kbb,
                                    'K_AA' = Kaa,
                                    'K_HF' = Khf)
          }
        }
      }
    }
  }
}

# Best set of scores
best_components
# Best Concordance Index
best_CI

#--------------------------------------------------------------------#
# Final Multivariate Functional Linear Cox Regression Model (MFLCRM) #
#--------------------------------------------------------------------#
cols = as.vector(best_components$baseline)
data_mflcrm = merge(survival_df[,..cols], fpca_ACE$PCscores[,1:(best_components$K_ACE+1)], by='id')
data_mflcrm = merge(data_mflcrm, fpca_BB$PCscores[,1:(best_components$K_BB+1)], by='id')
data_mflcrm = merge(data_mflcrm, fpca_AA$PCscores[,1:(best_components$K_AA+1)], by='id')
data_mflcrm = merge(data_mflcrm, fpca_HF$PCscores[,1:(best_components$K_HF+1)], by='id')

final_mflcrm <- coxph(Surv(death_time,death_status) ~ ., data=data_mflcrm[,-1])

#---------#
# Table 2 #
#---------#
summary(final_mflcrm)
tableHR <- function(model, digits=3){
  table = cbind('HR'=round(exp(coef(model)), digits),
                '2.5%CI'=round(exp(confint(model)[,1]), digits),
                '97.5%CI'=round(exp(confint(model)[,2]), digits) )
  return(table)
}
tableHR(final_mflcrm,digits=4)




