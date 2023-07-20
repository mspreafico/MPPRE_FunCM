
# Auxiliary functions for file 03_fpca_compensators.R
# Step 3: Summarize compensators through Functional Principal Component Analysis (FPCA)

library(tidyr)
library(latex2exp)

# Perform FPCA on compensators; return a list
fpca_prcomp <- function(cumulative_hazard, name=NULL, t_max=365){
  
  # Reformat cumulative hazard  in a matrix format
  patient_ids = unique(cumulative_hazard$id)
  np <- length(patient_ids) # number of rows: number of patients
  times <- seq(0,t_max,by=1)
  m <- length(times)  # number of columns: time
  evals_matrix = spread(cumulative_hazard, time, cumhaz)[,-1]
  rownames(evals_matrix) = as.character(patient_ids)
  
  # Get the mean curves
  mean_evals = colMeans(evals_matrix)
  
  # Perform FPCA using standard prcomp function as explained in 
  # Ramsay and Silverman (2005) Functional Data Analysis, Chaper 8, Section 8.4.1
  h <- (times[m] - times[1])/(m - 1)
  
  # FPCA: Eigenvalues and eigenfunctions
  pca <- prcomp(evals_matrix)
  efuncs <- pca$rotation*sqrt(1/h)
  evalues <- pca$sdev^2*h
  
  # Eigenfunctions
  # Sign chosen to have positive integral (Optional)
  for(i in 1:m) {
    tempfun <- efuncs[,i]
    tempsign <- sum(tempfun)
    efuncs[,i] <- ifelse(tempsign<0, -1,1) * tempfun
  }
  
  # Principal Components Scores
  Xt <- t(data.frame(evals_matrix))-mean_evals
  scores <- t(Xt) %*% efuncs*h
  if(!is.null(name)){
    colnames(scores) <- paste0(colnames(scores), '_', name)
  }
  scores <- cbind("id" = patient_ids, scores)
  rownames(scores) <- seq(1:np)
  
  # Create list
  fpca_list <- list("grid" = times,
                    "mean_hazard" = mean_evals,
                    "efuncs" = efuncs,
                    "evalues" = evalues,
                    "PVE" = evalues/sum(evalues),
                    "cumPVE" = cumsum(evalues/sum(evalues)),
                    "PCscores" = scores)
  
  return(fpca_list)
  
}


# Function to plot columns of Figures 7/8
component_plot <- function(list, k=1, h_type='h', cost,
                           xlim = NULL, ylim=NULL, ylim2=NULL,
                           xlab1=" ", ylab1=NULL, title1='Event h',
                           xlab2="time", ylab2=" ", title2=NULL, print.legend = T){
  
  efuncs = list$efuncs
  evalues = list$evalues
  mean_hazard = list$mean_hazard
  times = list$grid
  percent.var = list$percent.var
  m = length(times)
  
  if(is.null(xlim)){
    xlim <- range(times)
  }
  if(is.null(ylim)){
    ylim <- range(efuncs[,1:max(1,k)])
  }
  if(is.null(ylab1)){
    ylab1 <- paste0('Component ', k)
  }
  
  v1 <- efuncs[,k]
  plot(times, v1, type="l", xlim=xlim, ylab=ylab1, xlab=xlab1, ylim=ylim, lwd = 1.5,
       main=title1)
  text(2,-0.95, "Time [days]")
  legend(0,-0.05,
         legend =c(TeX(paste('$\\xi^{(',h_type,')}_',k,'(t)$',sep=''))), # text
         lty=1, lwd = 1.5,
         cex=1.25,y.intersp=0.85,box.lty=0
  )
  v2 <- v1 * cost * sqrt(evalues[k])
  widep <- ((1:m)%%10)==0
  plot(times, mean_hazard, xlim=xlim, type="l", ylim=ylim2, xlab=xlab2, ylab=ylab2, main=title2, lwd=1.5)
  points(times[widep],mean_hazard[widep] + v2[widep], lwd=3, pch="+", col='red')
  points(times[widep],mean_hazard[widep] - v2[widep], lwd=3, pch="-", col='blue')
  
  if(print.legend){
    if(cost==1){
      legend(0,7, # where
             legend=c(TeX(paste('$\\mu^{(',h_type,')}(t)$')), TeX(paste('$\\mu^{(',h_type,')}(t) + ','\\sqrt{\\nu^{(',h_type,')}_',k,'} \\xi^{(',h_type,')}_',k,'(t)$',sep = '')),
                      TeX(paste('$\\mu^{(',h_type,')}(t) - ','\\sqrt{\\nu^{(',h_type,')}_',k,'} \\xi^{(',h_type,')}_',k,'(t)$',sep = ''))), # text
             lty=c(1,NA,NA), lwd = c(1.5,1,1), pch = c(NA,'+','-'), col = c('black','red','blue'), # symbols
             cex=0.75,y.intersp=0.85,box.lty=0 # other options
      )
    }else{
      legend(0,7, # where
             legend=c(TeX(paste('$\\mu^{(',h_type,')}(t)$')), TeX(paste('$\\mu^{(',h_type,')}(t) + ',cost,'\\sqrt{\\nu^{(',h_type,')}_',k,'} \\xi^{(',h_type,')}_',k,'(t)$',sep = '')),
                      TeX(paste('$\\mu^{(',h_type,')}(t) - ',cost,'\\sqrt{\\nu^{(',h_type,')}_',k,'} \\xi^{(',h_type,')}_',k,'(t)$',sep = ''))), # text
             lty=c(1,NA,NA), lwd = c(1.5,1,1),pch = c(NA,'+','-'), col = c('black','red','blue'), # symbols
             cex=0.75,y.intersp=0.85,box.lty=0 # other options
      )
    }
  }
  
}