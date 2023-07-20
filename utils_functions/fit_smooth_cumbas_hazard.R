
# Auxiliary functions for file 02_modelling_compensators.R
# Step 2.2 Fit and smooth cumulative baseline hazard

library(cobs)
library(survival)

# Fit and smooth cumulative baseline hazard; return a list
fit_smooth_Lambda0 <- function(model){
  
  # Get estimated cumulative baseline hazard function
  bh = basehaz(model, centered = FALSE)
  t <- bh$time
  Lambda0 <- bh$hazard
  
  # Smooth version of Lambda0
  Lambda0S <- cobs(c(0,t), c(0,Lambda0), constraint=c("increase"), 
                   pointwise=matrix(c(0,0,0),nrow=1), nknots=20, lambda=0, toler.kn=0)
  
  return(list('times0' = t,
              'Lambda0' = Lambda0,
              'Lambda0S' = Lambda0S)
  )
}

# Function to plot Figure 4
plot_Lambda0 <- function(list_L0, h_type='h', title='Event h'){
  t <- list_L0$time
  Lambda0 <- list_L0$Lambda0
  Lambda0S <- list_L0$Lambda0S
  
  plot(Lambda0S$x, Lambda0S$fitted, type="l",
       main=title, ylab="Baseline cumulative hazard", xlab="Time [days]", col = 'red',
       cex.axis=1.5, cex.lab=1.5, cex.main=1.75)
  points(t,Lambda0,type="s",lty=2, col = 'blue')
  legend('bottomright',legend=c(TeX(paste('$\\tilde{\\Lambda}_0^{(',h_type,')}')),
                                TeX(paste('$\\hat{\\Lambda}_0^{(',h_type,')}'))),
         col=c("red", "blue"), lty=1:2, cex=1.5,y.intersp=2,box.lty=0)
  
}
