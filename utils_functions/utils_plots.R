
# Auxiliary functions for file figures.R

library(data.table)
library(latex2exp)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Turn recurrent events data (in the form of Table 3 - Appendix A.2) into a long format dataset
count_proc_long <- function(rec_data){
  count_proc = rec_data %>%
    mutate(time = map2(start, stop, seq, by=0.5)) %>%
    unnest(cols = time) %>%
    select(-start, -stop) %>%
    group_by(time)
  count_proc = data.table(count_proc)
  count_proc[time!=0.5, time := floor(time)]
  count_proc = count_proc[!duplicated(count_proc[,.(id,time)], fromLast=F)]
  count_proc = count_proc[time!=0.5]
  count_proc = count_proc[, .(id,time,enum)]
  return(count_proc)
}

# Compute martingale residuals from recurrent events data and compensators; return a long format data
residuals_long <- function(rec_data, cumulative_hazard){
  data_Nproc = count_proc_long(rec_data)
  residuals <- merge(data_Nproc, cumulative_hazard, by=c('id','time'))
  residuals[, M_res := cumhaz-enum]
  return(residuals)
}

# Function to plot Figure 3
plot_Ncount <- function(rec_data, h_type='h', title='Event h'){
  p<-ggplot(rec_data, aes(x= time, y=enum, group = factor(id), color=factor(id))) +
    geom_step() +
    xlab('Time [days]') + 
    ylab(TeX(paste('Counting proc $\\,$ $N_i^{($',h_type,')}(t)$'))) + 
    ggtitle(title) +
    theme(legend.position="none",
          axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.5)),
          plot.title = element_text(face="bold", size=rel(1.75), hjust = 0.5) )
  return(p)
}

# Function to plot Figure 5
plot_cumhaz <- function(cumulative_hazard, h_type='h', title='Event h'){
  p<-ggplot(cumulative_hazard, aes(x= time, y=cumhaz, group = factor(id), color=factor(id))) +
    geom_line() +
    xlab('Time [days]') + 
    ylab(TeX(paste('Compensators $\\,$ $\\hat{\\Lambda}_i^{($',h_type,')}(t)$'))) + 
    ggtitle(title) +
    theme(legend.position="none",
          axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.5)),
          plot.title = element_text(face="bold", size=rel(1.75), hjust = 0.5) )
  return(p)
}

# Function to plot Figure 6
plot_resid <- function(residuals, h_type='h', title='Event h'){
  mean_residuals = residuals[,list('mean_res' = mean(M_res)), by = 'time']
  p<-ggplot(data = mean_residuals, aes(x= time, y=mean_res)) +
    geom_line() +
    geom_line(data = residuals, aes(x= time, y=M_res, group = factor(id), color=factor(id))) +
    geom_line(data = mean_residuals, aes(x= time, y=mean_res), size=1.2) +
    xlab('Time [days]') + 
    ylab(TeX(paste('Residuals $\\,$ $\\hat{M}_i^{($',h_type,')}(t)$'))) + 
    ggtitle(title) +
    theme(legend.position="none",
          axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.5)),
          plot.title = element_text(face="bold", size=rel(1.75), hjust = 0.5) )
  return(p)
}
