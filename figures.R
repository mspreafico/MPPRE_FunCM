graphics.off()
rm(list = ls())

# Set directory
# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/Code" ## change to your working directory
setwd(directory_path)

#----------#
# Figure 3 #
#----------#
load('data/fake_recurrent.RData')
source('utils_functions/utils_plots.R')

ACE_Ncount = count_proc_long(ACE_df)
pACE = plot_Ncount(ACE_Ncount, h_type='ACE', title = 'Purchase of ACE')

BB_Ncount = count_proc_long(BB_df)
pBB = plot_Ncount(BB_Ncount, h_type='BB', title = 'Purchase of BB')

AA_Ncount = count_proc_long(AA_df)
pAA = plot_Ncount(AA_Ncount, h_type='AA', title = 'Purchase of AA')

HF_Ncount = count_proc_long(HF_df)
pHF = plot_Ncount(HF_Ncount, h_type='HF hosp', title = 'HF hospitalization')

x11()
ggarrange(pACE, pBB, pAA, pHF, ncol=2, nrow=2)


#----------#
# Figure 4 #
#----------#
load('results/baseline_haz.Rdata')
source('utils_functions/fit_smooth_cumbas_hazard.R')

x11()
par(mfrow=c(2,2))
plot_Lambda0(ACE_bashaz, h_type='ACE', title = 'Purchase of ACE')
plot_Lambda0(BB_bashaz, h_type='BB', title = 'Purchase of BB')
plot_Lambda0(AA_bashaz, h_type='AA', title = 'Purchase of AA')
plot_Lambda0(HF_bashaz, h_type='HF hosp', title = 'HF hospitalization')

#----------#
# Figure 5 #
#----------#
load('results/compensators.Rdata')
source('utils_functions/utils_plots.R')

pACE = plot_cumhaz(ACE_cumhaz, h_type='ACE', title = 'Purchase of ACE')
pBB = plot_cumhaz(BB_cumhaz, h_type='BB', title = 'Purchase of BB')
pAA = plot_cumhaz(AA_cumhaz, h_type='AA', title = 'Purchase of AA')
pHF = plot_cumhaz(HF_cumhaz, h_type='HF hosp', title = 'HF hospitalization')

x11()
ggarrange(pACE, pBB, pAA, pHF, ncol=2, nrow=2)


#----------#
# Figure 6 #
#----------#
load('data/fake_recurrent.Rdata')
load('results/compensators.Rdata')
source('utils_functions/utils_plots.R')

ACE_resid = residuals_long(ACE_df, ACE_cumhaz)
pACE = plot_resid(ACE_resid, h_type='ACE', title = 'Purchase of ACE')

BB_resid = residuals_long(BB_df, BB_cumhaz)
pBB = plot_resid(BB_resid, h_type='BB', title = 'Purchase of BB')

AA_resid = residuals_long(AA_df, AA_cumhaz)
pAA = plot_resid(AA_resid, h_type='AA', title = 'Purchase of AA')

HF_resid = residuals_long(HF_df, HF_cumhaz)
pHF = plot_resid(HF_resid, h_type='HF hosp', title = 'HF hospitalization')

x11()
ggarrange(pACE, pBB, pAA, pHF, ncol=2, nrow=2)


## FIGURES 7 and 8
load('results/fpca_compensators.Rdata')
source('utils_functions/fpca_prcomp.R')

# Range for FPCA plots
ylim <- range(c(fpca_ACE$efuncs[,1],
                fpca_BB$efuncs[,1],
                fpca_AA$efuncs[,1],
                fpca_HF$efuncs[,1],
                fpca_ACE$efuncs[,2],
                fpca_BB$efuncs[,2],
                fpca_AA$efuncs[,2],
                fpca_HF$efuncs[,2]))

#----------#
# Figure 7 #
#----------#
cost=1
k=1

x11()
par(mfcol=c(2,4), mar=c(3,2.5,1.5,0.1), mgp=c(1.5,0.5,0))
component_plot(fpca_ACE, k=k, h_type='ACE', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="First Components", title1='Purchase of ACE',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)

component_plot(fpca_BB, k=k, h_type='BB', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="", title1='Purchase of BB',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)
legend(-150,8.75, 
       legend=c(TeX(paste('$\\mu^{(h)}(t)$',' with ','$h \\in H$')), 
                TeX(paste('$\\mu^{(h)}(t) + \\sqrt{\\nu^{(h)}_',k,'} \\xi^{(h)}_',k,'(t)$  ',sep = '')),
                TeX(paste('$\\mu^{(h)}(t) - \\sqrt{\\nu^{(h)}_',k,'} \\xi^{(h)}_',k,'(t)$',sep = ''))),
       lty=c(1,1,1), lwd = c(1.5,1,1), pch = c(NA,'+','-'), col = c('black','red','blue'),
       cex=1.5,box.lty=0, xpd=NA, horiz=T, x.intersp=0.2
)
component_plot(fpca_AA, k=k, h_type='AA', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="", title1='Purchase of AA',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)
component_plot(fpca_HF, k=k, h_type='HF hosp', cost=cost,  ylim2=c(0,2), ylim=ylim,
               xlab1="", ylab1="", title1='HF hospitalization',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)

#----------#
# Figure 8 #
#----------#
cost=3
k=2

x11()
par(mfcol=c(2,4), mar=c(3,2.5,1.5,0.1), mgp=c(1.5,0.5,0))
component_plot(fpca_ACE, k=k, h_type='ACE', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="Second Components", title1='Purchase of ACE',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)
component_plot(fpca_BB, k=k, h_type='BB', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="", title1='Purchase of BB',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)
legend(-150,8.75, 
       legend=c(TeX(paste('$\\mu^{(h)}(t)$',' with ','$h \\in H$')), 
                TeX(paste('$\\mu^{(h)}(t) + ',cost,'\\sqrt{\\nu^{(h)}_',k,'} \\xi^{(h)}_',k,'(t)$  ',sep = '')),
                TeX(paste('$\\mu^{(h)}(t) - ',cost,'\\sqrt{\\nu^{(h)}_',k,'} \\xi^{(h)}_',k,'(t)$',sep = ''))),
       lty=c(1,1,1), lwd = c(1.5,1,1), pch = c(NA,'+','-'), col = c('black','red','blue'),
       cex=1.5,box.lty=0, xpd=NA, horiz=T, x.intersp=0.2, y.intersp=0
)
component_plot(fpca_AA, k=k, h_type='AA', cost=cost,  ylim2=c(0,7), ylim=ylim,
               xlab1="", ylab1="", title1='Purchase of AA',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)
component_plot(fpca_HF, k=k, h_type='HF hosp', cost=cost,  ylim2=c(0,2), ylim=ylim,
               xlab1="", ylab1="", title1='HF hospitalization',
               xlab2="Time [days]", ylab2="", title2=NULL, print.legend = F)



