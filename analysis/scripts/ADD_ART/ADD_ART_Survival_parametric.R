#-------------------- parametric survival models --------------------

rm(list=ls())
# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(survival)
library(survminer)
library(ldatools)
library(frailtypack)


path_data = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART/figures"

# Load data
setwd(path_data)

dat_recurrent = read.csv("ADDART_recurrent_14May2024.csv")

dat_long_eam = read.csv("ADDART_long_eam_14May2024.csv")

dat_long_tfv = read.csv("ADDART_long_tfv_14May2024.csv")

dat_long_joint = read.csv("ADDART_long_joint_14May2024.csv")

#------------------- data wrangling -------------------
dat_long_eam$Visit = as.factor(dat_long_eam$Visit)
dat_long_eam$`log2(VL+0.001)` = log2(dat_long_eam$VL + 0.001)

dat_long_eam = dat_long_eam %>%
  mutate(wp_nonadherence_prop_30day = n_heartbeats_30days/wp_days_30,
         sra_nonadherence_prop = 1 - sra1_adherence_prop)

###--------------data wrangling ----------------

# time to first viral load > 1000
dat_long_eam$Visit = as.numeric(dat_long_eam$Visit)
dat_long_eam$VBstatus = ifelse(dat_long_eam$VL > 1000, 1, 0)

Surv = with(dat_long_eam, Surv(time = Visit, event = VBstatus))
km_fit1a = survfit(Surv ~ 1, dat = dat_long_eam)
ggsurvplot(km_fit1a, dat = dat_long_eam, conf.type = "log-log")


#time to first non-suppression event
dat_long_eam$VL_supp_status = ifelse(dat_long_eam$VL > 50, 1, 0)
Surv2 = with(dat_long_eam, Surv(time = time_yrs, event = VL_supp_status))
km_fit2a = survfit(Surv2 ~ 1, data = dat_long_eam)
ggsurvplot(km_fit2a, data = dat_long_eam,
           risk.table = T,
           cumevents = T,
           conf.type = "log-log")

#recurrent
#create visit time start for interest
dat_recurrent = dat_recurrent %>%
  group_by(PID) %>%
  mutate(visit_start = lag(Visit, 1, default = 0)) %>%
  relocate(visit_start, .before = Visit)

#create log10TFV
dat_recurrent$log10_TFV = log10(dat_recurrent$TFV)



#---------------------- the cox-PH models -----------------------------
fit_ph_frail = coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) + frailty(PID),
                     data = dat_recurrent)

summary(fit_ph_frail)
cox.zph(fit_ph_frail)
plot(cox.zph(fit_ph_frail))

gg_coxsnell(fit_ph_frail, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)


fit_ph_marg = coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) + cluster(PID),
                    data = dat_recurrent)
summary(fit_ph_marg)

#---------------------- Weibull -------------------------------------

##--------------------------frailtypack
#note: in frailtypack, cluster is used to indicate the level for random effects (not a marginal model)

## calendar timescale
fit_wb_cal = frailtyPenal(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) + cluster(PID),
                         data = dat_recurrent,
                         hazard = "Weibull", 
                         recurrentAG = TRUE)
summary(fit_wb_cal)
fit_wb_cal

## gap timescale
fit_wb_gap = frailtyPenal(Surv(start, time_yrs, VL > 50) ~ timedep(n_heartbeats_30days) + Vmed9_everhadTB + sqrt(VL_baseline) + cluster(PID),
                         data = dat_recurrent,
                         hazard = "Weibull",
                         recurrentAG = FALSE,
                         RandDist = "Gamma")
summary(fit_wb_gap)
fit_wb_gap

plot(fit_wb_gap, data = dat_recurrent)

#that is all the functionality provided by frailtypack
# will have to do some fanangling to get residuals out, etc

