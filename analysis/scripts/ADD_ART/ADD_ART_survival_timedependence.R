

rm(list=ls())
# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(survival)
library(survminer)
library(ldatools)


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


#---------------------- CONDITIONAL MODEL ----------------
#define data for modelling
temp = dat_recurrent 
#------------- Voilations of PH
#fit which violates PH
fit_0 <- coxph(Surv(start, time_yrs, VL >= 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) +  frailty(PID),
               data = temp)
summary(fit_0)
sd(fit_0$frail)
cox.zph(fit_0)
par(mfrow=c(3,1))
plot(cox.zph(fit_0))

gg_coxsnell(fit_0, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)

#----------- removing observations and refitting
temp = dat_recurrent %>%
  filter(!(VL > 50 & time_yrs > 0.99))

fit_1a <- coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) + frailty(PID),
               data = temp)
summary(fit_1a)
cox.zph(fit_1a)
plot(cox.zph(fit_1a))

gg_coxsnell(fit_1a, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)
#---------- step function for time-varing coefficients
surv = Surv(dat_recurrent$start, dat_recurrent$time_yrs, dat_recurrent$VL >= 50)
survsplt = survSplit(surv ~ ., dat_recurrent, cut = c(0.4), episode = "wot")

fit_1 <- coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days:strata(wot) + Vmed9_everhadTB + sqrt(VL_baseline) + frailty(PID),
               data = survsplt)
summary(fit_1)

cox.zph(fit_1)
plot(cox.zph(fit_1))

gg_coxsnell(fit_1, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2) 



#------------- continuous function for time-varying
fit_2 <- coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + time_yrs*n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) + frailty(PID),
               data = temp
               )
summary(fit_2)

cox.zph(fit_2)

gg_coxsnell(fit_2, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2) 



tmp_fit <- coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days  + frailty(PID),
                 data = dat_recurrent
)
summary(tmp_fit)

cox.zph(tmp_fit)
plot(cox.zph(tmp_fit))

gg_coxsnell(tmp_fit, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)
