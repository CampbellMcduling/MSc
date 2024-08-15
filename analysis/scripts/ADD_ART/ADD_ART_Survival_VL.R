#ADD-ART Survival Modelling
# 02 April 2024

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
###---------------- supporting stats ---------------
n_events_per_id = dat_long_eam  %>% group_by(id) %>% summarise(n_events = sum(VL>50))
print(n_events_per_id, n=1000)
mean(n_events_per_id$n_events, na.rm = TRUE)
n_events_per_id %>% filter(n_events >= 1 & !is.na(n_events)) %>%
  summarise(mean = mean(n_events))

#create quantile variables for numerical covariates
colnames(dat_recurrent)
dat_recurrent %>% 
  filter(PID %in% unique(PID)) %>%
  summarise(median_TFV_baseline = quantile(TFV_baseline, na.rm = TRUE),
            median_vsra1_baseline = quantile((1-sra_1_baseline)*30,  na.rm = T),
            median_VL_baseline = quantile((VL_baseline), na.rm = TRUE))

dat_recurrent = dat_recurrent %>%
  ungroup() %>%
  mutate(poor_adherence_TFV_baseline = ifelse(TFV_baseline < 1060.4, 1, 0),
         poor_adherence_SR_baseline = ifelse((1-sra_1_baseline)*30 > 2, 1, 0),
         unsuppressed_VL_baseline = ifelse(VL_baseline > 49, 1, 0))

#aggregate adherence data for each individual
dat_recurrent_adherence = dat_recurrent %>%
  group_by(PID) %>%
  summarise(mean_n_heartbeats_30days = mean(n_heartbeats_30days, na.rm = TRUE),
            mean_vsra1 = mean(Vsra1_days_nonadherent, na.rm = TRUE),
            mean_TFV = mean(TFV, na.rm = TRUE))
dat_recurrent_adherence

#join 
dat_recurrent = dat_recurrent_adherence %>%
  left_join(dat_recurrent,
            by = "PID")

# global means 
tmp = dat_recurrent %>%
  summarise(mean(mean_vsra1, na.rm = TRUE),
            mean(mean_TFV, na.rm = TRUE),
            mean(mean_n_heartbeats_30days, na.rm = TRUE))

dat_recurrent = dat_recurrent %>%
  mutate(poor_adherence_SR = ifelse(mean_vsra1 > mean(mean_vsra1, na.rm = TRUE), 1, 0),
         poor_adherence_TFV = ifelse(mean_TFV < mean(mean_TFV, na.rm = TRUE), 1, 0),
         poor_adherence_eam = ifelse(mean_n_heartbeats_30days > mean(mean_n_heartbeats_30days, na.rm = TRUE), 1, 0))

dat_recurrent = dat_recurrent %>%
  mutate(BMI_cat = case_when(
    BMIcat == "< 18.5" ~ "< 25",
    BMIcat == "18.5-24.9" ~ "< 25",
    BMIcat == "25-29.9" ~ "25 - 29.9",
    BMIcat == "30 or +" ~ "> 30"
  ))


#time from visit 0
Surv3b = with(dat_recurrent, 
              Surv(time = start, time2 = time_yrs,
                   event = VL> 50, type = "counting"))




#--------------------------------- analysis ---------------------------------
#----------------------------- RECURRENT EVENTS -----------------------------
# forward selection procedure informed by theory and EDA



##------------------- MARGINAL MODELS (Anderson-Gill) ---------------

fit_ag0 = coxph(Surv3b ~ 1,
                cluster = PID,
                data = dat_recurrent) 
summary(fit_ag0)

###----------------summarised adherence models ----------------
#---- self report
fit_ag_sr1 = update(fit_ag0,
                   .~. + mean_vsra1)
summary(fit_ag_sr1)
#---- eam
fit_ag_eam1 = update(fit_ag0,
                   .~. + mean_n_heartbeats_30days)
summary(fit_ag_eam1)

#---- TFV
fit_ag_tfv1 = update(fit_ag0,
                   .~. + mean_TFV)
summary(fit_ag_tfv1)

#---- log10TFV
fit_ag_tfv1 = update(fit_ag0,
                   .~. + log10(mean_TFV))
summary(fit_ag_tfv1) #log10 TFV is more meaningful given the sensitivity of the response over the large range
#confint
exp(-1*coefficients(fit_ag_tfv1) + c(-1,1) * 1.96*sqrt(vcov(fit_ag_tfv1)[1,1]))

###----------------time-varying adherence models ----------------
fit_ag_sr1a = update(fit_ag0,
                     .~. + Vsra1_days_nonadherent)
summary(fit_ag_sr1a)

fit_ag_eam1a = update(fit_ag0,
                     .~. + n_heartbeats_30days)
summary(fit_ag_eam1a)


fit_ag_tfv1a = update(fit_ag0,
                     .~. + log10_TFV)
summary(fit_ag_tfv1a)
#confint
exp(-1*coefficients(fit_ag_tfv1a) + c(-1,1) * 1.96*sqrt(vcov(fit_ag_tfv1a)[1,1]))



#----------- adjusted models ----------------
#---- self report
fit_ag_sr2 = update(fit_ag_sr1a,
                   .~. + BMI_cat)
summary(fit_ag_sr2)

fit_ag_sr3 = update(fit_ag_sr1a,
                 .~. + Vmed5_rate_health)
summary(fit_ag_sr3)

fit_ag_sr4 = update(fit_ag_sr1a,
                 .~. + Vmed9_everhadTB)
summary(fit_ag_sr4)

fit_ag_sr5 = update(fit_ag_sr4,
                 .~. + sqrt(VL_baseline))
summary(fit_ag_sr5)

#----- eam
fit_ag_eam2 = update(fit_ag_eam1a,
                   .~. + BMI_cat)
summary(fit_ag_eam2)

fit_ag_eam3 = update(fit_ag_eam1a,
                 .~. + Vmed5_rate_health)
summary(fit_ag_eam3)

fit_ag_eam4 = update(fit_ag_eam1a,
                 .~. + Vmed9_everhadTB)
summary(fit_ag_eam4)

fit_ag_eam5 = update(fit_ag_eam4,
                 .~. + sqrt(VL_baseline))
summary(fit_ag_eam5)


#----- TFV
fit_ag_tfv2 = update(fit_ag_tfv1a,
                   .~. + BMI_cat)
summary(fit_ag_tfv2)

fit_ag_tfv3 = update(fit_ag_tfv1a,
                 .~. + Vmed5_rate_health)
summary(fit_ag_tfv3)

fit_ag_tfv4 = update(fit_ag_tfv1a,
                 .~. + Vmed9_everhadTB)
summary(fit_ag_tfv4)

fit_ag_tfv5 = update(fit_ag_tfv4,
                 .~. + sqrt(VL_baseline))
summary(fit_ag_tfv5)
#confint
exp(-1*coefficients(fit_ag_tfv5) -  (sqrt(diag(vcov(fit_ag_tfv5)))*1.96))
exp(-1*coefficients(fit_ag_tfv5) +  (sqrt(diag(vcov(fit_ag_tfv5)))*1.96))



#------------fit final models but frailties instead of marginal
fit_sr5_frail = coxph(Surv3b ~ Vsra1_days_nonadherent + Vmed9_everhadTB + sqrt(VL_baseline) + frailty(PID),
                      data = dat_recurrent)
fit_eam5_frail = coxph(Surv(start, time_yrs, VL > 50) ~ n_heartbeats_30days + Vmed9_everhadTB + sqrt(VL_baseline) +  frailty(PID),
                      data = dat_recurrent)
fit_tfv5_frail = coxph(Surv3b ~ log10_TFV + Vmed9_everhadTB + sqrt(VL_baseline) +  frailty(PID),
                      data = dat_recurrent)
summary(fit_sr5_frail)
summary(fit_eam5_frail)
summary(fit_tfv5_frail)


#----------------- DIAGNOSTICS --------------------------------------
model = fit_ag_tfv5
adher = dat_recurrent$n_heartbeats_30days
adhername = "EAM"
##----------------- AG model -------------
#check proportional hazards assumption
cox.zph(model)
par(mfrow=c(3,1))
plot(cox.zph(model), hr=F)

#check overall fit
gg_coxsnell(model, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2) #overall fit a bit suss

#check model specification
mg1 = residuals(coxph(Surv3b ~ Vmed9_everhadTB + sqrt(VL_baseline), #excludes adherence
                      data = dat_recurrent), type="deviance")
mg2 = residuals(coxph(Surv3b ~ adher + Vmed9_everhadTB, #excludes VL_baseline
                      data = dat_recurrent), type="deviance")

plot(mg1~adher , type="p",
     xlab = adhername,
     ylab = "Deviance residuals") 


plot(mg2~dat_recurrent$VL_baseline, type="p")
plot(mg2~sqrt(dat_recurrent$VL_baseline), type="p",
     xlab = "sqrt(VL_baseline)",
     ylab = "Deviance residuals")



#check outliers
ggcoxdiagnostics(model, type="deviance")     #strange p[atterns in these resids, a strip at resids = 1

ggcoxdiagnostics(model, type="dfbeta") #check for influentials


#check for influentials
ggcoxdiagnostics(model, type="score", ox.scale = "time")

#further investigate influential obs
model = fit_ag_tfv5
ggcoxdiagnostics(model, type="score", ox.scale = "time")
thresh = 2000
colnames(residuals(model, "score"))
indices = which(abs(residuals(model, "score"))[,1] > thresh)
dat_recurrent[indices,] %>%
  select(PID, Visit, time_yrs, VL,
         Vsra1_days_nonadherent, n_heartbeats_30days, n_intakes_30days,
         Vmed9_everhadTB, VL_baseline, TFV) %>% View()



#investigate the concentration of resids at -1
model_frail = coxph(Surv3b ~ TFV + sqrt(VL_baseline) + Vmed9_everhadTB + frailty(PID), data=dat_recurrent)
model_clust = coxph(Surv3b ~ TFV + sqrt(VL_baseline) + Vmed9_everhadTB + cluster(PID), data=dat_recurrent)
hist(residuals(model_clust, "deviance"))
hist(residuals(model_frail, "deviance"))
#scatterplot, color points by whether they are censored
plot(residuals(model_clust, "deviance"))
plot(residuals(model_clust, "deviance"), 
     col=ifelse(dat_recurrent$Visit == 11 & dat_recurrent$VL <=50, "red", "blue")) 

ggcoxdiagnostics(model_clust, type="deviance")
ggcoxdiagnostics(model_frail, type="deviance")


tst_frail = coxph(Surv3b ~  sqrt(VL_baseline) + Vmed9_everhadTB + frailty(PID), data=dat_recurrent)
plot(residuals(model_clust, "deviance"), 
     col=ifelse(dat_recurrent$Visit == 11 & dat_recurrent$VL <=50, "red", "blue")) 
hist(residuals(tst_frail, "deviance"))

#---------------- write objects to file ------------------
setwd(paste0(path_out, "/Models"))
saveRDS(VS_fit_1m, "ADD-ART_coxph_naive.rds")
saveRDS(VS_fit_robust, "ADD-ART_coxph_naive_robust.rds")
saveRDS(fit_ag4, "ADD-ART_AG.rds")

