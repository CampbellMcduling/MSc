#ADD-ART Joint modelling
# 14 March 2024
rm(list=ls())
#---- dependencies
library(tidyverse)
library(DHARMa)
library(lme4)
library(glmmTMB)
library(ggeffects)
library(JMbayes2)
library(GLMMadaptive)

path_models = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART/Models"
path_data = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"

#------- load data and model objects ---------------
setwd(path_models)
fit_glmm = readRDS("ADD-ART_glmm.rds")
fit_coxph_naive = readRDS("ADD-ART_coxph_naive.rds")
fit_coxph_naive_robust = readRDS("ADD-ART_coxph_naive_robust.rds")
fit_coxph_ag = readRDS("ADD-ART_AG.rds")

##-----------long data -----------------
#load data
setwd(path_data)
dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)

dat$id = dat$PID



dat = dat %>%
  group_by(id) %>%
  mutate(time = as.numeric(difftime(VISITdate, VISITdate[1], units = "days")),
         time_yrs = time/365 ) %>%
  relocate(., .after = Visit)



#format longitudinal data
dat_long = dat %>%
  select(id, Visit, time, time_yrs, missed_visit, VISITdate, n_intakes, n_heartbeats, n_none, wp_adherence_prop, wp_days_covered,
         TFV, HCT, BMI, BMIcat, bmi_1st, AGE, GENDER, TFV_baseline, HCT_baseline, VL_baseline, VL) 

#to deal with missing EAM data at time 0, we will need to impute these values
#for now, lets take the first observation carried backwards (Visit 1)
dat_long = dat_long %>%
  group_by(id) %>%
  mutate(n_heartbeats = case_when(Visit == 0 ~ n_heartbeats[2],
                                  TRUE ~ n_heartbeats),
         n_intakes = case_when(Visit == 0 ~ n_intakes[2],
                               TRUE ~ n_intakes),
         n_none = case_when(Visit == 0 ~ n_none[2],
                           TRUE ~ n_none))



##-----------recurrent data -----------------
#get data in counting process format
#counting process format: id is represented by rows as study entry, event 1 time, event 2 time, event m time, final follow-up
#remove all missing viral load events
dat_recurrent = dat %>%
  filter(!is.na(VL))

all_VS = dat_recurrent %>% filter(missed_visit == 0) %>% #all observations at visit 12 follow up or with VL_ns event
  filter(VL >= 50) 
last_obs = dat_recurrent %>% group_by(id) %>% #all censored observations at final follow up 
  filter(missed_visit == 0) %>% 
  slice(n())  %>%
  filter(VL < 50) %>%
  ungroup()

dat_recurrent = bind_rows(all_VS, last_obs) %>%
  arrange(id, Visit)
#create interval times
dat_recurrent = dat_recurrent %>% group_by(id) %>%
  mutate(start = lag(time_yrs, n = 1, default = 0)) %>%
  relocate(c(start, time_yrs), .after = Visit)


#remove baseline
dat_recurrent = dat_recurrent %>% filter(Visit > 0) %>% relocate(VL, .after = time_yrs)


#finally create enum variable: cumulative number of events for each id
dat_recurrent = dat_recurrent %>%
  group_by(id) %>%
  mutate(event_number = cumsum(VL >= 50)) %>%
  relocate(event_number, .after = VL)


##-------------------- missing values ---------------------
n_distinct(dat_recurrent$id)
n_distinct(dat_long$id)
dat_long %>% filter( id %in% setdiff(dat_long$id, dat_recurrent$id)) %>% View() #these guys only have baseline
dat_long = dat_long %>%
  filter(id %in% dat_recurrent$id)


#filter out all rows with missing values for covariates and outcomes
dat_recurrent = dat_recurrent %>%
  filter(!is.na(Visit), !is.na(VL), !is.na(VL_baseline), !is.na(Vmed9_everhadTB),
         !is.na(n_heartbeats), !is.na(AGE), !is.na(GENDER))
dat_long = dat_long %>%
  filter(!is.na(n_heartbeats)) %>%
  filter(id %in% dat_recurrent$id)

n_distinct(dat_recurrent$id)
n_distinct(dat_long$id)



#------------------------- joint modelling -----------------------
##======== Fit joint model ======
#refit long model with glmmadaptive
summary(fit_glmm)
dat_long$log_wp_days = log(dat_long$n_heartbeats + dat_long$n_intakes + 1)
off = dat_long$log_wp_days
summary(off)
fit_glmmadap = mixed_model(
  fixed = n_heartbeats ~ time_yrs + I(time_yrs^2) + offset(off),
  random = ~ time_yrs | id,
  family = negative.binomial(),
  data = dat_long,
  na.action = na.omit
  
)
summary(fit_glmmadap)
plot(simulateResiduals(fit_glmmadap))



#extend to recurrent events model
#refit with offset in dataframe
fit_glmmadap = mixed_model(
  fixed = n_heartbeats ~ time_yrs + I(time_yrs^2) + offset(log_wp_days),
  random = ~ time_yrs | id,
  family = negative.binomial(),
  data = dat_long,
  na.action = na.omit
  
)
summary(fit_coxph_ag)
fit_surv_recurr = coxph(Surv(start, time_yrs, VL >= 50) ~ log(VL_baseline + 1) + Vmed9_everhadTB + cluster(id), 
                        data = dat_recurrent)
which(dat_recurrent$start >= dat_recurrent$time_yrs)
summary(fit_surv_recurr)
jm_fit_2 = jm(fit_surv_recurr, 
              fit_glmmadap, 
              time_var = "time",
              id_var = "id",
              recurrent = "gap",
              functional_forms = ~ value(n_heartbeats))
summary(jm_fit_2)

traceplot(jm_fit_2)




##======== try gmvjoint ======
library(gmvjoint)

long_form = list(
  n_heartbeats ~ time + I(time^2) + I(time^3) + offset(log_visit) + (time|id)
)
family = list("negbin")

surv_form = list(
  Surv(time, VL_ns) ~ log(VL_baseline + 1) + Vmed9_everhadTB + cluster(id)
)

jm_fit = joint(
  long_form,
  surv_form,
  data=dat_jm,
  family = family
)
