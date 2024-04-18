#ADD-ART Survival Modelling
# 02 April 2024
rm(list=ls())


#load dependencies
library(tidyverse)
library(survival)
library(survminer)
library(KMsurv)
library(ldatools)


path_in = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"



#---- load data
setwd(path_in)
dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)

#----------- data wrangling ----------
##------ model agnostic data wrangling ------
dat$AGE_centered = scale(dat$AGE, scale = F, center = T) #demean age for interpretable intercepts
#number of days between arv initiation and visit 1
dat$days_since_initiation = as.numeric(difftime(dat$BLdate, dat$arv_initiation_date, units = "days"))
dat$days_since_initiation_centered = scale(dat$days_since_initiation, scale = F, center = T) 

#calculate number of days between baseline visit each subsequent visit
dat = dat %>% group_by(PID) %>%
  mutate(days_since_visit_0 = as.numeric(difftime(VISITdate, VISITdate[1], units = "days"))) %>%
  relocate(days_since_visit_0, .after = Visit)


##---------------- time to first viral non-suppression --------------
#create survival indicator for first viral non-suppression: ignore subsequent events
dat = dat %>%
  mutate(VL_ns = ifelse(VL >= 50, 1, 0)) %>%
  relocate(VL_ns, .after = VL)

#find time to first VL_ns event for each PID
first_VS = dat %>% filter(VL_ns == TRUE) %>%
  group_by(PID) %>% slice(1) %>% ungroup() 

#if PID not in above, then they are censored and last attended visit is used
last_VS = dat %>% group_by(PID) %>% filter(missed_visit == 0) %>%
  slice(n()) %>% ungroup() %>% filter(!PID %in% first_VS$PID)
#combine first and last
dat_singleevent = bind_rows(first_VS, last_VS) %>% arrange(PID)

surv_VS1st = Surv(time = dat_singleevent$days_since_visit_0, event = dat_singleevent$VL_ns)



#create survival object for lost to care: 3 consecutive missed visits
#find when PID has 3 consecutive missed visits == 1
dat$lost_to_care = 0
dat = dat %>% group_by(PID) %>%
  mutate(lost_to_care = ifelse(zoo::rollapply(missed_visit, 3, sum, fill = NA, align = "right", partial = TRUE) == 3, 1, 0)) %>%
  relocate(lost_to_care, .after = missed_visit)

surv_LTC = Surv(time = dat$days_since_visit_0, event = dat$lost_to_care)
head(cbind(dat$PID, dat$Visit, dat$missed_visit, dat$lost_to_care, surv_LTC), 20)



##----- recurrent events data wrangling ------
#get data in counting process format
#counting process format: PID is represented by rows as study entry, event 1 time, event 2 time, event m time, final follow-up
#remove all missing viral load events
dat = dat %>%
  filter(!is.na(VL))

all_VS = dat %>% filter(missed_visit == 0) %>% #all observations at visit 12 follow up or with VL_ns event
  filter(VL_ns == TRUE | Visit == 12) 
last_VS = dat %>% group_by(PID) %>% #all observations at final follow up when it is not visit 12
  filter(missed_visit == 0) %>% 
  slice(n()) %>% filter(Visit != 12) %>% ungroup()

#create interval times
all_VS = all_VS %>% group_by(PID) %>%
  mutate(start = lag(days_since_visit_0, n = 1, default = 0)) %>%
  relocate(start, .after = Visit)

last_VS = last_VS %>% mutate(start = 0) %>%
  relocate(start, .after = Visit)

#merge
dat_recurrent = bind_rows(all_VS, last_VS) %>%
  arrange(PID, Visit) %>%
  relocate(VL_ns, .after = days_since_visit_0)

#remove baseline
dat_recurrent = dat_recurrent %>% filter(Visit > 0)




#finally create enum variable: cumulative number of events for each PID
dat_recurrent = dat_recurrent %>%
  group_by(PID) %>%
  mutate(event_number = cumsum(VL_ns)) %>%
  relocate(event_number, .after = VL_ns)

Surv_cp = Surv(time = dat_recurrent$start, 
               time2 = dat_recurrent$days_since_visit_0,
               dat_recurrent$VL_ns, 
               type = "counting")
ggsurvplot(survfit(Surv_cp ~ 1), data = dat_recurrent)



#--------------------------------- analysis ---------------------------------
#---------------- EDA ---------------
##--- Viral nonsuppression ---
km_fit = survfit(surv_VS1st ~ 1, data = dat_singleevent)
ggsurvplot(km_fit,
           data=dat)






#----------------- Cox PH Model:  ---------------
##---- Time to first Viral nonsuppression event -------
VSfit_null = coxph(surv_VS1st ~ 1, 
                   data = dat_singleevent)
summary(VSfit_null)

VSfit_null_2 = coxph(Surv(days_since_visit_0, VL_ns) ~ 1,
                     data = dat, id=PID)
summary(VSfit_null_2)
#----- demographic covariates
VSfit_1a = update(VSfit_null,
                  . ~ . + AGE_centered)
VSfit_1b = update(VSfit_1a,
                  . ~ . + I(AGE_centered^2))
anova(VSfit_null, VSfit_1a, VSfit_1b) #no age terms are close to significant
summary(VSfit_1a)

VSfit_1c = coxph(surv_VS1st ~ GENDER, data = dat_singleevent)
anova(VSfit_null, VSfit_1c) #gender is borderline significant
summary(VSfit_1c)

VSfit_1d = coxph(surv_VS1st ~ GENDER * AGE_centered, data = dat_singleevent)
summary(VSfit_1d) #interaction is not significant but improves gender effect
anova(VSfit_1c, VSfit_1d)



#----- medical history covariates: baseline VL, TB diagnosis, anemia, HIV illness
VSfit_1e = update(VSfit_1d,
                  . ~ . +  VL_baseline)
summary(VSfit_1e) #baseline VL is significant

VS_fit_1f = update(VSfit_1e,
                   . ~ . +  Vmed9_everhadTB)
summary(VS_fit_1f) #TB diagnosis is significant

VS_fit_1g = update(VS_fit_1f,
                   . ~ . + Vmed9_everhadTB:VL_baseline)
summary(VS_fit_1g) #interaction is not significant
anova(VS_fit_1f, VS_fit_1g)

VS_fit_1h = update(VS_fit_1f,
                   .~. + Vmed10_anemia_ever)
summary(VS_fit_1h)

VS_fit_1i = update(VS_fit_1f,
                   .~. + Vmed7_HIVillness_ever)
summary(VS_fit_1i)

VS_fit_j = update(VS_fit_1f,
                  .~. + days_since_initiation_centered)
summary(VS_fit_j)

#----- clinical covariates: baseline hematocrit, BMI
VS_fit_1k = update(VS_fit_1f,
                  . ~ . + HCT_baseline)
summary(VS_fit_1k)

VS_fit_1l = update(VS_fit_1f,
                  . ~ . + bmi_1st)
summary(VS_fit_1l)

#----- reintroduce covariates that are included in literature (Jennings et al 2022)
# baseline HCT, BMI, age, gender
# excluding BMI because many missings
VS_fit_1m = update(VS_fit_1f,
                   . ~ . + HCT_baseline + AGE_centered * GENDER)
summary(VS_fit_1m)


#----------------------------- RECURRENT EVENTS -----------------------------
##----------------- Cox PH Model w unordered failure times ---------------
#fit above model with robust se
VS_fit_robust = coxph(surv_VS1st ~ VL_baseline + Vmed9_everhadTB + HCT_baseline + AGE_centered * GENDER,
                       cluster = PID,
                       data= dat_singleevent)
summary(VS_fit_robust)

#frailty model: random effect for PID
VS_fit_frailty_1 = coxph(surv_VS1st ~ VL_baseline + Vmed9_everhadTB + HCT_baseline + AGE_centered + GENDER + frailty(PID, dist="gamma"),
                         data = dat_singleevent)
summary(VS_fit_frailty_1) #this has convergence error

VS_fit_frailty_2 = coxph(surv_VS1st ~ VL_baseline + Vmed9_everhadTB + HCT_baseline + AGE_centered + GENDER + frailty(PID, dist="gauss"),
                         data = dat_singleevent)
summary(VS_fit_frailty_2) #convergence issues


###------------------- Anderson-Gill Model ---------------
fit_ag1 = coxph(Surv(start, days_since_visit_0, VL_ns) ~ VL_baseline + Vmed9_everhadTB + HCT_baseline + AGE_centered*GENDER +
                  cluster(PID), data = dat_recurrent)
summary(fit_ag1)

fit_ag2 = coxph(Surv(start, days_since_visit_0, VL_ns) ~ VL_baseline + Vmed9_everhadTB + HCT_baseline + AGE_centered + GENDER +
                  cluster(PID), data = dat_recurrent)
summary(fit_ag2)

fit_ag3 = coxph(Surv(start, days_since_visit_0, VL_ns) ~ VL_baseline + Vmed9_everhadTB + HCT_baseline +
                  cluster(PID), data = dat_recurrent)
summary(fit_ag3)

fit_ag4 = coxph(Surv(start, days_since_visit_0, VL_ns) ~ log(VL_baseline+1) + Vmed9_everhadTB +
                  cluster(PID), data = dat_recurrent)
summary(fit_ag4) #magnitude of significant coefficients barely change when excluding other covariates


#----------------- DIAGNOSTICS --------------------------------------
##----------------- naive model -------------
#check proportional hazards assumption
cox.zph(VS_fit_1m)
plot(cox.zph(VS_fit_1m)) #lookin good, good lookin

#check overall fit
gg_coxsnell(VS_fit_1m, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)

#check model specification
plot(lowess(residuals(VS_fit_1m, type="martingale")~tmp$VL_baseline), type="l",
     xlab="age.centered", ylab="Martingale residuals")

#check outliers
ggcoxdiagnostics(VS_fit_1m, type="deviance") 

#check for influentials
ggcoxdiagnostics(VS_fit_1m, type="score", ox.scale = "time")

##----------------- robust model -------------
cox.zph(VS_fit_robust)
plot(cox.zph(VS_fit_robust)) #not lookin good

gg_coxsnell(VS_fit_robust, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2)

plot(lowess(residuals(VS_fit_robust, type="martingale")~dat$VL_baseline), type="l",
     xlab="age.centered", ylab="Martingale residuals")

ggcoxdiagnostics(VS_fit_robust, type="deviance")

ggcoxdiagnostics(VS_fit_robust, type="score", ox.scale = "time")

##----------------- AG model -------------
#check proportional hazards assumption
cox.zph(fit_ag4)
plot(cox.zph(fit_ag4)) #lookin good, good lookin

#check overall fit
gg_coxsnell(fit_ag4, type="cumu") +
  geom_abline(intercept=0, slope=1, col=2) #overall fit a bit suss

#check model specification
mg1 = residuals(coxph(Surv(start, days_since_visit_0, VL_ns) ~ VL_baseline, data = dat_recurrent), type="martingale")
mg2 = residuals(coxph(Surv(start, days_since_visit_0, VL_ns) ~ Vmed9_everhadTB, data = dat_recurrent), type="martingale")

plot(mg1~dat_recurrent$Vmed9_everhadTB, type="p")

plot(mg2~log(dat_recurrent$VL_baseline[which(!is.na(dat_recurrent$Vmed9_everhadTB))]+1), type="p"
     )
#check outliers
ggcoxdiagnostics(fit_ag4, type="deviance")  

#check for influentials
ggcoxdiagnostics(fit_ag4, type="score", ox.scale = "time")



#---------------- write objects to file ------------------
setwd(paste0(path_out, "/Models"))
saveRDS(VS_fit_1m, "ADD-ART_coxph_naive.rds")
saveRDS(VS_fit_robust, "ADD-ART_coxph_naive_robust.rds")
saveRDS(fit_ag4, "ADD-ART_AG.rds")

