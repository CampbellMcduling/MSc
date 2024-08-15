#ADD-ART Joint modelling
# 14 March 2024
rm(list=ls())
#---- dependencies
library(tidyverse)
library(DHARMa)
library(lme4)
library(ggeffects)
library(JMbayes2)
library(GLMMadaptive)

path_models = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART/models"
path_data = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"

#------- load data and model objects ---------------
setwd(path_models)
fit_eam_pois_quad = readRDS("EAM30day_pois_quad_final.rds")
fit_eam_pois_cub = readRDS("EAM30day_pois_cubic_final.rds")
fit_eam_nb_quad = readRDS("EAM30day_nb_quad_final.rds")
fit_eam_nb_cub = readRDS("EAM30day_nb_cubic_final.rds")
fit_coxph_naive = readRDS("ADD-ART_coxph_naive.rds")
fit_coxph_naive_robust = readRDS("ADD-ART_coxph_naive_robust.rds")
fit_coxph_ag = readRDS("ADD-ART_AG.rds")

##-----------long data -----------------
#load data
setwd(path_data)
dat_recurrent = read.csv("ADDART_recurrent_14May2024.csv")

dat_long_eam = read.csv("ADDART_long_eam_14May2024.csv")

dat_long_tfv = read.csv("ADDART_long_tfv_14May2024.csv")

dat_long_joint = read.csv("ADDART_long_joint_14May2024.csv")


#------------------------- joint modelling: EAM data -----------------------
##-----------------longitudinal model -----------------
off = dat_long_eam$log_wp_days
summary(off)


#for joint model fit: refit each model with offset in dataframe
fit_eam_pois_quad = update(fit_eam_pois_quad, 
                           fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + 
                             offset(log_wp_days_30))
fit_eam_pois_cub = update(fit_eam_pois_cub,
                          fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + I(time_yrs^3) + 
                            offset(log_wp_days_30))
fit_eam_nb_quad = update(fit_eam_nb_quad,
                         fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + 
                           offset(log_wp_days_30),
                         family = GLMMadaptive:: negative.binomial)
fit_eam_nb_cub = update(fit_eam_nb_cub,
                        fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + I(time_yrs^3) + 
                          offset(log_wp_days_30),
                        family = GLMMadaptive:: negative.binomial)


##----------------- cox model: recurrent events -----------------
#frailty model
fit_frailty = coxph(Surv(start, time_yrs, VL >= 50) ~
                      sqrt(VL_baseline) + Vmed9_everhadTB + n_heartbeats_30days + frailty(id),
                    data = dat_recurrent)
summary(fit_frailty)
tst = gg_coxsnell(fit_frailty)  + geom_abline(intercept = 0, slope =1)
tst
cs = tst$data$coxsnell
hist(cs, breaks = 20, col = "lightblue", border = "black", main = "Cox-Snell residuals")
hist(rexp(length(cs), rate = 1), breaks = 20, col = "lightblue", border = "black", main = "Exponential distribution")
ggcoxdiagnostics(fit_frailty, type = "deviance")
ggcoxdiagnostics(fit_frailty, type = "score")
cox.zph(fit_frailty)
plot(cox.zph(fit_frailty))

summary(fit_coxph_ag)
fit_surv_recurr = coxph(Surv(start, time_yrs, VL >= 50) ~
                          sqrt(VL_baseline) + Vmed9_everhadTB,
                        data = dat_recurrent)
summary(fit_surv_recurr)


##----------------- fit joint model -----------------
#----poisson quadratic
set.seed(2024)
jm_fit_pois_quad = jm(fit_surv_recurr, 
              fit_eam_pois_quad, 
              time_var = "time_yrs",
              id_var = "id",
              recurrent = "gap",
              functional_forms = ~ value(n_heartbeats_30),
              control = list(n_iter = 10000, n_chains = 6))
summary(jm_fit_pois_quad)

traceplot(jm_fit_pois_quad) #nice
densplot(jm_fit_pois_quad) # nice

set.seed(2024)
jm_fit_pois_cub = jm(fit_surv_recurr, 
                      fit_eam_pois_cub, 
                      time_var = "time_yrs",
                      id_var = "id",
                      recurrent = "gap",
                      functional_forms = ~ value(n_heartbeats_30),
                      control = list(n_iter = 10000, n_chains = 6))
summary(jm_fit_pois_cub)

traceplot(jm_fit_pois_cub)
densplot(jm_fit_pois_cub)


#----negative binomial quadratic
set.seed(2024)
jm_fit_nb_quad = jm(fit_surv_recurr, 
              fit_eam_nb_quad, 
              time_var = "time_yrs",
              id_var = "id",
              recurrent = "gap",
              functional_forms = ~ value(n_heartbeats_30),
              control = list(n_iter = 10000, n_chains = 6))
summary(jm_fit_nb_quad)

densplot(jm_fit_nb_quad)
traceplot(jm_fit_nb_quad)

#----negative binomial cubic
set.seed(2024)
jm_fit_nb_cub = jm(fit_surv_recurr, 
              fit_eam_nb_cub, 
              time_var = "time_yrs",
              id_var = "id",
              recurrent = "gap",
              functional_forms = ~ value(n_heartbeats_30),
              control = list(n_iter = 10000, n_chains = 6, n_burnin = 500))
summary(jm_fit_nb_cub)

#after comparing model fits, we settle on the negative binomial cubic model
#refit with more iterations, larger burnin
set.seed(2024)
jm_fit_nb_cub = jm(fit_surv_recurr, 
                   fit_eam_nb_cub, 
                   time_var = "time_yrs",
                   id_var = "id",
                   recurrent = "gap",
                   functional_forms = ~ value(n_heartbeats_30),
                   control = list(n_iter = 50000, n_chains = 6, n_burnin = 10000))
summary(jm_fit_nb_cub)


#save and load model
saveRDS(jm_fit_nb_cub, paste0(path_models, "/jm_fit_nb_cub.rds"))
jm_fit_nb_cub = readRDS(paste0(path_models, "/jm_fit_nb_cub.rds"))

traceplot(jm_fit_nb_cub)
densplot(jm_fit_nb_cub)

par(mfrow = c(2,2))
traceplot(jm_fit_nb_cub, c("betas"))
par(mfrow = c(2,2))
traceplot(jm_fit_nb_cub, c("D"))
traceplot(jm_fit_nb_cub, c("sigmas"))


par(mfrow = c(1, 2))
traceplot(jm_fit_nb_cub, c("gammas"))
par(mfrow = c(1, 1))
traceplot(jm_fit_nb_cub, c("alphas"))

par(mfrow = c(2, 2))
traceplot(jm_fit_nb_cub, c("bs_gammas"))
traceplot(jm_fit_nb_cub, c("tau_bs_gammas"))

#need custom code to investigate frailty convergence
frailties = jm_fit_nb_cub$mcmc$frailty

sig_frail = jm_fit_nb_cub$mcmc$sigmaF
#convert sig_frail list to df
sig_frail_df = as.data.frame(matrix(unlist(sig_frail), ncol = length(sig_frail), byrow = FALSE))
#pivot to long format
sig_frail_df_long = pivot_longer(sig_frail_df, cols = everything(), names_to = "chain", values_to = "sigmaF")

ggplot(sig_frail_df_long, aes(x = sigmaF, y =, group=chain, fill = chain)) +
  geom_density(alpha = 0.5) + theme_minimal()

plot(1:nrow(sig_frail_df), sig_frail_df$V1, 
     type = "l", col = "black", lwd = 3,
     xlab = "Iterations", ylab = "sigmaF")
lines(1:nrow(sig_frail_df), sig_frail_df$V2, col = "lightgreen")
lines(1:nrow(sig_frail_df), sig_frail_df$V3, col = "darkblue")
lines(1:nrow(sig_frail_df), sig_frail_df$V4, col = "purple")
lines(1:nrow(sig_frail_df), sig_frail_df$V5, col = "lightblue", lwd = 3)
lines(1:nrow(sig_frail_df), sig_frail_df$V6, col = "hotpink", lwd =3)

#--------------------------------- model visualisations --------------
### lets try visualise the longitudinal submodel and compare to the long model when fit independently of the survival
dat_recurrent$log_wp_days_30 = log(dat_recurrent$n_heartbeats_30days + dat_recurrent$n_intakes_30days + 1)
newdat = list(newdataL = dat_long_eam, newdataE = dat_recurrent)
pred_longjm = predict(jm_fit_nb_cub,
                     process = "longitudinal",
                     type = "subject_specific",
                     type_pred = "response",
                     newdata = newdat,
                     return_newdata = TRUE)
plot(pred_longjm)

pred_long_solo = predict(fit_eam_nb_cub,
                        newdata = dat_long_eam,
                        type = "subject_specific",
                        type_pred = "response",
                        return_newdata = TRUE,
                        se.fit = TRUE)
str(pred_long_solo)


dat_long_compare = dat_long_eam %>%
  mutate(pred_jm = pred_longjm$pred_n_heartbeats_30days,
         upper_jm = pred_longjm$upp_n_heartbeats_30days,
         lower_jm = pred_longjm$low_n_heartbeats_30days,
         pred_longsolo = pred_long_solo$pred,
         upper_longsolo = pred_long_solo$upp,
         lower_longsolo = pred_long_solo$low)


#visualise
#select subset of PIDs
set.seed(420)
id_samp = sample(unique(dat_long_eam$id), 9)

ggplot(dat_long_compare %>% filter(id %in% id_samp), aes(x = time, y = n_heartbeats_30days, group = id)) +
  geom_line(aes(y = pred_jm), col = "red") +
  geom_ribbon(aes(ymin = lower_jm, ymax = upper_jm), fill = "red", alpha = 0.2) +
  geom_line(aes(y = pred_longsolo), col = "blue") +
  geom_ribbon(aes(ymin = lower_longsolo, ymax = upper_longsolo), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = n_heartbeats_30days), col = "black") +
  facet_wrap(~id) +
  theme_minimal()

#amend above plot by adding legend for observed, JM predicted and long solo predicted lines
ggplot(dat_long_compare %>% filter(id %in% id_samp), aes(x = time, y = n_heartbeats_30days, group = id)) +
  geom_line(aes(y = pred_jm), col = "red", linetype = "solid") +
  geom_ribbon(aes(ymin = lower_jm, ymax = upper_jm), fill = "red", alpha = 0.2) +
  geom_line(aes(y = pred_longsolo), col = "blue", linetype = "solid") +
  geom_ribbon(aes(ymin = lower_longsolo, ymax = upper_longsolo), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = n_heartbeats_30days), col = "black", linetype = "dashed") +
  facet_wrap(~id) +
  theme_minimal() +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  labs(col = "Line type")

#------------- try some mean level predictions -----------
pred_longjm_mean = predict(jm_fit_nb_cub,
                          process = "longitudinal",
                          type = "mean",
                          type_pred = "response",
                          newdata = newdat,
                          return_newdata = TRUE)
str(pred_longjm_mean)

pred_long_solo_mean = predict(fit_eam_nb_cub,
                             newdata = dat_long_eam,
                             type = "mean",
                             type_pred = "response",
                             return_newdata = TRUE,
                             se.fit = TRUE)
str(pred_long_solo_mean)


dat_long_compare_mean = dat_long_eam %>%
  mutate(pred_jm_mean = pred_longjm_mean$pred_n_heartbeats_30days,
         pred_longsolo_mean = pred_long_solo_mean$pred)

ggplot(dat_long_compare_mean %>% filter(id %in% id_samp[2]), aes(x = time, y = n_heartbeats_30days, group = id)) +
  geom_line(aes(y = pred_jm_mean), col = "red") +
  geom_line(aes(y = pred_longsolo_mean), col = "blue") +
  theme_minimal() +
  #custom legend (dashed black lines for observed, solid red for JM, solid blue for long solo)
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  labs(col = "Line type")

ggplot(dat_long_compare_mean %>% filter(id %in% id_samp), aes(x = time, y = n_heartbeats_30days, group = id)) +
  geom_line(data = dat_long_compare_mean %>% filter(id == id_samp[2]),aes(y = pred_jm_mean, color = "Joint Sub-model Fit", linetype = "Joint Sub-model Fit")) +
  geom_line(data = dat_long_compare_mean %>% filter(id == id_samp[2]),aes(y = pred_longsolo_mean, color = "Longitudinal Model Fit", linetype = "Longitudinal Model Fit")) +
  geom_line(aes(y = n_heartbeats_30days, color = "Observed", linetype = "Observed")) +
  scale_color_manual(values = c("Observed" = "black", "Joint Sub-model Fit" = "red", "Longitudinal Model Fit" = "blue")) +
  scale_linetype_manual(values = c("Observed" = "dashed", "Joint Sub-model Fit" = "solid", "Longitudinal Model Fit" = "solid")) +
  labs(color = "Legend", linetype = "Legend") +
  theme_minimal()


#------------- some survival and longitudinal predictions -----------
# following https://drizopoulos.github.io/JMbayes2/articles/Dynamic_Predictions.html
#follow up time
t0 <- 1
t1 <- 2
ND_long <- dat_long_eam %>%
  filter(id %in% id_samp) %>%
  filter(time_yrs < t0)
ND_surv <- dat_recurrent %>%
  filter(id %in% id_samp) 
tmp = dat_long_eam %>%
  filter(id %in% id_samp) %>%
  filter(time_yrs <= t0)
#get tmp in counting process format
tmp = tmp %>%
  group_by(id) %>%
  mutate(start = lag(time_yrs, default = 0)) %>%
  relocate(start, .before = time_yrs) %>%
  filter(time_yrs != 0) %>%
  left_join(dat_recurrent %>% filter(id %in% id_samp) %>% select(id, Vmed9_everhadTB) %>% distinct(), by = "id")

ND = list(newdataL = ND_long, newdataE = tmp)

predLong1 <- predict(jm_fit_nb_cub,
                    newdata = ND,
                    return_newdata = TRUE)
plot(predLong1, outcomes = 1)
plot(predLong1, subject = "AD122", outcomes = 1)


predLong2 <- predict(jm_fit_nb_cub,
                    newdata = ND,
                    return_newdata = TRUE,
                    times = seq(t0, t1, length.out = 51))
plot(predLong2, outcomes = 1, subject = "AD261")

predSurv <- predict(jm_fit_nb_cub,
                    newdata = ND,
                    return_newdata = TRUE,
                    process = "event",
                    times = seq(t0, t1, length.out = 51))

unique(id_samp)
plot(predSurv, subject = "AD261")

plot(predLong2, predSurv, subject = "AD261")
plot(predLong2, predSurv, subject = "AD279")
plot(predLong2, predSurv, subject = "AD122")
plot(predLong2, predSurv, subject = "AD316")





#------------------------- joint modelling: DBS data -----------------------
##-----------------longitudinal model -----------------
#center continuous covariates for interpretation
dat_long_tfv$AGE_center = scale(dat_long_tfv$AGE, scale = FALSE)
dat_long_tfv$HCT_baseline_center = scale(dat_long_tfv$HCT_baseline, scale = FALSE)

#------functional form of time
fit_dbs_i = lme(sqrt(TFV) ~ time_yrs,
              random = ~ time_yrs|id,
              data = dat_long_tfv,
              method = "ML")
summary(fit_dbs_i)
plot(fit_dbs_i)

fit_dbs_ii = update(fit_dbs_i,
                    . ~ . + I(time_yrs^2))
summary(fit_dbs_ii)
anova(fit_dbs_i, fit_dbs_ii) #the quadratic term is necessary
plot(fit_dbs_ii)

#------- adjust for covariates
#AGE_center
fit_dbs_iii = lme(sqrt(TFV) ~ time_yrs + I(time_yrs^2) + TFV_baseline + AGE_center,
                  random = ~ time_yrs|id,
                  data = dat_long_tfv,
                  method = "ML")
summary(fit_dbs_iii)


fit_dbs_iv = update(fit_dbs_iii,
                    sqrt(TFV) ~ time_yrs * AGE_center)
summary(fit_dbs_iv) #interaction not necessary

#GENDER
fit_dbs_v = lme(sqrt(TFV) ~ time_yrs + TFV_baseline+ AGE_center + GENDER,
                random = ~ time_yrs|id,
                data = dat_long_tfv,
                method = "ML")
summary(fit_dbs_v) 

fit_dbs_vi = update(fit_dbs_v,
                    .~. + AGE_center*GENDER)
summary(fit_dbs_vi)
anova(fit_dbs_i, fit_dbs_iii, fit_dbs_v, fit_dbs_vi) 

# baseline haematocrit
fit_dbs_vii = update(fit_dbs_v,
                     .~. + HCT_baseline_center)
summary(fit_dbs_vii)

#random effect compare
fit_dbs_re1 =update(fit_dbs_vii,
                    method = "REML",
                    random = ~ 1|id)
fit_dbs_re2 =update(fit_dbs_re1,
                    random = ~ time_yrs|id)
anova(fit_dbs_re1, fit_dbs_re2) #slope trumps only intercept

fit_dbs = fit_dbs_re2
ggemmeans(tst, terms = 'time_yrs [all]') %>% plot()

###---------- some diagnostics -------------
plot(fit_dbs) #some heteroskedasticity - may need to model variance
qqnorm(resid(fit_dbs))
qqnorm(unlist(random.effects(fit_dbs)[1]), breaks = 20) #intercept RE's are not symmetric
qqnorm(unlist(random.effects(fit_dbs)[2]), breaks = 20) #slope RE's look nice n normal


##----------------- fit joint model -----------------
jm_dbs_fit_1 = jm(fit_surv_recurr, 
              fit_dbs, 
              time_var = "time_yrs",
              id_var = "id",
              recurrent = "gap",
              functional_forms = ~ value(sqrt(TFV)),
              control = list(n_iter = 10000, n_chains = 6))

summary(jm_dbs_fit_1)

traceplot(jm_dbs_fit_1)
densplot(jm_dbs_fit_1)





#------------------------- joint modelling: multivariate -----------------------
###-----------refit EAM model to full DBS dataset ------
fit_glmmadap_zinb_2 = update(fit_glmmadap_zinb,
                             data = dat_long_tfv)
summary(fit_glmmadap_zinb_2)
summary(fit_glmmadap_zinb) #this does change the estimates unfortunately

##----------------- fit joint model -----------------

###----------parameterisation: current value-------

long_models = list(fit_glmmadap_zinb_2, fit_dbs)
fForms_val <- list("n_heartbeats" = ~ value(n_heartbeats),
               "sqrt(TFV)" = ~ value(sqrt(TFV)))
jm_mvfit_1 = jm(fit_surv_recurr, 
              long_models, 
              time_var = "time_yrs",
              id_var = "id",
              recurrent = "gap",
              functional_forms = fForms_val)
summary(jm_mvfit_1)
traceplot(jm_mvfit_1)


###----------parameterisation: area -------
fForms_area <- list("n_heartbeats" = ~ area(n_heartbeats),
               "sqrt(TFV)" = ~ area(sqrt(TFV)))
jm_mvfit_2 = jm(fit_surv_recurr,
                long_models,
                time_var = "time_yrs",
                id_var = "id",
                recurrent = "gap",
                functional_forms = fForms_area)
summary(jm_mvfit_2)
traceplot(jm_mvfit_2)



###----------parameterisation: slope -------



###----------parameterisation: slope + value -------





###----------parameterisation: area -------



###----------parameterisation: value + area -------

