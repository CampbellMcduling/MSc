#ADD-ART Longitudinal modelling
# 14 March 2024
rm(list=ls())
#---- dependencies
library(tidyverse)
library(DHARMa)
library(lme4)
library(glmmTMB)
library(ggeffects)

#---------------------- custom functions ----------------------
resids = function(model){

  fixed = fixef(model)
  ranef = ranef(model)
  resid = getResiduals(model)
  
  disp_test = testDispersion(model)
  zero_test = testZeroInflation(model)
  print(disp_test)
  print(zero_test)
  
  fit_v_resid = plot(fitted(model), resid)
  
  #extract covariates and plot against residuals
  dharma = simulateResiduals(model)
  covariates = model.matrix(model)

  resid_v_pred = list()
  fit_v_cov = list()
  for(i in 2:ncol(covariates)){
    resid_v_pred[i] = plotResiduals(dharma, form = covariates[,i], xlab = colnames(covariates)[i], quantreg = T)
    #fit_v_cov[i] = plot(fitted(model) ~ covariates[,i], xlab = colnames(covariates)[i])
  }
  
  dharma = plot(dharma, quantreg=T)
  

  
  return(list(fit_v_resid = fit_v_resid, resid_v_pred = resid_v_pred, dharma = dharma))
  
}

compare_modelfits = function(models){
  #models is a list of models, return common model fit criteria for each of model in model in a table
  metrics = data.frame(AIC = NA, BIC = NA, loglik=NA)
  for(i in 1:length(models)){
    model = models[[i]]
    AIC = AIC(model)
    BIC = BIC(model)
    loglik = logLik(model)
    metrics[i,] = c(AIC, BIC, loglik)
  }
  rownames(metrics) = names(models)
  metrics = metrics %>% arrange(BIC)
  
  return(metrics)
}


path_in = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"

#---- load data
setwd(path_in)
dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)

dat$wp_days = dat$n_intakes + dat$n_heartbeats
dat$log_wp_days = log(dat$wp_days+1)
#create non-adherence proportion
dat = dat %>%
  mutate(wp_nonadherence_prop = n_heartbeats/wp_days) %>%
  relocate(wp_nonadherence_prop, .after = wp_adherence_prop)
dat=dat %>% relocate(wp_days, .after = n_heartbeats)

dat = dat %>% filter(Visit != 0)
#view missing n_intakes
# dat %>% filter(is.na(n_intakes)) %>% View()

#censor EAM data: if missed visits, then EAM data is missing
dat_censored <- dat %>%
  mutate(across(
    c("n_intakes", "n_heartbeats", "n_none", "wp_days", "wp_adherence_prop", "wp_nonadherence_prop"),
    ~ case_when(
      missed_visit == 1 ~ NA_real_,
      TRUE ~ .
    )
  ))

dat_eam = read.csv("ADDART_long_eam_14May2024.csv")
dat_eam$wp_days = dat_eam$n_intakes + dat_eam$n_heartbeats
dat_eam$log_wp_days = log(dat_eam$wp_days+1)
offset = dat_eam$log_wp_days


#--------------------- beta GLMM w/ censored data -----------------------
#this is not an appropriate way distribution to model this data with, 
# since it assumes continuous proportions, not proportions of the type k/n
dat_censored2 <- dat_censored %>%
  mutate(wp_adherence_prop = case_when(
    wp_adherence_prop == 0 ~ 0.0001,
    wp_adherence_prop == 1 ~ 0.9999,
    TRUE ~ wp_adherence_prop
  ),
  wp_nonadherence_prop = case_when(
    wp_nonadherence_prop == 0 ~ 0.0001,
    wp_nonadherence_prop == 1 ~ 0.9999,
    TRUE ~ wp_nonadherence_prop
  ))
fit_beta_1a = glmmTMB(
  formula = wp_adherence_prop ~ Visit + I(Visit^2) + (1|PID),
  family = beta_family,
  data = dat_censored2
)
summary(fit_beta_1a)
resids(fit_beta_1a)

fit_beta_1b = glmmTMB(
  formula = wp_nonadherence_prop ~ Visit + I(Visit^2) + (1|PID),
  family = beta_family,
  data = dat_censored2
)
summary(fit_beta_1b)
resids(fit_beta_1b)

#--------------------- binomial GLMM w/ censored data (logistic) -----------------------
#binom
fit_bin_1a = glmmTMB(
  formula = cbind(n_heartbeats, wp_days) ~ Visit + I(Visit^2) + (1|PID),
  family = binomial,
  data = dat_censored
)
summary(fit_bin_1a)
resids(fit_bin_1a) #zero-inflated and heteroskedastic, possibly underdispersed

#beta binom
fit_betabin_1a = glmmTMB(
  formula = cbind(n_heartbeats, wp_days) ~ Visit + I(Visit^2) + (1|PID),
  family = betabinomial,
  data = dat_censored
)
summary(fit_betabin_1a) #fixed effect parameter estimates are larger in magnitude and uncertainty with betabin, random effects slightly smaller
resids(fit_betabin_1a) #dharma residuals seem much reasonable in betabin; still heteroskedastic with unexpected dispersion in upper tail
compare_modelfits(list(fit_bin_1a = fit_bin_1a, fit_betabin_1a = fit_betabin_1a)) #model fit indicates betabin is better

#model dispersion as a function of fixed effect
fit_betabin_1b = update(fit_betabin_1a,
                       dispformula = ~ Visit + I(Visit^2))
summary(fit_betabin_1b)
resids(fit_betabin_1b)

#try model non-adherence with zero inflation, without modeling dispersion. This should improve both underdispersion and zero-inflation if correct
fit_betabin_1c = update(fit_betabin_1a,
                        ziformula = ~ 1)
summary(fit_betabin_1c)
compare_modelfits(list(fit_betabin_1a = fit_betabin_1a, 
                       fit_betabin_1b = fit_betabin_1b, 
                       fit_betabin_1c = fit_betabin_1c))
resids(fit_betabin_1c) #no observable improvement, modelfits indicate its worse

fit_betabin_1c_ii = update(fit_betabin_1c,
                           ziformula = ~ Visit)
summary(fit_betabin_1c_ii)
resids(fit_betabin_1c_ii)
compare_modelfits(list(fit_betabin_1a = fit_betabin_1a, 
                       fit_betabin_1b = fit_betabin_1b, 
                       fit_betabin_1c = fit_betabin_1c,
                       fit_betabin_1c_ii = fit_betabin_1c_ii)) #this is closer to the dispersion one

#control dispersion
fit_betabin_1d = update(fit_betabin_1c,
                        dispformula = ~ Visit + I(Visit^2))
summary(fit_betabin_1d)
resids(fit_betabin_1d)
compare_modelfits(list(fit_betabin_1a = fit_betabin_1a, 
                       fit_betabin_1b = fit_betabin_1b, 
                       fit_betabin_1c = fit_betabin_1c,
                       fit_betabin_1c_ii = fit_betabin_1c_ii,
                       fit_betabin_1d = fit_betabin_1d))

#--------------------- POISSON GLMM w/ censored data -----------------------
##---------- determine fixed effects
# functional form of time considered up to cubic
fit_1_null = glmmTMB(
  formula = n_heartbeats ~ 1 + offset(offset) ,
  family = poisson(link = "log"),
  data = dat_eam
)
summary(fit_1_null)

fit_1a = glmmTMB(
  formula = n_heartbeats ~ Visit + offset(log_wp_days) + (1|id),
  family = poisson(link = "log"),
  data = dat_eam
)
summary(fit_1a)
resids(fit_1a)

fit_1b = update(fit_1a, 
                formula = n_heartbeats ~ Visit + I(Visit^2) + offset(log_wp_days) + (1|id))
summary(fit_1b)
resids(fit_1b)

fit_1c = update(fit_1b,
                formula = n_heartbeats ~ Visit + I(Visit^2) + I(Visit^3) + offset(log_wp_days) + (1|id))
summary(fit_1c)

models <- list(fit_1_null = fit_1_null, fit_1a = fit_1a, fit_1b = fit_1b, fit_1c = fit_1c)
compare_modelfits(models)
resids(fit_1c) #these models all have significant zero-inflation and are underdispersed with heteroskedasticity


#----------- adjusting for zero inflation
fit_1d_i = update(fit_1c,
                ziformula = ~ 1 )
fit_1d_ii = update(fit_1c,
                  ziformula = ~ Visit)
summary(fit_1d_i)
summary(fit_1d_ii)

resids(fit_1d_i)
resids(fit_1d_ii)


#----------  try account for dispersion with negative binomial model
fit_1e_i <- update(fit_1d_ii, family = nbinom1(link = "log"))
summary(fit_1e_i)
resids(fit_1e_i)

fit_1e_ii <- update(fit_1e_i, dispformula = ~ Visit )
summary(fit_1e_ii)
resids(fit_1e_ii)

compare_modelfits(list(
  fit_1a = fit_1a,
  fit_1b = fit_1b,
  fit_1c = fit_1c,
  fit_1d_i = fit_1d_i,
  fit_1d_ii = fit_1d_ii,
  fit_1e_i = fit_1e_i,
  fit_1e_ii = fit_1e_ii
)
)

#--------- try best model with random slope
fit_1e_iii = update(fit_1e_ii,
                    formula = n_heartbeats ~ Visit + I(Visit^2) + I(Visit^3)  + offset(log_wp_days) + (Visit|id),
                    ziformula = ~ 1,
                    dispformula = ~ Visit)

summary(fit_1e_iii)
resids(fit_1e_iii)

models <- list(fit_1e_i = fit_1e_i,
               fit_1e_ii = fit_1e_ii, 
               fit_1e_iii = fit_1e_iii)
compare_modelfits(models)
#--------- trying generalized poisson
#trying with generalized poisson - with the caveat that I have no idea what that is
tst = update(fit_1c,
             family = genpois)
summary(tst)
resids(tst)

#try control dispersion of gen poisson
tst_1 = update(tst,
               ziformula = ~ Visit)
summary(tst_1)
resids(tst_1)

#---------- compare residuals of each of the best models under each parametric assumption
models = list(
  fit_betabin_1d = fit_betabin_1d,
  fit_1e_iii = fit_1e_iii,
  tst_1 = tst_1
)
compare_modelfits(models)

#-------------------refit model with glmmadaptive
summary(fit_1e_iii)
fit_1e_iii_adapt = mixed_model(
  fixed = n_heartbeats ~ Visit + I(Visit^2) + I(Visit^3) + offset(log_wp_days),
  random = ~ Visit | id,
  family = "poisson",
  data = dat_eam
  
)

#--------------- write objects to file --------------
setwd(paste0(path_out, "/Models"))
saveRDS(fit_1e_iii, "ADD-ART_glmm.rds")
