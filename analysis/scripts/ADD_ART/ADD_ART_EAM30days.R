#ADD-ART EAM Modelling: 30 day WP counts
#in this script, we use a GLMM to model the number of missed WP doses in the past 30 days, over time
rm(list=ls())
#dependencies
library(tidyverse)
library(GLMMadaptive)
library(ggeffects)
source("~/R_custom_functions/resids_mixed_models.R")
source("~/R_custom_functions/modelfit_metrics.R")
#args
path_models = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART/models"
path_data = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"


#load data
setwd(path_data)
dat_long_eam = read.csv("ADDART_long_eam_14May2024.csv")


dat_long_eam = dat_long_eam %>%
  mutate(log_wp_days_30 = log(wp_days_30 + 1)) %>%
  relocate(log_wp_days_30, .after = n_none_30days)
off = dat_long_eam$log_wp_days_30


#------------------------- EDA -----------------------
#quick EDA to look at distribution of n_heartbeats, mean trajectory over time, variance profile over time
dat_long_eam %>%
  ggplot(aes(x=n_heartbeats_30days)) +
  geom_histogram(bins = 30)

dat_long_eam %>%
  ggplot(aes(x=n_heartbeats_30days)) +
  geom_histogram() +
  facet_wrap(~Visit)

dat_long_eam %>%
  ggplot(aes(x = Visit, y = n_heartbeats_30days, group=id)) +
  geom_line()

#mean profiles
dat_long_eam %>%
  ggplot(aes(x = time_yrs, y = n_heartbeats_30days)) +
  geom_point() +
  geom_smooth(method = "lm")

dat_long_eam %>%
  ggplot(aes(x = Visit, y = n_heartbeats_30days)) +
  stat_summary(fun.y = "mean", geom = "line")


dat_long_eam %>%
  ggplot(aes(x = Visit, y = n_heartbeats_30days/wp_days_30)) +
  stat_summary(fun.y = "mean", geom = "line")


#variance profiles
dat_long_eam %>%
  ggplot(aes(x = Visit, y = n_heartbeats_30days)) +
  stat_summary(fun.y = "sd", geom="line")




#----------------- Modelling -------------------
##---- POISSON --------

###---- fixed effects -------
fit_pois_null = mixed_model(
  fixed = n_heartbeats_30days ~ 1 + offset(off),
  random = ~ 1|id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_null)

#add time
fit_pois_1 = mixed_model(
  fixed = n_heartbeats_30days ~  offset(off) + time_yrs,
  random = ~ 1|id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_1)

fit_pois_2 = mixed_model(
  fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + offset(off),
  random = ~ 1 | id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_2)
ggemmeans(fit_pois_2, terms = "time_yrs[all]") %>% plot()
plot(simulateResiduals(fit_pois_2))

fit_pois_3 = mixed_model(
  fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + I(time_yrs^3) + offset(off),
  random = ~ 1 | id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_3)
ggemmeans(fit_pois_3, terms = "time_yrs[all]") %>% plot()
plot(simulateResiduals(fit_pois_3))


compare_modelfits(list("null"=fit_pois_null, "linear"=fit_pois_1, "quad"=fit_pois_2, "cube"=fit_pois_3))
anova(fit_pois_null, fit_pois_1)
anova(fit_pois_1, fit_pois_2)
anova(fit_pois_2, fit_pois_3)


###----- random effects (independent) ------
fit_pois_2a = mixed_model(
  fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + offset(off),
  random = ~ time_yrs || id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_2a)
anova(fit_pois_2, fit_pois_2a) # big improvement in fit
plot(simulateResiduals(fit_pois_2))
plot(simulateResiduals(fit_pois_2a)) # big improvement in resids, still some overdispersion


fit_pois_3a = mixed_model(
  fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + I(time_yrs^3) + offset(off),
  random = ~ time_yrs || id,
  family = poisson,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_pois_3a)
anova(fit_pois_3, fit_pois_3a)
plot(simulateResiduals(fit_pois_3))
plot(simulateResiduals(fit_pois_3a))

#### NOTE: we need random slopes term, whether to include cubic time is less clear. 
# will work with quadratic FE for now. cubic may be overparameterized and create issues with JM estimation down the line


### relax independence of REs
fit_pois_2b = update(fit_pois_2a,
                     random = ~ time_yrs | id)
anova(fit_pois_2, fit_pois_2a)
anova(fit_pois_2a, fit_pois_2b) #null rejected, so we need to relax independence

fit_pois_3b = update(fit_pois_3a,
                     random = ~ time_yrs | id)
anova(fit_pois_3, fit_pois_3a)
anova(fit_pois_3a, fit_pois_3b) #null rejected, so we need to relax independence

# settle on final poisson model
fit_pois_final = fit_pois_2b
summary(fit_pois_final)
ggemmeans(fit_pois_final, terms = "time_yrs[all]") %>% plot()
ggemmeans(fit_pois_3b, terms = "time_yrs[all]") %>% plot()

plot(simulateResiduals(fit_pois_final), quantreg=T) #lots of outliers, some heteroskedasticity
plot(simulateResiduals(fit_pois_3b), quantreg=T)

##---- NEG BINOM --------
# we have some overdispersion in the residuals, so we will try a negative binomial model
#fit model with same specs but family = negbin
fit_nb_1 = update(fit_pois_final, 
                  family = GLMMadaptive::negative.binomial)
summary(fit_nb_1)
ggemmeans(fit_nb_1, terms = "time_yrs[all]") %>% plot()
plot(simulateResiduals(fit_nb_1), quantreg=T) #looks to be a little underdispersion now, not too concerning

#check independent REs
fit_nb_1a = update(fit_nb_1,
                   random = ~ time_yrs || id)
anova(fit_nb_1a, fit_nb_1) #rejected, relax independence

#check the same model with cubic time
fit_nb_2 = mixed_model(
  fixed = n_heartbeats_30days ~ time_yrs + I(time_yrs^2) + I(time_yrs^3) + offset(off),
  random = ~ time_yrs | id,
  family = GLMMadaptive::negative.binomial,
  data = dat_long_eam,
  na.action = na.omit
)
summary(fit_nb_2)
ggemmeans(fit_nb_2, terms = "time_yrs[all]") %>% plot()
plot(simulateResiduals(fit_nb_2), quantreg=T) #looks to be a little underdispersion now, not too concerning

# looks a bit better with cubic time ito heteroskedasticity, but not a huge difference
#there is not much else to do but save the final models

setwd(path_models)
saveRDS(fit_pois_final, "EAM30day_pois_quad_final.rds")
saveRDS(fit_pois_3b, "EAM30day_pois_cubic_final.rds")
saveRDS(fit_nb_1, "EAM30day_nb_quad_final.rds")
saveRDS(fit_nb_2, "EAM30day_nb_cubic_final.rds")


