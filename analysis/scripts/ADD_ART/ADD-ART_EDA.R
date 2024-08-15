# ADD-ART Exploratory Analysis
# 05 March 2024
rm(list=ls())
# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(survival)
library(survminer)



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



#------------------ desc stats ------------------------





##----------------- outcomes --------------
###-------- adherence ---------
dat_long_eam %>% 
  summarise(EAM = quantile(n_heartbeats),
            SR = quantile(Vsra1_days_nonadherent, na.rm = T),
            TFV = quantile(TFV, na.rm = T))

#------------------- data visualization -------------------
setwd(path_out)
##--------------- cross-sectional ------
###------ univariate
#----adherence
#EAM missed
dat_long_eam %>%
  ggplot(aes(x = n_heartbeats_30days)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(x = "Number of EAM missed doses in last 30 days")
ggsave("WP_heartbeats_hist.pdf",
      height = 6, width = 8)
dat_long_eam %>%
  ggplot(aes(x = n_heartbeats_30days)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Number of EAM missed doses in last 30 days")
ggsave("WP_heartbeats_by_visit_hist.pdf",
      height = 6, width = 8)

#days between visits
# dat_long_eam %>% filter(Visit != 0) %>%
#   ggplot(aes(x = days_between_visits)) +
#   geom_histogram() + facet_wrap(~Visit) +
#   theme_bw() +
#   labs(x = "Days between visits")
# ggsave("days_between_visits_by_visit_hist.pdf",
#       height = 6, width = 8)
dat_long_eam %>%
  ggplot(aes(x = n_heartbeats_30days/wp_days_30)) +
  geom_histogram() + 
  theme_bw() +
  labs(x = "Proportion of EAM missed doses in last 30 days")
ggsave("WP_nonadherence_prop.pdf",
       height = 6, width = 8)

dat_long_eam %>%
  ggplot(aes(x = n_heartbeats_30days/wp_days_30)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Proportion of EAM missed doses in last 30 days")
ggsave("WP_nonadherence_prop_by_visit_hist.pdf",
      height = 6, width = 8)

# #days covered by wp
# dat %>% filter(Visit != 0) %>%
#   ggplot(aes(x = wp_days_covered)) +
#   geom_histogram() + facet_wrap(~Visit) +
#   theme_bw() +
#   labs(x = "Proportion of days between visits covered by EAM signals")
# ggsave("WP_days_covered_by_visit_hist.pdf",
#       height = 6, width = 8)

#---- TFV-DP
dat_long_eam %>%
  ggplot(aes(x = TFV)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "TFV-DP (fmol/punch)")
ggsave("TFV_hist.pdf",
       height = 6, width = 8)


dat_long_eam %>%
  ggplot(aes(x = TFV)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "TFV-DP (fmol/punch)")
ggsave("TFV_by_visit_hist.pdf",
      height = 6, width = 8)

#---- sqrt TFV
dat_long_eam %>%
  ggplot(aes(x = sqrt(TFV))) +
  geom_histogram() +
  theme_bw() +
  labs(x = "sqrt(TFV-DP) (fmol/punch)")
ggsave("sqrt_TFV_hist.pdf",
       height = 6, width = 8)

dat_long_eam %>% 
  ggplot(aes(x = sqrt(TFV))) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "sqrt(TFV-DP) (fmol/punch)")
ggsave("sqrt_TFV_by_visit_hist.pdf",
      height = 6, width = 8)

#---- FTC
# dat_long_eam %>%
#   ggplot(aes(x = FTC)) +
#   geom_histogram() + facet_wrap(~Visit) +
#   theme_bw() +
#   labs(x = "FTC-TP")
# ggsave("FTC_by_visit_hist.pdf",
#       height = 6, width = 8)

#---- HTC
dat_long_eam %>% 
  ggplot(aes(x = HCT)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Hematocrit (%)")
ggsave("HCT_by_visit_hist.pdf",
      height = 6, width = 8)


#---- self-reported adherence
dat_long_eam %>%
  ggplot(aes(x = Vsra1_days_nonadherent)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Self-reported number of missed doses in last 30 days")
ggsave("VSRA1_hist.pdf",
       height = 6, width = 8)

dat_long_eam %>%
  ggplot(aes(x = Vsra1_days_nonadherent)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported number of missed doses in last 30 days")
ggsave("VSRA1_by_visit_hist.pdf",
       height = 6, width = 8)


dat_long_eam %>%
  ggplot(aes(x = sra1_adherence_prop)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported proportion non-adherence in last 30 days")
ggsave("SRA1_by_visit_hist.pdf",
      height = 6, width = 8)

dat_long_eam %>%
  ggplot(aes(x = sra2_goodjob_prop)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported adherence 'good job' in last 30 days")
ggsave("SRA2_by_visit_hist.pdf",
      height = 6, width = 8)

dat_long_eam %>%
  ggplot(aes(x = sra3_way_supposed_to_prop)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported adherence 'manner' in last 30 days")
ggsave("SRA3_by_visit_hist.pdf",
      height = 6, width = 8)


#---- viral load
dat_long_eam %>%
  ggplot(aes(x = VL)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Viral load (copies/mL)")
ggsave("VL_by_visit_hist.pdf",
      height = 6, width = 8)
dat_long_eam %>%
  ggplot(aes(x = `log2(VL+0.001)`)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "log2(Viral load + 0.001)")
ggsave("log2_VL_by_visit_hist.pdf",
      height = 6, width = 8)

dat_long_eam %>%
  ggplot(aes(x = log2(VL))) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "log2(Viral load)")
ggsave("log2_VL_by_visit_hist_nozeros.pdf",
      height = 6, width = 8)


###------ bivariate -----
#adherence vs adherence pairs at baseline
dat_long_eam %>% 
  ggpairs(columns = c("wp_nonadherence_prop_30day", "sra_nonadherence_prop", "TFV"),
          title = "") + theme_bw()
ggsave("adherence_pairs_baseline.pdf",
      height = 8, width = 8)

#adherence vs adherence pairs at visit 1
dat_long_eam %>% filter(Visit == 0) %>%
  ggpairs(columns = c("wp_nonadherence_prop_30day", "sra_nonadherence_prop", "TFV"),
          title = "Adherence measures at visit 0")
ggsave("adherence_pairs_visit0.pdf",
      height = 8, width = 8)

#adherence vs vl pairs at baseline
dat_long_eam %>% filter(Visit == 0) %>%
  ggpairs(columns = c("wp_nonadherence_prop_30day", "sra_nonadherence_prop", "TFV", "log2(VL+0.001)"),
          title = "Adherence measures at visit 1")
ggsave("adherence_vl_pairs_visit0.pdf",
      height = 8, width = 8)


##--------------- longitudinal ------
###----- spaghetti plots -----
# EAM counts
dat_long_eam %>% 
  ggplot(aes(x = Visit,
             y = n_heartbeats_30days,
             group = id,
             color = id)) +
  geom_line() + 
  theme_bw() +
  theme(legend.position = "none")  +
  labs(y = "number of EAM missed doses in last 30 days")

ggsave("heartbeat_counts_spaghetti.pdf",
       height = 6, width = 8)

# (num WP intakes / num WP days) vs visit
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = wp_nonadherence_prop_30day,
             group = id,
             color = id)) +
         geom_line() + 
        theme_bw() +
         theme(legend.position = "none") +
  labs(y = "Proportion of EAM missed doses in last 30 days")
ggsave("wp_nonadherence_prop_spaghetti.pdf",
       height = 6, width = 8)


# TFV-DP vs visit
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = TFV,
             group = id,
             color = id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "TFV-DP (fmol/punch)")
ggsave("TFV_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#SRA1 vs visit
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = Vsra1_days_nonadherent,
             group = id,
             color = id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Self-reported missed doses in last 30 days")
ggsave("SRA1_by_visit_spaghetti.pdf",
       height = 6, width = 8)


# # FTC vs visit
# dat_long_eam %>%
#   ggplot(aes(x = Visit,
#              y = FTC,
#              group = id,
#              color = id)) +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(y = "FTC-TP (fmol/punch)")
# ggsave("FTC_by_visit_spaghetti.pdf",
#        height = 6, width = 8)
# 
# #FTC vs visit (no 12.5)
# dat_long_eam %>% filter(FTC != 12.5) %>%
#   ggplot(aes(x = Visit,
#              y = FTC,
#              group = id,
#              color = id)) +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(y = "FTC-TP (fmol/punch)")
# ggsave("FTC_by_visit_spaghetti_no12.5.pdf",
#        height = 6, width = 8)

# HCT vs visit
dat_long_eam %>% filter(Visit != 0 & Visit != 2 & Visit != 4 & Visit != 6 & Visit != 8 & Visit != 10) %>%
  ggplot(aes(x = Visit,
             y = HCT,
             group = id,
             color = id)) +
  geom_line() + 
  theme_bw() +
  theme(legend.position = "none")  +
  labs(y = "Hematocrit (%)")
ggsave("HCT_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#VL vs visit
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = VL,
             group = id,
             color = id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Viral load (copies/mL)")
ggsave("VL_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#log2(VL+0.00001)
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = `log2(VL+0.001)`,
             group = id,
             color = id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "log2(Viral load + 0.001)")
ggsave("log2_VL_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#log2(VL) - no undetectable VLs
dat_long_eam %>%
  ggplot(aes(x = Visit,
             y = log2(VL),
             group = id,
             color = id)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "log2(Viral load)")
ggsave("log2_VL_by_visit_spaghetti_nozeros.pdf",
       height = 6, width = 8)



###----- summarized profiles ----
dat_adherence_summarised = dat_long_eam %>%
  group_by(Visit) %>%
  summarise(
    mean_n_heartbeats_30days = mean(n_heartbeats_30days, na.rm = TRUE),
    sd_n_heartbeats_30days = sd(n_heartbeats_30days, na.rm = TRUE),
    mean_vsra1 = mean(Vsra1_days_nonadherent, na.rm = TRUE),
    sd_vsra1 = sd(Vsra1_days_nonadherent, na.rm = TRUE),
    mean_TFV = mean(TFV, na.rm = TRUE),
    sd_TFV = sd(TFV, na.rm = TRUE)
  )

dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_n_heartbeats_30days,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmax(0, mean_n_heartbeats_30days - sd_n_heartbeats_30days),
                  ymax = pmax(0, mean_n_heartbeats_30days + sd_n_heartbeats_30days)),
              alpha = 0.3) +
  theme_bw() +
  labs(y = "Mean EAM missed doses in last 30 days")
ggsave("meanprof_n_heartbeats_30days.pdf",
       height = 6, width = 8)

dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_vsra1,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmax(0, mean_vsra1 - sd_vsra1),
                  ymax = mean_vsra1 + sd_vsra1),
              alpha = 0.3) +
  theme_bw() +
  labs(y = "Mean self-reported missed doses in last 30 days")
ggsave("meanprof_vsra1.pdf",
       height = 6, width = 8)

dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_TFV,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmax(0, mean_TFV - sd_TFV),
                  ymax = mean_TFV + sd_TFV),
              alpha = 0.3) +
  theme_bw() +
  labs(y = "Mean TFV-DP (fmol/punch)")
ggsave("meanprof_TFV.pdf",
       height = 6, width = 8)

#variance profiles
dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = sd_n_heartbeats_30days,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(y = "Sample standard deviation of EAM missed doses in last 30 days")
ggsave("sdprof_n_heartbeats_30days.pdf",
       height = 6, width = 8)

dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = sd_vsra1,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(y = "Sample standard deviation of self-reported missed doses in last 30 days")
ggsave("sdprof_vsra1.pdf",
       height = 6, width = 8)

dat_adherence_summarised %>%
  ggplot(aes(x = Visit,
             y = sd_TFV,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(y = "Sample standard deviation of TFV-DP (fmol/punch)")
ggsave("sdprof_TFV.pdf",
       height = 6, width = 8)

##---------------- Survival profiles -----------


###----------------- kaplan meier ------------------
setwd(path_out)
#visit time
Surv3a = with(dat_recurrent, 
             Surv(time = visit_start, time2 = Visit,
                  event = VL> 50, type = "counting"))
km_fit3a = survfit(Surv3a ~ 1, data = dat_recurrent)
ggsurvplot(km_fit3a, data = dat_recurrent, 
           risk.table = T,
           cumevents = T,
           conf.type = "log-log")


#time from visit 0
Surv3b = with(dat_recurrent, 
              Surv(time = start, time2 = time_yrs,
                   event = VL> 50, type = "counting"))
km_fit3b = survfit(Surv3b ~ 1, data = dat_recurrent)
ggsurvplot(km_fit3b, data = dat_recurrent, 
           risk.table = T,
           cumevents = T,
           conf.type = "log-log")




###----------------- time to first non-suppression, stratified by baseline covariates ---

#----- adherence baseline
km_fit_ad1 = survfit(Surv3b ~ poor_adherence_TFV_baseline, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_ad2 = survfit(Surv3b ~ poor_adherence_SR_baseline, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad2, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

#vl baseline
km_fit_ad3 = survfit(Surv3b ~ unsuppressed_VL_baseline, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad3, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

#----- adherence summaries
km_fit_ad4 = survfit(Surv3b ~ poor_adherence_SR, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad4, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_ad5 = survfit(Surv3b ~ poor_adherence_TFV, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad5, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_ad6 = survfit(Surv3b ~ poor_adherence_eam, 
                     data = dat_recurrent)
ggsurvplot(km_fit_ad6, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")


#----- other baseline covariates
colnames(dat_recurrent)
km_fit_cov1 = survfit(Surv3b ~ Vmed9b_TBstatus,
                      data = dat_recurrent)

km_fit_cov1 = survfit(Surv3b ~ Vmed9_everhadTB,
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_cov1 = survfit(Surv3b ~ Vmed7_HIVillness_ever,
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_cov1 = survfit(Surv3b ~ Vmed5_rate_health >3,
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_cov1 = survfit(Surv3b ~ GENDER,
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_cov1 = survfit(Surv3b ~ AGE > mean(AGE),
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")

km_fit_cov1 = survfit(Surv3b ~ HCT_baseline > mean(HCT_baseline, na.rm = TRUE),
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")


dat_recurrent = dat_recurrent %>%
  mutate(BMI_cat = case_when(
    BMIcat == "< 18.5" ~ "< 25",
    BMIcat == "18.5-24.9" ~ "< 25",
    BMIcat == "25-29.9" ~ "25 - 29.9",
    BMIcat == "30 or +" ~ "> 30"
  ))
table(dat_recurrent$BMIcat)
table(dat_recurrent$BMI_cat)
dat_recurrent$BMI_cat = factor(dat_recurrent$BMI_cat,
                                 levels = c("< 25", "25 - 29.9", "> 30"))

km_fit_cov1 = survfit(Surv3b ~ BMI_cat,
                      data = dat_recurrent)
ggsurvplot(km_fit_cov1, data = dat_recurrent,
           conf.type = "log-log",
           risk.table = "nrisk_cumevents")
