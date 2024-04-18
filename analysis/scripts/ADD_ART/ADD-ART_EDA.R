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



path_in = "~/Downloads/MSc/Dissertation/Analysis/Input/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/Analysis/Output/ADD-ART/Figures"

# Load data
setwd(path_in)

dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)

#------------------- data wrangling -------------------
dat$Visit = as.factor(dat$Visit)
dat$`log2(VL+0.001)` = log2(dat$VL + 0.001)





#------------------- data visualization -------------------
setwd(path_out)
##--------------- cross-sectional ------
###------ univariate
#----adherence
#WP counts
dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = n_intakes)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Number of EAM intakes between visits")
ggsave("WP_intakes_by_visit_hist.pdf",
      height = 6, width = 8)

#days between visits
dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = days_between_visits)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Days between visits")
ggsave("days_between_visits_by_visit_hist.pdf",
      height = 6, width = 8)


dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = wp_adherence_prop)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Proportion of possible EAM intakes between visits")
ggsave("WP_adherence_prop_by_visit_hist.pdf",
      height = 6, width = 8)

#days covered by wp
dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = wp_days_covered)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Proportion of days between visits covered by EAM signals")
ggsave("WP_days_covered_by_visit_hist.pdf",
      height = 6, width = 8)

#---- TFV-DP
dat %>%
  ggplot(aes(x = TFV)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "TFV-DP (fmol/punch)")
ggsave("TFV_by_visit_hist.pdf",
      height = 6, width = 8)

#---- FTC
dat %>%
  ggplot(aes(x = FTC)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "FTC-TP")
ggsave("FTC_by_visit_hist.pdf",
      height = 6, width = 8)

#---- HTC
dat %>% filter(Visit != 1 & Visit != 3 & Visit != 5 & Visit != 7 & Visit != 9 & Visit != 11) %>% # remove visits with no HTC
  ggplot(aes(x = HCT)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Hematocrit (%)")
ggsave("HCT_by_visit_hist.pdf",
      height = 6, width = 8)


#---- self-reported adherence
dat %>%
  ggplot(aes(x = sra1_adherence_prop)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported adherence in last 30 days")
ggsave("SRA1_by_visit_hist.pdf",
      height = 6, width = 8)

dat %>%
  ggplot(aes(x = Vsra2_good_job)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported adherence 'good job' in last 30 days")
ggsave("SRA2_by_visit_hist.pdf",
      height = 6, width = 8)

dat %>%
  ggplot(aes(x = Vsra3_way_supposed_to)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Self-reported adherence 'manner' in last 30 days")
ggsave("SRA3_by_visit_hist.pdf",
      height = 6, width = 8)


#---- viral load
dat %>%
  ggplot(aes(x = VL)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "Viral load (copies/mL)")
ggsave("VL_by_visit_hist.pdf",
      height = 6, width = 8)
dat %>%
  ggplot(aes(x = `log2(VL+0.001)`)) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "log2(Viral load + 0.001)")
ggsave("log2_VL_by_visit_hist.pdf",
      height = 6, width = 8)

dat %>%
  ggplot(aes(x = log2(VL))) +
  geom_histogram() + facet_wrap(~Visit) +
  theme_bw() +
  labs(x = "log2(Viral load)")
ggsave("log2_VL_by_visit_hist_nozeros.pdf",
      height = 6, width = 8)


###------ bivariate
#adherence vs adherence pairs at baseline
dat %>% filter(Visit == 0 & FTC != 12.5) %>%
  ggpairs(columns = c("TFV", "FTC", "sra1_adherence_prop"),
          title = "Adherence measures at baseline")
ggsave("adherence_pairs_baseline_noFTC12.pdf",
      height = 6, width = 8)

#adherence vs adherence pairs at visit 1
dat %>% filter(Visit == 1 & FTC != 12.5) %>%
  ggpairs(columns = c("TFV", "FTC",  "wp_adherence_prop","sra1_adherence_prop"),
          title = "Adherence measures at visit 1")
ggsave("adherence_pairs_visit1_noFTC12.pdf",
      height = 6, width = 8)

# adherence vs viral load vs HTC at baseline
dat %>% filter(Visit == 0 & FTC != 12.5) %>%
  ggpairs(columns = c("TFV", "FTC", "log2(VL+0.001)", "HCT"),
          title = "Adherence, viral load, and HTC at baseline")

# adherence vs viral load vs HTC at visit 2
dat %>% filter(Visit == 2 & FTC != 12.5) %>%
  ggpairs(columns = c("TFV", "FTC", "wp_adherence_prop", "VL", "HCT"),
          title = "Adherence, viral load, and HTC at visit 2")




##--------------- longitudinal ------
###----- spaghetti plots -----
# intake counts
dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = Visit,
             y = n_intakes,
             group = PID,
             color = PID)) +
  geom_line() + 
  theme_bw() +
  theme(legend.position = "none")  +
  labs(y = "number of EAM intakes between visits")

ggsave("intake_counts_by_visit_spaghetti.pdf",
       height = 6, width = 8)

# (num WP intakes / num WP days) vs visit
dat %>% filter(Visit != 0) %>%
  ggplot(aes(x = Visit,
             y = wp_adherence_prop,
             group = PID,
             color = PID)) +
         geom_line() + 
        theme_bw() +
         theme(legend.position = "none") +
  labs(y = "Proportion of possible EAM intakes between visits")
ggsave("wp_adherence_by_visit_spaghetti.pdf",
       height = 6, width = 8)


# TFV-DP vs visit
dat %>%
  ggplot(aes(x = Visit,
             y = TFV,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "TFV-DP (fmol/punch)")
ggsave("TFV_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#SRA1 vs visit
dat %>%
  ggplot(aes(x = Visit,
             y = sra1_adherence_prop,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Self-reported adherence (proportion)")
ggsave("SRA1_by_visit_spaghetti.pdf",
       height = 6, width = 8)


# FTC vs visit
dat %>%
  ggplot(aes(x = Visit,
             y = FTC,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "FTC-TP (fmol/punch)")
ggsave("FTC_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#FTC vs visit (no 12.5)
dat %>% filter(FTC != 12.5) %>%
  ggplot(aes(x = Visit,
             y = FTC,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "FTC-TP (fmol/punch)")
ggsave("FTC_by_visit_spaghetti_no12.5.pdf",
       height = 6, width = 8)

# HCT vs visit
dat %>% filter(Visit != 1 & Visit != 3 & Visit != 5 & Visit != 7 & Visit != 9 & Visit != 11) %>%
  ggplot(aes(x = Visit,
             y = HCT,
             group = PID,
             color = PID)) +
  geom_line() + 
  theme_bw() +
  theme(legend.position = "none")  +
  labs(y = "Hematocrit (%)")
ggsave("HCT_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#VL vs visit
dat %>%
  ggplot(aes(x = Visit,
             y = VL,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Viral load (copies/mL)")
ggsave("VL_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#log2(VL+0.00001)
dat %>%
  ggplot(aes(x = Visit,
             y = `log2(VL+0.001)`,
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "log2(Viral load + 0.001)")
ggsave("log2_VL_by_visit_spaghetti.pdf",
       height = 6, width = 8)

#log2(VL) - no undetectable VLs
dat %>%
  ggplot(aes(x = Visit,
             y = log2(VL),
             group = PID,
             color = PID)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "log2(Viral load)")
ggsave("log2_VL_by_visit_spaghetti_nozeros.pdf",
       height = 6, width = 8)



###----- mean longitudinal profiles ----
# summarise data by visit, note FTC and VL use measures robust to the large 12.5 values
dat_long_summarised = dat %>%
  group_by(Visit) %>%
  summarise(mean_wp_adherence_prop = mean(wp_adherence_prop, na.rm = T),
            var_wp_adherence_prop = var(wp_adherence_prop, na.rm = T),
            mean_TFV = mean(TFV, na.rm = T),
            var_TFV = var(TFV, na.rm = T),
            median_FTC = median(FTC, na.rm = T),
            mad_FTC = mad(FTC, na.rm = T),
            mean_SRA1 = mean(sra1_adherence_prop, na.rm = T),
            var_SRA1 = var(sra1_adherence_prop, na.rm = T),
            mean_HCT = mean(HCT, na.rm = T),
            var_HCT = var(HCT, na.rm = T))
dat_long_summarised = as.data.frame(dat_long_summarised)
#replace NaN with NA
dat_long_summarised <- dat_long_summarised %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .)))

# WP
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_wp_adherence_prop,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_wp_adherence_prop - sqrt(var_wp_adherence_prop),
                  ymax = mean_wp_adherence_prop + sqrt(var_wp_adherence_prop)),
              alpha = 0.2) +
  theme_bw() +
  labs(x = "Visit",
       y = "Mean EAM adherence")
ggsave("wp_adherence_by_visit_meanprof.pdf",
       height = 6, width = 8)

#TFV
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_TFV,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_TFV - sqrt(var_TFV),
                  ymax = mean_TFV + sqrt(var_TFV)),
              alpha = 0.2) +
  theme_bw() +
  labs(x = "Visit",
       y = "Mean TFV-DP (fmol/punch)")
ggsave("TFV_by_visit_meanprof.pdf",
       height = 6, width = 8)

#FTC
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = median_FTC,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = median_FTC - mad_FTC,
                  ymax = median_FTC + mad_FTC),
              alpha = 0.2) +
  theme_bw() +
  labs(x = "Visit",
       y = "Median FTC-DP")
ggsave("FTC_by_visit_meanprof.pdf",
       height = 6, width = 8)

#SRA1
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = mean_SRA1,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_SRA1 - sqrt(var_SRA1),
                  ymax = mean_SRA1 + sqrt(var_SRA1)),
              alpha = 0.2) +
  theme_bw() +
  labs(x = "Visit",
       y = "Mean Self-Reported adherence proportion")
ggsave("SRA1_by_visit_meanprof.pdf",
       height = 6, width = 8)

#HCT
dat_long_summarised %>% filter(Visit != 1 & Visit != 3 & 
                                 Visit != 5 & Visit != 7 & Visit != 9 & Visit != 11) %>%
  ggplot(aes(x = Visit,
             y = mean_HCT,
             group = 1)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_HCT - sqrt(var_HCT),
                  ymax = mean_HCT + sqrt(var_HCT)),
              alpha = 0.2) +
  theme_bw() +
  labs(x = "Visit",
       y = "Mean HCT (%)")
ggsave("HCT_by_visit_meanprof.pdf",
       height = 6, width = 8)

### ----- variance profiles -----
# WP
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = var_wp_adherence_prop,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Visit",
       y = "WP adherence Sample Variance")
ggsave("wp_adherence_by_visit_varprof.pdf",
       height = 6, width = 8)


# TFV
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = var_TFV,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Visit",
       y = "TFV-DP Sample Variance")
ggsave("TFV_by_visit_varprof.pdf",
       height = 6, width = 8)

# FTC
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = mad_FTC,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Visit",
       y = "MAD of FTC-DP")
ggsave("FTC_by_visit_varprof.pdf",
       height = 6, width = 8)

#SRA1
dat_long_summarised %>%
  ggplot(aes(x = Visit,
             y = var_SRA1,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Visit",
       y = "Self-Reported Adherence Sample Variance")
ggsave("SRA1_by_visit_varprof.pdf",
       height = 6, width = 8)

#HCT
dat_long_summarised %>% filter(Visit != 1 & Visit != 3 & 
                                 Visit != 5 & Visit != 7 & Visit != 9 & Visit != 11) %>%
  ggplot(aes(x = Visit,
             y = var_HCT,
             group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Visit",
       y = "HCT Sample Variance")
ggsave("HCT_by_visit_varprof.pdf",
       height = 6, width = 8)
##---------------- Survival profiles -----------
# time to first viral load > 1000
dat$Visit = as.numeric(dat$Visit)
dat$VBstatus = ifelse(dat$VL > 1000, 1, 0)

Surv = with(dat, Surv(time = Visit, event = VBstatus))
km_fit1a = survfit(Surv ~ 1, data = dat)
ggsurvplot(km_fit1a, data = dat, conf.type = "log-log")


#time to first non-suppression event
dat$VL_supp_status = ifelse(dat$VL > 50, 1, 0)
Surv2 = with(dat, Surv(time = Visit, event = VL_supp_status))
km_fit2a = survfit(Surv2 ~ 1, data = dat)
ggsurvplot(km_fit2a, data = dat, conf.type = "log-log")



#time to first missed visit
Surv3 = with(dat, Surv(time=Visit, event = missed_visit))
km_fit3a = survfit(Surv3 ~ 1, data = dat)
ggsurvplot(km_fit3a, data = dat, conf.type = "log-log")
