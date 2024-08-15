# ADD-ART Missing Data Management
rm(list=ls())
library(tidyverse)

path_in = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
##-----------long data -----------------
#load data
setwd(path_in)
dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)

dat$id = dat$PID



dat = dat %>%
  group_by(id) %>%
  mutate(time = as.numeric(difftime(VISITdate, VISITdate[2], units = "days")),
         time_yrs = time/365 ) %>%
  relocate(., .after = Visit)


##-----------missing data 1 -----------------
#-----------------EAM data: missing at baseline, so we remove baseline visits and rescale visits to start at 0
dat = dat %>%
  filter(Visit > 0) %>%
  mutate(Visit = Visit - 1)

#---------------- missed visits
#see who attended less than 3 visits
visit_missers = dat %>% group_by(id) %>% 
  summarise(n_missed_visits = sum(missed_visit)) %>%
  filter(n_missed_visits >= 10) %>% select(id)
n_distinct(visit_missers$id)
dat %>% filter(id %in% visit_missers$id) %>% print(n=Inf) #this is early dropout. they dont' have much EAM data either so will be excluded entirely

#remove early dropouts
dat = dat %>% filter(!(id %in% visit_missers$id))
n_distinct(dat$id)


#-----------------missing covariates
#filter out all rows with missing values for covariates and outcomes
tmp = dat %>%
  filter(!is.na(VL_baseline), !is.na(Vmed9_everhadTB),
         !is.na(AGE), !is.na(GENDER), !is.na(HCT_baseline))
dat %>% filter(id %in% setdiff(dat$id, tmp$id)) %>%
  select(id, Visit, VL, VL_baseline, Vmed9_everhadTB, n_heartbeats, VL, VL_baseline, AGE, GENDER, HCT_baseline) %>%
  print(n=Inf) #these are the pids with missing covariates

# check missing BMI
dat %>% ungroup() %>% filter(Visit == 0) %>%
  summarise(n_miss = sum(is.na(bmi_1st))) #too many NA's, will have to exclude from analysis

dat = tmp

##----------- missing outcome data
##--- EAM
missing_eam = dat %>% group_by(Visit, id) %>%
  summarise(n_miss = sum(is.na(n_intakes))) %>%
  mutate(id_visit = paste0(id, Visit)) %>%
  filter(n_miss > 0)
missing_eam

##--- TFV
#below missing many TFV
dat %>% group_by(id) %>%
  summarise(n_miss = sum(is.na(TFV))) %>%
  filter(n_miss >= 8)

missing_tfv = dat %>% group_by(Visit, id) %>%
  summarise(n_miss = sum(is.na(TFV))) %>%
  mutate(id_visit = paste0(id, Visit)) %>%
  filter(n_miss > 0) 
missing_tfv
n_distinct(missing_tfv$id)
missing_tfv %>% group_by(id) %>% filter(n_miss > 8)

#below have missing outcome data but apparently attended visit
dat %>% filter(missed_visit == 0 & is.na(TFV)) %>% print(n=Inf) 

n_distinct(missing_tfv$id)

##--- VL
dat %>% group_by(id) %>%
  summarise(n_miss = sum(is.na(VL))) %>%
  filter(n_miss > 6)

##-----------recurrent data -----------------
#get data in counting process format
#counting process format: id is represented by rows as study entry, event 1 time, event 2 time, event m time, final follow-up
dat_recurrent = dat 

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
  relocate(c(start, time_yrs, VL), .after = Visit)


#remove events at baseline
n_distinct(dat_recurrent$id)
tmp = dat_recurrent %>%
  filter(Visit > 0) 
n_distinct(tmp$id)

dat_recurrent = tmp

##------------ long data ----------------
#format longitudinal data
dat_long = dat %>%
  select(id, Visit, time, time_yrs, missed_visit, VISITdate, n_intakes_30days, n_heartbeats_30days, n_none_30days, n_intakes, n_heartbeats, n_none, wp_adherence_prop, wp_days_covered,
         TFV, HCT, BMI, BMIcat, bmi_1st, AGE, GENDER, TFV_baseline, HCT_baseline, VL_baseline, VL, 
         sra1_adherence_prop, sra2_goodjob_prop, sra3_way_supposed_to_prop, Vsra1_days_nonadherent) 
dat_long$id_visit = paste0(dat_long$id, dat_long$Visit)


#offset variables
dat_long = dat_long %>%
  mutate(wp_days = n_heartbeats + n_intakes,
         log_wp_days = log(wp_days+1),
         wp_days_30 = n_heartbeats_30days + n_intakes_30days,
         log_wp_days_30 = log(wp_days_30+1))

#remove observations matching id and visit in missing_eam
dat_long_eam = dat_long %>% filter(!(id_visit %in% missing_eam$id_visit))

#remove observations matching id and visit in missing_tfv
dat_long_tfv = dat_long %>% filter(!(id_visit %in% missing_tfv$id_visit))

dat_long_joint = dat_long %>% filter(!(id_visit %in% c(missing_eam$id_visit, missing_tfv$id_visit)))



#-------------- write data out
setwd(path_out)
write.csv(dat_recurrent,
          "ADDART_recurrent_14May2024.csv",
          row.names = FALSE)

write.csv(dat_long_eam,
          "ADDART_long_eam_14May2024.csv",
          row.names = FALSE)

write.csv(dat_long_tfv,
          "ADDART_long_tfv_14May2024.csv",
          row.names = FALSE)

write.csv(dat_long_joint,
          "ADDART_long_joint_14May2024.csv",
          row.names = FALSE)


