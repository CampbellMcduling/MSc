#ADD-ART Descriptive Statistics
rm(list=ls())
#16 Jan 2024

#declare dependencies
library(tidyverse)
library(readxl)
library(psych)

#declare paths
path_in = "~/Downloads/MSc/Dissertation/data_management/output/ADD_ART"
path_out = "~/Downloads/MSc/Dissertation/analysis/output/ADD_ART"

#read in data
setwd(path_in)
dat = read.csv("ADDART_visitwise_07March2024.csv",
               row.names = 1)
#read in set used for analysis to identify excluded participants
dat_analysis = read.csv("~/Downloads/MSc/Dissertation/data_management/output/ADD_ART/ADDART_recurrent_14May2024.csv")

dat = dat %>%
  filter(Visit == 0)

dat_analysis = dat_analysis %>%
  filter(!duplicated(PID))

dat_excluded = dat %>%
  filter(!PID %in% dat_analysis$PID)

#subset by variables of interest
dat_baseline = dat %>%
  select(PID, VISITdate, missed_visit,
         VL, HCT, TFV, FTC, contains("Vsra"),
         AGE, GENDER, IfWORK, EDUC, FoodINSEC,
         AnyFoodIns, SAMISS, STIGMA, DISCLOSE, BMI, bmi_1st,
         demo2b_birthcountry, demo4_homelang, demo7_HHincome, demo11_HHsize,
         demo12_children, demo12a_numchildren, demo16_race, 
         Vmed1_yourclinic, Vmed3_testpos_year, Vmed5_rate_health,
         Vmed7_HIVillness_ever, Vmed9_everhadTB,Vmed10_anemia_ever,Vmed11_cd4_evercounted, Vmed11a_cd4_count_latest,
         Vmed12_VL_test_ever, AreyoucurrentlyPregnant)

#replace missing value codes with NA
dat_baseline[dat_baseline == 777] = NA
dat_baseline[dat_baseline == 888] = NA
dat_baseline[dat_baseline == 999] = NA

#investigate number of missings
colSums(is.na(dat_baseline))
plot(rowSums(is.na(dat_baseline)))



#---------------------------------------- Descriptive stats -------------------------------------------
describe(dat_baseline, 
         skew = FALSE)
summary(dat_baseline$VL <50)
summary(dat_baseline$VL)

dat_baseline %>% select(VL) %>%
  filter(VL > 50) %>%
  summary(VL)

summary(dat_baseline$FTC)
hist(dat_baseline$FTC, 30)

summary(dat_baseline$Vmed3_testpos_year)
hist(dat_baseline$Vmed3_testpos_year, 30)

