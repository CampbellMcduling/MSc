# ADD-ART Definitive Data Management
# 05 March 2024
# Simplify the data management process for ADD-ART data
#Wp data: adjust curiosity dosing, tally up intakes and total signals
# combining WP and clinic data: count intake and number of signals between
# we want: long df with time-varying covariates and repeated baseline data
#time-varying: adherence:counts, self-reported, tfv-dp, viral load
#repeated baseline: socio-demographic data

rm(list=ls())
# Load libraries
library(tidyverse)
library(foreign)
library(lubridate)
library(readxl)

pss2date <- function(x) {as.Date(x/86400, origin = "1582-10-14")}

path_in = "/Users/campbell.mcduling/Downloads/MSc/Dissertation/data_management/input/ADD_ART"
path_out = "/Users/campbell.mcduling/Downloads/MSc/Dissertation/data_management/output/ADD_ART"

#-------------- read in data ---------------
# import wp data
setwd(path_in)
dat_wp = read_xlsx("WP_merged_14DEC2020.xlsx", sheet=1)


#import socio-behavioural and demographic data
dat_clinic <- read.spss("ADDART_21Oct2022.sav", to.data.frame=TRUE)



#------------------- WISEPILL CLEANING -----------------------
#------- variable formats
dat_wp = dat_wp %>% rename(PID=pid)
dat_wp$Date = as.Date(dat_wp$Date, format = "%Y-%m-%d")
dat_wp$`Date&Time` = as.POSIXct(dat_wp$`Date&Time`, format = "%Y-%m-%d %H:%M:%S")
#isolate time from date&time
dat_wp$Time = format(dat_wp$`Date&Time`, format = "%H:%M:%S")

#------- curiosity dosing
#look at distribution of number of intakes per day
dat_wp %>%
  group_by(PID, Date) %>%
  summarise(n_intakes = sum(Type == "INTAKE")) %>%
  ggplot(aes(n_intakes)) + geom_histogram() + theme(legend.position = "none")

#first, remove exact duplicates with same PID, Date, Type, Time: only keep one observation per duplicate
dat_wp2 = dat_wp[!duplicated(dat_wp[, c("PID", "Date&Time", "Type")]),]

# for each PID, only one heartbeat should be recorded per day
# if there are multiple heartbeats, keep the one with the latest time
  



#rules: only one of each Type can be recorded per PID per day, so duplicates of PID, Date, Type should be removed
dat_wp2 = dat_wp2[!duplicated(dat_wp2[, c("PID", "Date", "Type")]),]


# if there are multiple intakes, keep the one with the latest time
tmp = dat_wp2 %>%
  group_by(PID, Date) %>%
  summarise(n_intakes = sum(Type == "INTAKE")) 
curiosities = tmp[which(tmp$n_intakes > 1),]
curiosities$PID_Date = paste(curiosities$PID, curiosities$Date, sep = "_")

#view all the entries that are on the same day
dat_wp2$PID_Date = paste(dat_wp2$PID, dat_wp2$Date, sep = "_")
doubles = dat_wp2 %>%
  filter(PID_Date %in% curiosities$PID_Date) %>% 
  arrange(PID, Date, `Date&Time`) 

#if dat_wp$PID_Date is in curiosities$PID_Date, keep the entry with the latest time, and if there are multiple latest entries with the same time, keep 1
doubles %>%
  group_by(PID, Date) %>%
  filter(`Date&Time` != max(`Date&Time`))-> duplicates

#remove duplicates from dat_wp
dat_wp2 = dat_wp2 %>%
  anti_join(duplicates, by = c("PID", "Date", "Type", "Time"))


#check that there is one intake per day
dat_wp2 %>%
  group_by(PID, Date) %>%
  summarise(n_intakes = sum(Type == "INTAKE")) %>%
  filter(n_intakes > 1)

#---------- if a heartbeat is recorded on the same day as an intake, remove the heartbeat
# if a heartbeat is recorded on the same day as an intake, remove the heartbeat
dat_wp2 %>%
  group_by(PID, Date) %>%
  filter(Type == "INTAKE" | Type == "HEARTBEAT") %>%
  filter(n() > 1) %>%
  filter(Type == "HEARTBEAT") %>%
  arrange(PID, Date, `Date&Time`) %>%
  slice(1) -> heartbeats

#remove heartbeats from dat_wp2
dat_wp2 = dat_wp2 %>%
  anti_join(heartbeats, by = c("PID", "Date", "Type", "Time"))


#------------------- CLINIC DATA CLEANING -----------------------
#fully populate PID column 
for(i in 2:nrow(dat_clinic)){
  if(dat_clinic$SID[i] == dat_clinic$SID[i-1]){
    dat_clinic$PID[i] = dat_clinic$PID[i-1]
  }
}

#------- visit dates
#apply pss2date to all columns with "date" in name
dat_clinic = dat_clinic %>% mutate(across(contains("date"), pss2date))

#check date variables
sum(na.omit(dat_clinic$VLdate != dat_clinic$DBSdate)) #VL and DBS dates are all equal, except for missings below
sum(is.na(dat_clinic$VLdate))
sum(is.na(dat_clinic$DBSdate))
dat_clinic[which(is.na(dat_clinic$VLdate)),]
dat_clinic[which(is.na(dat_clinic$DBSdate)),]


#create new variable - visit date - using the VL, DBS and V12 dates
dat_clinic$VISITdate = dat_clinic$DBSdate

head(dat_clinic$VISITdate)

for( i in 1:nrow(dat_clinic)){
  if(is.na(dat_clinic$VISITdate[i])){dat_clinic$VISITdate[i] = dat_clinic$VLdate[i]}
}
sum(is.na(dat_clinic$VISITdate))
dat_clinic[which(is.na(dat_clinic$VISITdate)),]  #can replace this missing with V12 date
dat_clinic[which(is.na(dat_clinic$VISITdate)),]$VISITdate =  dat_clinic[which(is.na(dat_clinic$VISITdate)),]$V12date

#make sure date classes are the same
dat_wp$Date = as.Date(dat_wp$Date)
dat_clinic$VISITdate = as.Date(dat_clinic$VISITdate)

dat_clinic = dat_clinic %>% relocate(VISITdate, .after = 3)

#manually correcting obvious data entry errors
#one date was input as January but supposed to be february, for everything else keep the original date
dat_clinic = dat_clinic %>%
  mutate(VISITdate = case_when(
    PID == "AD133" & VISIT == 10 ~ as.Date("2018-02-12"),
    TRUE ~ VISITdate
  ))

#------ MISSING VISITS
#create template df with 12 visits for each PID
template_df = data.frame(PID = rep(unique(dat_clinic$PID), each = 13),
                         VISIT = rep(0:12, times = length(unique(dat_clinic$PID))))


#merge template with dat_clinic
dat_clinic2 = merge(template_df, dat_clinic, by = c("PID", "VISIT"), all.x = TRUE)


#create new variable, missed visit which == 1 if visit date is missing
dat_clinic2$missed_visit = ifelse(is.na(dat_clinic2$VISITdate), 1, 0)
dat_clinic2 = dat_clinic2 %>% relocate(missed_visit, .after = 4)

#if visitdate is missing, replace with previous visit date + 30 days
for(i in 1:nrow(dat_clinic2)){
  if(is.na(dat_clinic2$VISITdate[i])){
    dat_clinic2$VISITdate[i] = dat_clinic2$VISITdate[i-1] + 30
  }
}

##------------------------- WP COUNTS: between visit dates -------------------------
#sort both by PID and visitdate
dat_wp2 <- dat_wp2[order(dat_wp2$PID, dat_wp2$Date), ]
dat_clinic2 <- dat_clinic2[order(dat_clinic2$PID, dat_clinic2$VISITdate), ]
#is every PID in dat_wp also in dat_clinic?
n_distinct(dat_wp2$PID)
n_distinct(dat_clinic2$PID)
#in wp data, create new variable, visit, which is the visit number for each PID 
dat_wp2$VISIT = 0
for(i in 1:nrow(dat_clinic2)){
  dat_wp2$VISIT[which(dat_wp2$PID == dat_clinic2$PID[i] & dat_wp2$Date >= dat_clinic2$VISITdate[i])] = dat_clinic2$VISIT[i] + 1
}
dat_wp2 = dat_wp2 %>% relocate(VISIT, .after = 3)

#tally the number of intakes, heartbeats and signals for each PID and visit
tmp = dat_wp2 %>%
  group_by(PID, VISIT) %>%
  summarise(n_intakes = sum(Type == "INTAKE"),
            n_heartbeats = sum(Type == "HEARTBEAT"),
            n_none = sum(Type == "NONE"))

#merge with dat_clinic2
dat_clinic2 = dat_clinic2 %>%
  left_join(tmp, by = c("PID", "VISIT"))

#for each PID, find days between consecutive VISITdates
dat_clinic2 = dat_clinic2 %>%
  group_by(PID) %>%
  mutate(days_between_visits = as.numeric(round(VISITdate - lag(VISITdate, default = first(VISITdate)), 0))) %>%
  relocate(days_between_visits, .after = 5)

#create proportion variables
dat_clinic2 = dat_clinic2 %>%
  mutate(wp_adherence_prop = n_intakes / (n_intakes + n_heartbeats),
         wp_days_covered = (n_intakes + n_heartbeats)/(days_between_visits)
         )

dat_clinic2 = dat_clinic2 %>%
  relocate(c(n_intakes, n_heartbeats, n_none, 
             wp_adherence_prop, wp_days_covered), .after = 6)

## ----------------------- WP COUNTS: 30 DAY BEFORE VISIT COUNTS ------------------------
##---------------- FOR EACH PID-VISIT: tally counts 30-days prior to visit date
#sort both by PID and visitdate
dat_wp2 <- dat_wp2[order(dat_wp2$PID, dat_wp2$Date), ]
dat_clinic2 <- dat_clinic2[order(dat_clinic2$PID, dat_clinic2$VISITdate), ]
#is every PID in dat_wp also in dat_clinic?
n_distinct(dat_wp2$PID)
n_distinct(dat_clinic2$PID)

dat_clinic2$PIDVisit = paste0(dat_clinic2$PID, dat_clinic2$VISIT)

dat_clinic2 = dat_clinic2 %>%
  mutate(n_intakes_30days = NA,
         n_heartbeats_30days = NA,
         n_none_30days = NA) %>%
  relocate(n_intakes_30days, n_heartbeats_30days, n_none_30days, .after = 5)

for(i in 1:nrow(dat_clinic2)){
  tmp = dat_wp2 %>% filter(PID == dat_clinic2$PID[i] & Date >= (dat_clinic2$VISITdate[i] - 30) & Date < dat_clinic2$VISITdate[i])
  #if tmp is empty, skip to next i
  if(nrow(tmp) == 0){next}
  dat_clinic2$n_intakes_30days[i] = sum(tmp$Type == "INTAKE")
  dat_clinic2$n_heartbeats_30days[i] = sum(tmp$Type == "HEARTBEAT")
  dat_clinic2$n_none_30days[i] = sum(tmp$Type == "NONE")
}

#---------------------------------- repeating baseline values -----------------
#isolate baseline data and remove columns with all NA values
dat_baseline = dat_clinic2 %>% filter(VISIT == 0)
dat_baseline = dat_baseline[, colSums(is.na(dat_baseline)) < nrow(dat_baseline)]
  

#isolate time-varying data
colnames(dat_clinic2)
dat_timevarying = dat_clinic2 %>%
  select(1:14,
         VL, HCT, TFV, FTC,
         Vsra1, Vsra2, Vsra3, Vsra4, 
         sra_1, sra_2, sra_3, sra_percent, sra_sd
         )
#exclude columns from baseline that are in time-varying
dat_baseline$VL_baseline = dat_baseline$VL
dat_baseline$HCT_baseline = dat_baseline$HCT
dat_baseline$TFV_baseline = dat_baseline$TFV
dat_baseline$FTC_baseline = dat_baseline$FTC
dat_baseline$sra_1_baseline = dat_baseline$sra_1
dat_baseline = dat_baseline %>%
  select(PID, colnames(dat_baseline)[!colnames(dat_baseline) %in% colnames(dat_timevarying)])

#merge time-varying data with baseline data
dat_long = dat_timevarying %>%
  left_join(dat_baseline, by = "PID", suffix = c("", "_baseline"))

# which cols have we lost
colnames(dat_clinic2)[!colnames(dat_clinic2) %in% colnames(dat_long)] #dont care about these, #goodriddance



#------------------rename some columns for interpretability ---------------
colnames(dat_long)

dat_long = dat_long %>%
  rename(Visit = VISIT,
         Vsra1_days_nonadherent = Vsra1,
         Vsra2_good_job = Vsra2,
         Vsra3_way_supposed_to = Vsra3,
         sra1_adherence_prop = sra_1,
         sra2_goodjob_prop = sra_2,
         sra3_way_supposed_to_prop = sra_3,
         sra_prop_mean = sra_percent,
         demo1_age = demo1,
         demo2a_birthprovince = demo2a,
         demo2b_birthcountry = demo2b,
         demo3_gender = demo3,
         demo4_homelang = demo4,
         demo4a_homelang_other = demo4a,
         demo5_education = demo5,
         demo6_working = demo6,
         demo6a_timesincework = demo6a,
         demo6b_employmentstatus = demo6b,
         demo6c_worktype = demo6c,
         demo6d_nightshifts = demo6d,
         demo6f_toattendclinic_timeoffwork = demo6f,
         demo7_HHincome = demo7,
         demo8_monthsinCPT = demo8months,
         demo8_yearsinCPT = demo8years,
         demo9_relationship = demo9,
         demo10_dwelling = demo10,
         demo11_HHsize = demo11,
         demo12_children = demo12,
         demo12a_numchildren = demo12a,
         demo12b_bybirth = demo12b,
         demo12c_toattendclinic_paysitter = demo12c,
         demo14a_foodmoney_own = demo14a,
         demo14b_foodmoney_otherinHH = demo14b,
         demo14c_foodmoney_otheroutHH = demo14c,
         demo15_foodsecurity_howoftennofood = demo15,
         demo16_race = demo16,
         demo18_citizenship = demo18,
         demo19_anyoflastyear_outsideCPT = demo19,
         demo20a_grant_disability_ever = demo20a,
         demo20b_grant_pension_ever = demo20b,
         demo20c_grant_child_ever = demo20c,
         demo20d_grant_other_ever = demo20d,
         demo21_transport_toclinic = demo21,
         demo21a_transport_toclinic_other = demo21a,
         demo21b_trasnport_toclinic_cost = demo21b,
         demo_22a_talk_worker = demo22a,
         demo_22b_talk_psychologist = demo22b,
         demo_22c_talk_doc_or_psychiatrist = demo22c,
         demo22d_talk_nurse = demo22d,
         demo22e_talk_pastor = demo22e,
         demo22f_talk_healer = demo22f,
         Vmed1_yourclinic = Vmed1,
         Vmed2_yourclinic_times = Vmed2,
         Vmed3_testpos_year = Vmed3,
         Vmed4_testpos_where = Vmed4,
         Vmed4a_testpos_other = Vmed4a,
         Vmed5_rate_health = Vmed5,
         Vmed6_everspentnight = Vmed6,
         Vmed6a_hospitalstays_past2yrs = Vmed6a,
         Vmed7_HIVillness_ever = Vmed7,
         Vmed9_everhadTB = Vmed9,
         Vmed9b_TBstatus = Vmed9b,
         Vmed9c_TB_receivingtreatment = Vmed9c,
         Vmed10_anemia_ever = Vmed10,
         Vmed10a_anemia_cause = Vmed10a,
         Vmed11_cd4_evercounted = Vmed11,
         Vmed11a_cd4_count_latest = Vmed11a,
         Vmed11b_cd4_count_lowest = Vmed11b,
         Vmed12_VL_test_ever = Vmed12,
         Vmed12a_VL_test_latest = Vmed12a,
         Vmed12b_VL_test_highest = Vmed12b,
         
  )


#------------------- recode missing values ----------------
#replace all 888 values with NA, new method
dat_long[dat_long == 888] = NA
colSums(dat_long == 999, na.rm = T)
dat_long[dat_long == 999] = NA


#------------------- write to file ------------------------
#write to file
setwd(path_out)
write.csv(dat_long, "ADDART_visitwise_07March2024.csv")

