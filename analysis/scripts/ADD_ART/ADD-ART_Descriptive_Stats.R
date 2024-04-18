#ADD-ART Descriptive Statistics
#16 Jan 2024

#declare dependencies
library(tidyverse)
library(readxl)
library(psych)

#declare paths
path_in = "~/OneDrive - University of Cape Town/2024/Dissertation/Analysis/Input/ADD_ART"
path_out = "~/OneDrive - University of Cape Town/2024/Dissertation/Analysis/Output/ADD_ART"

#read in data
setwd(path_in)
dat = read_csv("ADDART_21Oct2022_with_monthly_proportions.csv")
dat = dat[,-1]
dat_baseline = dat %>%
  filter(VISIT == 0)

head(dat_baseline)
colnames(dat_baseline)

demo_colnames = paste0(rep("demo", 13), c(1,3,5,6, "6d", 7, 9,12, "12c", 15, 19, 21, "21b"))

#subset by variables of interest
dat_baseline = dat_baseline %>%
  select(PID, VISITdate, AGE, GENDER, IfWORK, EDUC, FoodINSEC,
         AnyFoodIns, demo_colnames,
         TFV, TFVmn, VLdate, VL, HCT, FTC)

#replace missing value codes with NA
dat_baseline[dat_baseline == 888] = NA
dat_baseline[dat_baseline == 999] = NA

#investigate number of missings
colSums(is.na(dat_baseline))
plot(rowSums(is.na(dat_baseline)))



#---------------------------------------- Descriptive stats -------------------------------------------
describe(dat_baseline)


