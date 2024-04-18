##=========================== RETAIN - INTERIM RESULTS ============
rm(list=ls())
setwd("~/Downloads/MSc/Dissertation/Analysis/Input/RETAIN")

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

dat = read_xlsx("RETAIN V1 Summary Report, 24-Feb-2024.xlsx",
                sheet = 1)
dat = dat[-c(171:180),]
str(dat)
colnames(dat)
dat = dat %>% rename("SR adherence %"="...6", "WP adherence %" = "...11")
n_distinct(dat$PID)

#remove rows with missing Flagged by
dat = dat %>% filter(!is.na(`Flagged by`))

dat$`Flagged by` = factor(dat$`Flagged by`, levels = c("WP", "PR", "VL"))
dat$logVL = log(dat$VL)

setwd("~/Downloads/MSc/Dissertation/Analysis/Output/RETAIN/Figures")
#FLAGGING
#Method
ggplot(dat, aes(x=`Flagged by`, fill=`Flagged by`))+ geom_bar() + theme_pubclean()
#days from start to flag
ggplot(dat, aes(x=`Days`,fill=`Flagged by`)) + geom_density(alpha=0.8)  + theme_pubclean() +
  labs(x="Days until flagging", title = "Days until awareness of poor adherence") + theme(legend.position = 'none')
ggsave("days_vs_flagtype_density.pdf", width=6, height=4, units="in")

ggplot(dat, aes(x=`Days`,fill=`Flagged by`)) + geom_histogram(alpha=0.8)  + theme_pubclean()
ggplot(dat, aes(y=`Days`,x=`Flagged by`, fill=`Flagged by`)) + geom_boxplot()  + theme_pubclean() +
  labs(y='Days until flagging', title = "Days until awareness of poor adherence") + theme(legend.position = 'none')
ggsave("days_vs_flagtype_boxplot.pdf", width=6, height=4, units="in")

#ADHERENCE
#self report
summary(dat[,3:6])
ggplot(dat, aes(x=`SR number missed /30`)) + geom_bar() + theme_pubclean()
ggplot(dat, aes(x=`SR adherence %`)) + geom_bar() + theme_pubclean() + labs(x="Self-report Adherence (%/30days)")
ggplot(dat, aes(x=`SR good job`)) + geom_bar() + theme_pubclean() 

#WP
summary(dat$`WP adherence %`)
ggplot(dat, aes(x=WP, fill=`Flagged by`)) + geom_histogram()
ggplot(dat, aes(x=`WP adherence %`)) + geom_histogram() + labs(x="Proportion of WP doses taken in last 30 days")

#WP vs SR
dat2 = melt(dat, id.vars = c("PID")) 
dat2 = dat2 %>% filter(variable == "SR adherence %" | variable == "WP adherence %")
dat2$value = as.numeric(dat2$value)
dat2$variable = ifelse(dat2$variable=="SR adherence %", "Self-reported", "Wisepill")
ggplot(dat2 ) + geom_boxplot(aes(x=variable, y=value*100, fill=variable))  + theme_pubclean() + theme(legend.position = 'none') + labs(x="", y="Adherence over last 30 days (%)")

#VIRAL LOAD
summary(dat$logVL)
ggplot(dat, aes(x=logVL)) + geom_histogram() 

ggplot(dat, aes(x=logVL)) + geom_histogram((aes(y=..count../sum(..count..)))) +
  geom_vline(xintercept=3.92, linetype="dotted", color="red") + labs(y="Proportion")
sum(na.omit(log(dat$logVL))<3.92)/sum(!is.na(dat$logVL)) #proportion unsuppressed

ggplot(dat, aes(x=logVL, fill=`Flagged by`)) + geom_density(alpha=0.7) 

ggplot(dat, aes(y=logVL)) + geom_boxplot(aes(x=`Flagged by`, fill=`Flagged by`))+ theme_pubclean() + theme(legend.position = 'none')

#UTRA
summary(dat$UTRA)
ggplot(dat, aes(x=UTRA, fill=`Flagged by`)) + geom_bar()  + theme_pubclean()

#PL TFV
summary(dat$`pl TFV (ng/ml)`)
ggplot(dat, aes(x=`pl TFV (ng/ml)`)) + geom_histogram((aes(y=..count../sum(..count..)))) + labs(y="Proportion")
ggplot(dat, aes(x=log(`pl TFV (ng/ml)`))) + geom_histogram()

ggplot(dat, aes(y=`pl TFV (ng/ml)`)) + geom_boxplot(aes(x=`Flagged by`, fill=`Flagged by`))+ theme_pubclean() + theme(legend.position = 'none')


#TDF-DP DBS
summary(dat$`TDF-DP (fmol/punch)`)
ggplot(dat, aes(x=`TDF-DP (fmol/punch)`)) + geom_histogram((aes(y=..count../sum(..count..)))) + labs(y="Proportion")
ggplot(dat, aes(x=log(`TDF-DP (fmol/punch)`))) + geom_density()

ggplot(dat, aes(y=`TDF-DP (fmol/punch)`)) + geom_boxplot(aes(x=`Flagged by`, fill=`Flagged by`))+ geom_hline(yintercept = 800, linetype=2, color='red') +
  theme_pubclean() + theme(legend.position = 'none') 

#add custom legend to above plot that labels the red dotted line
ggplot(dat, aes(y=`TDF-DP (fmol/punch)`)) + geom_boxplot(aes(x=`Flagged by`, fill=`Flagged by`)) + geom_hline(yintercept = 800, linetype=2, color='red') +
  theme_pubclean() + theme(legend.position = 'none') +
  labs(title="TFV-DP concentrations at Visit 1", y="TFV-DP (fmol/punch)", x="Flagged by") +
  annotate("text", x=3, y=1100, label="800 fmol/punch", color="red")

ggsave("TFV-DP_boxplot.pdf", width=6, height=4, units="in")

ggplot(dat, aes(x=`TDF-DP (fmol/punch)`)) + geom_density(aes( fill=`Flagged by`), alpha=0.7) +
  theme_pubclean() + theme(legend.position = 'none') 


kruskal.test(data=dat, logVL~`Flagged by`)
kruskal.test(data=dat, `TDF-DP (fmol/punch)`~`Flagged by`)
kruskal.test(data=dat, `pl TFV (ng/ml)`~`Flagged by`)
kruskal.test(data=dat, `Days`~`Flagged by`)

