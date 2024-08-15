# ADD-ART Group-based Trajectory Modeling - with summary adherence on study MONTHLY BASIS
# C McDuling
# 27 June 2024
rm(list=ls())

library(tidyverse)
library(crimCV)
library(ggpubr)
library(survival)
library(survminer)
library(gridExtra)

#================ SOME DIAGNOSTICS - assess fit of underlying (cubic) model ======
#these functions come from (https://andrewpwheeler.com/2015/09/29/some-plots-to-go-with-group-based-trajectory-models-in-r/)
long_traj <- function(model,data){
  df <- data.frame(data)
  vars <- names(df)
  prob <- model['gwt'] #posterior probabilities
  df$GMax <- apply(prob$gwt,1,which.max) #which group # is the max
  df$PMax <- apply(prob$gwt,1,max)       #probability in max group
  df$Ord <- 1:dim(df)[1]                 #Order of the original data
  prob <- data.frame(prob$gwt)
  names(prob) <- paste0("G",1:dim(prob)[2]) #Group probabilities are G1, G2, etc.
  longD <- reshape(data.frame(df,prob), varying = vars, v.names = "y", 
                   timevar = "x", times = 1:length(vars), 
                   direction = "long") #Reshape to long format, time is x, y is original count data
  return(longD)                        #GMax is the classified group, PMax is the probability in that group
}

weighted_means <- function(model,long_data){
  G_names <- paste0("G",1:model$ng)
  long_data=na.omit(long_data)
  G <- long_data[,G_names]
  W <- G*long_data$y                                    #Multiple weights by original count var
  Agg <- aggregate(W,by=list(x=long_data$x),FUN="sum")  #then sum those products
  mass <- colSums(model$gwt)                            #to get average divide by total mass of the weight
  for (i in 1:model$ng){
    Agg[,i+1] <- Agg[,i+1]/mass[i]
  }
  long_weight <- reshape(Agg, varying=G_names, v.names="w_mean",
                         timevar = "Group", times = 1:model$ng, 
                         direction = "long")           #reshape to long
  return(long_weight)
}

pred_means <- function(model){
  prob <- model$prob               #these are the model predicted means
  Xb <- model$X %*% model$beta     #see getAnywhere(plot.dmZIPt), near copy
  lambda <- exp(Xb)                #just returns data frame in long format
  p <- exp(-model$tau * t(Xb))
  p <- t(p)
  p <- p/(1 + p)
  mu <- (1 - p) * lambda
  t <- 1:nrow(mu)
  myDF <- data.frame(x=t,mu)
  long_pred <- reshape(myDF, varying=paste0("X",1:model$ng), v.names="pred_mean",
                       timevar = "Group", times = 1:model$ng, direction = "long")
  return(long_pred)
}

pred_means_se <- function(model){
  prob <- model$prob               #these are the model predicted means
  Xb <- model$X %*% model$beta     #see getAnywhere(plot.dmZIPt), near copy
  lambda <- exp(Xb)                #just returns data frame in long format
  p <- exp(-model$tau * t(Xb))
  p <- t(p)
  p <- p/(1 + p)
  mu <- (1 - p) * lambda
  se <- sqrt(mu/nrow(model$param)) # se according to poisson, normal approx
  t <- 1:nrow(mu)
  myDF <- data.frame(x=t,mu = mu,se = se)
  long_pred <- myDF %>%
    pivot_longer(cols = -x, 
                 names_to = c(".value", "Group"), 
                 names_pattern = "(\\D+)(\\d+)") %>%
    arrange(Group, x)
  return(long_pred)
}



occ <- function(long_data){
  subdata <- subset(long_data,x==1)
  agg <- aggregate(subdata$PMax,by=list(group=subdata$GMax),FUN="mean")
  names(agg)[2] <- "AvePP" #average posterior probabilites
  agg$Freq <- as.data.frame(table(subdata$GMax))[,2]
  n <- agg$AvePP/(1 - agg$AvePP)
  p <- agg$Freq/sum(agg$Freq)
  d <- p/(1-p)
  agg$OCC <- n/d #odds of correct classification
  agg$ClassProp <- p #observed classification proportion
  #predicted classification proportion
  agg$PredProp <- colSums(as.matrix(subdata[,grep("^[G][0-9]", names(subdata), value=TRUE)]))/sum(agg$Freq) 
  #Jeff Ward said I should be using PredProb instead of Class prop for OCC
  agg$occ_pp <- n/ (agg$PredProp/(1-agg$PredProp))
  return(agg)
}

#fn to return sample means and standard dev for each group
obs_means = function(longdata){
  
  summarized <- longD %>%
    group_by(GMax,x) %>%
    summarize(obs_mean=mean(y, na.rm=TRUE),
              obs_sd = sd(y, na.rm=TRUE))
  
  summarized = summarized %>%
    rename(Group = GMax, Visit = x)
  
  return(summarized)
}

#=============== data management ================
setwd("~/Downloads/MSc/Dissertation/data_management/output/ADD_ART")
dat = read.csv("ADDART_long_eam_14May2024.csv")

dat_recurrent = read.csv("ADDART_recurrent_14May2024.csv")

# #drop any observations past month 12
# dat = dat %>% filter(Visit > 0)

#create matrix for intakes per month, wide format
dat.intakes= dat %>% select(id, Visit, n_heartbeats_30days) %>%
  pivot_wider(names_from = Visit,
              values_from = n_heartbeats_30days)


#create matrix of time-at-risk offset: t_ij = fraction of month days that WP device was active
dat.offset.wide = dat %>% select(id, Visit, wp_days_30) %>%
  pivot_wider(names_from = Visit,
              values_from = wp_days_30)

#drop id column from both matrices
rownames(dat.intakes) = dat.intakes$id
dat.intakes = dat.intakes[,-c(1)]
dat.intakes = as.matrix(dat.intakes)

rownames(dat.offset.wide) = dat.offset.wide$id
dat.offset.wide = dat.offset.wide[,-c(1)]
dat.offset.wide = as.matrix(dat.offset.wide)
#convert to proportion of 30
dat.offset.wide = dat.offset.wide/30


#replace missings with negative value
dat.intakes[is.na(dat.intakes)] = -99
dat.offset.wide[is.na(dat.offset.wide)] = -99

#plots all trajectories
ggplot(dat, aes(as.factor(Visit), n_heartbeats_30days, group=id)) + geom_line(aes(col=id, group=id))+
  theme(legend.position = "none")

#=============== modeling ================
#======= USING CRIMCV PACKAGE
N = nrow(dat.intakes); p = ncol(dat.intakes)

#risk matrix needs to have no zero values and no -99s
risk.dummy = dat.offset.wide
risk.dummy[risk.dummy==0] = 1e-14
risk.dummy[risk.dummy==-99] = 1



#====== QUADRATIC - increasing # groups
set.seed(2023)
gbtm2.2 = crimCV(dat.intakes, ng=2, dpolyp = 2, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm2.3 = crimCV(dat.intakes, ng=3, dpolyp = 2, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm2.4 = crimCV(dat.intakes, ng=4, dpolyp = 2, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm2.5 = crimCV(dat.intakes, ng=5, dpolyp = 2, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm2.6 = crimCV(dat.intakes, ng=6, dpolyp = 2, rcv=TRUE, Risk=risk.dummy)

#compare models
gbtm2.compare = data.frame(groups=c(2:6))
gbtm2.compare$loglik = c(gbtm2.2$llike, gbtm2.3$llike, gbtm2.4$llike, gbtm2.5$llike, gbtm2.6$llike)
gbtm2.compare$AIC = c(gbtm2.2$AIC, gbtm2.3$AIC, gbtm2.4$AIC, gbtm2.5$AIC, gbtm2.6$AIC)
gbtm2.compare$BIC = c(gbtm2.2$BIC, gbtm2.3$BIC, gbtm2.4$BIC, gbtm2.5$BIC, gbtm2.6$BIC)
gbtm2.compare$CVE = c(gbtm2.2$cv, gbtm2.3$cv, gbtm2.4$cv, gbtm2.5$cv, gbtm2.6$cv)
gbtm2.compare # groups along all criteria

plot(gbtm2.2)
plot(gbtm2.3)
plot(gbtm2.4)
plot(gbtm2.5)
plot(gbtm2.6)

#======= CUBIC - increasing # groups 
set.seed(2023)
gbtm3.2 = crimCV(dat.intakes, ng=2, dpolyp = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.3 = crimCV(dat.intakes, ng=3, dpolyp = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.4 = crimCV(dat.intakes, ng=4, dpolyp = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.5 = crimCV(dat.intakes, ng=5, dpolyp = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.6 = crimCV(dat.intakes, ng=6, dpolyp = 3, rcv=TRUE, Risk=risk.dummy)

#compare # of groups fit
gbtm3.compare = data.frame(groups=c(2:6))
gbtm3.compare$loglik = c(gbtm3.2$llike, gbtm3.3$llike, gbtm3.4$llike, gbtm3.5$llike, gbtm3.6$llike)
gbtm3.compare$AIC = c(gbtm3.2$AIC, gbtm3.3$AIC, gbtm3.4$AIC, gbtm3.5$AIC, gbtm3.6$AIC)
gbtm3.compare$BIC = c(gbtm3.2$BIC, gbtm3.3$BIC, gbtm3.4$BIC, gbtm3.5$BIC, gbtm3.6$BIC)
gbtm3.compare$CVE = c(gbtm3.2$cv, gbtm3.3$cv, gbtm3.4$cv, gbtm3.5$cv, gbtm3.6$cv)
gbtm3.compare # groups along all criteria

plot(gbtm3.2)
plot(gbtm3.3)
plot(gbtm3.4)
plot(gbtm3.5)
plot(gbtm3.6)
gbtm3_sum = summary(gbtm3.3)
gbtm4_sum = summary(gbtm3.4)
gbtm5_sum = summary(gbtm3.5)
gbtm6_sum = summary(gbtm3.6)

colSums(gbtm3_sum)
colSums(gbtm4_sum)
colSums(gbtm5_sum)
colSums(gbtm6_sum) # 6 groups are always rejected because of the small numbers making up the extra group

#create variable for group allocation and assign based off of max probability
dat.intakes2 = as.data.frame(dat.intakes)
dat.intakes2$GBTM3_group = rep(NA, nrow(dat.intakes2))
dat.intakes2$GBTM4_group = rep(NA, nrow(dat.intakes2))
dat.intakes2$GBTM5_group = rep(NA, nrow(dat.intakes2))
for(i in 1:nrow(gbtm3_sum)){
  allocation = which.max(gbtm3_sum[i,])
  dat.intakes2$GBTM3_group[i] = allocation
}
for(i in 1:nrow(gbtm4_sum)){
  allocation = which.max(gbtm4_sum[i,])
  dat.intakes2$GBTM4_group[i] = allocation
}
for(i in 1:nrow(gbtm5_sum)){
  allocation = which.max(gbtm5_sum[i,])
  dat.intakes2$GBTM5_group[i] = allocation
}


#====================== DIAGNOSTICS ===================
#choose best model now
gbtm = gbtm3.5

group_members = dat.intakes2 %>% 
  group_by(GBTM5_group) %>%
  summarise(n = n())
group_members$GBTM5_group = as.factor(group_members$GBTM5_group)


dat.intakes.nas = dat.intakes
dat.intakes.nas[dat.intakes.nas==-99] = NA

longD <- long_traj(model=gbtm,data=dat.intakes.nas)
x <- weighted_means(model=gbtm,long_data=longD)
pred <- pred_means(model=gbtm)
pred$Group = as.factor(pred$Group)
#estimate se
pred2 <- pred %>%
  left_join(group_members, by = c("Group" = "GBTM5_group")) %>%
  mutate(se = sqrt(pred_mean/n),
         lower = pred_mean - 1.96*se,
         upper = pred_mean + 1.96*se)



#plots individual observed trajectories
labels = NULL
for(i in 1:nrow(group_members)){
  labels = append(labels, 
                  paste0("Group ", i, " (n = ", group_members[i,2], ")"))
}
labels
custom_labels <- setNames(labels, 1:nrow(group_members))
longD$GMax = as.factor(longD$GMax)
p <- ggplot(data=longD, aes(x=x,y=y,group=Ord)) + geom_line(aes(colour = GMax),alpha = 0.3) +
  labs(x = "Visits", y = "Observed EAM missed doses") +
  facet_wrap(~GMax, labeller = labeller(GMax = custom_labels)) +
  theme_minimal() + theme(text = element_text(size = 16), legend.position = "none") 
p

#plot mean trajectories of each group
summ_dat = obs_means(longD)
summ_dat
summ_dat$Group = as.factor(summ_dat$Group)
p1 <- ggplot() +
  geom_line(data = summ_dat, 
            aes(x = Visit,
                y = obs_mean,
                col=Group)) +
  geom_errorbar(data = summ_dat,
                aes(x = Visit,
                    ymin = obs_mean - obs_sd,
                    ymax = obs_mean + obs_sd,
                    col=Group),
                width = 0.1) +
  labs(x = "Visits", y = "Observed mean EAM missed doses") +
  facet_wrap(~Group, labeller = labeller(Group = custom_labels)) + 
  theme_minimal() + theme(text = element_text(size = 16), legend.position = "none")
p1


#plots predicted trajectories overlayed with observed weighted means
p2 = ggplot() + geom_line(data=pred2, aes(x=x,y=pred_mean,col=Group)) + 
  geom_line(data=x, aes(x=x,y=w_mean,col=as.factor(Group)), alpha = 0.8, lty = 2) + 
  geom_point(data=x, aes(x=x,y=w_mean,col=as.factor(Group)), alpha = 0.6) + 
  geom_errorbar(data=pred2, 
                aes(x=x,y=pred_mean,ymin=lower,ymax=upper,col=Group), width = 0.1) +
  xlab("Visits") + ylab("Predicted and weighted means") +
  theme_minimal() + theme(text = element_text(size = 16)) 
p2

p2 <- ggplot() + 
  # Line plot for predicted means
  geom_line(data = pred2, aes(x = x, y = pred_mean, col = Group)) + 
  # Shaded area for confidence intervals
  geom_ribbon(data = pred2, aes(x = x, ymin = lower, ymax = upper, fill = Group), alpha = 0.2) +
  # Dashed line plot for weighted means
  geom_line(data = x, aes(x = x, y = w_mean, col = as.factor(Group)), alpha = 0.8, lty = 2) + 
  # Points for weighted means
  geom_point(data = x, aes(x = x, y = w_mean, col = as.factor(Group)), alpha = 0.6) + 
  # Labels and theme
  xlab("Clinic Visits") + ylab("Predicted and weighted mean missed doses") +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_minimal() + 
  theme(text = element_text(size = 16)) 
p2 

#add the poisson confidence interval (normal approximation)
p2 <- ggplot() + 
  # Line plot for predicted means
  geom_line(data = pred2, aes(x = x, y = pred_mean, col = Group)) + 
  # Shaded area for confidence intervals
  geom_ribbon(data = pred2, aes(x = x, ymin = lower, ymax = upper, fill = Group), alpha = 0.2) +
  # Dashed line plot for weighted means
  geom_line(data = x, aes(x = x, y = w_mean, col = as.factor(Group)), alpha = 0.8, lty = 2) + 
  # Points for weighted means
  geom_point(data = x, aes(x = x, y = w_mean, col = as.factor(Group)), alpha = 0.6) + 
  # Labels and theme
  xlab("Clinic Visits") + ylab("Predicted and observed mean missed doses") +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_minimal() + 
  theme(text = element_text(size = 16))
p2

#just the predicted trajectories
p3 = ggplot() + geom_line(data=pred, aes(x=x,y=pred_mean,col=Group)) + 
  xlab("Visits") + ylab("Predicted means") +
  theme_minimal() + theme(text = element_text(size = 16))
p3

#distribution of posterior probabilities
gbtm_sum = data.frame(summary(gbtm))
gbtm_sum$id = rownames(gbtm_sum)
gbtm_sum = gbtm_sum %>% gather(key="group", value="prob", -id)

gbtm_sum %>% ggplot() + 
  geom_histogram(aes(x = prob, fill = group)) +
  labs(x = "Posterior probability", y = "Frequency") +
  facet_wrap(~group, labeller = labeller(group = custom_labels)) +
  theme_minimal() + theme(text = element_text(size = 16), legend.position = "none") 


##====================== more stuff ===================
#can we extract standard errors to construct confidence intervals?
?crimCV

str(gbtm)


#black and white plot for abstract submission
p2 <- ggplot() + 
  # Line plot for predicted means
  geom_line(data = pred2, aes(x = x, y = pred_mean, group = Group)) + 
  # Shaded area for confidence intervals
  geom_ribbon(data = pred2, aes(x = x, ymin = lower, ymax = upper, group = Group), alpha = 0.2) +
  # Dashed line plot for weighted means
  geom_line(data = x, aes(x = x, y = w_mean, group = as.factor(Group)), alpha = 0.8, lty = 2) + 
  # Points for weighted means
  geom_point(data = x, aes(x = x, y = w_mean, group = as.factor(Group)), alpha = 0.6) + 
  # Labels and theme
  xlab("Clinic Visits") + ylab("Predicted and observed mean missed doses") +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  theme_minimal() + 
  theme(text = element_text(size = 16), legend.position = "none")
p2



# we explore time to viral non-suppression for each group via kaplan-meier curves
gbtm5_assignment = cbind(unique(dat$id),dat.intakes2$GBTM5_group)
gbtm5_assignment = data.frame(id = gbtm5_assignment[,1],
                              Group = gbtm5_assignment[,2])

dat_recurrent = dat_recurrent %>%
  left_join(gbtm5_assignment, by = c("id" = "id"))

Surv3a = with(dat_recurrent, 
              Surv(time = start, time2 = time_yrs,
                   event = VL> 50, type = "counting"))
km_fit3a = survfit(Surv3a ~ Group, data = dat_recurrent)
ggsurvplot(km_fit3a, data = dat_recurrent, 
           risk.table = T,
           cumevents = T)

#--------------------- test baseline characteristics between groups-----
#filter for 1 row per unique id
dat_recurrent = dat_recurrent %>%
  group_by(id) %>%
  slice(1)

dat_recurrent = dat_recurrent %>%
  mutate(SAMISStotal = sum(RecSamiss1, RecSamiss2, RecSamiss3, RecSamiss4, RecSamiss5, RecSamiss6, RecSamiss7,
                           RecSamiss8, RecSamiss9, RecSamiss10, RecSamiss11, RecSamiss12, RecSamiss13, RecSamiss14,
                           RecSamiss15, RecSamiss16
                           ))
  

cols = c("STIGMA", "DISCLOSE", "LEQtotal",
         "K10total", "SAMISStotal",
         "MSStotal", "HTKtotal", "BAMtotal", "PCRtotal", "HIBStotal")

plt = list()
for(col in cols){
  plt[[col]] = dat_recurrent %>%
    ggplot(aes(y = !!sym(col), x = Group, fill = Group)) +
    geom_boxplot() +
    # Add Kruskal-Wallis test
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    theme_minimal() + theme(text = element_text(size = 16)) +
    theme(legend.position = "none") 
}


do.call("grid.arrange", c(plt[1:5], ncol = 3))
do.call("grid.arrange", c(plt[6:10], ncol = 3))

