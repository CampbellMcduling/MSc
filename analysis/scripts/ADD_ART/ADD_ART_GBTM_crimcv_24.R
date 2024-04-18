# ADD-ART Group-based Trajectory Modeling - with summary adherence on study MONTHLY BASIS
# C McDuling
# 26 March 2024
rm(list=ls())

library(tidyverse)
library(crimCV)
library(ggpubr)

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
dat = read.csv("ADDART_visitwise_07March2024.csv", row.names = 1)


#drop any observations past month 12
dat = dat %>% filter(Visit > 0)

#create matrix for intakes per month, wide format
dat.intakes= dat %>% select(PID, Visit, n_intakes) %>%
  pivot_wider(names_from = Visit,
              values_from = n_intakes)


#create matrix of time-at-risk offset: t_ij = fraction of month days that WP device was active
dat.offset.wide = dat %>% select(PID, Visit, wp_days_covered) %>%
  pivot_wider(names_from = Visit,
              values_from = wp_days_covered)

#drop PID column from both matrices
rownames(dat.intakes) = dat.intakes$PID
dat.intakes = dat.intakes[,-c(1)]
dat.intakes = as.matrix(dat.intakes)

rownames(dat.offset.wide) = dat.offset.wide$PID
dat.offset.wide = dat.offset.wide[,-c(1)]
dat.offset.wide = as.matrix(dat.offset.wide)

#replace missings with negative value
dat.intakes[is.na(dat.intakes)] = -99
dat.offset.wide[is.na(dat.offset.wide)] = -99

#plots all trajectories
ggplot(dat, aes(as.factor(Visit), n_intakes, group=PID)) + geom_line(aes(col=PID, group=PID))+
  theme(legend.position = "none")

#=============== modeling ================
#======= USING CRIMCV PACKAGE
N = nrow(dat.intakes); p = ncol(dat.intakes)

#risk matrix needs to have no zero values and no -99s
risk.dummy = dat.offset.wide
risk.dummy[risk.dummy==0] = 1e-14
risk.dummy[risk.dummy==-99] = 1

#======= CUBIC - increasing # groups 
set.seed(2023)
gbtm3.0 = crimCV(dat.intakes, ng=2, dpolyp = 3, dpolyl = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.1 = crimCV(dat.intakes, ng=3, dpolyp = 3, dpolyl = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.2 = crimCV(dat.intakes, ng=4, dpolyp = 3, dpolyl = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.3 = crimCV(dat.intakes, ng=5, dpolyp = 3, dpolyl = 3, rcv=TRUE, Risk=risk.dummy)
set.seed(2023)
gbtm3.4 = crimCV(dat.intakes, ng=6, dpolyp = 3, dpolyl = 3, rcv=TRUE, Risk=risk.dummy)

#compare # of groups fit
gbtm3.compare = data.frame(groups=c(2:6))
gbtm3.compare$loglik = c(gbtm3.0$llike, gbtm3.1$llike, gbtm3.2$llike, gbtm3.3$llike, gbtm3.4$llike)
gbtm3.compare$AIC = c(gbtm3.0$AIC, gbtm3.1$AIC, gbtm3.2$AIC, gbtm3.3$AIC, gbtm3.4$AIC)
gbtm3.compare$BIC = c(gbtm3.0$BIC, gbtm3.1$BIC, gbtm3.2$BIC, gbtm3.3$BIC, gbtm3.4$BIC)
gbtm3.compare$CVE = c(gbtm3.0$cv, gbtm3.1$cv, gbtm3.2$cv, gbtm3.3$cv, gbtm3.4$cv)
gbtm3.compare # groups along all criteria

plot(gbtm3.0)
plot(gbtm3.1)
plot(gbtm3.2)
plot(gbtm3.3)
plot(gbtm3.4)
gbtm3_sum = summary(gbtm3.1)
gbtm4_sum = summary(gbtm3.2)
gbtm5_sum = summary(gbtm3.3)
gbtm6_sum = summary(gbtm3.4)

colSums(summary(gbtm3.2)) #approx number per group
colSums(summary(gbtm3.3)) #approx number per group
colSums(summary(gbtm3.4)) #approx number per group


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

dat.intakes.nas = dat.intakes
dat.intakes.nas[dat.intakes.nas==-99] = NA

longD <- long_traj(model=gbtm3.3,data=dat.intakes.nas)
x <- weighted_means(model=gbtm3.3,long_data=longD)
pred <- pred_means(model=gbtm3.3)
pred$Group = as.factor(pred$Group)

#plots individual observed trajectories
p <- ggplot(data=longD, aes(x=x,y=y,group=Ord)) + geom_line(alpha = 0.1) + facet_wrap(~GMax)
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
                width = 0.1) + facet_wrap(~Group)
p1


#plots predicted trajectories overlayed with observed weighted means
p2 = ggplot() + geom_line(data=pred, aes(x=x,y=pred_mean,col=Group)) + 
  geom_line(data=x, aes(x=x,y=w_mean,col=as.factor(Group))) + 
  geom_point(data=x, aes(x=x,y=w_mean,col=as.factor(Group))) + xlab("Months since study entrance") +theme_pubclean()
p2


