# Master Thesis KU Leuven
# Adam Mohammad
# Title: A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
# Obtained from Joint Models for Longitudinal and Survival Data
# August 2017
# Code for simulating data from the joint models

rm(list=ls())
#### Simulation parameters Longitudinal ####
N=1000
K=15
t_max=25
Bound_kn = c(0, 25)
Int_kn = c(2.5, 9)
betas = c("Group0" = 3.4554,
           "Group1" = 2.9470,
           "Covariate"=2.1240,
           "Group0:Time1" = 1.0027,
           "Group1:Time1" = 0.9709, 
           "Group0:Time2" = 4.1290,
           "Group1:Time2" = 4.0893,
           "Group0:Time3" = 6.2182,
           "Group1:Time3" = 6.6909)

D=matrix(c(0.5686193, 0.2126076, 0.1547322, 0.4354939,
           0.2126076, 1.6721086, 2.3299235, 2.1926166,
           0.1547322, 2.3299235, 5.0230656, 2.8873934,
           0.4354939, 2.1926166, 2.8873934, 4.0286104), 4, 4)
D=(D + t(D)) / 2
sigma_small=0.564
sigma_high=5*sigma_small

#### Simulation Parameters survival ####
gammas = c("(Intercept)" = -5.7296, "Group" = 0.48,"Covariate"=1.25) # coefficients for baseline covariates
alpha0 = 0.4672 # association parameter  current value
alpha1 = 0.4044 # association parameter current value and slope 1
alpha2 = 1.3616 # association parameter current value and slope 2
alpha3 = 0.0365 # association parameter cumulative
phi = 0.9518    # shape for the Weibull baseline hazard
meanCens0 = 10  # mean of the uniform censoring distribution for group 0
meanCens1 = 14  # mean of the uniform censoring distribution for group 1

####Packages ####
#devtools::install_github("drizopoulos/JMbayes")

library(splines)
library(MASS)
library(JMbayes)  

####Data Generation ####
sigma=sigma_small
alldatasets=list()
set.seed(0478542)
for (sc in 1:3){
  data=list()
  for(it in 1:100){
    out=list()
    source(paste("scenario",sc,".R",sep = ""))
    out=list(train_data=train_data, test_data=test_data,train_data2=train_data[!duplicated(train_data$id),])
    data[[it]]=out
  }
  alldatasets[[sc]]=data
}

datasmall=alldatasets
save(datasmall, file="datasmall.RData")

sigma=sigma_high
alldatasets=list()
set.seed(0478542)
for (sc in 1:3){
  data=list()
  for(it in 1:250){
    out=list()
    source(paste("scenario",sc,".R",sep = ""))
    out=list(train_data=train_data, test_data=test_data,train_data2=train_data[!duplicated(train_data$id),])
    data[[it]]=out
  }
  alldatasets[[sc]]=data
}

datahigh=alldatasets


save(datahigh, file="datahigh.RData")
save.image("resultsdatageneration.RData")
save.image()