# Master Thesis KU Leuven
# Adam Mohammad
# Title: A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
# Obtained from Joint Models for Longitudinal and Survival Data
# August 2017
# Code for fitting joint models and calculating their predictive performance

rm(list = ls())
load("alldata.RData")
source("functions.R")

library(doParallel)
library(nlme)

#### PE, AUC Calculation####
n.burnin=10000L
n.iter=200000L
n.adapt=10000L
M=200L
time_points=c(5.5,7.5,9.5,11.5,13.5,15.5)

lmcont=lmeControl(opt="optim", maxIter = 1000L,msMaxIter = 1000L)

####
set.seed(0478542)
#cn=detectCores()

dForm = list(fixed = ~ 0 + group:dns(time, knots=c(2.5, 9), B = c(0, 25)), 
             indFixed = 4:9,random = ~ 0 + dns(time, knots=c(2.5, 9), B = c(0, 25)), 
             indRandom = 2:4)

iForm = list(fixed = ~ 0 + group:time +Covariate:time + group:ins(time, knots=c(2.5, 9), B = c(0, 25)),
             indFixed = 1:9, random = ~ 0 + time + ins(time, knots=c(2.5, 9), B = c(0, 25)), 
             indRandom = 1:4)


data1=alldata

myresults=list()
#registerDoParallel(cores=cn)
registerDoParallel(cores=7L)
myresults=foreach(r=1:1050)%dopar%{
  library(nlme)
  library(survival)
  library(JMbayes)
  library(splines)
  
  dat=data1[[r]]
  
  lm=lme(y ~ 0 + group + Covariate+ group:ns(time, knots=c(2.5, 9), B = c(0, 25)), 
         random=list(id=pdDiag(~1+ns(time, knots=c(2.5, 9), B = c(0, 25)))),
         data=dat$train_data, method = "REML",control =lmcont)
  
  sm=coxph(Surv(TimeEv, event) ~ group+Covariate, data = dat$train_data2, x = TRUE, model=TRUE)
  
  joint1=jointModelBayes(lmeObject = lm,survObject = sm,timeVar = "time",
                         control =list(n.iter=n.iter,n.burnin=n.burnin,n.adapt=n.adapt,keepRE=TRUE) )
  
  joint2=jointModelBayes(lmeObject = lm,survObject = sm,timeVar = "time",
                         param = "td-both",extraForm = dForm,
                         control =list(n.iter=n.iter,n.burnin=n.burnin,n.adapt=n.adapt,keepRE=TRUE) )
  
  joint3=jointModelBayes(lmeObject = lm,survObject = sm,timeVar = "time",
                         param = "td-extra", extraForm = iForm,
                         control =list(n.iter=n.iter,n.burnin=n.burnin,n.adapt=n.adapt,keepRE=TRUE) )
  
  result_model1=PeAuc(object=joint1,newdata = dat$test_data,timepoints = time_points,idVar = "id", M=M )
  result_model2=PeAuc(object=joint2,newdata = dat$test_data,timepoints = time_points,idVar = "id", M=M )
  result_model3=PeAuc(object=joint3,newdata = dat$test_data,timepoints = time_points,idVar = "id", M=M )
  result_landmarking=PeAucCox(sm,dat$test_data,time_points)
  
  result=list(result_model1=result_model1,result_model2=result_model2,result_model3=result_model3,
              result_landmarking=result_landmarking,time_model1=joint1$time[3],time_model2=joint2$time[3],
              time_model3=joint3$time[3],name=dat$name)
  out=result
  
}

save(myresults,file="results.RData")

save.image("results_simulation.RData")
save.image()


