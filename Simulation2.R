# Master Thesis KU Leuven
# Adam Mohammad
# Title: A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
# Obtained from Joint Models for Longitudinal and Survival Data
# August 2017
# Code for combining simulated datasets into one list

rm(list=ls())
load("datasmall.Rdata")
load("datahigh.Rdata")

datasmall2=datasmall
datahigh2=datahigh


for (k in 1:3){
  for (l in 1:100){
    datasmall2[[k]][[l]]$name=paste("Small",k,sep="")
  }
}

for (k in 1:3){
  for (l in 1:250){
    datahigh2[[k]][[l]]$name=paste("High",k,sep="")
  }
}

data=list()
data1=list()
for (j in 1:3){
  data=append(data,datasmall2[[j]])
  data1=append(data1,datahigh2[[j]])
}

alldata=append(data,data1)
save(alldata,file="Alldata.RData")

save.image("resultsdatageneration2.RData")
