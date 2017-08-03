# Master Thesis KU Leuven
# Adam Mohammad
# Title: A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
# Obtained from Joint Models for Longitudinal and Survival Data
# August 2017
# Code for functions calculating Prediction error and AUC in simulation studies

#### Function AUC PE JM####
PeAuc=function(object, newdata, timepoints,idVar="id", M){
  PE=matrix(nrow = length(timepoints),ncol=length(timepoints))
  colnames(PE)=as.character(timepoints)
  rownames(PE)=as.character(timepoints)
  AUC=PE
  PEsim=PE
  AUCsim=PE
  NUM=vector(length = length(timepoints))
  t=length(timepoints)
  
  for (j in 1:I(t-1)){
    for (k in I(j+1):t){
      pes=list()
      pes=prederrJM(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,simulate = TRUE, M=M)
      aucs=list()
      aucs=aucJM(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,simulate = TRUE, M=M)
      PEsim[j,k]=pes$prederr
      AUCsim[j,k]=aucs$auc
      pe=list()
      pe=prederrJM(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,simulate = FALSE)
      auc=list()
      auc=aucJM(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,simulate = FALSE)
      PE[j,k]=pe$prederr
      AUC[j,k]=auc$auc
      
    }
    NUM[j]=pe$nr
  }
  result=list(PE=PE, AUC=AUC, PEsim=PEsim,AUCsim=AUCsim,NUM=NUM)
  return(result)
}

#### Function AUC PE Landmarking ####

aucCox=JMbayes:::aucJM.coxph
prederrCox=JMbayes:::prederrJM.coxph

PeAucCox=function(object, newdata, timepoints, idVar="id", respVar = "y", timeVar = "time", evTimeVar = "TimeEv"){
  PE_Cox=matrix(nrow = length(timepoints),ncol=length(timepoints))
  colnames(PE_Cox)=as.character(timepoints)
  rownames(PE_Cox)=as.character(timepoints)
  AUC_Cox=PE_Cox
  NUM=vector(length = length(timepoints))
  t=length(timepoints)
  
  for (j in 1:I(t-1)){
    for (k in I(j+1):t){
      pe=list()
      pe=prederrCox(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,respVar = respVar, 
                    timeVar = timeVar, evTimeVar = evTimeVar, summary="value",lossFun="square",interval = FALSE)
      auc=list()
      auc=aucCox(object,newdata,Tstart=timepoints[j],Thoriz = timepoints[k],idVar = idVar,respVar = respVar, 
                 timeVar = timeVar, evTimeVar = evTimeVar, summary="value")
      PE_Cox[j,k]=pe$prederr
      AUC_Cox[j,k]=auc$auc
    }
    NUM[j]=pe$nr
  }
  result=list(PE_Cox=PE_Cox, AUC_Cox=AUC_Cox, NUM=NUM)
  return(result)
}