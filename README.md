# Master Thesis KU Leuven
This repository provides the R code for the simulation in my Master Thesis at KU Leuven (Master in Statistics Programme). 
The title of the thesis was _A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
Obtained from Joint Models for Longitudinal and Survival Data_.

Before running the code, `devtools` and `JMbayes` packages should be installed. 
The latter should be installed from the corresponding GitHub repository: `devtools::install_github("drizopoulos/JMbayes")`.

## Simulating Datasets 
`Simulation1.R` contains the code for simulating datasets from joint models with three association structures used in the thesis
(current value, current value and slope, cumulative effects). 
Model parameters are defined in the first part of the code and the rest of the code simulates datasets.

The code for simulating a dataset and spliting it into training and test dataset can be found in:
* `scenario1.R` for joint model with the current value association structure
* `scenario2.R` for joint model with the current value and slope association structure
* `scenario3.R` for joint model with the cumulative effects association structure

For low error variance setting 100 datasets are simulated from each scenario, while more datasets are simulated for high variance scenario.

`Simulation2.R` combines all simulated datasets into one list. 

## Fitting Joint Models and Assessing Model Performance
`Simulation3.R` contains the code for fitting joint models with the three association structures 
on a training dataset and calculating predictive accuaracy measures (AUC and Prediction Error) by using a corresponding test dataset. 

The time points used for calculating the predictive accuracy measures and JMbayes settings can be modified at the beginning of the code. 
R code for calculating Prediction Errors and AUCs can be found in `functions.R`.

If a HPC Cluster (a "Supercomputer") is used for calculations, `Simulation3_HPC.r` should be used. 
Due to the fact that in high variance setting, some joint models do not converge, the use of HPC is preferred.

## Additional Resources
* `Thesis_Simulation.Rproj` is a R project file for RStudio. 
* `Literature Review.pdf` gives an overview of the available literature on the topic.

