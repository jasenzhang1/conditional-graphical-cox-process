# GMpp

### This repository contains the R codes for the paper "Graphical Modeling of Multivariate Cox Process"

The following R codes can be used to obtain the simulation results in the paper. 

run_Model1.R: run the simulation setting "Model 1" for one replication.

run_Model2.R: run the simulation setting "Model 2" for one replication.

run_Model3.R: run the simulation setting "Model 3" for one replication.

run_Model4.R: run the simulation setting "Model 4" for one replication.

The results returned by these R codes include the F1 score, the true positive rate (TPR), the false positive rate (FPR), the true negative rate (TNR), and the false negative rate (FNR) for each method.




### Make sure to install all required packages (see below).

## To install all the other required packages, run the following codes in R 

require_packages = c("data.table", "Matrix", "Rcpp", "methods",
                      "parallel","pracma", "PtProcess", "MASS")

new_packages = lapply(require_packages, require, character.only = TRUE)

new_packages = require_packages[!unlist(new_packages)]

if(length(new_packages)) {install.packages(new_packages)}
