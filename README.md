# SAM

This is a fork of [https://github.com/fishfollower/SAM](SAM), where we model the age-dependent
parameters `logFpar` and `logSdLogObs` using smoothing splines instead of the standard manual
grouping method. The model extension is described in the paper "Adding smoothing splines to the SAM
model improves stock assessment".

The package can be installed by calling `make` from your terminal, when located inside this directory.

All code used to create the results in the paper is available in the folder `scripts/`. This
folder consists of three scripts. The first script, `model_fitting.R`, is used for fitting four
competing versions of the SAM model to various different fish stock data sets. Cross-validation and
forward-validation studies are also performed in this script. The second script,
`examination_of_results.R`, is used to examine the results of all the model fitting performed in
`model_fitting.R`. Finally, the script `helper_functions.R` consists of a set of helper functions
that are used in the two previously mentioned scripts.

