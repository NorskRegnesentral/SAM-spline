# SAM

This is a fork of [SAM](https://github.com/fishfollower/SAM), where we model the age-dependent
parameters `logFpar` and `logSdLogObs` using smoothing splines instead of the standard manual
grouping method. The model extension is described in the paper "Adding smoothing splines to the SAM
model improves stock assessment", available on [arXiv](https://arxiv.org/abs/2502.20788).

The package can be installed by calling `make` from your terminal, when located inside this directory.

All code used to create the results in the paper is available in the folder `scripts/`. This
folder consists of three scripts. The first script, `model_fitting.R`, is used for fitting four
competing versions of the SAM model to various different fish stock data sets. Cross-validation and
forward-validation studies are also performed in this script. The second script,
`examination_of_results.R`, is used to examine the results of all the model fitting performed in
`model_fitting.R`. Finally, the script `helper_functions.R` consists of a set of helper functions
that are used in the two previously mentioned scripts.

Recent changes were made to the SAM package while we were developing our spline models. These recent
changes have been merged into this fork, and minor changes have been made to the code so that the
spline models will be usable with the newest version of SAM. However, the actual version of SAM that
was used for creating all the results in the paper, can be found by checking out tags starting with
"SAM-spline-paper", e.g., the tag "SAM-spline-paper-arxiv-v1", which points to the commit that was
used for creating all the results in the first published version of the paper on github. To checkout
this commit you can write `git checkout tags/SAM-spline-paper-arxiv-v1`.
