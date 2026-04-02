### file for the common paths files etc for all the sub models
### to be called by each submodel setup_NNN.R

### libraries
source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
# library(Rcpp, lib="/home/users/simon.brown/extremes/R/packages")
# library(evgam, lib="/home/users/simon.brown/extremes/R/packages")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
# library(ncdf4)
# library(PCICt)
library(evd)
library(float)
source("../libs/fn_JointMakeStationary.R")
source("../libs/fn_HotDays.R")
# source("../libs/fn_JointSimHD.R")
# source("../libs/lib_HotDay.R")

library(data.table)
library(PCICt)
library(gratia)
library(glue)
source("../science/copilot_fn.R")

