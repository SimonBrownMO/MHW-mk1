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
source("../../libs/fn_MakeStationary.R")
source("../../libs/fn_HotDays.R")
# source("../libs/fn_JointSimHD.R")
# source("../libs/lib_HotDay.R")

library(data.table)
library(PCICt)
library(gratia)
library(glue)
source("../../science/copilot_fn.R")


############################################################################################################
### general ################################################################################################
############################################################################################################
# INDIR_O     <- '/home/users/simon.brown/extremes/heatwaves/mhw/DATA/'
# MSSAVEDIR0  <- paste(INDIR_O,'RRR/VVV/',sep='')  # RRR=region, VVV=version

# # random constants
# cr          <- '\n'
# deg0C       <- 273.16    # conversion to celcius
# st_version  <- "v1"         # 'v2' 
# datestamp   <- "2025-10-15" # "2025-10-29" # paste(format(Sys.time(), "%Y-%m-%d"),sep='_') 

# ### pre-proc data ###############################################################
# DODIAGPRE         <- TRUE
# # set names to be used
# st_infile_o       <- 'ostia_cdr_nrt_regions.RData'

# ### global temperatures
# st_obs_gmst <- "/home/users/simon.brown/extremes/R/general/read_global_annual_temp.R"
# gmst_ref_period <- 1981:2000  
