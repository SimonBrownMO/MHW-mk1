# source("hd_test.R")

# library(Rcpp, lib="/home/h03/hadsx/extremes/R/packages")
# library(evgam, lib="/home/h03/hadsx/extremes/R/packages")
# library(mgcv)
# library(qgam)
# library(ncdf4)
# library(PCICt)
# library(ismev)

source("/home/h03/hadsx/extremes/heatwaves/code/joint_mk1/fn_JointMakeStationary.R")
source("/home/h03/hadsx/extremes/heatwaves/code/joint_mk1/fn_JointHotDays.R")
source("/home/h03/hadsx/extremes/heatwaves/code/joint_mk1/lib_HotDay.R")
# source("/home/h03/hadsx/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")

# ### mk1
# DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/"
# reload_MS_ref(paste0(DATA,"MakeStationary/ukgd_cpm85_5k_x146y44.MSref.2023-06-20-160454.RData"))
# reload_HD_ref(paste0(DATA,"HotDays/vHD1/ukgd_cpm85_5k_x146y44.HDref.2023-06-21-171445.RData"))
# load(st_preproc, verb=T)
# load(st_msdata01, verb=T)
# load(st_events, verb=T)


### 2024.03
DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/"
load(paste0(DATA,"MakeStationary/ukgd_cpm85_5k_x146y44.MSref.2023-06-20-160454.RData", verb=TRUE))
> str1(msqgam)
List of 16
 $ coord            :List of 4
 $ DOQGAM           : logi TRUE
 $ ALLOWGPDSHAPETIME: logi FALSE
 $ DOMSPLOT         : logi TRUE
 $ TRIMTIME         : logi FALSE
 $ FUTURE_PROJECTION: logi TRUE
 $ do.ptiles        : num [1:99] 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 ...
 $ bs_cc            : chr "ds"
 $ ms.k             :List of 3
 $ fmla.MSqgam      :List of 2
 $ fmla.MSgpd       :List of 10
 $ chosen.MSgpd     : chr "SBDTi.G0"
 $ ms.thgpd.u       : num 0.94
 $ nmax_qgam        : num 20000
 $ gmst_ref_period  : int [1:20] 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 ...
 $ files            :List of 18
    $ INDIR_O    : chr "/data/users/haduk/uk_climate_data/supported/haduk-grid/v1.2.0.0/data/ceda/HadUK-Grid/v1.2.0.ceda/5km/tasmax/day/v20230328/"
    $ st_infile_o: chr "tasmax_hadukgrid_uk_5km_day*"
    $ INDIR_M    : chr "/project/ukcp/land-cpm/uk/5km/rcp85/01/tasmax/day/v20210615/"
    $ st_infile_m: chr "tasmax_rcp85_land-cpm_uk_5km_01_day*"
    $ st_version : chr "v_j2"
    $ MSSAVEDIR  : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/"
    $ st_base    : chr "ukgd_cpm85_5k_x146y44.RData"
    $ datestamp  : chr "2023-06-20-160454"
    $ stmetadata : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/HotDaySettings_2023-06-20-160454.R"
    $ st_msref   : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44.MSref.2023-0"| __truncated__
    $ st_preproc : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_preproc.RData"
    $ st_qgam    : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MSqgam.RData"
    $ st_ms      : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MS.RData"
    $ st_q2p     : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MSq2p.RData"
    $ st_msgpd   : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MSgpd.RData"
    $ st_msdata01: chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MSdata01.RData"
    $ st_cpm_gmst: chr "/home/h03/hadsx/extremes/tawn/jordan/code/cpm/ukcp18_cpm_global_tempearture_covariates.RData"
    $ st_obs_gmst: chr "/home/h03/hadsx/extremes/R/general/read_global_annual_temp.R"

load(paste0(DATA,"HotDays/vHD1/ukgd_cpm85_5k_x146y44.HDref.2023-06-21-171445.RData"), verb=TRUE)
> str(HDreference,2)
    List of 8
    $ file_names  :List of 12
    ..$ msfiles     :List of 18
    ..$ datestamp   : chr "2023-06-21-171445"
    ..$ stmetadata  : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/HotDays/vHD1//HotDaysSettings_2023-06-21-171445.R"
    ..$ st_hdref    : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/HotDays/vHD1//ukgd_cpm85_5k_x146y44.HDref.2023-06"| __truncated__
    ..$ st_base     : chr "ukgd_cpm85_5k_x146y44.RData"
    ..$ st_events   : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/HotDays/vHD1//ukgd_cpm85_5k_x146y44_EventsTh0.94.RData"
    ..$ st_InitEvent: chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MS.RData"
    ..$ st_InitValue: chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MS.RData"
    ..$ st_Ht       : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MS.RData"
    ..$ st_Term     : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/ukgd_cpm85_5k_x146y44_MS.RData"
    ..$ st_HwSim    : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/sim_laplace/ukgd_cpm85_5k_x146y44_MS.RData"
    ..$ st_t0       : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/sim_t0/ukgd_cpm85_5k_x146y44_MS.RData"
    $ hdinfo      :List of 30
    ..$ MSSAVEDIR        : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/MakeStationary/"
    ..$ HDSAVEDIR        : chr "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/HotDays/vHD1//"
    ..$ st_version       : chr "v_j2"
    ..$ FUTURE_PROJECTION: logi TRUE
    ..$ imodel.initE     : chr "initE.dc"
    ..$ imodel.initV     : chr "NOTUSED"
    ..$ HTCOV            : chr "NOTUSED"
    ..$ imodel.ht1       : chr "ht1.top"
    ..$ imodel.htN       : chr "ht.top"
    ..$ imodel.term1     : chr "term.day1.tdc"
    ..$ imodel.termN     : chr "term.tdac"
    ..$ USESPICESIM      : logi TRUE
    ..$ DOHTREFILL       : logi TRUE
    ..$ DOTERM           : logi TRUE
    ..$ CCINITDAY1       : logi TRUE
    ..$ CCINITDAYN       : logi TRUE
    ..$ CCGPD            : logi TRUE
    ..$ CCHTDAY1         : logi TRUE
    ..$ CCHTDAYN         : logi TRUE
    ..$ CCTERMDAY1       : logi TRUE
    ..$ CCTERMDAYN       : logi TRUE
    ..$ LIMITMINTERMPROB : logi TRUE
    ..$ mintermprob      : num [1:40] 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ...
    ..$ DOSIMPLOT        : logi TRUE
    ..$ DOT0PLOT         : logi TRUE
    ..$ DOSIMBYDAYPLOTS  : logi TRUE
    ..$ DEBUGSIM         : logi FALSE
    ..$ DOPRINTOBSTERM   : logi TRUE
    ..$ cr               : chr "\n"
    ..$ deg0C            : num 273
    $ events_info :List of 6
    ..$ DOEVENTSPLOT : logi TRUE
    ..$ hw.lev.u     : num 0.94
    ..$ event.length : num 40
    ..$ USEHOTSEASON : logi TRUE
    ..$ FINDHOTSEASON: logi TRUE
    ..$ thHotSeason  : num 6
    $ fitinitevent:List of 4
    ..$ SAVEINITEVENTDIAG: logi TRUE
    ..$ DOINITEVENTPLOT  : logi TRUE
    ..$ event.th.u       : num 0.94
    ..$ ievent.k         :List of 2
    $ fitinitvalue:List of 4
    ..$ SAVEINITVALUEDIAG: logi TRUE
    ..$ DOINITVALUEPLOT  : logi TRUE
    ..$ value.th.u       : num 0.94
    ..$ ivalue.k         :List of 2
    $ fitht       :List of 14
    ..$ SAVEHTDIAG               : logi TRUE
    ..$ DOHTPLOT                 : logi TRUE
    ..$ PARAMETERPLOT            : logi TRUE
    ..$ PARPLTNEW                : logi FALSE
    ..$ ht.th.u                  : num 0.94
    ..$ cr.th.u                  : num 0.94
    ..$ ht.cov.sig.th            : num 0.9
    ..$ lam.pen0                 : num 0
    ..$ FORCEALPHASTRUCTURE      : logi FALSE
    ..$ ALLOW_MU_SIGMA_DEPENDENCE: logi FALSE
    ..$ HTmaxorder               : num 1
    ..$ lagnumber                : num 1
    ..$ nsim1                    : num 20000
    ..$ ht.exclude.l             : num -2
    $ fitterm     :List of 10
    ..$ SAVETERMDIAG      : logi TRUE
    ..$ DOTERMPLOT        : logi TRUE
    ..$ FORCETERMPROBFALL : logi FALSE
    ..$ term.th.u         : num 0.94
    ..$ term.prob.fall.tol: num 0
    ..$ minprob.day1      : num 0.23
    ..$ minprob.dayN      : num 0.056
    ..$ topout.tmp.l      : num 12
    ..$ iterm.k           :List of 4
    ..$ mintermprob       : num [1:40] 0.27 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ...
    $ sim         :List of 13
    ..$ hw.simthresh.u  : num 0.94
    ..$ hw.simthresh.l  : num 2.12
    ..$ iorder          : num 1
    ..$ durationmax     : num 40
    ..$ nyears2sim      : num 1e+06
    ..$ nyears2sim_total: num 1e+06
    ..$ simcc           : num -0.334
    ..$ simyear         : chr "1850"
    ..$ USESPICESIM     : logi TRUE
    ..$ spicetmp        : chr "/scratch/hadsx/heatwave/tmp_sim"
    ..$ nspicejobs      : num 20
    ..$ spice.time      : num 300
    ..$ spice.memory    : num 10000
#

load(HDreference$file_names$msfile$st_msdata01, verb=TRUE)
    Loading objects:
    data01
    data01.std.param
    > str1(data01)
    'data.frame':	59011 obs. of  10 variables:
    $ x    : num  11.58 9.05 11.06 12.44 9.1 ...
    $ time : num  1960 1960 1960 1960 1960 ...
    $ stime: num  -0.5 -0.5 -0.5 -0.5 -0.5 ...
    $ doy  : int  1 2 3 4 5 6 7 8 9 10 ...
    $ sdoy : num  0.00273 0.00546 0.0082 0.01093 0.01366 ...
    $ isobs: num  1 1 1 1 1 1 1 1 1 1 ...
    $ class: Factor w/ 2 levels "obs","mod": 1 1 1 1 1 1 1 1 1 1 ...
    $ gmst : num  -0.341 -0.341 -0.341 -0.341 -0.341 ...
    $ uqgam: num  0.89 0.655 0.852 0.943 0.665 ...
    $ u    : num  0.89 0.655 0.852 0.945 0.665 ...

load(HDreference$file_names$st_events,verb=TRUE)
    Loading objects:
    events01
    > str1(events01)
    List of 4
    $ obs          :List of 15
    $ mod          :List of 15
    $ events01_ally:List of 2
    $ info         :List of 6

ht.th.u   <- HDreference$fitht$ht.th.u
cr.th.u   <- HDreference$fitht$cr.th.u
scale.doy <- list(obs=data01.std.param$doy$o_max, mod=data01.std.param$doy$m_max)
stcovs    <- c('Index', 'CC', 'DOY', 'Age')

htJ <- fn_fitJointHt(events01, ht.th.u, cr.th.u, scale.doy, IsJoint=TRUE)

## not needed
# SjbPlotHtParametersCov(htJ$htN,new=TRUE)

## not needed
#  lag01 <- fn_lag_events01(events01, scale.doy)
#     list2env(lag01$om,env=parent.frame())
# SjbPlotHtModelPoints(htJ$htN, om12[om.1E1,])

fn_htJointPlot(htJ, stplot="test-ht.pdf")

#
