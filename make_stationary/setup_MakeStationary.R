# source("setup_MakeStationary.R")

### setup varaibles to perform MakeStationary

### libraries
source("../setup_all.R")

############################################################################################################
### general ################################################################################################
############################################################################################################
INDIR_O     <- '/home/users/simon.brown/extremes/heatwaves/mhw/DATA/'
MSSAVEDIR0  <- paste(INDIR_O,'RRR/VVV/',sep='')  # RRR=region, VVV=version
INDIR_M     <- '/data/users/simon.brown/extremes/mhw/projections/'

# random constants
cr          <- '\n'
deg0C       <- 273.16    # conversion to celcius

if(st_version=='v_Hobday') {
    # v_Hobday: match Hobday as best we can, event threshold 0.90, climC term linear with time
    #           uses doExtractEvents_Hobday.R
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ stime  +s(sdoy, bs='cc',k=ms.k$doy)  , ~ s(sdoy))
    chosen.MSgpd.name <- "Hobday" #  best match to hobday
    clim.year         <- 1996  # middle year of climatology period
} else if(st_version=='v1') {
    # v1: as v_Hobday but allowing linear trend with time but fixed annual cycle
    #     uses doExtractEvents_Hobday_wtime.R
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ stime  +s(sdoy, bs='cc',k=ms.k$doy)  , ~ s(sdoy))
    chosen.MSgpd.name <- "Hobday" #  best match to hobday
    clim.year         <- NULL  # no fixed climatology
} else if(st_version=='v2') {
    # v2: as for v1 but allowing linear trend with GMST but fixed annual cycle
    #     uses doExtractEvents_Hobday_wtime.R
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ gmst  +s(sdoy, bs='cc',k=ms.k$doy)  , ~ s(sdoy))
    chosen.MSgpd.name <- "DLG" # doy + linear gmst
    clim.year         <- NULL  # no fixed climatology
} else if(st_version=='v3') {
    # standard LST model for reference - not currently advocating it
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ s(sdoy, bs="cc", k=ms.k$doy) + s(gmst, bs='tp', k=ms.k$gmst) +    ti(sdoy, gmst, bs=c("cc", "tp"))  , ~ s(sdoy))
    chosen.MSgpd.name <- "DLG" # doy + linear gmst
    clim.year         <- NULL  # no fixed climatology
} else if(st_version=='v4') {
    # standard LST model for reference - not currently advocating it
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ NA  , ~ NA)
    chosen.MSgpd.name <- "NA" # doy + linear gmst
    clim.year         <- NULL  # no fixed climatology
} else if(st_version=='v5') {
    # NB now adding in climae model data
    # multistep approach - remove annual cycle and climate change term from the mean, then QGAM, then EVGAM
        ms.k              <- list(doy=12, gmst=4)
        fmla.MSqgam       <- list(x ~ NA  , ~ NA)
        chosen.MSgpd.name <- "NA" # doy + linear gmst
        clim.year         <- NULL  # no fixed climatology

    DOSTEP1 <- TRUE
    STEP1   <- "ONESTEP" # "TWOSTEP" # or 
    # if(STEP1 == "TWOSTEP") {
        fmla1A  <- as.formula(paste('x      ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs="ts", k=7, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))',sep=''))
        fmla1B  <- as.formula(paste('resid1A ~ s(stime, bs="ts", k=2*ny.om, by=ftype)',sep=''))
    # } else if(STEP1 == "ONESTEP") {
        fmla0   <- as.formula(paste('x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs="ts", k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYearOM, bs="sz")',sep=''))
        fmla1AB <- as.formula(paste('x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs="ts", k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYearOM, bs="sz", sp=sp.sz)',sep=''))
    # }   
    fmla2   <- NULL #as.formula(paste('resid1 ~ s(stime, bs="ts", k=2*ny.om, by=ftype)',sep=''))
    fmla3   <- NULL #as.formula(paste('resid2 ~ s(stime, bs="ts", k=2*ny.om, by=ftype)',sep=''))

} else {
    stop("st_version not recognised")
}

# # v3: very smooth spline on gmst covariate with interaction term with sdoy
#  NOT CURRENTLY USED
# x ~ s(sdoy, bs = "cc", k = ms.k$doy) + s(gmst, bs = bs.cc, k = ms.k$gmst) +    ti(sdoy, gmst, bs = c("cc", "tp"))
# # seems to be just fitting to interannual noise
# ms.k              <- list(doy=12, gmst=4)
# fmla.MSqgam       <- list(x ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))
# chosen.MSgpd.name <- "SDTi.G0" #  match qgam
# clim.year         <- NULL  # no fixed climatology


############################################################################################################
### pre-proc data ##########################################################################################
DODIAGPRE    <- TRUE
# set names to be used
st_infile_o  <- paste(INDIR_O,'ostia_cdr_nrt_regions.RData', sep='')
st_infile_m  <- paste(INDIR_M,'1980-2100-Region_mean_timeseries_Daily_r001i1p00000.RData', sep='')
st_base      <- basename(st_infile_o)


### global temperatures
st_obs_gmst     <- "/home/users/simon.brown/extremes/R/general/read_global_annual_temp.R"
st_mod_gmst     <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/cpm_gmst.RData"  # cpm_gmst$m001$gmst
gmst_ref_period <- 1981:2000  


### model and location specifics ################################################################################
MSSAVEDIR         <- sub('RRR',do.region, sub('VVV',st_version,    MSSAVEDIR0))

subdir            <- c('Preproc/plots','step1/plots','Qgam/plots','Evgam/plots','Q2p/plots','P2q/plots','MSdata01/plots', 'MSref')
for(d1 in subdir) system(paste("mkdir -p",paste(MSSAVEDIR,d1,sep='') ) )

### create required names and directories (used later in this file so needs to be here)
stmetadata      <- paste(MSSAVEDIR,paste(paste('HotDaySettings',datestamp,sep='_'),sep=''),'.R',sep='')
st_msref        <- sub('DDD',datestamp, sub('regions.RData','regions.MSref.DDD.RData',paste(MSSAVEDIR,'MSref/',  st_base,sep='')))
st_preproc      <-                      sub('regions.RData','regions.preproc.RData',  paste(MSSAVEDIR,'Preproc/',st_base,sep=''))

############################################################################################################
### doMakeStationary     ###################################################################################
DOQGAM            <- TRUE  # multithresholds is slow.  If previously run set this to false as they were saved
DODIAGQGAM        <- TRUE
DODIAGMSGPD       <- TRUE
DODIAGQ2P         <- TRUE
do.ptiles         <- seq(from=1, to=99, by=1)/100.0     

###  qgam ###################################################################################
do.ptiles  <- seq(from=1, to=99, by=2)/100.0     #   seq(from=1, to=99, by=1)/100.0
bs.cc      <- 'ts' # form of smoother for climate change term: tp, ts, ds, cr, cs  see smooth.terms in mgcv
ms.k       <- list(doy=9, gmst=4)
fmla.qgam  <- list(residStep1 ~ ftype +s(sdoy, bs='cc', k=ms.k$doy,by=ftype) +s(gmst,  bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst,  bs=c('cc',bs.cc)), ~ s(sdoy) +ftype)
nmax_qgam  <- 20000   # 20000. limit the number of points that qgam tries to fit quantiles to

###  makeHTdata #############################################################################
st_htdata_qgam  <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_qgam.RData',  st_base),sep='' )
st_htdata_evgam <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_evgam.RData', st_base),sep='' )
st_htdata_q2p   <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_q2p.RData',   st_base),sep='' )
st_htdata_p2q   <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_p2q.RData',   st_base),sep='' )

### evgam ###################################################################################
ms.thgpd.u             <- 0.90
do.ptiles              <- sort(unique(c(do.ptiles,0.50,ms.thgpd.u)))
fmla.MSgpd             <- list()
fmla.MSgpd$SDTi.GDTi <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) , ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) )
fmla.MSgpd$SDTi.GDT  <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) , ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst)  )
fmla.MSgpd$SDTi.GD   <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) , ~ s(sdoy, bs='cc',k=ms.k$doy)   )
fmla.MSgpd$SDTi.G0   <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) , ~ 1 )
fmla.MSgpd$SDT.G0    <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst)                                   , ~ 1 )
fmla.MSgpd$SD.G0     <- list(excess ~ s(sdoy, bs='cc',k=ms.k$doy)                                                                   , ~ 1 )
fmla.MSgpd$ST.G0     <- list(excess ~                             +s(gmst, bs=bs.cc, k=ms.k$gmst)                                   , ~ 1 )
fmla.MSgpd$S0.G0     <- list(excess ~ 1                                                                                             , ~ 1 )
fmla.MSgpd$Hobday    <- list(excess ~ stime  +s(sdoy, bs='cc',k=ms.k$doy)                                                           , ~ 1 )
fmla.MSgpd$DLG       <- list(excess ~ gmst   +s(sdoy, bs='cc',k=ms.k$doy)                                                           , ~ 1 )
fmla.MSgpd$S0DTi.G0  <- list(excess ~ ftype +s(sdoy, bs='cc',k=ms.k$doy,by=ftype) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc',bs.cc)) , ~ 1 )

chosen.MSgpd.name <- "S0DTi.G0" #  OR NULL if want to dynamically select best fit.
                          # better to decide what form MSgpd should take offline - see: explore_MS_evgam.R

st_step1    <- paste(MSSAVEDIR,'step1/',    sub('.RData','_MSstep1.RData',  st_base),sep='' )
st_qgam     <- paste(MSSAVEDIR,'Qgam/',    sub('.RData','_MSqgam.RData',  st_base),sep='' )
st_q2p      <- paste(MSSAVEDIR,'Q2p/',     sub('.RData','_MSq2p.RData',   st_base),sep='' )
st_p2q      <- paste(MSSAVEDIR,'P2q/',     sub('.RData','_MSp2q.RData',   st_base),sep='' )
st_msgpd    <- paste(MSSAVEDIR,'Evgam/',   sub('.RData','_MSgpd.RData',   st_base),sep='' )
st_msdata01 <- paste(MSSAVEDIR,'MSdata01/',sub('.RData','_MSdata.RData',  st_base),sep='' )

MSconfig                   <- list()
MSconfig$DOQGAM            <- DOQGAM
MSconfig$DODIAGPRE         <- DODIAGPRE
MSconfig$DODIAGQGAM        <- DODIAGQGAM
MSconfig$DODIAGMSGPD       <- DODIAGMSGPD
MSconfig$DODIAGQ2P         <- DODIAGQ2P

MSconfig$do.ptiles         <- do.ptiles
MSconfig$bs.cc             <- bs.cc
MSconfig$ms.k              <- ms.k
MSconfig$fmla.MSqgam       <- fmla.MSqgam
MSconfig$fmla.MSgpd        <- fmla.MSgpd
MSconfig$chosen.MSgpd.name <- chosen.MSgpd.name
MSconfig$ms.thgpd.u        <- ms.thgpd.u
MSconfig$nmax_qgam         <- nmax_qgam
MSconfig$gmst_ref_period   <- gmst_ref_period

MSconfig$st_version      <- st_version     
MSconfig$datestamp       <- datestamp      
MSconfig$do.region       <- do.region   

MSconfig$files <- list(

    INDIR_O     = INDIR_O,
    st_infile_o = st_infile_o,
    st_infile_m = st_infile_m,
    st_version  = st_version,
    MSSAVEDIR   = MSSAVEDIR,
    st_base     = st_base,
    datestamp   = datestamp,
    stmetadata  = stmetadata,
    st_msref    = st_msref,
    st_preproc  = st_preproc,
    st_step1    = st_step1,
    st_qgam     = st_qgam,
    st_q2p      = st_q2p,
    st_p2q      = st_p2q,
    st_msgpd    = st_msgpd,
    st_msdata01 = st_msdata01,
    st_obs_gmst = st_obs_gmst,
    st_htdata_qgam  = st_htdata_qgam,
    st_htdata_evgam = st_htdata_evgam,
    st_htdata_q2p   = st_htdata_q2p,
    st_htdata_p2q   = st_htdata_p2q
)

cat("File definitions are:",cr)
print(str(MSconfig$files))

save(file=MSconfig$files$st_msref, MSconfig)
cat("MSconfig saved to :", MSconfig$files$st_msref, cr)

cat(cr)

cat("###############################################################################",cr)
sysdate   <- Sys.time()
stsysdate <- print(sysdate)
cat("All Done setup_MakeStationary.R",cr)
cat("###############################################################################",cr,cr)
