# source("setup_MakeStationary.R")

### setup varaibles to perform MakeStationary

### libraries
source("../setup_all.R")

############################################################################################################
### general ################################################################################################
############################################################################################################
INDIR_O     <- '/home/users/simon.brown/extremes/heatwaves/mhw/DATA/'
MSSAVEDIR0  <- paste(INDIR_O,'RRR/VVV/',sep='')  # RRR=region, VVV=version

# random constants
cr          <- '\n'
deg0C       <- 273.16    # conversion to celcius

st_version  <- 'v1' # "v_Hobday"         # 
    ## v_Hobday: match Hobday as best we can, event threshold 0.90, climC term linear with time
    #ms.k              <- list(doy=12, gmst=4)
    #fmla.MSqgam       <- list(x ~ stime  +s(sdoy, bs='cc',k=ms.k$doy)  , ~ s(sdoy))
    #chosen.MSgpd.name <- "Hobday" #  best match to hobday
    #clim.year         <- 1996  # middle year of climatology period

    # v1: as v_Hobday but allowing linear trend with time in the annual cycle
    ms.k              <- list(doy=12, gmst=4)
    fmla.MSqgam       <- list(x ~ stime  +s(sdoy, bs='cc',k=ms.k$doy)  , ~ s(sdoy))
    chosen.MSgpd.name <- "Hobday" #  best match to hobday
    clim.year         <- NULL  # no fixed climatology

    # v3: very smooth gmst covariate
    # ms.k         <- list(doy=12, gmst=4)
    # fmla.MSqgam  <- list(x ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))

datestamp   <- "2025-11-07" # "2025-10-29" # paste(format(Sys.time(), "%Y-%m-%d"),sep='_') 

### pre-proc data ###############################################################
DODIAGPRE         <- TRUE
# set names to be used
st_infile_o       <- 'ostia_cdr_nrt_regions.RData'

### global temperatures
st_obs_gmst <- "/home/users/simon.brown/extremes/R/general/read_global_annual_temp.R"
gmst_ref_period <- 1981:2000  


### model and location specifics ################################################################################
MSSAVEDIR         <- sub('RRR',do.region, sub('VVV',st_version,    MSSAVEDIR0))

subdir            <- c('Preproc/plots','Qgam/plots','Evgam/plots','Q2p/plots','P2q/plots','MSdata01/plots', 'MSref')
for(d1 in subdir) system(paste("mkdir -p",paste(MSSAVEDIR,d1,sep='') ) )
st_base           <- paste(INDIR_O,st_infile_o, sep='')

### create required names and directories (used later in this file so needs to be here)
stmetadata      <- paste(MSSAVEDIR,paste(paste('HotDaySettings',datestamp,sep='_'),sep=''),'.R',sep='')
st_msref        <- sub('DDD',datestamp, sub('regions.RData','regions.MSref.DDD.RData',paste(MSSAVEDIR,'MSref/',  basename(st_base),sep='')))
st_preproc      <-                      sub('regions.RData','regions.preproc.RData',  paste(MSSAVEDIR,'Preproc/',basename(st_base),sep=''))

############################################################################################################
### MakeStationaryQgam.R ###################################################################################
DOQGAM            <- TRUE  # multithresholds is slow.  If previously run set this to false as they were saved
DODIAGQGAM        <- TRUE
DODIAGMSGPD       <- TRUE
DODIAGQ2P         <- TRUE
do.ptiles         <- seq(from=1, to=99, by=1)/100.0     

###  qgam ###################################################################################
bs.cc             <- 'tp' # form of smoother for climate change term: tp, ts, ds, cr, cs  see smooth.terms in mgcv
                          # ds seems to be the best but slowest, then ts, then cs, last tp & cr equal poorest
                          # ds requires ms_cc_smooth>=4 but ms_cc_smooth=3 seems best for obs
#top ms.k         <- list(doy=12, gmst=4)  #  list(doy=6, gmst=5)  #  -1 means use defaults (see ?mgcv::s)
                          #  for doy use k=10 (5-15 ok)  for gmst use k=4 (3-7 ok)
                          #  for interaction term ti(sdoy,gmst) use k=c(10,4)
                          #  if using bs='ds' then ms_cc_smooth must be >=4
                          #  if using bs='ts' then ms_cc_smooth must be >=3
                          #  if using bs='cs' then ms_cc_smooth must be >=2
                          #  if using bs='tp' or 'cr' then ms_cc_smooth must be >=1
### CHECK THIS BIT VERY CAREFULLY IF CHANGING
#top fmla.MSqgam  <- list(x ~ s(sdoy, bs='cc',k=ms.k$doy) +s(gmst, bs=bs.cc, k=ms.k$gmst) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))
nmax_qgam    <- 20000   # 20000. limit the number of points that qgam tries to fit quantiles to

###  makeHTdata #############################################################################
st_htdata_qgam  <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_qgam.RData',  basename(st_base)),sep='' )
st_htdata_evgam <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_evgam.RData', basename(st_base)),sep='' )
st_htdata_q2p   <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_q2p.RData',   basename(st_base)),sep='' )
st_htdata_p2q   <- paste(MSSAVEDIR,'HTdata/',    sub('.RData','_p2q.RData',   basename(st_base)),sep='' )

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


# chosen.MSgpd.name <- "SDTi.G0" #  OR NULL if want to dynamically select best fit.
                          # better to decide what form MSgpd should take offline - see: explore_MS_evgam.R

st_qgam     <- paste(MSSAVEDIR,'Qgam/',    sub('.RData','_MSqgam.RData',  basename(st_base)),sep='' )
st_q2p      <- paste(MSSAVEDIR,'Q2p/',     sub('.RData','_MSq2p.RData',   basename(st_base)),sep='' )
st_p2q      <- paste(MSSAVEDIR,'P2q/',     sub('.RData','_MSp2q.RData',   basename(st_base)),sep='' )
st_msgpd    <- paste(MSSAVEDIR,'Evgam/',   sub('.RData','_MSgpd.RData',   basename(st_base)),sep='' )
st_msdata01 <- paste(MSSAVEDIR,'MSdata01/',sub('.RData','_MSdata.RData',basename(st_base)),sep='' )

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
    st_infile_o = st_base,
    st_version  = st_version,
    MSSAVEDIR   = MSSAVEDIR,
    st_base     = st_base,
    datestamp   = datestamp,
    stmetadata  = stmetadata,
    st_msref    = st_msref,
    st_preproc  = st_preproc,
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
