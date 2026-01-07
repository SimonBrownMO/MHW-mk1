# source("main_ExtractEvents.R")

st.pwd    <- system("pwd", intern=TRUE)

ensure you do not want to run main_ExtractEvents_Hobday.R with
    doExtractEvents_Hobday.R
    doExtractEvents_Hobday_wtime.R

### setup variables #####################################
source(paste(st.pwd,"/../setup_all.R",sep=''))
event.th.u    <- 0.90    # definition of when a heatwave event begins/ends
event.length  <- 180      # maximum allowed length of heatwave.  For wrt-peak this is fwd and bwd so events can be up to 2x longer than this (days)
USEHOTSEASON  <- TRUE
FINDHOTSEASON <- TRUE    # TRUE== determine hot season based on annual cycle greater than thHotSeason
thHotSeason   <- 5.0     # annual range degC must exceed this if a hot season is to be defined
                        # for regions with a strong annual cycle it is a good idea to exclude weather types that are not associated with hot days in an absolute sense
facHotSeas    <- 0.25   # fraction of the year that defines the hot season. 
                        # E.g. 0.5 = half the year, 0.25 = quarter of the year eg ~JJA
eventsDOY     <- NULL   # Specify here desired DOY for events to be extracted if find FINDHOTSEASON is FALSE
DOEVENTSPLOT  <- TRUE

events                <- list()
events$DOEVENTSPLOT   <- DOEVENTSPLOT
events$event.th.u     <- event.th.u
events$event.length   <- event.length
events$USEHOTSEASON   <- USEHOTSEASON
events$FINDHOTSEASON  <- FINDHOTSEASON
events$thHotSeason    <- thHotSeason

### currently missing from MSconfig
do.region <- "NWS"  

### v_Hobday
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"

### v1
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"

### v2
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-10-29.RData"

### v3
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v3/MSref/ostia_cdr_nrt_regions.MSref.2025-11-04.RData"


cat(cr,"MSref.file",cr)
cat(MSref.file,cr,cr)
load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

st_events <- sub('.RData','_EventsThMMM.RData',basename(st_base))
st_events <- paste(MSSAVEDIR,'Events/',sub('MMM',event.th.u,st_events),sep='')

# if(st_version=='v1') {
#     load(st_msdata01, verb=TRUE)
#     data01$doy <- data01$doy + 0.5
#     st_msdata01 <- sub('/MSdata01/','/v1/MSdata01/',st_msdata01)
#     st_qgam     <- sub('/Qgam/','/v1/Qgam/',         st_qgam)
#     MSSAVEDIR   <- sub('RRR',do.region, sub('VVV',  st_version,    MSSAVEDIR0))
# }

### setup HotDays variables #####################################
    # source(paste(st.pwd,"/setup_HotDays.R",sep=''))
    # cat("HotDays main: Completed setup_HotDays.R",cr,"#########################",cr,cr)

# if(st_version=='v1') {
#     st_msdata01 <- sub('/MSdata01/','/v1/MSdata01/',st_msdata01)
#     MSSAVEDIR   <- sub('RRR',do.region, sub('VVV',st_version,    MSSAVEDIR0))
# }


    load(st_msdata01, verb=TRUE)
    cat("Loaded msdata01 from ", st_msdata01, cr,cr)

    
### extract events #####################################
    source(paste(st.pwd,"/doExtractEvents.R",sep=''))

    load(st_events,   verb=TRUE)
    iobs <- which(data01$isobs==1)
    # if(DOEVENTSPLOT) fn_plotEvents(events01$obs, data01[iobs,], nplot=3, savefile=sub('Events/ukgd','Events/plots/ukgd',sub('.RData','.pdf',st_events)) )
    if(DOEVENTSPLOT) fn_plotEvents(events01$obs, data01[iobs,], nplot=3, savefile=NULL )
    cat("HotDays main: Completed doExtractEvents.R",cr,"#########################",cr,cr)
