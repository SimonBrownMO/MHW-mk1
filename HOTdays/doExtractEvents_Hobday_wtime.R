# source("doExtractEvents_Hobday_wtime.R")   

# library(qgam)

# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-05.RData"


# load(MSref.file, verb=TRUE)
# list2env(MSconfig , envir = .GlobalEnv)
# list2env(MSconfig$files , envir = .GlobalEnv)

# st_events <- sub('.RData','_EventsThMMM.RData',basename(st_base))
# st_events <- paste(MSSAVEDIR,'Events/',sub('MMM',event.th.u,st_events),sep='')
# st_events <- sub('.RData','_Hobday.RData', st_events)

load(st_qgam, verb=TRUE)
th90     <- qdo(fit.qgam, newdata=data01, type='response', qu=0.9, predict)
th50     <- qdo(fit.qgam, newdata=data01, type='response', qu=0.5, predict)
th5090 <- list(th50=th50, th90=th90)

# sanity plot
r1 <- range(data01$x)
 plot(data01$time, data01$x, pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
lines(data01$time, th90, col='red', lwd=1)
        readline("Continue?")

### extract events #####################################
# calculate anomalies
data01$anomaly <- data01$x - th90

 plot(data01$time, data01$anomaly, pch=20, cex=.3, main=sub('RRR', do.region, 'Hobday anomalies RRR') )    
grid()
        readline("Continue?")

### extract events on Hobday anomalies
i0     <- which(data01$isobs==1)
ch_hobwt <- fn_extractEventsGeneric(data01$anomaly[i0], 0.0, event.length)

st.pdf.hseas   <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.pdf',basename(st_base)),sep='/')
if(!dir.exists(dirname(st.pdf.hseas))) system(paste("mkdir -p",dirname(st.pdf.hseas) ) )

hotseas        <- fn_findHotSeasonSingle(data01[i0,], thHotSeason, facHotSeas=facHotSeas, st.qgam=st_qgam, st.pre=st_preproc, st.pdf=st.pdf.hseas, DOPLOT=TRUE)
ch_hobwt$hotseas <- list(IsHotSeason=hotseas$IsHotSeason, hotSdoy=hotseas$hotSdoy$o, ihotseas=hotseas$ihotseas$o)
ch_hobwt$time      <- data01$time[i0]


if(USEHOTSEASON) {
    events01_ally <- list(obs=ch_hobwt)
    events01     <- events01_ally

    # obs
    hs_doy <- events01$obs$hotseas$hotSdoy
    # hs_doy <- events01$obs$hotseas$ihotseas
    ihs <- NULL
    # for(i in 1:dim(events01$obs$ch.sev)[2]) if (length(which(hs_doy %in% events01$obs$ch.doy[,i]))>0) ihs <- c(ihs,i)
    for(i in seq_along(events01$obs$ch.sev)) if (length(which(hs_doy %in% data01$doy[ch_hobwt$ch.i[,i]]))>0) ihs <- c(ihs,i)
                                                # if any day of the event is within the hot season then keep the whole event
    events01$obs$ch.pk$chains       <- events01$obs$ch.pk$chains[, ihs]
    events01$obs$ch.pk$absichains     <- events01$obs$ch.pk$absichains[ihs]
    events01$obs$ch.pk$ichains      <- events01$obs$ch.pk$ichains[, ihs]
    events01$obs$ch.pk$indeciesByYear <- events01$obs$ch.pk$indeciesByYear[ihs]
    events01$obs$ch.pk$thresh.u       <- events01$obs$ch.pk$thresh.u
    events01$obs$ch.pk$chain.length   <- events01$obs$ch.pk$chain.length
    events01$obs$ch.st.l              <- events01$obs$ch.st.l[,ihs]
    events01$obs$ch.st.o              <- events01$obs$ch.st.o[,ihs]
    # events01$obs$ch.doy               <- events01$obs$ch.doy[,ihs]
    events01$obs$ch.age               <- events01$obs$ch.age[,ihs]
    events01$obs$ch.i                 <- events01$obs$ch.i[,ihs]
    events01$obs$ch.time              <- events01$obs$ch.time[ihs]
    events01$obs$ch.sev               <- events01$obs$ch.sev[ihs]
    events01$obs$ch.mean.sev          <- events01$obs$ch.mean.sev[ihs]
    events01$obs$ch.pkval             <- events01$obs$ch.pkval[ihs]
    events01$obs$ch.pkval.day         <- events01$obs$ch.pkval.day[ihs]
    events01$obs$ch.duration          <- events01$obs$ch.duration[ihs]
    # events01$obs$ch.gmst              <- events01$obs$ch.gmst[ihs]
    events01$obs$ihs                  <- events01$obs$ihs
   
    events01$events01_ally <- events01_ally

} else {
    events01               <- list(obs=ch_hobwt)
    events01$events01_ally <- NULL
}
events01$info       <- events
events01$hotseas    <- hotseas

# events01.hobwt <- events01
if(!dir.exists(dirname(st_events))) system(paste("mkdir -p",dirname(st_events) ) )
cat("Saved events to :", st_events, cr, cr )
save(file=st_events, events01, th5090)


#