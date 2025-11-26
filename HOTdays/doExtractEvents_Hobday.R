# source("doExtractEvents_Hobday.R")   

# library(qgam)

# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-05.RData"


# load(MSref.file, verb=TRUE)
# list2env(MSconfig , envir = .GlobalEnv)
# list2env(MSconfig$files , envir = .GlobalEnv)

# st_events <- sub('.RData','_EventsThMMM.RData',basename(st_base))
# st_events <- paste(MSSAVEDIR,'Events/',sub('MMM',event.th.u,st_events),sep='')
# st_events <- sub('.RData','_Hobday.RData', st_events)

# load(st_msdata01, verb=TRUE)

load(st_qgam, verb=TRUE)
th90    <- qdo(fit.qgam, newdata=data01, type='response', qu=0.9, predict)
th50    <- qdo(fit.qgam, newdata=data01, type='response', qu=0.5, predict)
th.qgam <- list(th50=th50, th90=th90)

# Hobday threshold calculation
# Hobday et al. (2016) define the threshold for a hot day as the 90th percentile
# of daily mean temperatures over a climatological period (here 1982-2011),
data_clim       <- data01
time_clim       <- clim.year + data01$time - trunc(data01$time) # middle of 1982:2111
stime_clim      <- (time_clim - mean(time_clim))/ ( max(trunc(data01$time)) - min(trunc(data01$time)) )
data_clim$stime <- stime_clim

th90_hob <- qdo(fit.qgam, newdata=data_clim, type='response', qu=0.9, predict)            
th50_hob <- qdo(fit.qgam, newdata=data_clim, type='response', qu=0.5, predict)            
th5090   <- list(th50=th50_hob, th90=th90_hob)

# sanity plot
r1 <- range(data01$x)
 plot(data01$time, data01$x, pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
lines(data01$time, th90_hob, col='red', lwd=1)
        readline("Continue?")

### extract events #####################################
# calculate anomalies
data01$anomaly <- data01$x - th90_hob

 plot(data01$time, data01$anomaly, pch=20, cex=.3, main=sub('RRR', do.region, 'Hobday anomalies RRR') )    
grid()

# compare my clim with segolens clim
load("../DATA/ostia_cdr_nrt_regions.RData", verb=T)
th90_clim_sego <-  l.mhw[[do.region]] - l.mhw_an[[do.region]]
i0 <- which(trunc(data01$time) >= 2000 & trunc(data01$time) <= 2004)
plot(data01$time[i0],th90_hob[i0]-th90_clim_sego[i0], pch=20, cex=.3)
grid()
        readline("Continue?")

r1 <- range(c(th90_clim_sego[i0], th90_hob[i0]), na.rm=TRUE)
plot(data01$time[i0],th90_clim_sego[i0], pch=20, cex=.3, ylim=r1, main='Hobday 90th percentile climatology')
lines(data01$time[i0],th90_hob[i0], col='red', lwd=2)
grid()
        readline("Continue?")

### SO the qgam 90th percentile clim is similar enough to Segolens climatology
### to allow comparison of event extraction, but my clim is smoother through
### the annual cycle and biased high wrt Segolens clim by 0.67  sd 0.16
### BUT dont have the 50th percentile climatology form Segolens to compare
### catagory definitions of Hobday


### extract events on Hobday anomalies
i0     <- which(data01$isobs==1)
ch_hob <- fn_extractEventsGeneric(data01$anomaly[i0], 0.0, event.length)

st.pdf.hseas   <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.pdf',basename(st_base)),sep='/')
if(!dir.exists(dirname(st.pdf.hseas))) system(paste("mkdir -p",dirname(st.pdf.hseas) ) )

hotseas        <- fn_findHotSeasonSingle(data01[i0,], thHotSeason, facHotSeas=facHotSeas, st.qgam=st_qgam, st.pre=st_preproc, st.pdf=st.pdf.hseas, DOPLOT=TRUE)
ch_hob$hotseas <- list(IsHotSeason=hotseas$IsHotSeason, hotSdoy=hotseas$hotSdoy$o, ihotseas=hotseas$ihotseas$o)
ch_hob$time      <- data01$time[i0]


if(USEHOTSEASON) {
    events01_ally <- list(obs=ch_hob)
    events01     <- events01_ally

    # obs
    hs_doy <- events01$obs$hotseas$hotSdoy
    # hs_doy <- events01$obs$hotseas$ihotseas
    ihs <- NULL
    # for(i in 1:dim(events01$obs$ch.sev)[2]) if (length(which(hs_doy %in% events01$obs$ch.doy[,i]))>0) ihs <- c(ihs,i)
    for(i in seq_along(events01$obs$ch.sev)) if (length(which(hs_doy %in% data01$doy[ch_hob$ch.i[,i]]))>0) ihs <- c(ihs,i)
                                                # if any day of the event is within the hot season then keep the whole event
    events01$obs$ch.pk$chains         <- events01$obs$ch.pk$chains[, ihs]
    events01$obs$ch.pk$absichains     <- events01$obs$ch.pk$absichains[ihs]
    events01$obs$ch.pk$ichains        <- events01$obs$ch.pk$ichains[, ihs]
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
    events01               <- list(obs=ch_hob)
    events01$events01_ally <- NULL
}
events01$info       <- events
events01$hotseas    <- hotseas

# events01.hob <- events01
if(!dir.exists(dirname(st_events))) system(paste("mkdir -p",dirname(st_events) ) )
cat("Saved events to :", st_events, cr, cr )
save(file=st_events, events01, th5090)


#