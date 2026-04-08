# source("doExtractEvents.R")

# > str(data01)
# 'data.frame':	59011 obs. of  7 variables:
#  $ x    : num  11.58 9.05 11.06 12.44 9.1 ...
#  $ time : num  -0.5 -0.5 -0.5 -0.5 -0.5 ...
#  $ doy  : num  0.00273 0.00546 0.0082 0.01093 0.01366 ...
#  $ isobs: num  1 1 1 1 1 1 1 1 1 1 ...
#  $ class: Factor w/ 2 levels "obs","cpm": 1 1 1 1 1 1 1 1 1 1 ...
#  $ uqgam: num  0.918 0.656 0.868 0.959 0.666 ...
#  $ u    : num  0.918 0.656 0.868 0.965 0.666 ...

# ? redundant
# if(!min(data01$u<ms.thgpd.u, na.rm=TRUE)){
#     cat("Fixing data01$u",cr)
#     i1 <- which(data01$u<=ms.thgpd.u)
#     data01$u[i1] <- data01$uqgam[i1]
# }

i0   <- which(data01$isobs==1)
yo1  <- trunc(median(unique(trunc(data01$time[i0]))))   # find a single yearin the middle  
ym1  <- trunc(median(unique(trunc(data01$time[-i0]))))  # find a single yearin the middle
iyo1 <- which(trunc(data01$time)==yo1 & data01$isobs==1)
iym1 <- which(trunc(data01$time)==ym1 & data01$isobs==0)

if(FINDHOTSEASON) {
    # st.pdf.hseas <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.pdf',basename(st_base)),sep='/')
    # hotseas <- fn_findHotSeasonSingle(data01, thHotSeason, fracHotSeas=fracHotSeas, st.qgam=st_qgam, st.pre=st_preproc, st.pdf=st.pdf.hseas, DOPLOT=TRUE)
    # hotseas <- fn_findHotSeason(data01, thHotSeason, st.qgam=st_qgam, st.pre=st_preproc, st.pdf='test.pdf', DOPLOT=TRUE)
    st.pdf.hseas <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.obs.pdf',basename(st_base)),sep='/')
    ihseas.o     <- fn_findHotSeasonSimple(data01$ft1A[iyo1], thHotSeason, fracHotSeas=0.5, DOPLOT=TRUE, st.pdf=st.pdf.hseas)
    hseas.doy.o  <- data01$doy[iyo1][ihseas.o]
    i1.o         <- which(data01$doy %in% hseas.doy.o & data01$isobs==1 )
    hotseas.o    <- list(IsHotSeason=TRUE, hotSdoy=hseas.doy.o, ihotseas=i1.o)

    st.pdf.hseas <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.mod.pdf',basename(st_base)),sep='/')
    ihseas.m     <- fn_findHotSeasonSimple(data01$ft1A[iym1], thHotSeason, fracHotSeas=0.5, DOPLOT=TRUE, st.pdf=st.pdf.hseas)
    hseas.doy.m  <- data01$doy[iym1][ihseas.m]
    i1.m         <- which(data01$doy %in% hseas.doy.m & data01$isobs==0 )
    hotseas.m    <- list(IsHotSeason=TRUE, hotSdoy=hseas.doy.m, ihotseas=i1.m)

    hotseas      <- list(IsHotSeason=TRUE, hotSdoy=list(o=hseas.doy.o, m=hseas.doy.m), ihotseas=list(o=i1.o, m=i1.m) )

} else hotseas <- list(IsHotSeason=FALSE)

ch_o         <- fn_extractEvents(data01[i0,],  event.th.u, event.length)
ch_o$hotseas <- hotseas.o
ch_o$time    <- data01$time[i0] # is this needed?

ch_m         <- fn_extractEvents(data01[-i0,], event.th.u, event.length)
# NEED TO FIX THE INDECIES OF THE MODEL EVENTS AS THEY ARE realtive to data01[-i0,] in the fn_extractEvents call
ch_m$ch.i    <- ch_m$ch.i + length(i0) 
ch_m$hotseas <- hotseas.m
ch_m$time    <- data01$time[-i0] # is this needed?

# x.l   <- qlaplace(data01$u)
# ch.pk <- ChainsFromLaplace(tx2.l, dtx2,  seas_length, chain.lev.u=event.th.u, chain.length=event.length, DOPLOT=FALSE)

if(USEHOTSEASON) {
    events01_ally <- list(obs=ch_o, mod=ch_m)
    events01      <- events01_ally

    # obs
    ihs.o <- NULL
    for(i in 1:length(ch_o$ch.sev)) if (length(which(hotseas$hotSdoy$o %in% ch_o$ch.doy[,i]))>0) ihs.o <- c(ihs.o,i)
                                                # if any day of the event is within the hot season then keep the whole event
    events01$obs$ch.pk$chains.l       <- events01$obs$ch.pk$chains.l[, ihs.o]
    events01$obs$ch.pk$absichains     <- events01$obs$ch.pk$absichains[ihs.o]
    events01$obs$ch.pk$ichains.l      <- events01$obs$ch.pk$ichains.l[, ihs.o]
    events01$obs$ch.pk$indeciesByYear <- events01$obs$ch.pk$indeciesByYear[ihs.o]
    events01$obs$ch.pk$thresh.u       <- events01$obs$ch.pk$thresh.u
    events01$obs$ch.pk$chain.length   <- events01$obs$ch.pk$chain.length
    events01$obs$ch.st.l              <- events01$obs$ch.st.l[,ihs.o]
    events01$obs$ch.st.o              <- events01$obs$ch.st.o[,ihs.o]
    events01$obs$ch.doy               <- events01$obs$ch.doy[,ihs.o]
    events01$obs$ch.sdoy              <- events01$obs$ch.sdoy[,ihs.o]
    events01$obs$ch.age               <- events01$obs$ch.age[,ihs.o]
    events01$obs$ch.i                 <- events01$obs$ch.i[,ihs.o]
    events01$obs$ch.time              <- events01$obs$ch.time[ihs.o]
    events01$obs$ch.sev               <- events01$obs$ch.sev[ihs.o]
    events01$obs$ch.mean.sev          <- events01$obs$ch.mean.sev[ihs.o]
    events01$obs$ch.pkval             <- events01$obs$ch.pkval[ihs.o]
    events01$obs$ch.pkval.day         <- events01$obs$ch.pkval.day[ihs.o]
    events01$obs$ch.duration          <- events01$obs$ch.duration[ihs.o]
    events01$obs$ch.gmst              <- events01$obs$ch.gmst[ihs.o]
    events01$obs$ihs                  <- events01$obs$ihs.o

    # mod
    ihs.m <- NULL
    for(i in 1:length(ch_m$ch.sev)) if (length(which(hotseas$hotSdoy$m %in% ch_m$ch.doy[,i]))>0) ihs.m <- c(ihs.m,i)
                                                # if any day of the event is within the hot season then keep the whole event
    events01$mod$ch.pk$chains.l       <- events01$mod$ch.pk$chains.l[, ihs.m]
    events01$mod$ch.pk$absichains     <- events01$mod$ch.pk$absichains[ihs.m]
    events01$mod$ch.pk$ichains.l      <- events01$mod$ch.pk$ichains.l[, ihs.m]
    events01$mod$ch.pk$indeciesByYear <- events01$mod$ch.pk$indeciesByYear[ihs.m]
    events01$mod$ch.pk$thresh.u       <- events01$mod$ch.pk$thresh.u
    events01$mod$ch.pk$chain.length   <- events01$mod$ch.pk$chain.length
    events01$mod$ch.st.l              <- events01$mod$ch.st.l[,ihs.m]
    events01$mod$ch.st.o              <- events01$mod$ch.st.o[,ihs.m]
    events01$mod$ch.doy               <- events01$mod$ch.doy[,ihs.m]
    events01$mod$ch.sdoy              <- events01$mod$ch.sdoy[,ihs.m]
    events01$mod$ch.age               <- events01$mod$ch.age[,ihs.m]
    events01$mod$ch.i                 <- events01$mod$ch.i[,ihs.m]
    events01$mod$ch.time              <- events01$mod$ch.time[ihs.m]
    events01$mod$ch.sev               <- events01$mod$ch.sev[ihs.m]
    events01$mod$ch.mean.sev          <- events01$mod$ch.mean.sev[ihs.m]
    events01$mod$ch.pkval             <- events01$mod$ch.pkval[ihs.m]
    events01$mod$ch.pkval.day         <- events01$mod$ch.pkval.day[ihs.m]
    events01$mod$ch.duration          <- events01$mod$ch.duration[ihs.m]
    events01$mod$ch.gmst              <- events01$mod$ch.gmst[ihs.m]
    events01$mod$ihs                  <- events01$mod$ihs.m

    events01$events01_ally <- events01_ally

} else {
    events01               <- list(obs=ch_o, mod=ch_m)
    events01$events01_ally <- NULL
}
events01$info       <- events
events01$hotseas    <- hotseas

if(!dir.exists(dirname(st_events))) system(paste("mkdir -p",dirname(st_events) ) )
cat("Saved events to :", st_events, cr )
save(file=st_events, events01)


#
