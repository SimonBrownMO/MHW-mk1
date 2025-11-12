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

if(!min(data01$u<ms.thgpd.u, na.rm=TRUE)){
    cat("Fixing data01$u",cr)
    i1 <- which(data01$u<=ms.thgpd.u)
    data01$u[i1] <- data01$uqgam[i1]
}

i0   <- which(data01$isobs==1)

ch_o <- fn_extractEvents(data01[i0,], event.th.u, event.length)

# ch_m <- fn_extractEvents(data01[-i0,], event.th.u, event.length)

# x.l   <- qlaplace(data01$u)
# ch.pk <- ChainsFromLaplace(tx2.l, dtx2,  seas_length, chain.lev.u=event.th.u, chain.length=event.length, DOPLOT=FALSE)


if(FINDHOTSEASON) {
    st.pdf.hseas <- paste(sub('/Events','/Hotseason',dirname(st_events)),sub('.RData','.hotseason.pdf',basename(st_base)),sep='/')
    hotseas <- fn_findHotSeasonSingle(data01, thHotSeason, facHotSeas=facHotSeas, st.qgam=st_qgam, st.pre=st_preproc, st.pdf=st.pdf.hseas, DOPLOT=TRUE)
} else hotseas <- list(IsHotSeason=FALSE)

ch_o$hotseas <- list(IsHotSeason=hotseas$IsHotSeason, hotSdoy=hotseas$hotSdoy$o, ihotseas=hotseas$ihotseas$o)

## add time to each
ch_o$time <- data01$time[i0]

if(USEHOTSEASON) {
    events01_ally <- list(obs=ch_o)
    events01     <- events01_ally

    # obs
    hs_doy <- events01$obs$hotseas$hotSdoy
    ihs <- NULL
    for(i in 1:dim(events01$obs$ch.doy)[2]) if (length(which(hs_doy %in% events01$obs$ch.doy[,i]))>0) ihs <- c(ihs,i)
                                                # if any day of the event is within the hot season then keep the whole event
    events01$obs$ch.pk$chains.l       <- events01$obs$ch.pk$chains.l[, ihs]
    events01$obs$ch.pk$absichains     <- events01$obs$ch.pk$absichains[ihs]
    events01$obs$ch.pk$ichains.l      <- events01$obs$ch.pk$ichains.l[, ihs]
    events01$obs$ch.pk$indeciesByYear <- events01$obs$ch.pk$indeciesByYear[ihs]
    events01$obs$ch.pk$thresh.u       <- events01$obs$ch.pk$thresh.u
    events01$obs$ch.pk$chain.length   <- events01$obs$ch.pk$chain.length
    events01$obs$ch.st.l              <- events01$obs$ch.st.l[,ihs]
    events01$obs$ch.st.o              <- events01$obs$ch.st.o[,ihs]
    events01$obs$ch.doy               <- events01$obs$ch.doy[,ihs]
    events01$obs$ch.age               <- events01$obs$ch.age[,ihs]
    events01$obs$ch.i                 <- events01$obs$ch.i[,ihs]
    events01$obs$ch.time              <- events01$obs$ch.time[ihs]
    events01$obs$ch.sev               <- events01$obs$ch.sev[ihs]
    events01$obs$ch.mean.sev          <- events01$obs$ch.mean.sev[ihs]
    events01$obs$ch.pkval             <- events01$obs$ch.pkval[ihs]
    events01$obs$ch.pkval.day         <- events01$obs$ch.pkval.day[ihs]
    events01$obs$ch.duration          <- events01$obs$ch.duration[ihs]
    events01$obs$ch.gmst              <- events01$obs$ch.gmst[ihs]
    events01$obs$ihs                  <- events01$obs$ihs

   
    events01$events01_ally <- events01_ally

} else {
    events01               <- list(obs=ch_o)
    events01$events01_ally <- NULL
}
events01$info       <- events
events01$hotseas    <- hotseas

if(!dir.exists(dirname(st_events))) system(paste("mkdir -p",dirname(st_events) ) )
cat("Saved events to :", st_events, cr )
save(file=st_events, events01)


#
