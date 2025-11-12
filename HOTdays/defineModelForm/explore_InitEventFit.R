# source("explore_InitEventFit.R")

### Use this file to determine the form of Initialising an event model

st.rootd <- "/home/users/simon.brown/extremes/heatwaves/mhw"
source("setup_dev.R")
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-10-15.RData"
reload_MS_ref(MSref.file)
# list2env(MSconfig , envir = .GlobalEnv)
# list2env(MSconfig$files , envir = .GlobalEnv)
load(st_msdata01, verb=TRUE)

st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.94.RData"
load(st_events,verb=TRUE)
list2env(events01$info , envir = .GlobalEnv)

############################################################################################################
### doFitInitEvent.R #########################################################################################
# source("../setup_HotDays.R")
SAVEINITEVENTDIAG    <- TRUE
ievent.th.u          <- 0.90 # the GAM.fit can fail if this is too low with errors like: Error in eigen(hess1, symmetric = TRUE) : 0 x 0 matrix
ievent.k             <- list(doy=-1, cc=-1)  # doy=6 or 12, cc=5 or 3
ievent.covaraite     <- 'gmst'
# if(FUTURE_PROJECTION) ievent.k$cc <- -1 else ievent.k$cc <- -1     # GAM smoother (k) for climate change effect on heatwave initiation probability
# st_InitEvent <- sub('.RData','_IEventThXXX.RData',basename(st_base))
# st_InitEvent <- paste(HDSAVEDIR, 'InitEvent/', sub('XXX',event.th.u,st_InitEvent),sep='/')

# initI                   <- list()
# initI$st_InitEvent      <- st_InitEvent
# initI$SAVEINITEVENTDIAG <- SAVEINITEVENTDIAG
# initI$DOINITEVENTPLOT   <- DOINITEVENTPLOT
# initI$ievent.th.u       <- ievent.th.u
# initI$ievent.k          <- ievent.k
# initI$ievent.covaraite  <- ievent.covaraite
# HDconfig$initI          <- initI
### doFitInitEvent.R #########################################################################################



# HS_InitEvent <- fn_fitInitEvent(events01, data01, initI, IsJoint=TRUE, SAVEINITEVENTDIAG=SAVEINITEVENTDIAG)

### obs
evs                 <- events01$obs
iobs                <- which(data01$isobs==1)
init2               <- rep(0, length(iobs)) # make an array for all data
init2[evs$ch.i[2,]] <- 1  # 2=first day of heatwave, absolute index for the original data
if(evs$hotseas$IsHotSeason) {
    form.doy  <- 'cs'
    init2     <- init2[evs$hotseas$ihotseas]
    evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[evs$hotseas$ihotseas], cc=data01[evs$hotseas$ihotseas,ievent.covaraite], class='obs' )
} else {
    form.doy  <- 'cc' # ensure smooth is cyclic 'cc' if whole year is being modelled
    evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[iobs], cc=data01[iobs,ievent.covaraite], class='obs' )
}
    plot(evinit.o$init,pch=3,cex=.3)
    plot(evinit.o$sdoy,pch=3,cex=.3)
    plot(evinit.o$cc,pch=3,cex=.3)
    plot(evinit.o$sdoy,evinit.o$init,pch=3,cex=.3)
    plot(evinit.o$cc,evinit.o$init,pch=3,cex=.3)
    i1<-which(evinit.o$init==1)
    hist(evinit.o$sdoy[i1],20)
    hist(evinit.o$cc[i1],20)
fit.init.D  <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=evinit.o, family="binomial")
fit.init.T  <- gam(init ~                                     s(cc,bs='cr',k=ievent.k$cc), data=evinit.o, family="binomial")
fit.init.DT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='cr',k=ievent.k$cc), data=evinit.o, family="binomial")
# predict.gam(fit.init.D, newdata=data.frame(sdoy=0.5),type="response")
# InitModel$obs <- list(D=fit.init.D, T=fit.init.T, DT=fit.init.DT, ch.th.u=evs$ch.th.u, covariate=ievent.covaraite)

### 2025.10.31
    # inspecting all these NWS initiation has no dependence on doy or gmst

    # So is it a simple logistic with fixed probability of initiation each day?

    # DO WE NEED TO EXTEND THE SUMMER HALF YEAR TO CAPTURE MORE EVENTS?