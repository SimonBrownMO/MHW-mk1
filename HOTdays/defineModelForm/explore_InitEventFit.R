# source("explore_InitEventFit.R")

### Use this file to determine the form of Initialising an event model
st.pwd <- system("pwd", intern=TRUE)
# source(paste(st.pwd,"/../../setup_all.R",sep=''))
source("setup_dev.R")


    MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/MSref/ostia_cdr_nrt_regions.MSref.2026-03-26.RData"
    load(MSref.file, verb=TRUE)

    list2env(MSconfig ,       envir = .GlobalEnv)
    list2env(MSconfig$files , envir = .GlobalEnv)

    load(st_msdata01, verb=TRUE)

    source(paste(st.pwd,"/../../make_stationary/setup_MakeStationary.R",sep=''))

    st_events  <- list.files(glue(dirname(MSref.file),"/../Events"), pattern='_EventsTh', full.names=TRUE)
    st_events  <- st_events[grep('.RData', st_events)]
    load(st_events, verb=TRUE)
    list2env(events01$info , envir = .GlobalEnv)

    source(paste(st.pwd,"/../setup_HotDays.R",sep=''))

readline("Stop1")

DOHOTSEASON <- FALSE
pdf(paste('results/explore_InitEventFit.ally.pdf',sep='/'), width=10, height=10)


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
covariate <- initI$ievent.covaraite
ievent.k  <- initI$ievent.k
InitModel <- list()

# ### obs
    #     iobs   <- which(data01$isobs==1)
    #     init2  <- rep(0, length(iobs)) # make an array for all data
    #     if(events01$hotseas$IsHotSeason) {
    #         evs                 <- events01$obs
    #         init2[evs$ch.i[2,]] <- 1  # 2=first day of heatwave, absolute index for the original data
    #         form.doy            <- 'cs'
    #         init2               <- init2[evs$hotseas$ihotseas]
    #         evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[evs$hotseas$ihotseas], cc=data.frame(data01)[evs$hotseas$ihotseas,covariate], class='obs' )
    #     } else {
    #         evs                 <- events01$events01_ally$obs
    #         init2[evs$ch.i[2,]] <- 1  # 2=first day of heatwave, absolute index for the original data
    #         form.doy  <- 'cc' # ensure smooth is cyclic 'cc' if whole year is being modelled
    #         evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[iobs], cc=data.frame(data01)[iobs,covariate], class='obs' )
    #     }

    #     up.2by2()
    #         plot(evinit.o$init,pch=3,cex=.3)
    #         plot(evinit.o$sdoy,pch=3,cex=.3)
    #         plot(evinit.o$cc,pch=3,cex=.3)
    #         plot(evinit.o$sdoy,evinit.o$init,pch=3,cex=.3)

    #     up.2by2()
    #         plot(evinit.o$cc,evinit.o$init,pch=3,cex=.3)
    #         i1<-which(evinit.o$init==1)
    #         hist(evinit.o$sdoy[i1],20)
    #         hist(evinit.o$cc[i1],20)
    #     fit.init.D  <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=evinit.o, family="binomial")
    #     fit.init.T  <- gam(init ~                                     s(cc,bs='ts',k=ievent.k$cc), data=evinit.o, family="binomial")
    #     fit.init.DT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='ts',k=ievent.k$cc), data=evinit.o, family="binomial")

    #     up.1()
    #     plot(fit.init.D,page=1)
    #     plot(fit.init.T,page=1)
    #     plot(fit.init.DT,page=1)

    #     # 2026.04.01 - For whole year DT seems the way to go
    #     # summary(fit.init.DT)
    #     # Formula:
    #     # init ~ s(sdoy, bs = form.doy, k = ievent.k$doy) + s(cc, bs = "ts", 
    #     #     k = ievent.k$cc)

    #     # Parametric coefficients:
    #     #             Estimate Std. Error z value Pr(>|z|)    
    #     # (Intercept) -4.19731    0.06417  -65.41   <2e-16 ***

    #     # Approximate significance of smooth terms:
    #     #           edf Ref.df Chi.sq p-value  
    #     # s(sdoy) 2.078      8  6.272  0.0271 *
    #     # s(cc)   0.408      9  0.684  0.1948  

    #     # R-sq.(adj) =  0.000402   Deviance explained = 0.349%
    #     # UBRE = -0.84478  Scale est. = 1         n = 16723
###

### Joint
iall   <- seq_along(data01$isobs) # use all data for joint model
iobs   <- which(data01$isobs==1)
imod   <- which(data01$isobs!=1)
# imod1  <- tail(iobs,1)+1
init2  <- c(rep(0, length(iobs)), rep(0, length(imod))) # make an array for all data

if(DOHOTSEASON) {
    evs.o               <- events01$obs
    evs.m               <- events01$mod
    form.doy            <- 'cs'
    id01                <- c(evs.o$hotseas$ihotseas, evs.m$hotseas$ihotseas)
} else {
    evs.o               <- events01$events01_ally$obs
    evs.m               <- events01$events01_ally$mod
    form.doy  <- 'cc' # ensure smooth is cyclic 'cc' if whole year is being modelled
    id01                <- seq_along(data01$isobs)
}
init2[c(evs.o$ch.i[2,], evs.m$ch.i[2,])] <- 1  # 2=first day of heatwave, absolute index for the original data

evinit.j  <- data.frame(init=init2[id01], sdoy=data01$sdoy[id01], cc=data.frame(data01)[id01,covariate], ftype=data01$ftype[id01] )

up.2by2()
    plot(evinit.j$init,pch=3,cex=.3)
    plot(evinit.j$sdoy,pch=3,cex=.3)
    plot(evinit.j$cc,pch=3,cex=.3)
    plot(evinit.j$sdoy,evinit.j$init,pch=3,cex=.3)

up.2by2()
    plot(evinit.j$cc,evinit.j$init,pch=3,cex=.3)
    i1<-which(evinit.j$init==1)
    hist(evinit.j$sdoy[i1],20)
    hist(evinit.j$cc[i1],20)
fit.init.D  <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=evinit.j, family="binomial")
fit.init.T  <- gam(init ~                                     s(cc,bs='ts',k=ievent.k$cc), data=evinit.j, family="binomial")
fit.init.DT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='ts',k=ievent.k$cc), data=evinit.j, family="binomial")

plot(fit.init.D,page=1)
plot(fit.init.T,page=1)
plot(fit.init.DT,page=1)

fit.init.IDT <- gam(init ~ ftype + s(sdoy,bs=form.doy,k=ievent.k$doy,by=ftype) +s(cc,bs='ts',k=ievent.k$cc), data=evinit.j, family="binomial")
plot(fit.init.IDT,page=1)
summary(fit.init.IDT)
    # 2026.04.02 Run with this I think
        # summary(fit.init.IDT)
        # Family: binomial 
        # Link function: logit 

        # Formula:
        # init ~ ftype + s(sdoy, bs = form.doy, k = ievent.k$doy, by = ftype) + 
        #     s(cc, bs = "ts", k = ievent.k$cc)

        # Parametric coefficients:
        #             Estimate Std. Error z value Pr(>|z|)    
        # (Intercept) -4.19625    0.06417  -65.39   <2e-16 ***
        # ftypemod     0.05550    0.07706    0.72    0.471    
        # ---
        # Approximate significance of smooth terms:
        #                      edf Ref.df Chi.sq p-value   
        # s(sdoy):ftypeobs 2.05232      8  6.193 0.02745 * 
        # s(sdoy):ftypemod 3.95411      8 15.505 0.00154 **
        # s(cc)            0.01163      9  0.002 0.67051   
        # ---
        # R-sq.(adj) =  0.000395   Deviance explained = 0.316%
        # UBRE = -0.83991  Scale est. = 1         n = 52693

fit.init.IbDT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy,by=ftype) +s(cc,bs='ts',k=ievent.k$cc), data=evinit.j, family="binomial")
plot(fit.init.IbDT,page=1)
summary(fit.init.IbDT)

dev.off()


###############################################################################
### try using the HOTday fn
cat(cr,cr,"##############################",cr,"Using fn_fitInitEvent",cr)
list2env(HDconfig$initI , envir = .GlobalEnv)

fmla1 <- formula(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy,by=ftype) +s(cc,bs='ts',k=ievent.k$cc))

HS_InitEvent.ally <- fn_fitInitEvent(fmla1, events01$events01_ally, data01, HDconfig$initI, DOHOTSEASON=FALSE, IsJoint=TRUE, SAVEINITEVENTDIAG=TRUE)
HS_InitEvent.seas <- fn_fitInitEvent(fmla1, events01,               data01, HDconfig$initI, DOHOTSEASON=TRUE,  IsJoint=TRUE, SAVEINITEVENTDIAG=TRUE)

fn_plotInitEvent(HS_InitEvent.ally$Model, main0='Init ~', st.pdf='test.initevent.ally.pdf')
fn_plotInitEvent(HS_InitEvent.seas$Model, main0='Init ~', st.pdf='test.initevent.hseas.pdf')