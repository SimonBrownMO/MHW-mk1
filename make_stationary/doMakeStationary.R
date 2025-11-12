# source("doMakeStationary.R")

### the Hobday approach uses a stationary  defined climatology thus need to select the appropriate year
### for the centre of this climatology to perform the HOTday transforms for qgam, q2p, EVgam and p2q 

sysdate   <- Sys.time()
cat("###############################################################################",cr)
stsysdate <- print(sysdate)
cat("Start doMakeStationary",cr)
cat("###############################################################################",cr,cr)

if(st_version != "v_Hobday" ) {

    ### Calculate annualcycle and climate change signal with QGAM
    # fit QGAM

    if(MSconfig$DOQGAM) {

        cat(cr,"Fitting QGAM",cr)

        if(length(data01$time)>MSconfig$nmax_qgam) {
            set.seed(0)
            i0 <- sort(sample(seq_along(data01$time), MSconfig$nmax_qgam, replace=FALSE))
        } else i0 <- seq_along(data01$time)


        fn_fit_qgam(MSconfig$fmla.MSqgam, data01[i0,], st.qgam=MSconfig$files$st_qgam, ptiles=MSconfig$do.ptiles)

        load(MSconfig$files$st_qgam, verb=TRUE) # if doing q2p+p2q only

        cat(cr)
        cat(cr,"Generating q2p",cr)
        fn_gen_q2p(data01, st.q2p=MSconfig$files$st_q2p, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam)
        cat(cr,"Generating p2q",cr)
        load(MSconfig$files$st_q2p)
        fn_gen_p2q(data01, qgam.q2p.fn, ms.thgpd.u, st.p2q=MSconfig$files$st_p2q, st.q2p=MSconfig$files$st_q2p)

        gc()

    } # if doing qgam+q2p+p2q
    # readline("stop2")

    cat(cr,"doMakeStationary: (re)loading pre-calculated QGAM",cr,MSconfig$files$st_qgam,cr)
    load(MSconfig$files$st_qgam, verb=TRUE)
    load(MSconfig$files$st_q2p,  verb=TRUE)

    cat(cr,"Applying q2p",cr)
    data01$uqgam <- fn_apply_q2p(data01, qgam.q2p.fn, st.q2p=NULL, SANITY=FALSE)

    # calc MSgpd
    cat(cr,"Fitting EVGAM",cr)
    fn_fit_MSgpd(data01, ms.thgpd.u, fmla.MSgpd, chosen.MSgpd.name=chosen.MSgpd.name, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam, st.msgpd=MSconfig$files$st_msgpd)
    load(MSconfig$files$st_msgpd, verb=TRUE) # if doing only apply MSgpd

    # apply MSgpd
    cat(cr,"Applying EVGAM",cr)
    u1 <- fn_apply_MSgpd(data01, ms.thgpd.u, fit.qgam=NULL, st.qgam=st_qgam, st.msgpd=st_msgpd)

    data01$u     <- u1

    ### general testing plots
        stin   <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('MSqgam.RData','general.pdf',basename(stin)),sep='/')
        fn_diag_general(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)


    # st.u0 <- paste(MSconfig$files$MSSAVEDIR,sub('.RData','_u0.RData',basename(MSconfig$files$st_base)),sep='' )
    save(file=st_msdata01, data01, data01.std.param)
    cat(cr,"Saved msdata",cr,st_msdata01,cr)


    if(DODIAGQGAM) {
    cat(cr,"DODIAGQGAM",cr)
        stin <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_qgam(fit.qgam, savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGMSGPD) {
    cat(cr,"DODIAGMSGPD",cr)
        stin <- MSconfig$files$st_msgpd
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_evgam(savefile=stdiag, st.msgpd=stin, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGQ2P) {
    cat(cr,"DODIAGQ2P",cr)
        stin <- MSconfig$files$st_q2p
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_q2p(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    tidy()

} else {  # "v_Hobday"==TRUE # try and mimic Hobday

    ### Calculate annualcycle and climate change signal with QGAM
    # fit QGAM

    if(MSconfig$DOQGAM) {

        cat(cr,"Fitting QGAM",cr)

        if(length(data01$time)>MSconfig$nmax_qgam) {
            set.seed(0)
            i0 <- sort(sample(seq_along(data01$time), MSconfig$nmax_qgam, replace=FALSE))
        } else i0 <- seq_along(data01$time)


        fn_fit_qgam(MSconfig$fmla.MSqgam, data01[i0,], st.qgam=MSconfig$files$st_qgam, ptiles=MSconfig$do.ptiles)

        load(MSconfig$files$st_qgam, verb=TRUE) # if doing q2p+p2q only

        cat(cr)
        cat(cr,"Generating q2p",cr)
        fn_gen_q2p(data01, st.q2p=MSconfig$files$st_q2p, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam)
        cat(cr,"Generating p2q",cr)
        load(MSconfig$files$st_q2p)
        fn_gen_p2q(data01, qgam.q2p.fn, ms.thgpd.u, st.p2q=MSconfig$files$st_p2q, st.q2p=MSconfig$files$st_q2p)

        gc()

    } # if doing qgam+q2p+p2q
    # readline("stop2")

    cat(cr,"doMakeStationary: (re)loading pre-calculated QGAM",cr,MSconfig$files$st_qgam,cr)
    load(MSconfig$files$st_qgam, verb=TRUE)
    load(MSconfig$files$st_q2p,  verb=TRUE)

    cat(cr,"Applying q2p",cr)
    iclim.y <- which(data01$time > clim.year & data01$time <= (clim.year+1) )
    qgam.q2p.fn.clim <- qgam.q2p.fn[iclim.y]
    xx <- rep(NA, nrow(data01))
    for(y1 in unique(trunc(data01$time))) {
        ix   <- which(data01$time > y1 & data01$time <= (y1+1) )
         sink("/dev/null"); idoy <- findall.nearest(data01$doy[ix], data01$doy[iclim.y]); sink()
        for(j in seq_along(idoy)) {
            xx[ix[j]] <- qgam.q2p.fn.clim[[idoy[j]]](data01$x[ix[j]])
        }
    }
    data01$uqgam <- xx
    
    # calc MSgpd
    cat(cr,"Fitting EVGAM",cr)
    fn_fit_MSgpd(data01, ms.thgpd.u, fmla.MSgpd, chosen.MSgpd.name=chosen.MSgpd.name, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam, st.msgpd=MSconfig$files$st_msgpd)
    th.01.gpd <- qdo(fit.qgam, ms.thgpd.u, predict, newdata=data01 )
    attr(th.01.gpd, "thresh.u") <- ms.thgpd.u

    # Hobday threshold calculation
    # Hobday et al. (2016) define the threshold for a hot day as the 90th percentile
    # of daily mean temperatures over a climatological period (here 1982-2011),
    data_clim       <- data01
    time_clim       <- clim.year + (data01$time - trunc(data01$time)) # middle of 1982:2111, add on the fraction of year
    stime_clim      <- (time_clim - mean(time_clim))/ ( max(trunc(data01$time)) - min(trunc(data01$time)) )
    data_clim$stime <- stime_clim

    th.01.gpd <- qdo(fit.qgam, newdata=data_clim, type='response', qu=ms.thgpd.u, predict)            
    attr(th.01.gpd, "thresh.u") <- ms.thgpd.u

    data01$excess     <- data01$x - th.01.gpd
    data01$th.01.gpd  <- th.01.gpd
    data01$gpd.ktime  <- 1 # evgam bug fix, values dont mater
    data.gpd          <- subset(data01, excess > 0)
    data.gpd$ms.k     <- 0

    ### do fit
        sink("/dev/null")
        chosen.MSgpd      <- evgam(fmla.MSgpd[[chosen.MSgpd.name]], data.gpd, family="gpd", trace=2)
        sink()

        cat('Chosen MS GPD:',chosen.MSgpd.name,cr)
        print(summary(chosen.MSgpd))

        chosen.MSgpd$ms.thgpd.u <- ms.thgpd.u
        date01 <- data01[,c('time','doy')]
        save(file=MSconfig$files$st_msgpd, chosen.MSgpd, chosen.MSgpd.name, th.01.gpd, date01)

    # apply MSgpd
    load(MSconfig$files$st_msgpd, verb=TRUE) # if doing only apply MSgpd
    cat(cr,"Applying EVGAM",cr)
    u1 <- fn_apply_MSgpd(data01, ms.thgpd.u, fit.qgam=NULL, st.qgam=st_qgam, st.msgpd=st_msgpd)

    data01$excess    <- data01$x - th.01.gpd
    data01$th.01.gpd <- th.01.gpd

    gpd01par <- predict(chosen.MSgpd, newdata=data_clim, type="response")
    u01      <- data01$uqgam
    excess   <- data01$excess
    dt1      <- trunc(length(excess)/20)
    for(i in seq_along(excess)) {
        # if(excess[i]>0 & u01>ms.thgpd.u) {
        if(excess[i]>0) {
            sh          <- gpd01par$shape[i]
            sc          <- gpd01par$scale[i]
            u01[i] <- 1-(1-ms.thgpd.u)*pmax(0,(1+(sh*((excess[i])/sc))))^(-1/sh)
            # cat(i,sc,sh,u01[i],cr)
        } else if(is.finite(data01$uqgam[i])) {
            # assume th.01.gpd is truth over q2p
            # needed as sometimes the interpolation of q2p gives uqgam>u.ev even when excess<=0
            if(excess[i]<=0 & data01$uqgam[i]>ms.thgpd.u) {  
                u01[i] <- ms.thgpd.u
            }
        }
    }
    data01$u     <- u01

    # up.1()
    # plot(qlaplace(data01$uqgam), qlaplace(u01),pch=20,cex=.3)
    # i1<- which(excess<=0 & data01$uqgam>ms.thgpd.u)
    # points(qlaplace(data01$uqgam[i1]), qlaplace(u01[i1]),pch=20,cex=.6, col=3)
    # abline(h=qlaplace(0.9),col=2)
    # abline(v=qlaplace(0.9),col=2)


    ### general testing plots
        stin   <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('MSqgam.RData','general.pdf',basename(stin)),sep='/')
        fn_diag_general(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)


    # st.u0 <- paste(MSconfig$files$MSSAVEDIR,sub('.RData','_u0.RData',basename(MSconfig$files$st_base)),sep='' )
    save(file=st_msdata01, data01, data01.std.param)
    cat(cr,"Saved msdata",cr,st_msdata01,cr)


    if(DODIAGQGAM) {
    cat(cr,"DODIAGQGAM",cr)
        stin <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_qgam(fit.qgam, savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGMSGPD) {
    cat(cr,"DODIAGMSGPD",cr)
        stin <- MSconfig$files$st_msgpd
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_evgam(savefile=stdiag, st.msgpd=stin, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGQ2P) {
    cat(cr,"DODIAGQ2P",cr)
        stin <- MSconfig$files$st_q2p
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_q2p(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    tidy()

}


cat(cr,cr)
cat("###############################################################################",cr)
sysdate   <- Sys.time()
stsysdate <- print(sysdate)
cat("All Done doMakeStationary",cr)
cat("###############################################################################",cr,cr)








#
