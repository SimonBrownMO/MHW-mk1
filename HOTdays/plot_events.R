# source("plot_events.R")

### libraries
        # source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
        # source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
        # # library(Rcpp, lib="/home/users/simon.brown/extremes/R/packages")
        # # library(evgam, lib="/home/users/simon.brown/extremes/R/packages")
        # library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
        # library(mgcv)
        # library(qgam)
        # # library(ncdf4)
        # # library(PCICt)
        # library(evd)
        # library(float)
        # source("../libs/fn_MakeStationary.R")
        # source("../libs/fn_HotDays.R")
        # # source("../libs/fn_JointSimHD.R")
        # # source("../libs/lib_HotDay.R")
source(paste(st.pwd,"/../setup_all.R",sep=''))

# # Hobday
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday"
# CHPKARELAPLACE <- FALSE

# Hobday with time
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ time"
# CHPKARELAPLACE <- FALSE

# Hobday with GMST
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ GMST"
# CHPKARELAPLACE <- FALSE

# v5
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/MSref/ostia_cdr_nrt_regions.MSref.2026-03-26.RData"
st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
tit1       <- "v5 UKV:"
CHPKARELAPLACE <- TRUE
SAVEPLOT       <- TRUE

load(st_events, verb=TRUE)
load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)
load(st_msdata01, verb=TRUE)

# plotdata1 <- data01$anomaly
plotdata1 <- data01$residStep1

st_pdf <- sub('.RData','.EventsPlot.pdf',st_events)

if(SAVEPLOT) pdf(file=st_pdf,height=7,width=10)

    # events <- events01$obs
    # ch.lev.l <- qlaplace(events$ch.th.u)

    # c21      <- distinct21colours()
    # isortage <- sort(events$ch.duration,dec=TRUE,index=TRUE)$ix
    # for( j in isortage[1:nplot]){
    #     up.1()
    #     y0  <- unique(unlist(events$ch.time[[j]]))
    #     i01 <- which(trunc(data$time)==trunc(y0)[1])
    #     ylim0 <- c(-0.5, max(qlaplace(data$u[i01]),na.rm=TRUE))
    #     plot(data$time[i01],qlaplace(data$u[i01]), pch=3, cex=.3, main=trunc(y0)[1], ylim=ylim0)
    #     y1  <- unlist(lapply(events$ch.time, "[", 1))
    #     i1  <- which(trunc(y1)==trunc(y0)[1])
    #     c2 <- sample(c21, length(i1), rep=FALSE)
    #     for(k in seq_along(i1)) {
    #         ch.t <- events$ch.time[[i1[k]]]
    #         points(ch.t, events$ch.st.l[,i1[k]][seq_along(ch.t)],col=c2[k])
    #          lines(ch.t, events$ch.st.l[,i1[k]][seq_along(ch.t)],col=c2[k])
    #     }
    #     abline(h=ch.lev.l, lty=2, col='grey70')
    #     if(DOPAUSE) if(!SAVEPLOT) readline("Continue?")
    # }
    # if(!is.null(savefile))  dev.off()


#

ch.lev.l <- qlaplace(events01$info$event.th.u)
c21      <- distinct21colours()
r.l      <- range(qlaplace(data01$u[finite(events01$obs$ch.i[])]),na.rm=T)
r.l2     <- c(ch.lev.l-.5, max(qlaplace(data01$u[finite(events01$obs$ch.i[])]),na.rm=T))
i0       <- which(data01$isobs==1)

###  wrt time  #####################################################
if(TRUE) {
    # all years
    # obs
    plot(data01$time[i0], qlaplace(data01$u[i0]), pch=3, cex=.3, main=paste(tit1,"All obs data with events overlayed") )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")
    # mod
    plot(data01$time[-i0], qlaplace(data01$u[-i0]), pch=3, cex=.3, main=paste(tit1,"All moddata with events overlayed") )    
    for(k in 1:length(events01$mod$ch.sev)) {
        ch.t <- events01$mod$ch.time[[k]]
        ch.l <- qlaplace(data01$u[-i0][finite(events01$mod$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")

    # all summers
    # obs
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o & data01$isobs==1)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main=paste(tit1,"Summer obs data with events overlayed") )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")

    # mod
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$m & data01$isobs==0)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main=paste(tit1,"Summer model data with events overlayed") )    
    for(k in 1:length(events01$mod$ch.sev)) {
        ch.t <- events01$mod$ch.time[[k]]
        ch.l <- qlaplace(data01$u[-i0][finite(events01$mod$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")


    #  2022-2024
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o & data01$time>=2022 & data01$time<2025 & data01$isobs==1)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main=paste(tit1,"2022-2024 obs with events overlayed") )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")


    # with time laplace and observed temperature
    up.2()
    # isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o & data01$isobs==1)
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o & data01$time>=2015 & data01$time<2025 & data01$isobs==1)
    r.x  <- range(data01$x[isum], na.rm=TRUE) # range(events01$obs$ch.st.o, na.rm=TRUE)
    plot(data01$time[isum], data01$x[isum], pch=3, cex=.3, main=paste(tit1,"Obs events by date"), ylim=r.x ,xlab="Year", ylab="Temperature (C))")    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.o <- data01$x[finite(events01$obs$ch.i[,k])]
        points(ch.t, ch.o,col=c21[k %% 21 +1])
         lines(ch.t, ch.o,col=c21[k %% 21 +1],lwd=3)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3,    main=NULL,ylim=r.l, xlab="year", ylab="Temperature (Laplace scale)" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l,col=c21[k %% 21 +1])
         lines(ch.t, ch.l, col=c21[k %% 21 +1],lwd=3)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()

}
#

####  wrt doy  #####################################################
if(TRUE) {
    # all years
    # obs
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$o)-10, to=max(events01$hotseas$hotSdoy$o)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l2, ty='n', main=paste(tit1,"All obs data with events overlayed"),xlab="Day of year", ylab="Temperature (Laplace scale)" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.d <-          data01$doy[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    if(!SAVEPLOT) readline("Continue?")
    # mod
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$m)-10, to=max(events01$hotseas$hotSdoy$m)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l2, ty='n', main=NULL,xlab="Day of year", ylab="Temperature (Laplace scale)" )    
    for(k in 1:length(events01$mod$ch.sev)) {
        ch.d <-        data01$doy[-i0][finite(events01$mod$ch.i[,k])]
        ch.l <- qlaplace(data01$u[-i0][finite(events01$mod$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    if(!SAVEPLOT) readline("Continue?")
}
#

####  top 10 plots  #####################################################
up.1()
if(TRUE) {

    # top 10 wrt peak temperature
    i2 <- sort(events01_ally$obs$ch.pkval, dec=TRUE, index=TRUE)$ix
    xdoy <- 1:365
    r.l3 <- c(1, max(events01_ally$obs$ch.pkval))
    plot(xdoy, xdoy, ylim=r.l3, ty='n', main=paste(tit1,"top 20 obs by peak temperature") ,xlab="Day of year", ylab="Temperature (Laplace scale)")    
    for(k in i2[1:20]) {
        ch.d <-        data01$doy[finite(events01_ally$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01_ally$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE),
            labels=paste0(trunc(events01_ally$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")


    # top 10 wrt severity
    i2 <- sort(events01_ally$obs$ch.sev, dec=TRUE, index=TRUE)$ix
    xdoy <- 1:365
    r.l3 <- c(1, max(events01_ally$obs$ch.pkval))
    plot(xdoy, xdoy, ylim=r.l3, ty='n', main=paste(tit1,"top 20 obs by Severity"),xlab="Day of year", ylab="Temperature (Laplace scale)")    
    for(k in i2[1:20]) {
        ch.d <-        data01$doy[finite(events01_ally$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01_ally$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE),
            labels=paste0(trunc(events01_ally$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")

    # top 10 wrt duration
    i2 <- sort(events01_ally$obs$ch.dur, dec=TRUE, index=TRUE)$ix
    xdoy <- 1:365
    plot(xdoy, xdoy, ylim=r.l3, ty='n', main=paste(tit1,"top 20 obs by Duration"),xlab="Day of year", ylab="Temperature (Laplace scale)")    
    for(k in i2[1:20]) {
        ch.d <-        data01$doy[finite(events01_ally$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01_ally$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE),
            labels=paste0(trunc(events01_ally$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='red4', lwd=2)     
    grid()
    if(!SAVEPLOT) readline("Continue?")

}
#

####  with respect to peak day  ####################################################
if(TRUE) {

    # anomalies
    e.dur <- events01$obs$ch.duration
    p.dur <- -30:30
    r.a   <- range(p.dur)
    r.x   <- c(-0.1,max(events01$obs$ch.pk$chains, na.rm=TRUE))
    r.x   <- c(-0.1,max(plotdata1, na.rm=TRUE))
    c2    <- sample(c21, size=length(events01$obs$ch.sev), replace=TRUE)
    plot(p.dur, p.dur, ty='n', xlim=r.a, ylim=r.x, main=paste(tit1,"Obs events centered on peak"), xlab="Days relative to peak", ylab="Anomaly (C)" )   
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.i <- which(events01$obs$ch.pk$chains[,k]>0)
        if(length(ch.i)>2){
        ch.a <- ch.i - which.max(events01$obs$ch.pk$chains[,k])
        ch.x <- events01$obs$ch.pk$chains[ch.i,k]  #  qlaplace(data01$x[finite(events01$obs$ch.i[,k])])
        ch.x <- plotdata1[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        # points(ch.a, ch.x, col=c21[k %% 21 +1])
        lines(ch.a, ch.x, col=c21[k %% 21 +1])
        }
    }
    # abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    if(!SAVEPLOT) readline("Continue?")

    headtail <- function(x,b=1,e=1) {
    if(length(x) <= 2) return(c())
    x[(b+1):(length(x)-e)]
    }

    # laplace
    r.l2  <- c(ch.lev.l-.5, max( c( qlaplace(data01$u[finite(events01$obs$ch.i[])]),
                                    qlaplace(data01$u[finite(events01$mod$ch.i[])])),na.rm=T))
    p.dur <- -30:30
    r.a   <- range(p.dur)

    # obs
    # e.dur <- events01$obs$ch.duration
    c2    <- sample(c21, size=length(events01$obs$ch.sev), replace=TRUE)
    plot(p.dur, p.dur, ty='n', xlim=r.a, ylim=r.l2, main=paste(tit1,"Obs events centered on peak"), xlab="Days relative to peak", ylab="Laplace" )   
    for(k in 1:length(events01$obs$ch.sev)) {
        # ch.t <- events01$obs$ch.time[[k]]
        if(CHPKARELAPLACE) ch.i <- which(events01$obs$ch.pk$chains[,k]>ch.lev.l) else ch.i <- which(events01$obs$ch.pk$chains[,k]>0)
        ch.a <- ch.i - which.max(events01$obs$ch.pk$chains[,k])
        # ch.x <- events01$obs$ch.pk$chains[ch.i,k]  #  qlaplace(data01$x[finite(events01$obs$ch.i[,k])])
        ch.l <- headtail(qlaplace(data01$u[finite(events01$obs$ch.i[,k])]))
        # points(ch.a, ch.x, col=c21[k %% 21 +1])
        lines(ch.a, ch.l, col=c21[k %% 21 +1])
    }
    # abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    if(!SAVEPLOT) readline("Continue?")

    # model
    # e.dur <- events01$mod$ch.duration
    c2    <- sample(c21, size=length(events01$mod$ch.sev), replace=TRUE)
    plot(p.dur, p.dur, ty='n', xlim=r.a, ylim=r.l2, main=paste(tit1,"Mod events centered on peak"), xlab="Days relative to peak", ylab="Laplace" )   
    for(k in 1:length(events01$mod$ch.sev)) {
        # ch.t <- events01$mod$ch.time[[k]]
        if(CHPKARELAPLACE) ch.i <- which(events01$mod$ch.pk$chains.l[,k]>ch.lev.l) else ch.i <- which(events01$mod$ch.pk$chains[,k]>0)
        ch.a <- ch.i - which.max(events01$mod$ch.pk$chains.l[,k])
        # ch.x <- events01$mod$ch.pk$chains.l[ch.i,k]  #  qlaplace(data01$x[finite(events01$mod$ch.i[,k])])
        ch.l <- headtail(qlaplace(data01$u[-i0][finite(events01$mod$ch.i[,k])]))        
        # points(ch.a, ch.x, col=c21[k %% 21 +1])
        lines(ch.a, ch.l, col=c21[k %% 21 +1])
    }
    # abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    if(!SAVEPLOT) readline("Continue?")

}

if(SAVEPLOT) dev.off()
#
