# source("plot_events.R")

### libraries
source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
# library(Rcpp, lib="/home/users/simon.brown/extremes/R/packages")
# library(evgam, lib="/home/users/simon.brown/extremes/R/packages")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
# library(ncdf4)
# library(PCICt)
library(evd)
library(float)
source("../libs/fn_MakeStationary.R")
source("../libs/fn_HotDays.R")
# source("../libs/fn_JointSimHD.R")
# source("../libs/lib_HotDay.R")

# Hobday
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
tit1       <- "Hobday"
CHPKARELAPLACE <- FALSE

# Hobday with time
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ time"
CHPKARELAPLACE <- FALSE

# Hobday with GMST
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ GMST"
CHPKARELAPLACE <- FALSE



load(st_events, verb=TRUE)

load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_msdata01, verb=TRUE)

st_pdf <- sub('.RData','.EventsPlot.pdf',st_events)
# pdf(file=st_pdf,height=7,width=7)

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
    #     if(DOPAUSE) readline("Continue?")
    # }
    # if(!is.null(savefile))  dev.off()


#

ch.lev.l <- qlaplace(events01$info$event.th.u)
c21      <- distinct21colours()
    r.l  <- range(qlaplace(data01$u[finite(events01$obs$ch.i[])]),na.rm=T)
    r.l2  <- c(ch.lev.l-.5, 9)

###  wrt time  #####################################################
if(FALSE) {
    # all years
    plot(data01$time, qlaplace(data01$u), pch=3, cex=.3, main="All data with events overlayed" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

    # all summers
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main="All data with events overlayed" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")


    #  2022-2024
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o & data01$time>=2022 & data01$time<2025)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main="2022-2024 with events overlayed" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l, col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

    # just events with time
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o)
    plot(data01$time[isum], qlaplace(data01$u)[isum], ty='n', main="Events by date", ylim=r.l, xlab="Year", ylab="Temperature (Laplace scale)" )   
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l[seq_along(ch.t)],col=c21[k %% 21 +1])
        lines(ch.t, ch.l[seq_along(ch.t)],col=c21[k %% 21 +1])
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")


    # with time laplace and observed temperature
    up.2()
    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o)
    r.x  <- range(data01$x[isum], na.rm=TRUE) # range(events01$obs$ch.st.o, na.rm=TRUE)
    plot(data01$time[isum], data01$x[isum], pch=3, cex=.3, main="Events by date", ylim=r.x ,xlab="Year", ylab="Temperature (C))")    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.o <- data01$x[finite(events01$obs$ch.i[,k])]
        points(ch.t, ch.o,col=c21[k %% 21 +1])
        lines(ch.t, ch.o,col=c21[k %% 21 +1],lwd=3)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

    isum <- which(data01$doy %in% events01$hotseas$hotSdoy$o)
    # r1 <- range(qlaplace(data01$u)[isum], na.rm=TRUE) # range(events01$obs$ch.st.l, na.rm=TRUE)
    plot(data01$time[isum], qlaplace(data01$u)[isum], pch=3, cex=.3, main="Events by date", 
            ylim=r.l, xlab="year", ylab="Temperature (Laplace scale)" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        points(ch.t, ch.l,col=c21[k %% 21 +1])
        lines(ch.t, ch.l, col=c21[k %% 21 +1],lwd=3)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()

}
#

####  wrt doy  #####################################################
if(FALSE) {
    # all years
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$o)-10, to=max(events01$hotseas$hotSdoy$o)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l2, pch=3, cex=.3, main="All data with events overlayed",xlab="Day of year", ylab="Temperature (Laplace scale)" )    
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.d <-          data01$doy[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    readline("Continue?")
}
#

####  top 10 plots  #####################################################
if(FALSE) {

    # top 10 wrt peak temperature
    i0 <- sort(events01$obs$ch.pkval, dec=TRUE, index=TRUE)$ix
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$o)-10, to=max(events01$hotseas$hotSdoy$o)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l, pch=3, cex=.3, main="top 10 peak temperature" ,xlab="Day of year", ylab="Temperature (Laplace scale)")    
    for(k in i0[1:10]) {
        ch.d <-        data01$doy[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE)+rnorm(1,0.1,.1),
            labels=paste0(trunc(events01$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")


    # top 10 wrt severity
    i0 <- sort(events01$obs$ch.sev, dec=TRUE, index=TRUE)$ix
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$o)-10, to=max(events01$hotseas$hotSdoy$o)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l, pch=3, cex=.3, main="top 10 by Severity" ,xlab="Day of year", ylab="Temperature (Laplace scale)")    
    for(k in i0[1:10]) {
        ch.d <-        data01$doy[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE)+rnorm(1,0.1,.1),
            labels=paste0(trunc(events01$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

    # top 10 wrt duration
    i0 <- sort(events01$obs$ch.dur, dec=TRUE, index=TRUE)$ix
    xdoy <- seq(from=min(events01$hotseas$hotSdoy$o)-10, to=max(events01$hotseas$hotSdoy$o)+10, by=1)
    plot(xdoy, xdoy, ylim=r.l, pch=3, cex=.3, main="top 10 by Duration" ,xlab="Day of year", ylab="Temperature (Laplace scale)")    

    for(k in i0[1:10]) {
        ch.d <-        data01$doy[finite(events01$obs$ch.i[,k])]
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        lines(ch.d, ch.l,col=c21[k %% 21 +1], lwd=2)
        text(ch.d[which.max(ch.l)],  max(ch.l, na.rm=TRUE)+rnorm(1,0.1,.1),
            labels=paste0(trunc(events01$obs$ch.time[[k]][1])), pos=3, cex=.8)
    }
    abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

}
#

####  with respect to peak day  ####################################################
if(TRUE) {

    # anomalies
    e.dur <- events01$obs$ch.duration
    p.dur <- -40:40
    r.a   <- range(p.dur)
    r.x   <- c(-0.1,max(data01$anomaly, na.rm=TRUE))
    c2    <- sample(c21, size=length(events01$obs$ch.sev), replace=TRUE)
    plot(p.dur, p.dur, ty='n', xlim=r.a, ylim=r.x, main="Events centered on peak", xlab="Days relative to peak", ylab="Anomaly (C)" )   
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        ch.i <- which(events01$obs$ch.pk$chains[,k]>0)
        ch.a <- ch.i - which.max(events01$obs$ch.pk$chains[,k])
        ch.x <- events01$obs$ch.pk$chains[ch.i,k]  #  qlaplace(data01$x[finite(events01$obs$ch.i[,k])])
        ch.l <- qlaplace(data01$u[finite(events01$obs$ch.i[,k])])
        # points(ch.a, ch.x, col=c21[k %% 21 +1])
        lines(ch.a, ch.x, col=c21[k %% 21 +1])
    }
    # abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")

headtail <- function(x,b=1,e=1) {
  if(length(x) <= 2) return(c())
  x[(b+1):(length(x)-e)]
}
    # laplace
    e.dur <- events01$obs$ch.duration
    p.dur <- -40:40
    r.a   <- range(p.dur)
    c2    <- sample(c21, size=length(events01$obs$ch.sev), replace=TRUE)
    plot(p.dur, p.dur, ty='n', xlim=r.a, ylim=r.l, main="Events centered on peak", xlab="Days relative to peak", ylab="Laplace" )   
    for(k in 1:length(events01$obs$ch.sev)) {
        ch.t <- events01$obs$ch.time[[k]]
        if(CHPKARELAPLACE) ch.i <- which(events01$obs$ch.pk$chains[,k]>ch.lev.l) else ch.i <- which(events01$obs$ch.pk$chains[,k]>0)
        ch.a <- ch.i - which.max(events01$obs$ch.pk$chains[,k])
        ch.x <- events01$obs$ch.pk$chains[ch.i,k]  #  qlaplace(data01$x[finite(events01$obs$ch.i[,k])])
        ch.l <- headtail(qlaplace(data01$u[finite(events01$obs$ch.i[,k])]))
        # points(ch.a, ch.x, col=c21[k %% 21 +1])
        lines(ch.a, ch.l, col=c21[k %% 21 +1])
    }
    # abline(h=ch.lev.l, lty=2, col='grey70')     
    grid()
    readline("Continue?")






}








#
