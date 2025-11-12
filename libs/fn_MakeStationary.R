cr <- '\n'

##############################################################################################
### random utilities

# Returns the number of days in a given year (handles leap years)
days_in_year <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || (year %% 400 == 0)) {
    return(366)
  } else {
    return(365)
  }
}

##############################################################################################
### Make stationary functions
##############################################################################################
reload_MS_ref<-function(freference=NULL)
{
    load(file=freference, verb=T)
    assign("MSconfig", MSconfig, envir = .GlobalEnv)
    # lapply(list(MSconfig), list2env ,envir = .GlobalEnv)
    list2env(MSconfig,       envir = .GlobalEnv)
    list2env(MSconfig$files, envir = .GlobalEnv)

}


##############################################################################
fn_fit_qgam <- function(fmla.qgam, data01, st.qgam, ptiles=(1:99)/100) {

    fit.qgam        <- mqgam(fmla.qgam, data=data01, qu=ptiles)
    fit.qgam$ptiles <- ptiles
    fit.qgam$fmla   <- fmla.qgam

    cat('Chosen MS QGAM:',cr)
    print(qdo(fit.qgam, 0.5, summary))

    save(file=st.qgam, fit.qgam, ptiles)
    cat("fn_fit_qgam: QGAM saved to",st.qgam,cr)

}


##############################################################################
fn_gen_q2p <- function(data01, st.q2p, fit.qgam=NULL, st.qgam=NULL) {

    if (is.null(fit.qgam)) {
        if(!file.exists(st.qgam)) {
            cat("fn_gen_q2p: ERROR",cr,"Must either supply fit.qgam object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_gen_q2p: loading qgam",cr,st.qgam,cr)
            load(st.qgam)
        }
    }

    ptiles <- fit.qgam$ptiles

    ### fit approx functions to all quantiles for given doy and selected times
    q.doy.time   <- array(0.0, dim=c(length(ptiles), length(data01$x)))

    cat("fn_gen_q2p: Calculating quantiles for all data",cr)
    for(p in seq_along(fit.qgam$ptiles)) {
        # q.doy.time[p,]   <- qdo(fit.qgam, do.ptiles[p], predict, newdata=data.frame(doy=data01$doy, time=data01$time) )
        q.doy.time[p,]   <- as.vector(qdo(fit.qgam, fit.qgam$ptiles[p], predict, newdata=data01) )
        if(p%%5==0) cat("fn_gen_q2p: Done q.doy.time ptile",p,cr)
    }
    cat(cr)
    save(file=sub('_MSq2p','_MSth',st.q2p), q.doy.time)

    qgam.q2p.fn <- vector(mode='list',length=length(data01$sdoy))
    for(d1 in seq_along(data01$sdoy)) {
        qgam.q2p.fn[[d1]] <- approxfun(q.doy.time[,d1], ptiles)
        if(d1%%5000==0) cat("fn_gen_q2p: Done qgam.q2p.fn ",d1,cr)
    }

    date01 <- data01[,c('time','doy')]
    save(file=st.q2p, qgam.q2p.fn, q.doy.time, ptiles, date01)

}


##############################################################################
fn_apply_q2p <- function(data01, q2p, st.q2p=NULL, SANITY=FALSE) {

    if (!exists("q2p")) {
        if(!file.exists(st.q2p)) {
            cat("fn_apply_q2p: ERROR",cr,"Must either supply q2p object or path to file",cr)
            return(NULL)
        } else {
            load(st.q2p)
        }
    }

    # transform data
    u0 <- double(length(data01$x))
    for(i in seq_along(u0))  {
        u0[i] <- q2p[[i]](data01$x[i])
    }

    # sanity
    if(SANITY) {
        par(mfcol=c(2,1))
        yo <- pred.o$time
        plot(yo, qo.from.qgam[50,],ty='n',ylim=range(qo.from.qgam),main='QGAM')
        i1 <- c(1, 1:7*10, 80:90, 99)
        c1 <- rainbow(100)
        for(i in i1) lines(yo, qo.from.qgam[i,], col=c1[i])

        x1 <- seq(min(qo.from.qgam), max(qo.from.qgam), length=50)
        y1 <- seq(min(pred.o$time), max(pred.o$time), length=20)
        plot(x1, qgam.q2p.o.fn[[1]](x1), ty='n', xlab='input', ylab='Prob',main='transform fn')
        for(y in y1) lines(x1, qgam.q2p.o.fn[[y-y1[1]+1]](x1), col=y)
    }

    return(u0)
}


##############################################################################
fn_gen_p2q <- function(data01, qgam.q2p.fn, ms.thgpd.u, st.p2q, n.q1=105, st.q2p=NULL, p2q.years=NULL) {

    if (is.null(qgam.q2p.fn)) {
        if(!file.exists(st.q2p)) {
            cat("fn_gen_p2q: ERROR",cr,"Must either supply q2p object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_gen_p2q: loading q2p",cr,st.q2p,cr)
            load(st.q2p)
        }
    }

    x  <- data01$x

    # q1 <- c( seq(min(x), quantile(x, ms.thgpd.u-0.1), by=0.5),  seq(quantile(x, ms.thgpd.u-0.09), max(x),  by=0.1), max(x))
    # q1 <- c( seq(min(x), max(x),  length=n.q1) )

    if(is.null(p2q.years)) ido <- seq_along(data01$x) else ido <- which(floor(data01$time) %in% p2q.years)

    qgam.p2q.fn <- vector(mode='list',length=length(ido))
    for(d1 in ido) {
        q1   <- seq( min(q.doy.time[,d1])-0.1, max(q.doy.time[,d1])+0.1, length=n.q1) # 2025.11.12
        n2.q1 <- n.q1
        p1    <- qgam.q2p.fn[[d1]](q1)
        iNa   <- which(is.na(p1))
        # cat(d1, length(p1), length(iNa), cr)
        if( length(iNa==0)) qgam.p2q.fn[[d1]] <- approxfun(p1,q1) else
        if( length(p1[-iNa]) > 33 ) { # need at least every 3rd percentile
            qgam.p2q.fn[[d1]] <- approxfun(p1[-iNa],q1[-iNa])
        } else {
            n2.q1 <- n2.q1 + n.q1
            q1    <- c( seq(min(q.doy.time[,d1]), max(q.doy.time[,d1]),  length=n2.q1) )
            p1    <- qgam.q2p.fn[[d1]](q1)
            iNa   <- which(is.na(p1))
            while (length(p1[-iNa]) <= 33 & n2.q1<10000) {
                n2.q1 <- n2.q1 + n.q1
                q1    <- c( seq(min(q.doy.time[,d1]), max(q.doy.time[,d1]),  length=n2.q1) )
                p1    <- qgam.q2p.fn[[d1]](q1)
                iNa   <- which(is.na(p1))
                cat('fn_gen_p2q: Increasing n.q1 to',n2.q1,'for day',d1,cr)
            }
            qgam.p2q.fn[[d1]] <- approxfun(p1[-iNa],q1[-iNa])
        }
        if(d1%%5000==0) cat('fn_gen_p2q: Done',d1,cr)
    }
    cat('fn_gen_p2q: Done',d1,cr)

    date01 <- data01[,c('time','doy')]
    save(file=st.p2q, qgam.p2q.fn, date01)
    cat("Done save qgam.p2q.fn\n",st.p2q,cr)

}

##############################################################################
fn_uback2obs_qgam <- function( datain.u, qgam.p2q.fn, gpd_vals, gpd_thu, gpd_th, QUIET=FALSE) {

        if( length(gpd_th)!=dim(gpd_vals)[1] | length(gpd_th)!=length(qgam.p2q.fn) ) {
            cat("fn_uback2obs_qgam: qgam.p2q.fn, gpd_vals, gpd_th all need to be the same length",cr)
            return(NA)
        }

        datain.o <-  double(length(datain.u))
        # data above the th.u
        iex      <- which(datain.u > gpd_thu)
        if(length(gpd_th)==1) { # single sample
            gpd_vals <- as.double(gpd_vals) # this stops datain.o being turned into a list!
            # undo ecdf
            for(i in seq_along(datain.u)[-iex]) datain.o[i] <- qgam.p2q.fn(datain.u[i])
            # undo GPD
            for(i in iex) {
                p1bb         <- 1-(1- datain.u[i])/(1-gpd_thu)
                datain.o[i]  <- gpd_th + evd::qgpd(p1bb, 0, gpd_vals[1], gpd_vals[2])
            }
        } else { # a set of parameters for each ob
            # undo ecdf
            for(i in seq_along(datain.u)[-iex]) datain.o[i] <- qgam.p2q.fn[[i]](datain.u[i])
            # undo GPD
            for(i in iex) {
                p1bb         <- 1-(1- datain.u[i])/(1-gpd_thu)
                datain.o[i]  <- gpd_th[i] + evd::qgpd(p1bb, 0, gpd_vals[i,1], gpd_vals[i,2])
            }
        }

        return(datain.o)
}


##############################################################################
f_pevg <- function(evgam1){
    xx<-summary(evgam1)
    cat(cr)
    cat('  locscale isobs\t',myround(xx[[1]]$logscale$Estimate[2],3),'\t\tPr',myround(xx[[1]]$logscale$'Pr(>|t|)'[2],3),cr)
    cat('  locscale s(time)\t\t\tPr',myround(xx[[2]]$logscale$'Pr(>|t|)',3),cr)
    cat('     shape isobs\t',myround(xx[[1]]$shape$Estimate[2],3),   '\t\tPr',myround(xx[[1]]$shape$'Pr(>|t|)'[2],3),cr)
    cat('     shape s(time)\t\t\tPr',myround(xx[[2]]$shape$'Pr(>|t|)',3),cr)
    cat(cr)
}


##############################################################################
fn_fit_MSgpd <- function(data01, ms.thgpd.u, fmla.MSgpd, chosen.MSgpd.name=NULL, fit.qgam=NULL, st.qgam=NULL, st.msgpd=NULL) {

    if (is.null(fit.qgam)) {
        if(!file.exists(st.qgam)) {
            cat("fn_fit_MSgpd: ERROR",cr,"Must either supply fit.qgam object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_fit_MSgpd: loading qgam",cr,st.qgam,cr)
            load(st.qgam)
        }
    }

    th.01.gpd <- qdo(fit.qgam, ms.thgpd.u, predict, newdata=data01 )
    attr(th.01.gpd, "thresh.u") <- ms.thgpd.u

    # # create time covariate based on 50 percentile
    # d2        <- data01[,c('x','time','doy','isobs')]
    # d2$doy    <- 180
    # q50       <- qdo(fit.qgam, 0.5, predict, newdata=d2 )

    #    plot(data01$time,       data01$x, type='n', main='Obs and CPM 90th quantile wrt time')
    # points(data01$time[iob],         data01$x[iob], pch=46)
    # points(data01$time[-iob],       data01$x[-iob], pch=46, col=2)
    # points(data01$time[iob],   th.01.gpd[iob], pch=46, col=3)
    # points(data01$time[-iob], th.01.gpd[-iob], pch=46, col=4)

    data01$excess     <- data01$x - th.01.gpd
    data01$th.01.gpd  <- th.01.gpd
    data01$gpd.ktime  <- 1 # evgam bug fix, values dont mater
    # data01$cov.T      <- q50

    data.gpd       <- subset(data01, excess > 0)
    data.gpd$ms.k  <- 0

    ### do fit
    if(is.null(chosen.MSgpd.name)) {

        fits.MSgpd        <- vector(mode='list',length=length(fmla.MSgpd))
        names(fits.MSgpd) <- names(fmla.MSgpd)

        for(i in seq_along(fmla.MSgpd)) {
            cat("Fitting",i,cr,paste(fmla.MSgpd[[i]],collapse=';  '),cr)
            sink("/dev/null")
            fits.MSgpd[[i]] <- evgam(fmla.MSgpd[[i]], data.gpd, family="gpd", trace=2)
            sink()
        }

        fits.MSgpd.aic        <- lapply(fits.MSgpd, AIC)
        names(fits.MSgpd.aic) <- names(fmla.MSgpd)
        cat("AIC",cr)
        for(i in seq_along(fmla.MSgpd)) cat(names(fmla.MSgpd)[i],fits.MSgpd.aic[[i]],cr)

        cat("Best MS GPD fit",cr)
        ix1 <- which.min(fits.MSgpd.aic)
        cat(names(fmla.MSgpd)[ix1],fits.MSgpd.aic[[ix1]],cr)
        chosen.MSgpd.name <- names(fmla.MSgpd)[ix1]
        chosen.MSgpd      <- fits.MSgpd[[ix1]]

        cat('Chosen MS GPD:',chosen.MSgpd.name,cr)
        print(summary(chosen.MSgpd))

        chosen.MSgpd$ms.thgpd.u <- ms.thgpd.u
        date01 <- data01[,c('time','doy')]
        save(file=st.msgpd, fits.MSgpd, chosen.MSgpd, chosen.MSgpd.name, th.01.gpd, date01)

    } else {
        cat("Fitting GPD",chosen.MSgpd.name,cr)
        sink("/dev/null")
        chosen.MSgpd      <- evgam(fmla.MSgpd[[chosen.MSgpd.name]], data.gpd, family="gpd", trace=2)
        sink()

        cat('Chosen MS GPD:',chosen.MSgpd.name,cr)
        print(summary(chosen.MSgpd))

        chosen.MSgpd$ms.thgpd.u <- ms.thgpd.u
        date01 <- data01[,c('time','doy')]
        save(file=st.msgpd, chosen.MSgpd, chosen.MSgpd.name, th.01.gpd, date01)
    }


    #    return(list(u.gpd=u.gpd, gpd01par=gpd01par, th=th.01.gpd, th.po=th.gpd.pred.o, th.u=ms.thgpd.u,
    #                     fgpd.d=fgpd.d, chosen.fgpd.d=chosen.fgpd.d))

}


##############################################################################
fn_apply_MSgpd <- function(data01, ms.thgpd.u, fit.qgam=NULL, st.qgam=NULL, st.msgpd=NULL, QUIET=TRUE) {

    if(!file.exists(st.msgpd)) {
        cat("fn_apply_MSgpd: ERROR",cr,"Must path to MSgpd file",cr)
        return(NULL)
    } else {
        cat("fn_apply_MSgpd: loading MSgpd",cr,st.msgpd,cr)
        load(st.msgpd, verb=TRUE)
    }

    if (is.null(fit.qgam)) {
        if(!file.exists(st.qgam)) {
            cat("fn_apply_MSgpd: ERROR",cr,"Must either supply fit.qgam object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_apply_MSgpd: loading qgam",cr,st.qgam,cr)
            load(st.qgam)
        }
    }

    # calc threshold
    data01$gpd.ktime <- 1 # evgam bug fix, values dont mater
    data01$ms.k      <- 1 # evgam bug fix, values dont mater
    cat("fn_apply_MSgpd: Calculating gpd threshold",cr)
    th.01.gpd        <- qdo(fit.qgam, ms.thgpd.u, predict, newdata=data01 ) # predict(gam.d.th, newdata=data01)

    data01$excess    <- data01$x - th.01.gpd
    data01$th.01.gpd <- th.01.gpd

    # data.gpd0       <- subset(data01, excess > 0)

    ### apply GPD
    cat("fn_apply_MSgpd: Making gpd parameters",cr)
    gpd01par <- predict(chosen.MSgpd, newdata=data01, type="response")
    u01      <- data01$uqgam
    excess   <- data01$excess
    dt1      <- trunc(length(excess)/20)
    cat("fn_apply_MSgpd: Completed 0 of",length(excess),cr)
    for(i in seq_along(excess)) {
        # if(excess[i]>0 & u01>ms.thgpd.u) {
        if(excess[i]>0) {
            sh          <- gpd01par$shape[i]
            sc          <- gpd01par$scale[i]
            u01[i] <- 1-(1-ms.thgpd.u)*pmax(0,(1+(sh*((excess[i])/sc))))^(-1/sh)
            # cat(i,sc,sh,u01[i],cr)
        } else if(is.finite(data01$uqgam[i])) {
            if(excess[i]<=0 & data01$uqgam[i]>ms.thgpd.u) {  # needed as sometimes q2p and evgam disagree
                u01[i] <- ms.thgpd.u
            }
        }
        if(i%%dt1==0 & !QUIET) cat("fn_apply_MSgpd: Completed ",i,"of",length(excess),cr)
    }

    return(u01)

}


##############################################################################
fn_diag_pre_proc <- function(savefile=NULL, width=10, height=7, DOPAUSE=FALSE){

    if(is.null(savefile)) {x11(); DOPAUSE=TRUE} else pdf(file=savefile, width=width, height=height)

    plot(om_data$obs$time, om_data$obs$data, pch=46, main='Obs')
    grid()
    if(DOPAUSE) readline("Continue?")

    plot(om_data$obs$info$doy, om_data$obs$data, pch=46, main='Obs')
    grid()
    if(DOPAUSE) readline("Continue?")

    plot(om_data$obs$info$sdoy, om_data$obs$data, pch=46, main='Obs')
    grid()
    if(DOPAUSE) readline("Continue?")

    plot(as.POSIXct(om_data$mod$time), om_data$mod$data, pch=46, main='CPM')
    grid()
    if(DOPAUSE) readline("Continue?")

    plot(om_data$mod$info$doy, om_data$mod$data, pch=46, main='CPM')
    grid()
    if(DOPAUSE) readline("Continue?")

    plot(om_data$mod$info$sdoy, om_data$mod$data, pch=46, main='CPM')
    grid()
    if(DOPAUSE) readline("Continue?")

    if(!is.null(savefile)) dev.off()
}


##############################################################################
fn_diag_qgam <- function(fit.qgam, savefile=NULL, width=10, height=7, DOPAUSE=FALSE, member=0){

    if(is.null(savefile)) x11() else pdf(file=savefile, width=width, height=height)

    ip50 <- which.min(abs(fit.qgam$ptiles-0.5))
    ip10 <- which.min(abs(fit.qgam$ptiles-0.1))
    ip90 <- which.min(abs(fit.qgam$ptiles-0.9))

    sink("/dev/null")
    qdo(fit.qgam, fit.qgam$ptiles[ip50], plot, main='q-50', page=1)
    qdo(fit.qgam, fit.qgam$ptiles[ip10], plot, main='q-10', page=1)
    qdo(fit.qgam, fit.qgam$ptiles[ip90], plot, main='q-90', page=1)
    sink()

    par(mfcol=c(1,1))
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    par(mgp=c(3,1,0))

    years     <- trunc(data01$time)
      plot(data01$doy, data01$sdoy,  pch=46, main='Check iob and doy scaling')
    grid()
    if(DOPAUSE) readline("Continue?")

      plot(data01$time, data01$x, pch=46, main='Obs and CPM data check wrt time')
    grid()
    if(DOPAUSE) readline("Continue?")

      plot(data01$doy, data01$x, pch=46, main='Obs and CPM data check wrt doy')
    grid()
    if(DOPAUSE) readline("Continue?")

      plot(data01$sdoy, data01$x, pch=46, main='Obs and CPM data check wrt SCALED doy')
    grid()
    if(DOPAUSE) readline("Continue?")

    ### quantiles
      plot(data01$stime,       data01$x, type='n', main='Obs and CPM 90th quantile wrt SCALED time')
    points(data01$stime,       data01$x,      pch=46)
    points(data01$stime, q.doy.time[ip90, ], pch=46, col=3)
    grid()
    if(DOPAUSE) readline("Continue?")

      plot(data01$sdoy,       data01$x, type='n', main='Obs and CPM 90th quantile wrt SCALED DOY')
    points(data01$sdoy,         data01$x,    pch=46)
    points(data01$sdoy, q.doy.time[ip90, ], pch=46, col=3)
    grid()
    if(DOPAUSE) readline("Continue?")


    iy1.o <- which(years==min(years)+1)
    iy2.o <- which(years==max(years)-1)


    # just obs annual signal
      plot(data01$sdoy,       data01$x, type='n', main='Obs 90 & 50 th quantile wrt scaled DOY')
    points(data01$sdoy,                data01$x                       , pch=46,cex=.3)
     lines(data01$sdoy[iy1.o[-1]],    q.doy.time[ip90,][iy1.o[-1]], lty=1, col=2, lwd=2)
     lines(data01$sdoy[iy2.o[-1]],    q.doy.time[ip90,][iy2.o[-1]], lty=2, col=2, lwd=3)
     lines(data01$sdoy[iy1.o[-1]],    q.doy.time[ip50,][iy1.o[-1]], lty=1, col=4, lwd=2)
     lines(data01$sdoy[iy2.o[-1]],    q.doy.time[ip50,][iy2.o[-1]], lty=2, col=4, lwd=3)
    grid()
    legend('topleft',c('Obs','90 %ile','50 %ile','start','end'),col=c(1,2,4,1,1),pch=c(3,NA,NA,NA,NA),lty=c(NA,1,1,1,2),bty='n')

   
    # start end
      plot(data01$sdoy,       data01$x, type='n', main='Obs and CPM 90th quantile wrt SCALED DOY')
     lines(data01$sdoy[iy1.o[-1]],    q.doy.time[ip90,][iy1.o[-1]], pch=46, col=3)
    points(data01$sdoy[iy2.o],                  data01$x[iy2.o],     pch=1,  cex=.6)
     lines(data01$sdoy[iy2.o[-1]],    q.doy.time[ip90,][iy2.o[-1]], pch=46, col=3)
    grid()
    legend('topleft',c('Obs','Obs 90%ile start/end'),col=c(1,3),pch=c(3,NA),lty=c(NA,1),bty='n')

    # all quantiles
    k0 <- seq(1,length(ptiles),by=4) # seq(2,98,by=4)
    par(mfcol=c(2,2))
    par(mar=c(3,3,2,1))
    par(mgp=c(2,1,0))

    plot(data01$sdoy,       data01$x, type='n', main='Obs start all %iles wrt DOY')
    for(k in k0) lines(data01$sdoy[iy1.o[-1]],    q.doy.time[k, ][iy1.o[-1]])
    grid()

    plot(data01$sdoy,       data01$x, type='n', main='Obs end all %iles wrt DOY')
    for(k in k0) lines(data01$sdoy[iy2.o[-1]],    q.doy.time[k, ][iy2.o[-1]])
    grid()



    if(!is.null(savefile)) dev.off()

}


##############################################################################
fn_diag_evgam <- function(savefile=NULL, st.msgpd=NULL, width=10, height=7, DOPAUSE=FALSE){

    if(!file.exists(st.msgpd)) {
        # cat("fn_apply_MSgpd: ERROR",cr,"Must path to MSgpd file",cr)
        # return(NULL)
        cat("fn_apply_MSgpd: NOT loading a new MSgpd",cr)
    } else {
        cat("fn_apply_MSgpd: loading MSgpd",cr,st.msgpd,cr)
        load(st.msgpd, verb=TRUE)
    }

    if(is.null(savefile)) x11() else pdf(file=savefile, width=width, height=height)

    # smooths
    par(oma=c(0,0,1,0))
    plot(chosen.MSgpd)
    mtext(chosen.MSgpd.name,outer=T,at=c(0.01),adj=0,col=1,line=-0.3) # left justified
    if(DOPAUSE) readline("Continue?")

    # print out model summary ?onto plot
    s1 <- capture.output(summary(chosen.MSgpd))
    plot(0:1,0:1,ty='n',main=paste('Summary stats for EV GAM',chosen.MSgpd.name),xlab=NA, ylab=NA, yaxt='n', xaxt='n')
    x1 <- 0.1
    y1 <- 1
    for(i in seq_along(s1)) {
        text(x1,y1,s1[i],adj=0,cex=0.5)
        y1 <- y1 - 1/length(s1)
    }
    if(DOPAUSE) readline("Continue?")

    if(!is.null(savefile)) dev.off()

}


##############################################################################
fn_diag_q2p <- function(savefile=NULL, width=10, height=7, DOPAUSE=FALSE){

    if(is.null(savefile)) x11() else pdf(file=savefile, width=width, height=height)

    ### refit gam to check stationarity
    x.l  <- qlaplace(data01$u)
    if(length(x.l) > 20000) i1 <- sort(sample(seq_along(x.l),2e4)) else i1 <- seq_along(x.l)

    ### using sdoy and gmst as covariate
    # test.data <- data.frame(x.l=x.l[i1], sdoy=data01$sdoy[i1], gmst=data01$gmst[i1], isobs=data01$isobs[i1])
    # fit.test  <- mqgam( x.l ~ s(sdoy, bs='cc') +s(gmst), data=test.data, qu = c(0.90, 0.98))
    # test.data <- data.frame(x=x.l[i1], sdoy=data01$sdoy[i1], gmst=data01$gmst[i1], class=data01$class[i1])

    test.data     <- data01[i1,]
    test.data$x.l <- x.l[i1]
    # fit.test  <- mqgam( fit.qgam$fmla, data=test.data, qu = c(0.90, 0.98))
    # q90        <- qdo(fit.test, 0.9, predict, newdata=test.data)
    # iob       <- which(data01$isobs==1)
    # iob2      <- which(test.data$class=='obs')

    fit.test  <- mqgam( x.l ~ s(sdoy, bs='cc') +s(stime), data=test.data, qu = c(0.5, 0.90, 0.99))
    q50t      <- qdo(fit.test, 0.5,  predict, newdata=test.data)
    q90t      <- qdo(fit.test, 0.9,  predict, newdata=test.data)
    q99t      <- qdo(fit.test, 0.99, predict, newdata=test.data)
    fit.test  <- mqgam( x.l ~ s(sdoy, bs='cc') +s(gmst), data=test.data, qu = c(0.5, 0.90, 0.99))
    q50g      <- qdo(fit.test, 0.5,  predict, newdata=test.data)
    q90g      <- qdo(fit.test, 0.9,  predict, newdata=test.data)
    q99g      <- qdo(fit.test, 0.99, predict, newdata=test.data)
    iob       <- which(data01$isobs==1)
    iob2      <- which(test.data$isobs==1)

    if(!is.null(test.data$gmst)) {
      plot(   data01$gmst,         x.l, type='n', main='Stationarity check wrt GMST')
    points(   data01$gmst[ iob],   x.l[  iob], pch=46)
    points(   data01$gmst[-iob],   x.l[ -iob], pch=46, col=2)
    points(test.data$gmst[ iob2], q50g[ iob2], pch=46, col=3)
    points(test.data$gmst[-iob2], q50g[-iob2], pch=46, col=4)
    points(test.data$gmst[ iob2], q90g[ iob2], pch=46, col=3)
    points(test.data$gmst[-iob2], q90g[-iob2], pch=46, col=4)
    points(test.data$gmst[ iob2], q99g[ iob2], pch=46, col=3)
    points(test.data$gmst[-iob2], q99g[-iob2], pch=46, col=4)
    abline(h=qlaplace(0.5)) # q50g[1])
    abline(h=qlaplace(0.9)) # q90g[1])
    abline(h=qlaplace(0.99)) # q99g[1])
    grid()
    if(DOPAUSE) readline("Continue?")

    # plot(test.data$gmst, 100*(q90-q90[1])/q90[1],pch=46, main="% nonstationarity")
    # if(DOPAUSE) readline("Continue?")
    }

    if(!is.null(data01$stime)) {
      plot(   data01$stime,         x.l, type='n', main='Stationarity check wrt time')
    points(   data01$stime[ iob],   x.l[  iob], pch=46)
    points(   data01$stime[-iob],   x.l[ -iob], pch=46, col=2)
    points(test.data$stime[ iob2], q50t[ iob2], pch=46, col=3)
    points(test.data$stime[-iob2], q50t[-iob2], pch=46, col=4)
    points(test.data$stime[ iob2], q90t[ iob2], pch=46, col=3)
    points(test.data$stime[-iob2], q90t[-iob2], pch=46, col=4)
    points(test.data$stime[ iob2], q99t[ iob2], pch=46, col=3)
    points(test.data$stime[-iob2], q99t[-iob2], pch=46, col=4)
    abline(h=qlaplace(0.5)) # q50g[1])
    abline(h=qlaplace(0.9)) # q90g[1])
    abline(h=qlaplace(0.99)) # q99g[1])
    grid()    
    if(DOPAUSE) readline("Continue?")

    # plot(test.data$stime, 100*(q90-q90[1])/q90[1],pch=46, main="% nonstationarity")
    # if(DOPAUSE) readline("Continue?")
    }

    if(!is.null(data01$sdoy)) {
      plot(   data01$sdoy,       x.l, type='n', main='Stationarity check wrt sdoy')
    points(   data01$sdoy[ iob], x.l[  iob], pch=46)
    points(   data01$sdoy[-iob], x.l[ -iob], pch=46, col=2)
    points(test.data$sdoy[ iob2], q90g[ iob2], pch=46, col=3)
    points(test.data$sdoy[-iob2], q90g[-iob2], pch=46, col=4)
    abline(h=1)
    grid()
    if(DOPAUSE) readline("Continue?")
    }

    if(!is.null(data01$stime) & !is.null(data01$sdoy)) {
    ### using sdoy and stime as covariate
    # test.data <- data.frame(x.l=x.l[i1], sdoy=data01$sdoy[i1], stime=data01$stime[i1], isobs=data01$isobs[i1])

    # test.data$x.l <- x.l[i1]
    # fit.test  <- mqgam( x.l ~ s(sdoy, bs='cc') +s(stime), data=test.data, qu = c(0.90, 0.98))
    # q90        <- qdo(fit.test, 0.9, predict, newdata=test.data)
    # iob       <- which(data01$isobs==1)
    # iob2      <- which(test.data$isobs==1)

    #   plot(   data01$stime,       x.l, type='n', main='Stationarity check wrt time')
    # points(   data01$stime[ iob], x.l[  iob], pch=46)
    # points(   data01$stime[-iob], x.l[ -iob], pch=46, col=2)
    # points(test.data$stime[ iob2], q90[ iob2],    pch=46, col=3)
    # points(test.data$stime[-iob2], q90[-iob2]+.1, pch=46, col=4)
    # points(test.data$stime[ iob2], q98[ iob2],    pch=46, col=3)
    # points(test.data$stime[-iob2], q98[-iob2]+.1, pch=46, col=4)
    # abline(h=q90[1])
    # abline(h=q98[1])
    # grid()
    # if(DOPAUSE) readline("Continue?")
    }

    # plot(data01$u, data01$uqgam,pch=46,xlim=c(ms.thgpd.u-0.01,1.0),ylim=c(ms.thgpd.u-0.01,1.0), main='Impact of evGAM')
        i1 <- which(data01$uqgam>ms.thgpd.u-0.01)
      plot(x.l[i1], data01$uqgam[i1],pch=46,ylim=c(ms.thgpd.u-0.01,1.0), main='Impact of evGAM')
    points(x.l[i1], data01$u[i1],    pch=46,col=2)
    if(DOPAUSE) readline("Continue?")

    if(!is.null(savefile)) dev.off()

}

##############################################################################
fn_diag_general <- function(savefile=NULL, width=10, height=7, dodoy=200, DOPAUSE=FALSE){

    if(is.null(savefile)) x11() else pdf(file=savefile, width=width, height=height)

    ### testing plots
    i1s    <- which(data01$isobs==1 & data01$doy==dodoy)
    test01 <- data01[i1s,]
    # test.u <- fn_apply_MSgpd(test01, ms.thgpd.u, fit.qgam=NULL, st.qgam=st_qgam, st.msgpd=st_msgpd)
    # plot(test01$time, test01$u, pch=20, cex=.3, main="Test u obs doy=200")
    # grid()
    # plot(test01$time, qlaplace(test01$u), pch=20, cex=.3, main="Test u obs doy=200")
    # grid()

    # qgam
    q50s      <- qdo(fit.qgam, 0.5,  predict, newdata=test01)
    q90s      <- qdo(fit.qgam, 0.9,  predict, newdata=test01)
    q99s      <- qdo(fit.qgam, 0.99,  predict, newdata=test01)
      plot(data01$gmst, data01$x, pch=46, main="qgam evgam")
    lines(test01$gmst, q50s, lwd=2, col=2)
    lines(test01$gmst, q90s, lwd=2, col=3)
    lines(test01$gmst, q99s, lwd=2, col=4)
    # evgam
    gpd01par <- predict(chosen.MSgpd, newdata=test01, type="response")
    evth     <- qdo(fit.qgam, ms.thgpd.u,  predict, newdata=test01)
    ev99     <- evth + (gpd01par$scale/ gpd01par$shape)* ( ( (1 - 0.99)/(1 - ms.thgpd.u) )^(-gpd01par$shape) -1 )
    lines(test01$gmst, ev99, lwd=2, lty=2, col=1)
    # add winter
    i1w    <- which(data01$isobs==1 & data01$doy==20)
    test01 <- data01[i1w,]
    q50w      <- qdo(fit.qgam, 0.5,  predict, newdata=test01)
    q90w      <- qdo(fit.qgam, 0.9,  predict, newdata=test01)
    q99w      <- qdo(fit.qgam, 0.99,  predict, newdata=test01)
    lines(test01$gmst, q50w, lwd=2, col=2)
    lines(test01$gmst, q90w, lwd=2, col=3)
    lines(test01$gmst, q99w, lwd=2, col=4)
    gpd01par <- predict(chosen.MSgpd, newdata=test01, type="response")
    evth     <- qdo(fit.qgam, ms.thgpd.u,  predict, newdata=test01)
    ev99     <- evth + (gpd01par$scale/ gpd01par$shape)* ( ( (1 - 0.99)/(1 - ms.thgpd.u) )^(-gpd01par$shape) -1 )
    lines(test01$gmst, ev99, lwd=2, lty=2, col=1)

    # q2p
    test01 <- data01[i1s,]
    test01$x <- q90s[1]
    p17uqgam <- fn_apply_q2p(test01, qgam.q2p.fn[i1s], st.q2p=NULL, SANITY=FALSE)
    plot(test01$gmst, p17uqgam, ty='l', pch=20, cex=.3, main="q2p sanity check q90s[1]",ylim=c(0,1))
    grid()

    # p2q
    load(MSconfig$files$st_p2q,  temp_env <- new.env(), verb=TRUE)
    l.vars <- as.list(temp_env)
    test01 <- data01[i1s,]
    for(i in seq_along(i1s)) test01$o[i] <- l.vars$qgam.p2q.fn[[i1s[i]]](0.9)
    plot(test01$gmst, test01$o, ty='l', pch=20, cex=.3, main="p2q sanity check q90s[1]")
    grid()

    p0 <- c(0.001*1:9,0.01*1:99,1-0.001*9:1)
    # changes
    q.start <- l.vars$qgam.p2q.fn[[i1s[1]]](p0)
    q.end   <- l.vars$qgam.p2q.fn[[i1s[length(i1s)]]](p0)
    j1 <- which(!is.na(q.start) & !is.na(q.end))
    plot(q.start[j1], q.end[j1], pch=20, cex=.3, main="qqplot start vs end", xlim=range(c(q.start[j1],q.end[j1]),na.rm=TRUE), ylim=range(c(q.start[j1],q.end[j1]),na.rm=TRUE) )
    abline(0,1,col=2)
    grid()

    p.start.at.end <- sapply(q.start, qgam.q2p.fn[[i1s[length(i1s)]]])
    plot((0.01*1:99)[j1], p.start.at.end[j1], pch=20, cex=.3, main="p at start quantiles evaluated at end", xlim=c(0,1), ylim=c(0,1))
    abline(0,1,col=2)
    grid()

    # cdf start vs end
    plot(q.start, p0, xlim=range(c(q.start,q.end),na.rm=TRUE), ylim=c(0,1), pch=20, cex=.3, main="cdf start vs end")
    points(q.end, p0)
    grid()

    # pdf start vs end
    k1 <- which(!is.na(q.start))
    k2 <- which(!is.na(q.end))

    q2p.s <- smooth.spline(q.start[k1], p0[k1], spar=0.9)
    x.s   <- seq(min(q.start[k1]), max(q.start[k1]), length.out=100)
    p2.s  <- predict(q2p.s,x.s)

    q2p.e <- smooth.spline(q.end[k2], p0[k2], spar=0.9)
    x.e   <- seq(min(q.end[k2]), max(q.end[k2]), length.out=100)
    p2.e  <- predict(q2p.e,x.e)
    
    plot(p2.s$x[-1],diff(p2.s$y), ylim=c(0,0.03), xlim=range(c(q.start,q.end),na.rm=TRUE), pch=20, cex=.3, main="pdf start vs end")
    points(p2.e$x[-1],diff(p2.e$y))

    # evgam
    i1s    <- which(data01$isobs==1 & data01$doy==dodoy)
    test01 <- data01[i1s,]
    gpd01par <- predict(chosen.MSgpd, newdata=test01, type="response")
    evth     <- qdo(fit.qgam, ms.thgpd.u,  predict, newdata=test01)
    ev99     <- evth + (gpd01par$scale/ gpd01par$shape)* ( ( (1 - 0.99) /(1 - ms.thgpd.u) )^(-gpd01par$shape) -1 )
    ev999    <- evth + (gpd01par$scale/ gpd01par$shape)* ( ( (1 - 0.999)/(1 - ms.thgpd.u) )^(-gpd01par$shape) -1 )
      plot(test01$gmst,  ev99, pch=20, cex=.3, main="evgam 99th and 99.9th percentiles doy=200", ylim=range(c(ev99,ev999),na.rm=TRUE))
    points(test01$gmst, ev999, pch=20, cex=.3, col=2)
    grid()  
    
      plot(test01$time,  ev99, pch=20, cex=.3, main="evgam 99th and 99.9th percentiles doy=200", ylim=range(c(ev99,ev999),na.rm=TRUE))
    points(test01$time, ev999, pch=20, cex=.3, col=2)
    grid()  

    if(!is.null(savefile)) dev.off()

}

##############################################################################
trim_qgam <- function(fit.qgam){

    i.fit         <- which(names(fit.qgam)=='fit')
    fit.small     <- fit.qgam[-i.fit]
    fit.small$fit <- list()

    l.2rm <- c('residuals', 'fitted.values', 'linear.predictors', 'weights', 'prior.weights', 'working.weights', 'y', 'z', 'dw.drho', 'hat', 'offset')

    n.fit  <- names(fit.qgam$fit[[1]])
    i.rm   <- which(n.fit %in% l.2rm)
    i.keep <- seq_along(n.fit)[-i.rm]

    for(i in seq_along(fit.qgam$fit)) {
        fit1               <- fit.qgam$fit[[i]][i.keep] # seems to produce the smallest
        attr(fit1,'class') <- attr(fit.qgam$fit[[i]],'class')
        fit.small$fit[[i]] <- fit1
    }

    names(fit.small$fit) <- names(fit.qgam$fit)
    return(fit.small)
}

# ##############################################################################
# xxx <- function(){

# }

#
