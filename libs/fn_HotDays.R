##############################################################################################
### Make HotDays functions
##############################################################################################

### NB NB in transition from fn_JointHotDays.R

source("/home/users/simon.brown/old-home/extremes/R/EVlibs/simon_gpd_reference.R")

reload_HD_ref<-function(freference=NULL) {
    load(file=freference, verb=T)
    assign("HDconfig", HDconfig, envir = .GlobalEnv)
    # lapply(list(MSconfig), list2env ,envir = .GlobalEnv)
    list2env(HDconfig,               envir = .GlobalEnv)
    list2env(HDconfig$general,       envir = .GlobalEnv)
    list2env(HDconfig$events,        envir = .GlobalEnv)
    list2env(HDconfig$files,         envir = .GlobalEnv)
    list2env(HDconfig$initI        , envir = .GlobalEnv)
    list2env(HDconfig$initV        , envir = .GlobalEnv)
    list2env(HDconfig$fitht        , envir = .GlobalEnv)
    list2env(HDconfig$term         , envir = .GlobalEnv)
}

findall <- function(find.these,in.this, firstoccurrenceonly=FALSE) {
    if(length(find.these)>0 & length(in.this)>0) {
	n1 <- length(find.these)
	l.all <- vector(mode='list',length=n1)
	names(l.all) <- rep('index',n1)
	for (i in 1:n1) {
		i1 <- which(find.these[i]==in.this)
		if(firstoccurrenceonly)	l.all[[i]] <- i1[1] else
		                        l.all[[i]] <- i1
	}
	 i.all<-unlist(sapply(l.all, "["))

	 i.all
    } else return(NULL)
}

fn_prob_of_event <- function(st.varx2, ref_event_year) {

    source("../sim.from.ht.Gnoise.R")
    source("/home/users/simon.brown/old-home/extremes/heatwaves/code/libs/fn_useful.R")
    library(evd)

    VariableNames <- load(st.varx2)
    varx2 <- get(VariableNames[grep('d.',VariableNames)])
    rm(VariableNames)

    th.gpd.u.cet    <- varx2$cet$gpdfit$th.u
    th.gpd.u.pnd    <- varx2$pnd$gpdfit$th.u

    # all.y        <- sort(unique(varx2$cet$d.pred.o$time))
    all.y        <- seq(from=min(varx2$pnd$data$time), to=max(varx2$pnd$data$time))

    p.isnob <- which(varx2$pnd$data$m1==1)
    c.isnob <- which(varx2$cet$data$m1==1)
    p.isob  <- which(varx2$pnd$data$isobs==1)
    c.isob  <- which(varx2$cet$data$isobs==1)

    i76c     <- c.isob[which(varx2$cet$data$time[c.isob]==ref_event_year)]
    i76p     <- p.isob[which(varx2$pnd$data$time[p.isob]==ref_event_year)]
    d.76     <- list()
    d.76$cet <- varx2$cet$data01[i76c,]
    d.76$pnd <- varx2$pnd$data01[i76p,]

    ### pnd #####################################################################
    # p.newd       <- data.frame(time=plt.y, isobs=1, ktime=1, m1=0, m2=0, m3=0, m4=0, m5=0, m6=0, m7=0, m8=0, m9=0, m10=0, m11=0, m12=0)
    x1           <- varx2$pnd$pred.o[1,]
    x1$ktime     <- 1
    p.newd       <- x1[rep(seq_len(nrow(x1)), length(all.y)), ]
    iy           <- which(varx2$pnd$pred.o$time %in% all.y)
    gpd.p.y      <- varx2$pnd$pred.o$time[iy]
    p.th.gpd     <- varx2$pnd$gpdfit$th.po[iy]
    p.newd$time  <- gpd.p.y
    p.gpdpar     <- predict(varx2$pnd$gpdfit$fit.gpd.d, p.newd, type="response")
    Pxgtx.pnd76  <- NULL
    for(i in seq_along(gpd.p.y)) Pxgtx.pnd76[i] <- pgpd(d.76$pnd$x, loc=p.th.gpd[i], scale=p.gpdpar$scale[i], shape=p.gpdpar$shape[i],lo=F) *(1-th.gpd.u.pnd)

    ###########################################################################################
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # check that the probabilities of 76 in 1976 matches between Pxgtx.pnd76 and varx2$pnd$data$u
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ###########################################################################################
    # 2023.02.010 probabilities tie up so all good
    #  plot(varx2$pnd$data$time[p.isob], varx2$pnd$data$x[p.isob])
    # lines(varx2$pnd$data$time[p.isob], varx2$pnd$gpdfit$th[p.isob])
    # lines(varx2$pnd$data$time[p.isob], varx2$pnd$gpdfit$th[p.isob])

    #    plot(varx2$pnd$data$time[p.isob], varx2$pnd$data$u0[p.isob])
    #  points(varx2$pnd$data$time[p.isob], varx2$pnd$data$u[p.isob],       pch=3, cex=.3, col=2)
    #  points(varx2$pnd$data$time[p.isob], varx2$pnd$gpdfit$u.gpd[p.isob], pch=3, cex=.3, col=3)

    #   plot(varx2$pnd$data$time[p.isob], qlaplace(varx2$pnd$data$u0[p.isob]),ylim=c(-4,6))
    # points(gpd.p.y,                     qlaplace(1-Pxgtx.pnd76),            pch=1, cex=.3, col=4)
    # points(varx2$pnd$data$time[p.isob], qlaplace(varx2$pnd$data$u[p.isob]), pch=3, cex=.4, col=2)


    # fix those points that are below the threshold.
    inc <- which((1- Pxgtx.pnd76)<=(th.gpd.u.pnd))
    for(i in inc) {
        iy1          <- which(varx2$pnd$pred.o$time == gpd.p.y[i])
        if(d.76$pnd$x > varx2$pnd$qo.from.qgam[1,iy1]) { # is it greater than the minimum in the ecdf
            Pxgtx.pnd76[i] <- 1-varx2$pnd$qgam.q2p.o.fn[[iy1]](d.76$pnd$x)
        } else Pxgtx.pnd76[i] <- 1
    }
    # extract idx of year of interest
    p.j76          <- which(gpd.p.y==ref_event_year)
    Pxgtx.pnd76.76 <- Pxgtx.pnd76[p.j76]


    ### cet ###################################################################################

    ## calculate independent probabilities of cet > cet76
    # c.newd <- data.frame(time=plt.y, isobs=1, ktime=1, m1=0, m2=0, m3=0, m4=0, m5=0, m6=0, m7=0, m8=0, m9=0, m10=0, m11=0, m12=0)
    x1           <- varx2$cet$pred.o[1,]
    x1$ktime     <- 1
    c.newd       <- x1[rep(seq_len(nrow(x1)), length(all.y)), ]
    iy           <- which(varx2$cet$pred.o$time %in% all.y)
    gpd.p.y      <- varx2$cet$pred.o$time[iy]
    c.th.gpd     <- varx2$cet$gpdfit$th.po[iy]
    c.newd$time  <- gpd.p.y
    c.gpdpar     <- predict(varx2$cet$gpdfit$fit.gpd.d, c.newd, type="response")
    Pxgtx.cet76    <- NULL
    for(i in seq_along(gpd.p.y)) Pxgtx.cet76[i] <- pgpd(d.76$cet$x, loc=c.th.gpd[i], scale=c.gpdpar$scale[i], shape=c.gpdpar$shape[i],lo=F) *(1-th.gpd.u.cet)

    # fix those points that are below the threshold.
    # load(varx2$cet$st.qgam)
    inc <- which((1- Pxgtx.cet76)<=(th.gpd.u.cet))
    for(i in inc) {
        iy1          <- which(varx2$pnd$pred.o$time == gpd.p.y[i])
        if(d.76$cet$x > varx2$cet$qo.from.qgam[1,iy1]) { # is it greater than the minimum in the ecdf
            Pxgtx.cet76[i] <- 1-varx2$cet$qgam.q2p.o.fn[[iy1]](d.76$cet$x)
        } else Pxgtx.cet76[i] <- 1
    }
    # extract idx of year of interest
    c.j76          <- which(gpd.p.y==ref_event_year)
    Pxgtx.cet76.76 <- Pxgtx.cet76[c.j76]

    ### cet | pnd
    ## NB Prob of exceeding pnd76 changes through time
    ## so need to calculate Pcet_gtcet76.cond.pnd76 using the time varying Pxgtx.pnd76
    # need timeseries of the pnd.l
    ts.pnd.l                <- qlaplace(1-Pxgtx.pnd76)
    Pcetgt.cet76.cond.pnd76 <- NULL
    ht.stat <- varx2$HT
    for(i in seq_along(ts.pnd.l)){
        if(ts.pnd.l[i]>=qlaplace(ht.stat$mod.lev.u)) {
            sim0                       <- sim.from.ht.Gnoise(ht.stat,2e6, crit.lev.l=ts.pnd.l[i], GNOISE=TRUE, SHRINK=FALSE, bw=0.4)
            iex76                      <- which(sim0$y > qlaplace(1-Pxgtx.cet76[i]) )
            Pcetgt.cet76.cond.pnd76[i] <- (1 - plaplace(ts.pnd.l[i]) ) * length(iex76)/length(sim0$y) # P(pnd>76) * P(cet>76)
            if(i %% 10 == 0) cat(i,Pcetgt.cet76.cond.pnd76[i],cr)
        } else cat(i, 'pnd < ht.model.threshold',cr)
    }

    rp.ex76.cond.p     <- 1/Pcetgt.cet76.cond.pnd76
    rp.ex76.cond.p.ind <- 1/(Pxgtx.pnd76*Pxgtx.cet76)

    dout <- list()
    dout$c.gpdpar                <- c.gpdpar
    dout$p.gpdpar                <- p.gpdpar
    dout$rp.ex76.cond.p          <- rp.ex76.cond.p
    dout$rp.ex76.cond.p.ind      <- rp.ex76.cond.p.ind
    dout$Pcetgt.cet76.cond.pnd76 <- Pcetgt.cet76.cond.pnd76
    dout$Pxgtx.cet76             <- Pxgtx.cet76
    dout$Pxgtx.pnd76             <- Pxgtx.pnd76
    dout$p.j76                   <- p.j76
    dout$c.j76                   <- c.j76
    dout$all.y                   <- all.y
    dout$gpd.p.y                 <- gpd.p.y

    return(dout)
}

ChainsFromLaplace <- function(xdata.l, xdate, dtime=NULL, chain.lev.u, chain.length=chain.length, DOPLOT=FALSE, DOPAUSE=FALSE) {

    # derived from /home/users/simon.brown/old-home/extremes/R/gitEVlibs/hugo/temporal-sjb/ChainsFromUObsSjb3.R

    # source("ChainsFromLaplace.R")
    # source('/home/users/simon.brown/old-home/extremes/R/gitEVlibs/hugo/temporal-sjb/ChainsFromLaplace.R')

    # third attempt to extract chains from the obs but in a way that mimics the way we simulate chains
    #
    #     - A identify all points that are above threshold of interest
    #     - B find the maximum of these
    #     - C go forwads and backward till event drops below threshold or hits end of season in both directions or delta-time is > 1day (indicating end of season)
    #     - D mask out the points of this chain
    #     - E find the max of the remaining
    #     - F back to C untill all points masked out
    #
    #  xdata.l         : data_on_laplace_scale
    #  xdate           : time/dates of obs. can be POSIX PCICT or anything really.
    #  dtime           : the maximum time interval allowed between two Tmax.  A POSIX object.  If not set will use min(diff(xdate))
    #                    this is used to detect gaps in timeseries
    #  chain.lev.u     : modelling level on unit scale
    #
    # e.g. > obs_chains <- ChainsFromLaplace(tx.l, dtx, dtime=dtime, hw.lev.u3, chain.length=chain.length, DOPLOT=FALSE)

    # 2023.06
    #       removed seas_length dependence
    #       now assumes whole year will be used and hotseason subselection happens later

    chain.lev.l  <- qlaplace(chain.lev.u)

    if(is.null(dtime)) dtime <- 0.98 * (2* median(diff(xdate)))  # slightly less than two time intervals  # min(diff(xdate))
    cat("ChainsFromLaplace: Using dtime",dtime,cr)

    DONOTICEF <- TRUE
    DONOTICEB <- TRUE

    # check if there are any NA in data
    iNA <- which(is.na(xdata.l))
    if(length(iNA)>0){
        cat("ChainsFromLaplace: NA data found in xdata.l.  Setting to 0",cr)
        xdata.l[iNA] <- 0
    }
    if(DOPLOT) {
        plot(xdata.l,pch=46)
        abline(h=chain.lev.l)
        points(seq_along(xdata.l)[iNA],rep(-2,length(iNA),col=2),pch=3,cex=.3,col=2)
    }
    
    ichainsYr   <- list()
    all_chains  <- list()
    all_ichains <- list()
    all_absichains <- NULL
    all_count   <- 0
    xdata.i     <- seq_along(xdata.l)

    ichains    <- list()
    s1         <- xdata.l
    iremaining <- which(s1>=chain.lev.l)
    count      <- 0

    while(length(iremaining) > 0) {
        count <- count+1
        i0    <- which.max(s1[iremaining])
        xi    <- iremaining[i0]
        fi    <- xi
        bi    <- xi

        GOFWD <- TRUE
        while ((fi-xi)<chain.length & GOFWD & fi<length(s1)) {
            # if(s1[fi]>=chain.lev.l & ( xdate[fi+1]-xdate[fi]==dtime | ((|=360 |=365 |=366)&(=1) )) fi <- fi + 1 else GOFWD <- FALSE
            if(s1[fi]>=chain.lev.l & (xdate[fi+1]-xdate[fi]<=dtime) ) fi <- fi + 1 else GOFWD <- FALSE
        }
            # >359 to cope with dec-31 to jan-1
        #if( fi==xi ) fi <- fi +1 # fwd length=0
        if( ((fi-xi)>chain.length) & DONOTICEF ) {cat("CAUTION: FWD chain length not long enough",cr); DONOTICEF=FALSE}

        GOBWD <- TRUE
        while ((xi-bi)<chain.length & GOBWD & bi>1) {
            # if(s1[bi]>=chain.lev.l & ( xdate[fi+1]-xdate[fi]==dtime | xdate[fi+1]-xdate[fi]>359 )) bi <- bi - 1 else GOBWD <- FALSE
            if(s1[bi]>=chain.lev.l & (xdate[bi]-xdate[bi-1]<=dtime) ) bi <- bi - 1 else GOBWD <- FALSE
        }
        #if( bi==xi ) bi <- bi -1 # bwd length=0
        if( ((xi-bi)>chain.length) & DONOTICEB ) {cat("CAUTION: BWD chain length not long enough",cr); DONOTICEB=FALSE}
        # NB this produces indecies that are one before and one after the indecies that are above chain.lev.l

        if(bi<1)          b_buffers <- TRUE
        if(fi>length(s1)) f_buffers <- TRUE

        ichains[[count]] <- c(bi,xi,fi)
        i2               <- c(max(c(1,bi)),xi,min(c(fi,length(s1))) ) # make sure indecies are within range of dim(xdata.l)[1]
        ###if(i2[3]-i2[1]>(2*chain.length-1)) cat(cr,"ERROR-ChainsFromUObsSjb3: chain length not long nough. Needs to be >",ceiling((i2[3]-i2[1])/2),cr,cr)

        allChainLength <- 2*chain.length+2
        achain                                                                <- rep(0,length=allChainLength)
        achain[ (chain.length-(i2[2]-i2[1])) : (chain.length+(i2[3]-i2[2])) ] <- xdata.l[ i2[1]:i2[3] ]
        absi                                                                  <- xdata.i[ i2[1]  ]

        all_count                 <- all_count +1
        all_chains[[all_count]]   <- achain
        all_absichains[all_count] <- absi
        all_ichains[[all_count]]  <- ichains[[count]]

        s1[ i2[1]:i2[3] ]  <- chain.lev.l - 0.1*abs(chain.lev.l)

        iremaining <- which(s1>=chain.lev.l)

        #if(xi==1 | xi==2 | xi==149 | xi==150) browser()

        if(DOPLOT) {
            points(ichains[[count]][2],                    xdata.l[ichains[[count]][2]],                       col=count,pch=3,cex=.3)
                lines(ichains[[count]][1]:ichains[[count]][3],xdata.l[ichains[[count]][1]:ichains[[count]][3]],col=count)
               points(ichains[[count]][1]:ichains[[count]][3],xdata.l[ichains[[count]][1]:ichains[[count]][3]],col=count,pch=3,cex=.5)
            # cat('Count', count,'chain',ichains[[count]],cr)
            # cat("Points remaining",length(iremaining),cr)
            # browser()
            if(DOPAUSE) readline("Continue?")
        }

        if(all_count%%1000==0)cat("ChainsFromLaplace: HW count and points remaining:",all_count,length(iremaining),cr)

    } #end while

    ichainsYr[[count]] <- ichains

    all_chains2      <- unlist(all_chains)
    dim(all_chains2) <- c(allChainLength, length(all_chains2)/(allChainLength))
    all_ichains2     <- unlist(all_ichains)
    dim(all_ichains2)      <- c(3, length(all_ichains2)/(3))
    dimnames(all_ichains2) <- list(c('day-before','peak-day','day-after'))
    names(all_absichains)  <- 'Absolute index of day before heatwave of -summer- data'

    # return(list(chains.l=all_chains2, absichains=all_absichains, ichains.l=all_ichains2, indeciesByYear=ichainsYr, thresh.u=chain.lev.u, months=NULL, seas.length=seas_length, chain.length=chain.length))
    return(list(chains.l=all_chains2, absichains=all_absichains, ichains.l=all_ichains2, indeciesByYear=ichainsYr, thresh.u=chain.lev.u, chain.length=chain.length))

}

ChainsFromGeneric <- function(x1,  th.x1, chain.length=chain.length, DOPLOT=FALSE, DOPAUSE=FALSE) {

    # derived from ChainsFromLaplace
    xdate <- seq_along(x1)

    dtime <- 0.98 * (2* median(diff(xdate)))  # slightly less than two time intervals  # min(diff(xdate))
    cat("ChainsFromGeneral: Using dtime",dtime,cr)

    DONOTICEF <- TRUE
    DONOTICEB <- TRUE

    # check if there are any NA in data
    iNA <- which(is.na(x1))
    if(length(iNA)>0){
        cat("ChainsFromGeneral: NA data found in x1.  Setting to 0",cr)
        x1[iNA] <- 0
    }
    if(DOPLOT) {
        plot(xdate, x1,pch=46)
        abline(h=th.x1)
        points(seq_along(x1)[iNA], rep(-2,length(iNA),col=2),pch=3,cex=.3,col=2)
    }
    
    ichainsYr   <- list()
    all_chains  <- list()
    all_ichains <- list()
    all_absichains <- NULL
    all_count   <- 0
    xdata.i     <- seq_along(x1)

    ichains    <- list()
    s1         <- x1
    iremaining <- which(s1>=th.x1)
    count      <- 0

    while(length(iremaining) > 0) {
        count <- count+1
        i0    <- which.max(s1[iremaining])
        xi    <- iremaining[i0]
        fi    <- xi
        bi    <- xi

        GOFWD <- TRUE
        while ((fi-xi)<chain.length & GOFWD & fi<length(s1)) {
            if(s1[fi]>=th.x1 & (xdate[fi+1]-xdate[fi]<=dtime) ) fi <- fi + 1 else GOFWD <- FALSE
        }
        if( ((fi-xi)>chain.length) & DONOTICEF ) {cat("CAUTION: FWD chain length not long enough",cr); DONOTICEF=FALSE}

        GOBWD <- TRUE
        while ((xi-bi)<chain.length & GOBWD & bi>1) {
            if(s1[bi]>=th.x1 & (xdate[bi]-xdate[bi-1]<=dtime) ) bi <- bi - 1 else GOBWD <- FALSE
        }
        if( ((xi-bi)>chain.length) & DONOTICEB ) {cat("CAUTION: BWD chain length not long enough",cr); DONOTICEB=FALSE}
        # NB this produces indecies that are one before and one after the indecies that are above th.x1

        if(bi<1)          b_buffers <- TRUE
        if(fi>length(s1)) f_buffers <- TRUE

        ichains[[count]] <- c(bi,xi,fi)
        i2               <- c(max(c(1,bi)),xi,min(c(fi,length(s1))) ) # make sure indecies are within range of dim(x1)[1]

        allChainLength <- 2*chain.length+2
        achain                                                                <- rep(0,length=allChainLength)
        achain[ (chain.length-(i2[2]-i2[1])) : (chain.length+(i2[3]-i2[2])) ] <- x1[ i2[1]:i2[3] ]
        absi                                                                  <- xdata.i[ i2[1]  ]

        all_count                 <- all_count +1
        all_chains[[all_count]]   <- achain
        all_absichains[all_count] <- absi
        all_ichains[[all_count]]  <- ichains[[count]]

        # s1[ i2[1]:i2[3] ]  <- th.x1 - 0.1*abs(th.x1)
        s1[ i2[1]:i2[3] ]  <- th.x1 - 0.1*abs(min(x1))

        iremaining <- which(s1>=th.x1)

        if(DOPLOT) {
               points(ichains[[count]][2],                    x1[ichains[[count]][2]],                    col=count,pch=3,cex=.3)
                lines(ichains[[count]][1]:ichains[[count]][3],x1[ichains[[count]][1]:ichains[[count]][3]],col=count)
               points(ichains[[count]][1]:ichains[[count]][3],x1[ichains[[count]][1]:ichains[[count]][3]],col=count,pch=3,cex=.5)
            # cat('Count', count,'chain',ichains[[count]],cr)
            # cat("Points remaining",length(iremaining),cr)
            # browser()
            if(DOPAUSE) readline("Continue?")
        }

        if(all_count%%1000==0)cat("ChainsFromGeneral: HW count and points remaining:",all_count,length(iremaining),cr)

    } #end while

    ichainsYr[[count]] <- ichains

    all_chains2      <- unlist(all_chains)
    dim(all_chains2) <- c(allChainLength, length(all_chains2)/(allChainLength))
    all_ichains2     <- unlist(all_ichains)
    dim(all_ichains2)      <- c(3, length(all_ichains2)/(3))
    dimnames(all_ichains2) <- list(c('day-before','peak-day','day-after'))
    names(all_absichains)  <- 'Absolute index of day before heatwave of -summer- data'

    return(list(chains=all_chains2, absichains=all_absichains, ichains=all_ichains2, indeciesByYear=ichainsYr, thresh=th.x1, chain.length=chain.length))

}

### 2023 code
##############################################################################################
###  HotDays functions
##############################################################################################

##############################################################################
fn_findHotSeason <- function(data01, thHotSeason, st.qgam=NULL, st.pre=NULL, st.pdf=NULL, DOPLOT=TRUE){

    if (!exists("om_data")) {
        if(!file.exists(st.pre)) {
            cat("fn_findHotSeason: ERROR",cr,"Must either supply om_data object or path to file",cr)
            return(NULL)
        } else {
            load(st.pre)
        }
    }

    if (!exists("fit.qgam")) {
        if(!file.exists(st.pre)) {
            cat("fn_findHotSeason: ERROR",cr,"Must either supply fit.qgam object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_findHotSeason: Loading qgam",cr,st.qgam,cr)
            load(st.qgam)
        }
    }

    cl0 <- unique(data01$class)

    q50o <- qdo(fit.qgam, 0.5, predict, newdata=data.frame(gmst=0, sdoy=1:data01.std.param$doy$o_max/data01.std.param$doy$o_max, class=cl0[1]))
    q50m <- qdo(fit.qgam, 0.5, predict, newdata=data.frame(gmst=0, sdoy=1:data01.std.param$doy$m_max/data01.std.param$doy$m_max, class=cl0[2]))

    annCyc.o     <- abs(diff(range(q50o)))
    annCyc.m     <- abs(diff(range(q50m)))
    if( annCyc.o > thHotSeason | annCyc.m > thHotSeason) { # there is a hot season
        IsHotSeason <- TRUE

        ix      <- which.max(q50o)
        th1     <- q50o[ix]
        while(length(which(q50o>th1)) < ceiling(data01.std.param$doy$o_max/2)) th1 <- th1 -0.1
        hseas.o <- which(q50o>th1)

        ix      <- which.max(q50m)
        th1     <- q50m[ix]
        while(length(which(q50m>th1)) < ceiling(data01.std.param$doy$m_max/2)) th1 <- th1 -0.1
        hseas.m <- which(q50m>th1)

        cat("Hot season found. Range:",annCyc.o, annCyc.m,cr)

    } else {

        hseas.o     <- c( 1 : data01.std.param$doy$o_max )
        hseas.m     <- c( 1 : data01.std.param$doy$m_max )

        cat("No hot season found. Range:",annCyc.o, annCyc.m,cr)
    }

    hotSdoy <- list(o=hseas.o, m=hseas.m)

    if(DOPLOT) {
        if(!is.null(st.pdf)) pdf(file=st.pdf, width=7, height=7) else x11()
        up.1()
        plot(  q50o,main='Annual cycle & hot season', ylim=range(c(q50o,q50m)))
        points(hotSdoy$o, q50o[hotSdoy$o],cex=.3,pch=3,col=2)
        points(q50m,                                   col=3)
        points(hotSdoy$m, q50m[hotSdoy$m],cex=.3,pch=3,col=4)
        if(!is.null(st.pdf)) dev.off()
    }

    ihotseas <- list(o=which(om_data$obs$info$doy %in% hotSdoy$o), m=which(om_data$mod$info$doy %in% hotSdoy$m))
    ihotseas <- list(o=which(om_data$obs$info$doy %in% hotSdoy$o), m=which(om_data$mod$info$doy %in% hotSdoy$m))
    hotseas <- list(IsHotSeason=IsHotSeason, hotSdoy=hotSdoy, ihotseas=ihotseas)

    return(hotseas)

}

##############################################################################
fn_findHotSeasonSingle <- function(data01, thHotSeason, facHotSeas=0.5, st.qgam=NULL, st.pre=NULL, st.pdf=NULL, DOPLOT=TRUE){
    
    # facHotSeas is the fraction of the year that has to be in the hot season. 
    # E.g. 0.5 = half the year, 0.25 = quarter of the year eg ~JJA

    if (!exists("fit.qgam")) {
        if(!file.exists(st.pre)) {
            cat("fn_findHotSeason: ERROR",cr,"Must either supply fit.qgam object or path to file",cr)
            return(NULL)
        } else {
            cat("fn_findHotSeason: Loading qgam",cr,st.qgam,cr)
            load(st.qgam)
        }
    }

    iy1        <- which(trunc(data01$time)==sort(unique(trunc(data01$time)))[2])
    ndata      <- data01[iy1,]
    ndata$gmst <- 0
    q50        <- qdo(fit.qgam, 0.5, predict, newdata=ndata)

    annCyc     <- abs(diff(range(q50)))
    if( annCyc > thHotSeason) { # there is a hot season
        IsHotSeason <- TRUE

        ix      <- which.max(q50)
        th1     <- q50[ix]
        while(length(which(q50>th1)) < ceiling(max(data01$doy)*facHotSeas)) th1 <- th1 -0.1
        hseas <- which(q50>th1)

        cat("Hot season found. Range:",annCyc,cr)

    } else {
        hseas     <- c( 1 : max(data01$doy) )
        cat("No hot season found. Range:",annCyc,cr)
    }

    hotSdoy <- list(o=hseas)

    if(DOPLOT) {
        if(!is.null(st.pdf)) {
            system(paste("mkdir -p",dirname(st.pdf) ) )
            pdf(file=st.pdf, width=7, height=7) 
        } else x11()
        up.1()
        plot(  q50,main='Annual cycle & hot season', ylim=range(c(q50)))
        points(hotSdoy$o, q50[hotSdoy$o],cex=.3,pch=3,col=2)
        if(!is.null(st.pdf)) dev.off()
    }

    ihotseas <- list(o=which(data01$doy %in% hotSdoy$o))
    hotseas <- list(IsHotSeason=IsHotSeason, hotSdoy=hotSdoy, ihotseas=ihotseas)

    return(hotseas)

}

##############################################################################
fn_extractEvents <- function(data01, ch.lev.u, event.length) {

    ### NB NB NB extracted chains/events now have the day before and day after included
    cat("fn_extractEvents: NB NB NB extracted chains/events now have the day before and day after included",cr)

    ## fix points where makestationary has produced NAs
    ina <- which(is.na(data01$u))
    if(length(ina)>0) data01$u[ina] <- 0.1

    x.l      <- qlaplace(data01$u)
    ch.lev.l <- qlaplace(ch.lev.u)
    # get events wrt peak values
    ch.pk <- ChainsFromLaplace(x.l, data01$time,  chain.lev.u=ch.lev.u, chain.length=event.length, DOPLOT=TRUE)

    ### make wrt start
    # ChainsFromLaplace returns days above threshold.  Add on day before and day after
    ch.idx  <- apply(ch.pk$chains.l,2,FUN=function(a,b) {return(which(a>=b))} , ch.lev.l)
    max.len <- max(sapply(ch.idx,length))
    ch.st.l <- array(NA,dim=c(max.len+2, length(ch.pk$absichains)))
    ch.doy  <- array(NA,dim=dim(ch.st.l))
    ch.age  <- array(NA,dim=dim(ch.st.l))
    ch.i    <- array(NA,dim=dim(ch.st.l))
    ch.st.o <- array(NA,dim=dim(ch.st.l))
    ch.time <- list()
    ch.gmst <- double(length(ch.st.l[1,]))

    i_fail <- NULL
    for(i in 1:dim(ch.st.l)[2]) {
        #ichange are now absolute # idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]  -ch.pk$ichains[1,i] +ch.pk$absichains[i]
        idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]   # day before to day after

        if(length(idays) >= 3){
            # ch.st.l[1:(j2-j1+1)  , i] <- x.l[j1:j2]
               ch.i[1:(length(idays)), i] <- idays
            ch.st.l[1:(length(idays)), i] <-        x.l[idays]
            ch.st.o[1:(length(idays)), i] <-   data01$x[idays]
             ch.doy[1:(length(idays)), i] <- data01$doy[idays]
             ch.age[1:(length(idays)), i] <- c(0,1:(length(idays)-2),0) # make day before and after =0
            ch.time[[i]]                  <- data01$time[idays]
             ch.gmst[i]                   <- data01$gmst[idays[2]]
        } else i_fail <- c(i_fail, i)
    }

    # remove fails
    if(length(i_fail)>0) {
        ch.pk$chains.l       <- ch.pk$chains.l[, -i_fail]
        ch.pk$absichains     <- ch.pk$absichains[-i_fail]
        ch.pk$ichains.l      <- ch.pk$ichains.l[, -i_fail]
        ch.pk$indeciesByYear <- ch.pk$indeciesByYear[-i_fail]

        ch.i    <-    ch.i[, -i_fail]
        ch.st.l <- ch.st.l[, -i_fail]
        ch.st.o <- ch.st.o[, -i_fail]
        ch.doy  <-  ch.doy[, -i_fail]
        ch.age  <-  ch.age[, -i_fail]
        ch.time <- ch.time[  -i_fail]
        ch.gmst <- ch.gmst[  -i_fail]
    }

    # calc summary metrics
    ch.sev       <- double(length(ch.st.l[1,]))
    ch.mean.sev  <- double(length(ch.st.l[1,]))
    ch.pkval     <- double(length(ch.st.l[1,]))
    ch.pkval.day <- double(length(ch.st.l[1,]))
    ch.duration  <- double(length(ch.st.l[1,]))
    for(i in seq_along(ch.st.l[1,])) {
        j0              <- which(ch.st.l[  ,i] >= ch.lev.l)
        ch.sev[i]       <-   sum(ch.st.l[j0,i]-ch.lev.l,na.rm=T)
        ch.mean.sev[i]  <-  mean(ch.st.l[j0,i]-ch.lev.l,na.rm=T)
        ch.duration[i]  <- length(j0)
        ch.pkval.day[i] <- which.max(ch.st.l[j0,i])
        ch.pkval[i]     <- ch.st.l[j0[ch.pkval.day[i]],i]
        if(i%%trunc(length(ch.st.l[1,])/10)==0) cat("Obs severity Completed",i,'of',length(ch.st.l[1,]),cr)
    }

    cat("fn_extractEvents: No of chains found",length(ch.pk$absichains),cr)
    x_chains <- list(ch.pk=ch.pk, ch.st.l=ch.st.l, ch.st.o=ch.st.o,
                            ch.doy=ch.doy, ch.age=ch.age, ch.i=ch.i, ch.time=ch.time,
                            ch.sev=ch.sev, ch.mean.sev=ch.mean.sev, ch.pkval=ch.pkval,
                            ch.pkval.day=ch.pkval.day, ch.duration=ch.duration, ch.gmst=ch.gmst, ch.th.u=ch.lev.u)

    # x_chains$reference           <- makeStat$reference
    # x_chains$reference$hw.lev.u  <- hw.lev.u
    # x_chains$reference$hw.lev.u3 <- hw.lev.u3
    # x_chains$reference$u3.par    <- u3.par

    return(x_chains)

}

##############################################################################
fn_extractEventsGeneric <- function(x1, th.x1, event.length) {

    ### input any timeseries and a threshold.  Output chains/events

    ### NB NB NB extracted chains/events now have the day before and day after included
    cat("fn_extractEventsGeneric: NB NB NB extracted chains/events now have the day before and day after included",cr)

    ## fix points where makestationary has produced NAs
    ina <- which(is.na(x1))
    if(length(ina)>0) x1[ina] <- min(x1,na.rm=TRUE)

    # get events wrt peak values
    ch.pk <- ChainsFromGeneric(x1, th.x1, chain.length=event.length, DOPLOT=TRUE, DOPAUSE=FALSE)

    ### make wrt start
    # ChainsFromLaplace returns days above threshold.  Add on day before and day after
    ch.idx  <- apply(ch.pk$chains,2,FUN=function(a,b) {return(which(a>=b))} , th.x1)
    max.len <- max(sapply(ch.idx,length))
    ch.st.o <- array(NA,dim=c(max.len+2, length(ch.pk$absichains)))
    ch.age  <- array(NA,dim=dim(ch.st.o))
    ch.i    <- array(NA,dim=dim(ch.st.o))
    ch.time <- list()

    i_fail <- NULL
    for(i in 1:dim(ch.st.o)[2]) {
        #ichange are now absolute # idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]  -ch.pk$ichains[1,i] +ch.pk$absichains[i]
        idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]   # day before to day after

        if(length(idays) >= 3){
            # ch.st.o[1:(j2-j1+1)  , i] <- x.l[j1:j2]
               ch.i[1:(length(idays)), i] <- idays
            ch.st.o[1:(length(idays)), i] <-        x1[idays]
             ch.age[1:(length(idays)), i] <- c(0,1:(length(idays)-2),0) # make day before and after =0
            ch.time[[i]]                  <- data01$time[idays]
        } else i_fail <- c(i_fail, i)
    }

    # remove fails
    if(length(i_fail)>0) {
        ch.pk$chains         <- ch.pk$chains[, -i_fail]
        ch.pk$absichains     <- ch.pk$absichains[-i_fail]
        ch.pk$ichains.l      <- ch.pk$ichains.l[, -i_fail]
        ch.pk$indeciesByYear <- ch.pk$indeciesByYear[-i_fail]

        ch.i    <-    ch.i[, -i_fail]
        ch.st.o <- ch.st.o[, -i_fail]
        ch.age  <-  ch.age[, -i_fail]
        ch.time <- ch.time[  -i_fail]
    }

    # calc summary metrics
    ch.sev       <- double(length(ch.st.o[1,]))
    ch.mean.sev  <- double(length(ch.st.o[1,]))
    ch.pkval     <- double(length(ch.st.o[1,]))
    ch.pkval.day <- double(length(ch.st.o[1,]))
    ch.duration  <- double(length(ch.st.o[1,]))
    for(i in seq_along(ch.st.o[1,])) {
        j0              <- which(ch.st.o[  ,i] >= th.x1)
        ch.sev[i]       <-   sum(ch.st.o[j0,i]   -th.x1,na.rm=T)
        ch.mean.sev[i]  <-  mean(ch.st.o[j0,i]   -th.x1,na.rm=T)
        ch.duration[i]  <- length(j0)
        ch.pkval.day[i] <- which.max(ch.st.o[j0,i])
        ch.pkval[i]     <- ch.st.o[j0[ch.pkval.day[i]],i]
        if(i%%trunc(length(ch.st.o[1,])/10)==0) cat("Obs severity Completed",i,'of',length(ch.st.o[1,]),cr)
    }

    cat("fn_extractEvents: No of chains found",length(ch.pk$absichains),cr)
    x_chains <- list(ch.pk=ch.pk, ch.st.o=ch.st.o, 
                            ch.age=ch.age, ch.i=ch.i, ch.time=ch.time,
                            ch.sev=ch.sev, ch.mean.sev=ch.mean.sev, ch.pkval=ch.pkval,
                            ch.pkval.day=ch.pkval.day, ch.duration=ch.duration, th.x1=th.x1)

    return(x_chains)

}

##############################################################################
fn_extractEventsAbs <- function(simevents, th.abs) {

    ### expects output from main_HDsim

    cat("fn_extractEventsAbs: NB NB NB extracting chains/events with respect to a fixed absolute threshold",cr)
    cat("fn_extractEventsAbs: NB NB NB extracted chains/events now have the day before and day after included",cr)

    ## fix points where NAs
    ina <- which(is.na(data01$u))
    if(length(ina)>0) data01$u[ina] <- 0.1
    x.l      <- qlaplace(data01$u)
    ch.lev.l <- qlaplace(ch.lev.u)
    # get events wrt peak values
    ch.pk <- ChainsFromLaplace(x.l, data01$time,  chain.lev.u=ch.lev.u, chain.length=event.length, DOPLOT=FALSE)

    ### make wrt start
    # ChainsFromLaplace returns days above threshold.  Add on day before and day after
    ch.idx  <- apply(ch.pk$chains.l,2,FUN=function(a,b) {return(which(a>=b))} , ch.lev.l)
    max.len <- max(sapply(ch.idx,length))
    ch.st.l <- array(NA,dim=c(max.len+2, length(ch.pk$absichains)))
    ch.doy  <- array(NA,dim=dim(ch.st.l))
    ch.age  <- array(NA,dim=dim(ch.st.l))
    ch.i    <- array(NA,dim=dim(ch.st.l))
    ch.st.o <- array(NA,dim=dim(ch.st.l))
    ch.time <- list()
    ch.gmst <- double(length(ch.st.l[1,]))

    i_fail <- NULL
    for(i in 1:dim(ch.st.l)[2]) {
        #ichange are now absolute # idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]  -ch.pk$ichains[1,i] +ch.pk$absichains[i]
        idays <- ch.pk$ichains[1,i]:ch.pk$ichains[3,i]   # day before to day after

        if(length(idays) >= 3){
            # ch.st.l[1:(j2-j1+1)  , i] <- x.l[j1:j2]
               ch.i[1:(length(idays)), i] <- idays
            ch.st.l[1:(length(idays)), i] <-        x.l[idays]
            ch.st.o[1:(length(idays)), i] <-   data01$x[idays]
             ch.doy[1:(length(idays)), i] <- data01$doy[idays]
             ch.age[1:(length(idays)), i] <- c(0,1:(length(idays)-2),0) # make day before and after =0
            ch.time[[i]]                  <- data01$time[idays]
             ch.gmst[i]                   <- data01$gmst[idays[2]]
        } else i_fail <- c(i_fail, i)
    }

    # remove fails
    if(length(i_fail)>0) {
        ch.pk$chains.l       <- ch.pk$chains.l[, -i_fail]
        ch.pk$absichains     <- ch.pk$absichains[-i_fail]
        ch.pk$ichains.l      <- ch.pk$ichains.l[, -i_fail]
        ch.pk$indeciesByYear <- ch.pk$indeciesByYear[-i_fail]

        ch.i    <-    ch.i[, -i_fail]
        ch.st.l <- ch.st.l[, -i_fail]
        ch.st.o <- ch.st.o[, -i_fail]
        ch.doy  <-  ch.doy[, -i_fail]
        ch.age  <-  ch.age[, -i_fail]
        ch.time <- ch.time[  -i_fail]
        ch.gmst <- ch.gmst[  -i_fail]
    }

    # calc summary metrics
    ch.sev       <- double(length(ch.st.l[1,]))
    ch.mean.sev  <- double(length(ch.st.l[1,]))
    ch.pkval     <- double(length(ch.st.l[1,]))
    ch.pkval.day <- double(length(ch.st.l[1,]))
    ch.duration  <- double(length(ch.st.l[1,]))
    for(i in seq_along(ch.st.l[1,])) {
        j0              <- which(ch.st.l[  ,i] >= ch.lev.l)
        ch.sev[i]       <-   sum(ch.st.l[j0,i]-ch.lev.l,na.rm=T)
        ch.mean.sev[i]  <-  mean(ch.st.l[j0,i]-ch.lev.l,na.rm=T)
        ch.duration[i]  <- length(j0)
        ch.pkval.day[i] <- which.max(ch.st.l[j0,i])
        ch.pkval[i]     <- ch.st.l[j0[ch.pkval.day[i]],i]
        if(i%%trunc(length(ch.st.l[1,])/10)==0) cat("Obs severity Completed",i,'of',length(ch.st.l[1,]),cr)
    }

    cat("fn_extractEvents: No of chains found",length(ch.pk$absichains),cr)
    x_chains <- list(ch.pk=ch.pk, ch.st.l=ch.st.l, ch.st.o=ch.st.o,
                            ch.doy=ch.doy, ch.age=ch.age, ch.i=ch.i, ch.time=ch.time,
                            ch.sev=ch.sev, ch.mean.sev=ch.mean.sev, ch.pkval=ch.pkval,
                            ch.pkval.day=ch.pkval.day, ch.duration=ch.duration, ch.gmst=ch.gmst, ch.th.u=ch.lev.u)

    # x_chains$reference           <- makeStat$reference
    # x_chains$reference$hw.lev.u  <- hw.lev.u
    # x_chains$reference$hw.lev.u3 <- hw.lev.u3
    # x_chains$reference$u3.par    <- u3.par

    return(x_chains)

}

##############################################################################
fn_plotEvents <- function(events, data, nplot=3, savefile=NULL){
    # stpdf_ch <- paste(dirname(stout_ch),'/plots/',sub('.RData','.pdf',basename(stout_ch)),sep='')
    if(!is.null(savefile)) {
        pdf(file=savefile,height=7,width=7)
        DOPAUSE <- FALSE
    } else DOPAUSE <- TRUE

    ch.lev.l <- qlaplace(events$ch.th.u)

    c21      <- distinct21colours()
    isortage <- sort(events$ch.duration,dec=TRUE,index=TRUE)$ix
    for( j in isortage[1:nplot]){
        up.1()
        y0  <- unique(unlist(events$ch.time[[j]]))
        i01 <- which(trunc(data$time)==trunc(y0)[1])
        ylim0 <- c(-0.5, max(qlaplace(data$u[i01]),na.rm=TRUE))
        plot(data$time[i01],qlaplace(data$u[i01]), pch=3, cex=.3, main=trunc(y0)[1], ylim=ylim0)
        y1  <- unlist(lapply(events$ch.time, "[", 1))
        i1  <- which(trunc(y1)==trunc(y0)[1])
        c2 <- sample(c21, length(i1), rep=FALSE)
        for(k in seq_along(i1)) {
            ch.t <- events$ch.time[[i1[k]]]
            points(ch.t, events$ch.st.l[,i1[k]][seq_along(ch.t)],col=c2[k])
             lines(ch.t, events$ch.st.l[,i1[k]][seq_along(ch.t)],col=c2[k])
        }
        abline(h=ch.lev.l, lty=2, col='grey70')
        if(DOPAUSE) readline("Continue?")
    }
    if(!is.null(savefile))  dev.off()


}

### Initialisation
##############################################################################
fn_fitInitEvent <- function(events01, data01, initI, IsJoint=TRUE, SAVEINITEVENTDIAG=TRUE){

    # events01  : extraced events for both model and obs
        #  $ ch.pk       : Named list()
        #  $ ch.st.l     : num [1:19, 1:333] 1.96 2.54 2.49 5.76 2.24 ...
        #  $ ch.st.o     : num [1:19, 1:333] 29.8 31 30.8 35.7 30.2 ...
        #  $ ch.doy      : int [1:19, 1:333] 215 216 217 218 219 220 221 222 223 224 ...
        #  $ ch.age      : num [1:19, 1:333] 0 1 2 3 4 5 6 7 8 9 ...
        #  $ ch.i        :     [1:19, 1:333] 15921 15922 15923 ... index wrt the full year int
        #  $ ch.time     :List of 333
        #  $ ch.sev      : num [1:333] 18.27 8.01 41.57 9.93 7.02 ...
        #  $ ch.mean.sev : num [1:333] 2.03 4 2.45 1.99 2.34 ...
        #  $ ch.pkval    : num [1:333] 8.42 7.3 7.12 6.67 6.63 ...
        #  $ ch.pkval.day: num [1:333] 7 2 5 2 3 2 4 1 3 4 ...
        #  $ ch.duration : num [1:333] 9 2 17 5 3 3 6 3 4 4 ...
        #  $ ch.th.u     : num 0.94

    # data01    : the original (all) obs + model data that events were extracted from
        # 'data.frame':	59011 obs. of  10 variables:
        #  $ x    : num  11.58 9.05 11.06 12.44 9.1 ...
        #  $ time : num  1960 1960 1960 1960 1960 ...
        #  $ stime: num  -0.5 -0.5 -0.5 -0.5 -0.5 ...
        #  $ doy  : int  1 2 3 4 5 6 7 8 9 10 ...
        #  $ sdoy : num  0.00273 0.00546 0.0082 0.01093 0.01366 ...
        #  $ isobs: num  1 1 1 1 1 1 1 1 1 1 ...
        #  $ class: Factor w/ 2 levels "obs","mod": 1 1 1 1 1 1 1 1 1 1 ...
        #  $ gmst : num  -0.341 -0.341 -0.341 -0.341 -0.341 ...
        #  $ uqgam: num  0.89 0.655 0.852 0.943 0.665 ...
        #  $ u    : num  0.89 0.655 0.852 0.945 0.665 ...

    # perhaps too much is done in this subroutine and it should be slimmed down

    covariate <- initI$ievent.covaraite
    ievent.k  <- initI$ievent.k

    InitModel <- list()

    if(IsJoint) {

       ### model
        evs                 <- events01$mod
        imod                <- which(data01$isobs==0)
        init2               <- rep(0, length(imod)) # make an array for all data
        init2[evs$ch.i[2,]] <- 1  # 2=first day of heatwave, absolute index for the original data
        nob                 <- length(which(data01$isobs==1))
        if(evs$hotseas$IsHotSeason) {
            form.doy  <- 'cs'
            init2     <- init2[evs$hotseas$ihotseas]
            evinit.m  <- data.frame(init=init2, sdoy=data01$sdoy[nob+evs$hotseas$ihotseas], cc=data01[nob+evs$hotseas$ihotseas,covariate], class='mod' )
        } else {
            form.doy  <- 'cc' # ensure smooth is cyclic 'cc' if whole year is being modelled
            evinit.m  <- data.frame(init=init2, sdoy=data01$sdoy[imod], cc=data01[imod,covariate], class='mod' )
        }
            # plot(evinit.m$init,pch=3,cex=.3)
            # plot(evinit.m$sdoy,pch=3,cex=.3)
            # plot(evinit.m$cc,pch=3,cex=.3)
            # plot(evinit.m$sdoy,evinit.m$init,pch=3,cex=.3)
            # plot(evinit.m$cc,  evinit.m$init,pch=3,cex=.3)
            # i1<-which(evinit.m$init==1)
            # hist(evinit.m$sdoy[i1],20)
            # hist(evinit.m$cc[i1],20)
        fit.init.D  <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=evinit.m, family="binomial")
        if(diff(range(evinit.m$cc))>0.1) {  # control runs have constant gmst
            fit.init.T  <- gam(init ~                                     s(cc,bs='cr',k=ievent.k$cc), data=evinit.m, family="binomial")
            fit.init.DT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='cr',k=ievent.k$cc), data=evinit.m, family="binomial")
        } else {
            fit.init.T  <- NULL
            fit.init.DT <- NULL
        }
        # predict.gam(fit.init.D, newdata=data.frame(sdoy=0.5),type="response")
        InitModel$mod <- list(D=fit.init.D, T=fit.init.T, DT=fit.init.DT, ch.th.u=evs$ch.th.u, covariate=covariate)

       ### obs
        evs                 <- events01$obs
        iobs                <- which(data01$isobs==1)
        init2               <- rep(0, length(iobs)) # make an array for all data
        init2[evs$ch.i[2,]] <- 1  # 2=first day of heatwave, absolute index for the original data
        if(evs$hotseas$IsHotSeason) {
            form.doy  <- 'cs'
            init2     <- init2[evs$hotseas$ihotseas]
            evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[evs$hotseas$ihotseas], cc=data01[evs$hotseas$ihotseas,covariate], class='obs' )
        } else {
            form.doy  <- 'cc' # ensure smooth is cyclic 'cc' if whole year is being modelled
            evinit.o  <- data.frame(init=init2, sdoy=data01$sdoy[iobs], cc=data01[iobs,covariate], class='obs' )
        }
            # plot(evinit.o$init,pch=3,cex=.3)
            # plot(evinit.o$sdoy,pch=3,cex=.3)
            # plot(evinit.o$cc,pch=3,cex=.3)
            # plot(evinit.o$sdoy,evinit.o$init,pch=3,cex=.3)
            # plot(evinit.o$cc,evinit.o$init,pch=3,cex=.3)
            # i1<-which(evinit.o$init==1)
            # hist(evinit.o$sdoy[i1],20)
            # hist(evinit.o$cc[i1],20)
        fit.init.D  <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=evinit.o, family="binomial")
        fit.init.T  <- gam(init ~                                     s(cc,bs='cr',k=ievent.k$cc), data=evinit.o, family="binomial")
        fit.init.DT <- gam(init ~ s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='cr',k=ievent.k$cc), data=evinit.o, family="binomial")
        # predict.gam(fit.init.D, newdata=data.frame(sdoy=0.5),type="response")
        InitModel$obs <- list(D=fit.init.D, T=fit.init.T, DT=fit.init.DT, ch.th.u=evs$ch.th.u, covariate=covariate)

       ### both obs and model
        # inherits the hotseason definition of the obs, which is good
        om <- rbind(evinit.o,evinit.m)
        # plot(om$cc,  pch=3,cex=.3)
        fit.init.B   <- gam(init ~ class                                                                 , data=om, family="binomial")
        fit.init.D   <- gam(init ~        s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=om, family="binomial")
        fit.init.T   <- gam(init ~                                            s(cc,bs='cr',k=ievent.k$cc), data=om, family="binomial")
        fit.init.BD  <- gam(init ~ class +s(sdoy,bs=form.doy,k=ievent.k$doy)                             , data=om, family="binomial")
        fit.init.BT  <- gam(init ~ class                                     +s(cc,bs='cr',k=ievent.k$cc), data=om, family="binomial")
        fit.init.DT  <- gam(init ~        s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='cr',k=ievent.k$cc), data=om, family="binomial")
        fit.init.BDT <- gam(init ~ class +s(sdoy,bs=form.doy,k=ievent.k$doy) +s(cc,bs='cr',k=ievent.k$cc), data=om, family="binomial")
        InitModel$joint <- list(B=fit.init.B, D=fit.init.D, T=fit.init.T,
                                BD=fit.init.BD, BT=fit.init.BT, DT=fit.init.DT,
                                BDT=fit.init.BDT, ch.th.u=evs$ch.th.u, covariate=covariate)

    } else {

        cat("fn_fitInitEvent: Non-joint fitting not implemented yet",cr)

    }

    if(SAVEINITEVENTDIAG) {
        stout_InitDiag     <- paste(dirname(st_InitEvent),'diag',sub('.RData','.diag.txt',basename(st_InitEvent)),sep='/')
        conDiag    <- file(stout_InitDiag,'w')
        writeLines(stout_InitDiag, con=conDiag)
        writeLines(st_InitEvent,   con=conDiag)
        writeLines('\n',   con=conDiag)

        st1 <- c('obs','mod')
        st2 <- c('D','T','DT')
        for(i in st1) for (j in st2) {
            if(!is.null(InitModel[[i]][[j]])) {
            writeLines("### summary.gam(fit) #################################",con=conDiag)
            writeLines(paste(i,j),con=conDiag)
            writeLines(capture.output(summary.gam(InitModel[[i]][[j]])),con=conDiag)
            writeLines("######################################################\n\n",con=conDiag)
            }
        }
        st1 <- c('joint')
        st2 <- c('B','D','T','BDT')
        for(i in st1) for (j in st2) {
            if(!is.null(InitModel[[i]][[j]])) {
            writeLines("### summary.gam(fit) #################################",con=conDiag)
            writeLines(paste(i,j),con=conDiag)
            writeLines(capture.output(summary.gam(InitModel[[i]][[j]])),con=conDiag)
            writeLines("######################################################\n\n",con=conDiag)
            }
        }
        close(conDiag)
    }

    return(list(InitModel=InitModel, initI=initI))
}

##############################################################################
fn_plotInitEvent <- function(InitModel, ylim=c(0,0.1), main0='Init ~', st.pdf='plotInitEvent.pdf'){

    # cccov  <- initI$ievent.covaraite
    # cov.n  <- attributes(hw.term$terms)$term.labels
    # r.sdoy <- range(InitModel$model$sdoy)
    # newo.d <- data.frame(sdoy=seq(r.sdoy[1],   r.sdoy[2],   length=100),
    #                      cc=median(data01[,covar]),
    #                      class='obs')
    # newo.c <- data.frame(sdoy=median(data01$sdoy),
    #                      cc=seq(min(data01[,covar]),max(data01[,covar]),length=100),
    #                      class='obs')
    # newm.d       <- newo.d
    # newm.d$class <- 'mod'
    # newm.c       <- newo.c
    # newm.c$class <- 'mod'

    # pro.d <- predict.gam(InitModel, newdata=newo.d,type="response")
    # prm.d <- predict.gam(InitModel, newdata=newm.d,type="response")
    # pro.c <- predict.gam(InitModel, newdata=newo.c,type="response")
    # prm.c <- predict.gam(InitModel, newdata=newm.c,type="response")

    # up.2()
    #   plot(newo.d$sdoy, pro.d, ylim=ylim, ylab='Probability', main='HW Initiation prob wrt DOY')
    # points(newo.d$sdoy, prm.d, col=2, cex=.3)

    #   plot(newo.c$cc, pro.c, ylim=ylim, ylab='Probability', main=paste('HW Initiation prob wrt',covar))
    # points(newo.c$cc, prm.c, col=2, cex=.3)

    ##################################

    if(class(InitModel)[1]!="gam") {
        cat("fn_termJointGAMPlot: must be a GAM object",cr)
        return(NULL)
    }
        if(!is.null(st.pdf)) pdf(file=st.pdf) else x11()

    # hw.cov <- hw.term$data
    hw.cov <- InitModel$model
    cov.n  <- attributes(InitModel$terms)$term.labels
    if(any(grepl('class',cov.n))) {
        ISCLASS <- TRUE
        cov.n  <- cov.n[cov.n!='class']
    } else ISCLASS <- FALSE

    cov.q  <- NULL
    cov.x  <- NULL
    for(i in seq_along(cov.n)) {
        cov.q      <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
        cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
    }

    df1   <- NULL
    for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
    df1   <- data.frame(df1)
    df0   <- df1
    up.2()
    if(ISCLASS) {
        df1$class <- 'obs'
        df0$class <- 'mod'
        for(i in seq_along(cov.n)) {
            df0b             <- df0
            df0b[[cov.n[i]]] <- cov.x[[i]]
            df1b             <- df1
            df1b[[cov.n[i]]] <- cov.x[[i]]
            newy0            <- predict(InitModel,df0b,type="response")
            newy1            <- predict(InitModel,df1b,type="response")
            plot(cov.x[[i]], newy0, ylim=ylim, ty='l', lwd=2, main=paste(main0,cov.n[i]), ylab='Prob Initiation', xlab="Covariate")
            lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
            grid()
        }
        legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
    } else {
        df1$class <- 'obs'
        for(i in seq_along(cov.n)) {
            df1b             <- df1
            df1b[[cov.n[i]]] <- cov.x[[i]]
            newy1            <- predict(InitModel,df1b,type="response")
            plot(cov.x[[i]], newy1, ylim=ylim, ty='l', lwd=2, main=paste(main0,cov.n[i]), ylab='Prob Initiation', xlab="Covariate")
            grid()
        }
        legend('bottomleft',c('Obs'),col=1, pch=NA, lty=1, lwd=2, bty='n')
    }
    if(!is.null(st.pdf)) dev.off()
}


##############################################################################
fn_fitInitVal <- function(events01, gpd.th.u, scale.doy, fmla.IVgpd, IsJoint=TRUE){

    # gpd.th.u  <- initV$th.u
    # scale.doy <- list(obs=data01.std.param$doy$o_max, mod=data01.std.param$doy$m_max)

    gpd.th.l <- qlaplace(gpd.th.u) # here we make the assumption that the MakeStationary has worked and
                                   # chains are truly stationary

    if(IsJoint) {

        ### Fit the Gam-Gpd

 # excess ~ s(sdoy, bs = "cc", k = ms.k$doy, by = class) + class +
 #     s(gmst, bs = "ds", k = ms.k$gmst) + ti(sdoy, gmst, bs = c("cc", "tp"))

        ### prep data
        d1.o <- events01$obs$ch.st.l[2,] # chains start the day before heatwave
        d1.m <- events01$mod$ch.st.l[2,] # chains start the day before heatwave

        sdoy.o <- events01$obs$ch.doy[2,]/scale.doy$obs
        sdoy.m <- events01$mod$ch.doy[2,]/scale.doy$mod

        gmst.o <- events01$obs$ch.gmst
        gmst.m <- events01$mod$ch.gmst

        class.o <- rep('obs',length(d1.o))
        class.m <- rep('mod',length(d1.m))
        class.om <- factor(c( rep(0,length(d1.o)), rep(1,length(d1.m))),labels=c('obs','mod'))

        # evdata01        <- data.frame(x=c(d1.o,d1.m), sdoy=c(sdoy.o,sdoy.m), gmst=c(gmst.o,gmst.m), class=c(class.o,class.m), iv.k=1, doy=1, gmst=1)
        evdata01        <- data.frame(x=c(d1.o,d1.m), sdoy=c(sdoy.o,sdoy.m), gmst=c(gmst.o,gmst.m), class=class.om, iv.k=1, gmst=1)
        evdata01$excess <- evdata01$x - gpd.th.l
        iexcess         <- which(evdata01$excess >0)

        evdata01.ex     <- subset(evdata01, excess > 0)

        ### fit one model
        sink("/dev/null")
            fit.IVgpd <- evgam(fmla.IVgpd, evdata01.ex, family="gpd", trace=2)
        sink()

        #
            # ### fit all models and select chosen prior to simulating
            # fits.IVgpd        <- vector(mode='list',length=length(fmla.IVgpd))
            # names(fits.IVgpd) <- names(fmla.IVgpd)
            # sink("/dev/null")
            #     for(i in seq_along(fmla.IVgpd)) fits.IVgpd[[i]] <- evgam(fmla.IVgpd[[i]], evdata01.ex, family="gpd", trace=2)
            # sink()

            ### not sure this is still needed
            # for(i in seq_along(fmla.IVgpd)) {
            #     which.fix <- which(sapply(fits.IVgpd[[i]], class) == "gamlist")
            #     for (j in which.fix) fits.IVgpd[[i]][[j]]$Vp <- as.matrix(fits.IVgpd[[i]][[j]]$Vp)
            # }

            ### no longer trying to choose best fit here
            # fits.IVgpd.aic        <- lapply(fits.IVgpd, AIC)
            # names(fits.IVgpd.aic) <- names(fmla.IVgpd)
            # cat("AIC",cr)
            # for(i in seq_along(fmla.IVgpd)) cat(names(fmla.IVgpd)[i],fits.IVgpd.aic[[i]],cr)

            # for(i in seq_along(fmla.IVgpd)) {
            #     cat(cr,'################',cr,names(fmla.IVgpd)[i],cr)
            #     print(summary(fits.IVgpd[[i]]))
            # }

            # # I dont trust the AIC as it does not pick out the failed fits
            # cat("Best IV GPD fit",cr)
            # ix1 <- which.min(fits.IVgpd.aic)
            # cat(names(fmla.IVgpd)[ix1],fits.IVgpd.aic[[ix1]],cr)
            # chosen.IVgpd.name <- names(fmla.IVgpd)[ix1]
            # chosen.IVgpd      <- fits.IVgpd[[ix1]]

            # pval7 <- summary(fits.IVgpd[[7]])[[2]]$logscale$'Pr(>|t|)'
        #

    } else {

        cat("fn_fitInitVal: Non-joint fitting not implemented yet",cr)

    }
    fit.IVgpd$gpd.th.u <- gpd.th.u
    fit.IVgpd$fit.data <- evdata01
    return(fit.IVgpd)
}

##############################################################################
fn_diagInitVal <- function(fit.IVgpd, st.GpdInitDiag='InitValues.stats', SAVE=TRUE, NEWFILE=TRUE){

    if(SAVE) {
        if (NEWFILE) {
            conDiag    <- file(st.GpdInitDiag,'w')
            writeLines(st.GpdInitDiag,                   con=conDiag)
            writeLines("\n###################",               con=conDiag)
            writeLines(capture.output(print(fit.IVgpd$name)), con=conDiag)
            writeLines("summary(fit.IVgpd)",                  con=conDiag)
            writeLines(capture.output(summary(fit.IVgpd)),    con=conDiag)
            writeLines("###################\n\n",             con=conDiag)
        } else  {
            conDiag    <- file(st.GpdInitDiag,'a')
            writeLines("\n###################",               con=conDiag)
            writeLines(capture.output(print(fit.IVgpd$name)), con=conDiag)
            writeLines("summary(fit.IVgpd)",                  con=conDiag)
            writeLines(capture.output(summary(fit.IVgpd)),    con=conDiag)
            writeLines("###################\n\n",             con=conDiag)
        }

        close(conDiag)
    } else {
        summary(fit.IVgpd)
    }

}

##############################################################################
fn_plotInitVal <- function(fit.IVgpd, newdata, xplot=1, retp=NULL, st.pdf=NULL){

    if(!is.null(st.pdf)) pdf(st.pdf, width=7, height=7) else x11()
    if(is.null(dim(newdata))) newdata <- cbind(newdata)
    if(is.null(retp)) retp <- 1/runif(dim(newdata)[1])

    gpd.th.l <- qlaplace(fit.IVgpd$gpd.th.u)

    nd1        <- newdata
    nd1$sdoy[] <- median(nd1$sdoy)
    param1 <- predict(fit.IVgpd, nd1, type="response")
    q1     <- NULL
    for(i in seq_along(param1$scale)) q1[i]   <- gpdq(c(param1$scale[i], param1$shape[i]),gpd.th.l,sample(1/retp,1))

    nd2        <- newdata
    nd2$gmst[] <- median(nd2$gmst)
    param2 <- predict(fit.IVgpd, nd2, type="response")
    q2     <- NULL
    for(i in seq_along(param2$scale)) q2[i]   <- gpdq(c(param2$scale[i], param2$shape[i]),gpd.th.l,sample(1/retp,1))

    q0 <- list(q1,q2)
    n0 <- list(nd2,nd1) # bodge to plot against the right x
    if(length(fit.IVgpd$gotsmooth)>0) {
        # plot smooths
        plot(fit.IVgpd)
        mtext(fit.IVgpd$name,outer=TRUE,at=c(0.02),adj=0,col=1,line=-1)
        if(is.null(st.pdf)) readline("Continue?")

        # plot values against covariates
        # param0 <- predict(fit.IVgpd, newdata, type="response")
        # gpd.th.l <- qlaplace(fit.IVgpd$gpd.th.u)
        # q3     <- NULL
        # for(i in seq_along(param0$scale)) q3[i]   <- gpdq(c(param0$scale[i], param0$shape[i]),gpd.th.l,1/qu)
        # for(j in seq_along(xplot)) {
        #     ylim0 <- c(min(q0[[j]]*0.95), max(q0[[j]]*1.05))
        #     plot(newdata[,xplot[j]], q0[[j]], pch=3,cex=.3,ylim=ylim0, main='HW initialisation temperature',xlab=names(newdata[xplot[j]]), ylab="quantile")
        #     grid()
        #     mtext(fit.IVgpd$name,outer=TRUE,at=c(0.02),adj=0,col=1,line=-1)
        # }

        nd1        <- newdata
        nd1$sdoy[] <- median(nd1$sdoy)
        param1 <- predict(fit.IVgpd, nd1, type="response")
        q1     <- NULL
        # for(i in seq_along(param1$scale)) q1[i]   <- gpdq(c(param1$scale[i], param1$shape[i]),gpd.th.l,sample(1/retp,1))
        for(i in seq_along(param1$scale)) q1[i]   <- gpdq(c(param1$scale[i], param1$shape[i]),gpd.th.l,1/10)
        ylim0 <- c(min(q1*0.95), max(q1*1.05))
        plot(newdata$gmst, q1, pch=3,cex=.3,ylim=ylim0, main='HW initialisation temperature. 10y RP',xlab='gmst', ylab="quantile")
        grid()
        mtext(fit.IVgpd$name,outer=TRUE,at=c(0.02),adj=0,col=1,line=-1)

        nd1        <- newdata
        nd1$gmst[] <- median(nd1$gmst)
        param1 <- predict(fit.IVgpd, nd1, type="response")
        q1     <- NULL
        # for(i in seq_along(param1$scale)) q1[i]   <- gpdq(c(param1$scale[i], param1$shape[i]),gpd.th.l,sample(1/retp,1))
        for(i in seq_along(param1$scale)) q1[i]   <- gpdq(c(param1$scale[i], param1$shape[i]),gpd.th.l,1/10)
        ylim0 <- c(min(q1*0.95), max(q1*1.05))
        plot(newdata$sdoy, q1, pch=3,cex=.3,ylim=ylim0, main='HW initialisation temperature. 10y RP',xlab='sdoy', ylab="quantile")
        grid()
        mtext(fit.IVgpd$name,outer=TRUE,at=c(0.02),adj=0,col=1,line=-1)

    } else {
        # plot a histogram for the stationary models
        nd1        <- newdata
        param1 <- predict(fit.IVgpd, nd1, type="response")
        q1     <- NULL
        for(i in 1:(5*length(param1$scale))) q1[i]   <- gpdq(c(sample(param1$scale,1), sample(param1$shape,1)), gpd.th.l, sample(1/retp,1))
        hist(q1,trunc(length(q1)/50),main='HW initialisation temperature (stationary model)')
        mtext(fit.IVgpd$name,outer=TRUE,at=c(0.02),adj=0,col=1,line=-1)

    }
    if(!is.null(st.pdf)) dev.off()

}

##############################################################################################
### Heffernan-Tawn related ###################################################################
##############################################################################################

##############################################################################
fn_fitJointHt <- function(events01, ht.th.u, cr.th.u, scale.doy, IsJoint=TRUE){

    ht.th.l <- qlaplace(ht.th.u)
    cr.th.l <- qlaplace(cr.th.u)

    cat(cr)
    if(IsJoint) {

        #     ### prep data  #######################################################
            #   ### simpler lagging. unwrap the 2d chain array to 1d
            #     ## obs
            #     ndata    <- length(events01$obs$ch.st.l)
            #     o12      <- array(NA,dim=c(ndata,2))
            #     o12.age  <- array(NA,dim=c(ndata,2))
            #     o12.doy  <- array(NA,dim=c(ndata,2))
            #     o12.cc   <- array(NA,dim=c(ndata,2))
            #     # day t=t data
            #     o12[ ,1]     <- events01$obs$ch.st.l[1:ndata]  # unwrap the 2d array to 1d
            #     o12.age[ ,1] <-  events01$obs$ch.age[1:ndata]
            #     o12.doy[ ,1] <-  events01$obs$ch.doy[1:ndata]
            #     o12.cc[ ,1]  <- rep(events01$obs$ch.gmst, each=dim(events01$obs$ch.st.l)[1]) # need to do something different for d1.cc
            #     # day t=t+1 data
            #         o12[1:(ndata-1) ,2] <- events01$obs$ch.st.l[2:ndata]
            #     o12.age[1:(ndata-1) ,2] <-  events01$obs$ch.age[2:ndata]
            #     o12.doy[1:(ndata-1) ,2] <-  events01$obs$ch.doy[2:ndata]
            #      o12.cc[            ,2] <- o12.cc[ ,1]
            #     # filter removing cluster of points near zero, NAs and day before and after which have ages of 0
            #     o.igood <- NULL
            #     o12.age[is.na(o12.age)] <- 0
            #     # for( i in 1:dim(o12)[1])  if(all(o12[i,]>ht.exclude.l & o12.age[i,]>0)) o.igood <- c(o.igood, i)
            #     # o.iday1E1   <- which(o12.age[o.igood,1]==1)
            #     # o.iday1G1  <- which(o12.age[o.igood,1]>1)
            #     o.iday1E1  <- which(o12.age[,1]==1)
            #     o.iday1G1  <- which(o12.age[,1]>1)
            #         #   plot(o12[1:100,1],ty='l',ylim=c(0,12))
            #         # points(o12[1:100,1],pch=3,cex=.3)
            #         #  lines(o12[1:100,2],col=2)
            #         # points(o12[1:100,2],pch=3,cex=.3)
            #         # points(o.iday1G1,o12[o.iday1G1,1],col=4,pch=1,cex=1)
            #         # points(o.iday1E1,o12[o.iday1E1,1],col=3,pch=1,cex=1)

            #     ## model
            #     ndata    <- length(events01$mod$ch.st.l)
            #     m12      <- array(NA,dim=c(ndata,2))
            #     m12.age  <- array(NA,dim=c(ndata,2))
            #     m12.doy  <- array(NA,dim=c(ndata,2))
            #     m12.cc   <- array(NA,dim=c(ndata,2))
            #     # day t=t data
            #     m12[ ,1]     <- events01$mod$ch.st.l[1:ndata]  # unwrap the 2d array to 1d
            #     m12.age[ ,1] <-  events01$mod$ch.age[1:ndata]
            #     m12.doy[ ,1] <-  events01$mod$ch.doy[1:ndata]
            #     m12.cc[ ,1]  <- rep(events01$mod$ch.gmst, each=dim(events01$mod$ch.st.l)[1]) # need to do something different for d1.cc
            #     # day t=t+1 data
            #         m12[1:(ndata-1) ,2] <- events01$mod$ch.st.l[2:ndata]
            #     m12.age[1:(ndata-1) ,2] <-  events01$mod$ch.age[2:ndata]
            #     m12.doy[1:(ndata-1) ,2] <-  events01$mod$ch.doy[2:ndata]
            #      m12.cc[            ,2] <- m12.cc[ ,1]
            #     # filter removing cluster of points near zero, NAs and day before and after which have ages of 0
            #     m.igood <- NULL
            #     m12.age[is.na(m12.age)] <- 0
            #     # for( i in 1:dim(m12)[1])  if(all(m12[i,]>ht.exclude.l & m12.age[i,]>0)) m.igood <- c(m.igood, i)
            #     # m.iday1E1   <- which(m12.age[m.igood,1]==1)
            #     # m.iday1G1  <- which(m12.age[m.igood,1]>1)
            #     m.iday1E1  <- which(m12.age[,1]==1)
            #     m.iday1G1  <- which(m12.age[,1]>1)
            #         #   plot(m12[1:100,1],ty='l',ylim=c(0,12))
            #         # points(m12[1:100,1],pch=3,cex=.3)
            #         #  lines(m12[1:100,2],col=2)
            #         # points(m12[1:100,2],pch=3,cex=.3)
            #         # points(m.iday1G1,m12[m.iday1G1,1],col=4,pch=1,cex=1)
            #         # points(m.iday1E1,m12[m.iday1E1,1],col=3,pch=1,cex=1)

            #     ## combine
            #     om12      <- rbind(o12,m12)
            #     om12.age  <- rbind(o12.age,m12.age)
            #     ## scale doy via cos((doy - midsummer_doy)/yearlength)
            #     tmp1.o    <- (o12.doy-median(events01$obs$hotseas$hotSdoy))/scale.doy$obs
            #     tmp1.m    <- (m12.doy-median(events01$mod$hotseas$hotSdoy))/scale.doy$mod
            #     om12.cosd  <- cos(2*pi*(rbind(tmp1.o,tmp1.m)))
            #     om12.cc   <- rbind(o12.cc ,m12.cc )
            #     om12.iobs <- rbind(array(1,dim=dim(o12.cc)), array(0,dim=dim(m12.cc)))
            #     # om.igood <- c(o.igood,m.igood)
            #     # om.igood1E1  <- c(o.igood[o.iday1E1], m.igood[m.iday1E1 ] +length(o12[,1]))
            #     # om.igood1G1 <- c(o.igood[o.iday1G1],m.igood[m.iday1G1] +length(o12[,1]))
            #     #         k1 <- 5500:7500
            #     #      plot(k1,om12[k1,1],ty='l',ylim=c(0,12))
            #     #     lines(k1,om12[k1,2],col=2)
            #     #     points(om.igood1E1,om12[om.igood1E1,2],col=3,pch=3,cex=.3)
            #     om.1E1 <- c(o.iday1E1, m.iday1E1 +length(o12[,1]))
            #     om.1G1 <- c(o.iday1G1, m.iday1G1 +length(o12[,1]))
            #         #     k1 <- 6500:6800
            #         #  plot(k1,om12[k1,1],ty='l',ylim=c(0,12))
            #         # lines(k1,om12[k1,2],col=2)
            #         # points(om.1E1,om12[om.1E1,1],col=3,pch=1,cex=1)
            #         # points(om.1G1,om12[om.1G1,1],col='grey70',pch=1,cex=1)
            #         # points(om.1E1,om12[om.1E1,2],col=4,pch=1,cex=1)
        #   ### end lagging

       lag01 <- fn_lag_events01(events01, scale.doy)
       list2env(lag01$om,env=parent.frame())

       #########################################################################################
       ### Day+1|Day+0
       #
        ### fit HT1 model
        cat("Fitting HT1 - day+1|day+0",cr)
        lx       <- om12[om.1E1,] # x12[igood[iday2E2],]
        lx.u     <- plaplace(lx)

        ## purge NAs
        lx.u[is.na(lx)]   <- 0
        lx[  is.na(lx)]   <- 0
        lx[  is.na(lx.u)] <- 0
        lx.u[is.na(lx.u)] <- 0

        # sort covariates
        # covdoy0  <- LaggedDoy/yearlength
        # covdoy0 <- cos(2*pi*(doyAB[igood[iday1E1],]-(which.max(HDconfig$hotseas$annCycAll)))/yearlength)
        cov.iobs <- om12.iobs[om.1E1,]
        cov.cc   <-   om12.cc[om.1E1,]
        cov.cosd <-  om12.cosd[om.1E1,]
        cov.age  <-  om12.age[om.1E1,]
        ## 2024.02 ## cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.doy[,1], cov.age[,1])
        cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.cosd[,1], cov.age[,1])
        st.c     <- c('I','C','D','A')
        colnames(cov.ICDA) <- st.c

        f1       <- function(a) {which(st.c==a)}

        # remove all below ht.exclude.l
        ix0      <- apply(lx[,1:2],1,FUN=function(x,b,c) {return(all(x[1]>=b & x[2]>=c))}, b=ht.th.l, c=ht.exclude.l)
        # chi.cov0 <- NULL # min(covcc0[ix0,1:2],na.rm=T)
        d0       <- lx.u[ix0,1:2]
        mlu      <- ht.th.u
        clu      <- cr.th.u
        c0       <- cov.ICDA[ix0,]

        ht1.a0.b0.m0.s0 <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, lam.pen=lam.pen0)
        ht1.a0.b0.m0.s0$nmod <- 'a0_b0'
        cat("ht1.a0.b0.m0.s0 nllh par", cr,ht1.a0.b0.m0.s0$nllh, myround(ht1.a0.b0.m0.s0$par,3), cr)
        # par0    <- ht1.a0.b0.m0.s0$par # all are on the linked scale

        ### fit each covariate to each parameter in turn and order priority based on nllh
        # impose structure so that alpha and beta covariate dependence takes precidence over mu and sigma
        # alpha and beta
        in0                  <- NULL
        ht1.singleAB         <- list(alpha=NULL, beta=NULL)
        ht1.singleAB$alpha$I <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('I'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        ht1.singleAB$alpha$C <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('C'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        ht1.singleAB$alpha$D <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('D'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        # ht1.singleAB$alpha$A <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('A'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        ht1.singleAB$beta$I  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('I'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        ht1.singleAB$beta$C  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('C'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        ht1.singleAB$beta$D  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('D'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        # ht1.singleAB$beta$A  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('A'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)

        # ht1.singleAB$alpha$C <- SjbBveHTDepCovHrc(data=d0, alpcovind=1, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # ht1.singleAB$alpha$D <- SjbBveHTDepCovHrc(data=d0, alpcovind=2, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # # ht1.singleAB$alpha$A <- SjbBveHTDepCovHrc(data=d0, alpcovind=3, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # ht1.singleAB$beta$C  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=1, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # ht1.singleAB$beta$D  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=2, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # # ht1.singleAB$beta$A  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=3, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)

        ht1.1AB.nllh <- sapply(sapply(ht1.singleAB,'['),'[[','nllh')
        ht1.1AB.par  <- sapply(sapply(ht1.singleAB,'['),'[[','par')
        ht1.1AB.pchisq  <- NULL
        ht1.0         <- ht1.a0.b0.m0.s0 # ht1.AB[[1]]
        l.temp        <- c(ht1.singleAB$alpha, ht1.singleAB$beta)
        for(i in seq_along(l.temp)) ht1.1AB.pchisq[i] <-  pchisq(2*(ht1.0$nllh - l.temp[[i]]$nllh), df=length(l.temp[[i]]$par)-length(ht1.0$par) )

        ht1.1AB.mod <- c('alpha_I','alpha_C','alpha_D','beta_I','beta_C','beta_D')
        st.c2    <- c('I','C','D','I','C','D')
        st.A     <- list('I','C','D',NULL,NULL,NULL)
        st.B     <- list(NULL,NULL,NULL,'I','C','D')
        seq.AB   <- sort(ht1.1AB.nllh, index=TRUE)$ix   # sort on nllh
        # cat(cr,'fn_fitJointHt HT1: Order of importance for alpha and beta:',cr,ht1.1AB.mod[seq.AB],cr,cr)
        cat(cr,'fn_fitJointHt HT1: Order of importance for alpha and beta: \nname pchi',cr)
        for(i in seq_along(seq.AB)) cat(ht1.1AB.mod[seq.AB[i]],myround(ht1.1AB.pchisq[seq.AB[i]],3),cr)
        cat(cr,cr)

        ### Are any covariates having an effect ?
        if( any(2*( ht1.a0.b0.m0.s0$nllh -ht1.1AB.nllh) > qchisq(.9,df=1))) {
            cat("Adding covariates that have impact on nllh",cr)

            ### add these hierarchy
            # alpha and beta first
            # mo.AB  <- list(aD=list(al=2, be=NULL),   aA=list(al=3, be=NULL),   bD=list(al=NULL, be=2),   bA=list(al=NULL, be=3))
            # mo.AB  <- list(aC=list(al=1, be=NULL),   aD=list(al=2, be=NULL),   bC=list(al=NULL, be=1),   bD=list(al=NULL, be=2))
            # ht1.AB <- list()
            # al.ci <- NULL
            # be.ci <- NULL
            # for(i in seq_along(seq.AB)) {
            #     j          <- seq.AB[i]
            #     al.ci      <- c(al.ci, mo.AB[[j]]$al)
            #     be.ci      <- c(be.ci, mo.AB[[j]]$be)
            #     cat('cov indx al', al.ci, 'be', be.ci, cr)
            #     ht1.AB[[i]] <- SjbBveHTDepCovHrc(data=d0, alpcovind=al.ci, betcovind=be.ci, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
            # }
            # mo.AB  <- list(aC=list(al=1, be=NULL),   aD=list(al=2, be=NULL),   bC=list(al=NULL, be=1),   bD=list(al=NULL, be=2))
            ht1.AB <- list()
            ht1.AB.mod  <- NULL
            al.ci  <- NULL
            be.ci  <- NULL
            for(i in seq_along(seq.AB)) {
                j          <- seq.AB[i]
                al.ci      <- c(al.ci, unlist(st.A[seq.AB[i]]))
                be.ci      <- c(be.ci, unlist(st.B[seq.AB[i]]))
                cat('cov indx al', al.ci, 'be', be.ci, cr)
                # ht1.AB[[i]] <- SjbBveHTDepCovHrc(data=d0, alpcovind=unlist(st.A[seq.AB[i]]), betcovind=unlist(st.B[seq.AB[i]]), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
                ht1.AB.mod[i] <- paste(paste('a', paste(al.ci, sep='', collapse=''),sep=''), paste('b', paste(be.ci, sep='', collapse=''),sep=''),sep='_')
                ht1.AB[[i]] <- SjbBveHTDepCovHrc(data=d0, alpcovind=al.ci, betcovind=be.ci, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
                ht1.AB[[i]]$nmod <- ht1.AB.mod[i]
            }
            cat("ht1.AB.mod",ht1.AB.mod,cr,cr)

            ht1.AB.nllh <- sapply(ht1.AB,'[[','nllh')
            ht1.AB.par  <- sapply(ht1.AB,'[[','par')
            ht1.AB.pchisq  <- NULL
            ht1.topAB   <- ht1.a0.b0.m0.s0 # ht1.AB[[1]]
            for(i in seq_along(seq.AB)) {
                ht1.AB.pchisq[i] <-  pchisq(2*(ht1.topAB$nllh - ht1.AB[[i]]$nllh), df=length(ht1.AB[[i]]$par)-length(ht1.topAB$par) )
                if(ht1.AB.pchisq[i] > ht.cov.sig.th) ht1.topAB <- ht1.AB[[i]]
            }
            ht1.top <- ht1.topAB

            # SjbPlotHtParametersCov(ht1.top, newwindow=PARPLTNEW)
            cat(cr,"ht1.top nllh par", ht1.top$nmod, ht1.top$nllh, myround(ht1.top$par,3), cr, cr)

            ht1.1AB  <- list(mod=ht1.1AB.mod, nllh=ht1.1AB.nllh, par=ht1.1AB.par, pchisq=ht1.1AB.pchisq)
            ht1.AB   <- list(mod=ht1.AB.mod,  nllh=ht1.AB.nllh,  par=ht1.AB.par,  pchisq=ht1.AB.pchisq)
            # sigtests <- list(ht1.1AB.pchisq=ht1.1AB.pchisq, ht1.AB.pchisq=ht1.AB.pchisq, ht1.1MS.pchisq=NULL, ht1.MS.pchisq=NULL)

            ### compare with other approaches to fitting
            cat("Compare with this uninitialised fit aICD.bICD:",cr)
            # e for free on alpha and beta, no initialising
            in0 <- NULL
            ht1.au.bu.m0.s0  <- SjbBveHTDepCovHrc(data=d0, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, alpcovind=st.c[1:3], betcovind=st.c[1:3], mucovind=NULL, sigcovind=NULL, lam.pen=lam.pen0, initvec0=in0)
            # SjbPlotHtParametersCov(ht1.au.bu.m0.s0, newwindow=PARPLTNEW)
            cat("nllh par", ht1.au.bu.m0.s0$nllh, myround(ht1.au.bu.m0.s0$par,3), cr,cr)
            # # e for free on all, no initialising
            # in0 <- NULL
            # ht1.au.bu.me.se  <- SjbBveHTDepCovHrc(data=d0, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, alpcovind=1:2, betcovind=1:2, mucovind=1:2, sigcovind=1:2, nsim=nsim1, trace=DOTRACE, lam.pen=lam.pen0, initvec0=in0, DOPLOT=PLTFIT)
            # SjbPlotHtParametersCov(ht1.au.bu.me.se, newwindow=PARPLTNEW)
            # cat("ht1.au.bu.me.se nllh par", ht1.au.bu.me.se$nllh, ht1.au.bu.me.se$par, cr)
        } else {
            ht1.top  <- ht1.a0.b0.m0.s0
            ht1.1AB  <- list(mod=ht1.1AB.mod, nllh=ht1.1AB.nllh, par=ht1.1AB.par, pchisq=ht1.1AB.pchisq)
            ht1.AB   <- list(mod=NULL,  nllh=NULL,  par=NULL,  pchisq=NULL)
            # sigtests <- list(ht1.1AB.pchisq=ht1.1AB.pchisq, ht1.AB.pchisq=NULL, ht1.1MS.pchisq=NULL, ht1.MS.pchisq=NULL)
            cat(cr,'fn_fitJointHt HT1: No covariate effect found for alpha and beta:',cr,cr)
        }
        ht1.top$om <- d0
        ht1.info <- list(ht1.1AB=ht1.1AB, ht1.AB=ht1.AB)
        cat("fn_fitJointHt Day+1|Day+0 done",cr,"#################",cr,cr)
       #
       ### end Day+1|Day+0

       #########################################################################################
       ### DayN+1|DayN
       #

        cat("fn_fitJointHt Fitting HTN - dayN+1|day+N",cr)
        lx       <- om12[om.1G1,] # x12[igood[iday2E2],]
        lx.u     <- plaplace(lx)
        ## purge NAs
        lx.u[is.na(lx)]   <- 0
        lx[  is.na(lx)]   <- 0
        lx[  is.na(lx.u)] <- 0
        lx.u[is.na(lx.u)] <- 0

        # sort covariates
        # covdoy0  <- LaggedDoy/yearlength
        # covdoy0 <- cos(2*pi*(doyAB[igood[iday1E1],]-(which.max(HDconfig$hotseas$annCycAll)))/yearlength)
        cov.iobs <- om12.iobs[om.1G1,]
        cov.cc   <-   om12.cc[om.1G1,]
        cov.cosd <- om12.cosd[om.1G1,]
        cov.age  <-  om12.age[om.1G1,]
        ## 2024.02 ## cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.doy[,1], cov.age[,1])
        cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.cosd[,1], cov.age[,1])

        st.c     <- c('I','C','D','A')
        colnames(cov.ICDA) <- st.c

        f1       <- function(a) {which(st.c==a)}

        # remove all below ht.exclude.l
        ix0      <- apply(lx[,1:2],1,FUN=function(x,b,c) {return(all(x[1]>=b & x[2]>=c))}, b=ht.th.l, c=ht.exclude.l)
        # chi.cov0 <- NULL # min(covcc0[ix0,1:2],na.rm=T)
        d0       <- lx.u[ix0,1:2]
        c0       <- cov.ICDA[ix0,]
        mlu      <- ht.th.u
        clu      <- cr.th.u

        htN.a0.b0.m0.s0 <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, lam.pen=lam.pen0)
        cat("htN.a0.b0.m0.s0 nllh par", cr,htN.a0.b0.m0.s0$nllh, htN.a0.b0.m0.s0$par, cr)
        # par0    <- htN.a0.b0.m0.s0$par # all are on the linked scale

        ### fit each covariate to each parameter in turn and order priority based on nllh
        # impose structure so that alpha and beta covariate dependence takes precidence over mu and sigma
        # alpha and beta
        in0                  <- NULL
        htN.singleAB         <- list(alpha=NULL, beta=NULL)
        htN.singleAB$alpha$I <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('I'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$alpha$C <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('C'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$alpha$D <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('D'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$alpha$A <- SjbBveHTDepCovHrc(data=d0, alpcovind=f1('A'), betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$beta$I  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('I'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$beta$C  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('C'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$beta$D  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('D'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
        htN.singleAB$beta$A  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=f1('A'), mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)

        # htN.singleAB$alpha$C <- SjbBveHTDepCovHrc(data=d0, alpcovind=1, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # htN.singleAB$alpha$D <- SjbBveHTDepCovHrc(data=d0, alpcovind=2, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # # htN.singleAB$alpha$A <- SjbBveHTDepCovHrc(data=d0, alpcovind=3, betcovind=NULL, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # htN.singleAB$beta$C  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=1, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # htN.singleAB$beta$D  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=2, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)
        # # htN.singleAB$beta$A  <- SjbBveHTDepCovHrc(data=d0, alpcovind=NULL, betcovind=3, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, nsim=nsim1, initvec0=in0, DOPLOT=F)

        htN.1AB.nllh    <- sapply(sapply(htN.singleAB,'['),'[[','nllh')
        htN.1AB.par     <- sapply(sapply(htN.singleAB,'['),'[[','par')
        htN.1AB.pchisq  <- NULL
        htN.0         <- htN.a0.b0.m0.s0 # htN.AB[[1]]
        l.temp        <- c(htN.singleAB$alpha, htN.singleAB$beta)
        for(i in seq_along(l.temp)) htN.1AB.pchisq[i] <-  pchisq(2*(htN.0$nllh - l.temp[[i]]$nllh), df=length(l.temp[[i]]$par)-length(htN.0$par) )
        # for(i in seq_along(l.temp)) l.temp[[i]]$nllh

        htN.1AB.mod <- c('alpha_I','alpha_C','alpha_D','alpha_A','beta_I','beta_C','beta_D','beta_A')
        st.c2    <-    c('I','C','D','A','I','C','D','A')
        st.A     <- list('I','C','D','A',NULL,NULL,NULL,NULL)
        st.B     <- list(NULL,NULL,NULL,NULL,'I','C','D','A')
        seq.AB   <- sort(htN.1AB.nllh, index=TRUE)$ix
        # cat(cr,'fn_fitJointHt HTn: Order of importance for alpha and beta: ',cr,htN.1AB.mod[seq.AB],cr,cr)
        cat(cr,'fn_fitJointHt HTn: Order of importance for alpha and beta: \nname pchi',cr)
        for(i in seq_along(seq.AB)) cat(htN.1AB.mod[seq.AB[i]],myround(htN.1AB.pchisq[seq.AB[i]],3),cr)
        cat(cr,cr)

        ### Are any covariates having an effect ?
        if( any(2*( htN.a0.b0.m0.s0$nllh -htN.1AB.nllh) > qchisq(.9,df=1))) {
            cat("Adding covariates that have impact on nllh",cr)

            ### add these hierarchy
            # alpha and beta first
            # mo.AB  <- list(aD=list(al=2, be=NULL),   aA=list(al=3, be=NULL),   bD=list(al=NULL, be=2),   bA=list(al=NULL, be=3))
            # mo.AB  <- list(aC=list(al=1, be=NULL),   aD=list(al=2, be=NULL),   bC=list(al=NULL, be=1),   bD=list(al=NULL, be=2))
            # htN.AB <- list()
            # al.ci <- NULL
            # be.ci <- NULL
            # for(i in seq_along(seq.AB)) {
            #     j          <- seq.AB[i]
            #     al.ci      <- c(al.ci, mo.AB[[j]]$al)
            #     be.ci      <- c(be.ci, mo.AB[[j]]$be)
            #     cat('cov indx al', al.ci, 'be', be.ci, cr)
            #     htN.AB[[i]] <- SjbBveHTDepCovHrc(data=d0, alpcovind=al.ci, betcovind=be.ci, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
            # }
            # mo.AB  <- list(aC=list(al=1, be=NULL),   aD=list(al=2, be=NULL),   bC=list(al=NULL, be=1),   bD=list(al=NULL, be=2))
            htN.AB <- list()
            htN.AB.mod <- NULL
            al.ci  <- NULL
            be.ci  <- NULL
            for(i in seq_along(seq.AB)) {
                j          <- seq.AB[i]
                al.ci      <- c(al.ci, unlist(st.A[seq.AB[i]]))
                be.ci      <- c(be.ci, unlist(st.B[seq.AB[i]]))
                cat('cov indx al', al.ci, 'be', be.ci, cr)
                htN.AB.mod[i] <- paste(paste('a', paste(al.ci, sep='', collapse=''),sep=''), paste('b', paste(be.ci, sep='', collapse=''),sep=''),sep='_')
                htN.AB[[i]] <- SjbBveHTDepCovHrc(data=d0, alpcovind=al.ci, betcovind=be.ci, mucovind=NULL, sigcovind=NULL, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, initvec0=in0)
                htN.AB[[i]]$nmod <- htN.AB.mod[i]
            }
            cat("htN.AB.mod",htN.AB.mod,cr,cr)

            htN.AB.nllh <- sapply(htN.AB,'[[','nllh')
            htN.AB.par  <- sapply(htN.AB,'[[','par')
            htN.AB.pchisq  <- NULL
            htN.topAB   <- htN.a0.b0.m0.s0 # htN.AB[[1]]
            for(i in seq_along(seq.AB)) {
                htN.AB.pchisq[i] <-  pchisq(2*(htN.topAB$nllh - htN.AB[[i]]$nllh), df=length(htN.AB[[i]]$par)-length(htN.topAB$par) )
                if(htN.AB.pchisq[i] > ht.cov.sig.th) htN.topAB <- htN.AB[[i]]
            }
            htN.top <- htN.topAB

            # SjbPlotHtParametersCov(htN.top, newwindow=PARPLTNEW)
            cat(cr,"htN.top nllh par", htN.top$nmod, htN.top$nllh, myround(htN.top$par,3), cr, cr)

            htN.1AB  <- list(mod=htN.1AB.mod, nllh=htN.1AB.nllh, par=htN.1AB.par, pchisq=htN.1AB.pchisq)
            htN.AB   <- list(mod=htN.AB.mod,  nllh=htN.AB.nllh,  par=htN.AB.par,  pchisq=htN.AB.pchisq)
            # sigtests       <- list(htN.1AB.pchisq=htN.1AB.pchisq, htN.AB.pchisq=htN.AB.pchisq, htN.1MS.pchisq=NULL, htN.MS.pchisq=NULL)

            # ### compare with other approaches to fitting
            # cat("Compare fits with these:",cr)
            # # e for free on alpha and beta, no initialising
            # in0 <- NULL
            # htN.au.bu.m0.s0  <- SjbBveHTDepCovHrc(data=d0, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, alpcovind=st.c[1:3], betcovind=st.c[1:3], mucovind=NULL, sigcovind=NULL, lam.pen=lam.pen0, initvec0=in0)
            # # SjbPlotHtParametersCov(htN.au.bu.m0.s0, newwindow=PARPLTNEW)
            # cat("htN.au.bu.m0.s0 nllh par", htN.au.bu.m0.s0$nllh, htN.au.bu.m0.s0$par, cr)
            # # # e for free on all, no initialising
            # # in0 <- NULL
            # # htN.au.bu.me.se  <- SjbBveHTDepCovHrc(data=d0, mod.lev.u=mlu, crit.lev.u=clu, cvte=c0, alpcovind=1:2, betcovind=1:2, mucovind=1:2, sigcovind=1:2, nsim=nsim1, trace=DOTRACE, lam.pen=lam.pen0, initvec0=in0, DOPLOT=PLTFIT)
            # # SjbPlotHtParametersCov(htN.au.bu.me.se, newwindow=PARPLTNEW)
            # # cat("htN.au.bu.me.se nllh par", htN.au.bu.me.se$nllh, htN.au.bu.me.se$par, cr)
        } else {
            htN.top  <- htN.a0.b0.m0.s0
            htN.1AB  <- list(mod=htN.1AB.mod, nllh=htN.1AB.nllh, par=htN.1AB.par, pchisq=htN.1AB.pchisq)
            htN.AB   <- list(mod=NULL,  nllh=NULL,  par=NULL,  pchisq=NULL)
            # sigtests <- list(htN.1AB.pchisq=htN.1AB.pchisq, htN.AB.pchisq=NULL, htN.1MS.pchisq=NULL, htN.MS.pchisq=NULL)
            cat(cr,'fn_fitJointHt HTn: No covariate effect found for alpha and beta:',cr,cr)
        }
        htN.top$om <- d0
        htN.info <- list(htN.1AB=htN.1AB, htN.AB=htN.AB)

        cat("fn_fitJointHt DayN+1|Day+N done",cr,"#################",cr,cr)
       #
       ### end DayN+1|Day+N

        cat("fn_fitHt: Completed",cr)
        htJoint <- list(ht1=ht1.top, htN=htN.top, ht1.info=ht1.info, htN.info=htN.info, ht1.0=ht1.a0.b0.m0.s0, htN.0=htN.a0.b0.m0.s0)
        return(htJoint)
    } else {
        cat("fn_fitHt: only joint obs+model implemented at this point",cr)
        return(NULL)
    }
}

##############################################################################
fn_htJointDiagnostics <- function(htJ, stout_Ht='htJointDiagnostics.RData') {

    # must pass in joint ht produced by fn_fitJointHt

    stout_HtDiag     <- sub('.RData','.diag.txt',stout_Ht)
    conDiag    <- file(stout_HtDiag,'w')
    writeLines(stout_HtDiag,con=conDiag)
    writeLines(stout_Ht,    con=conDiag)

    ### to console first
    ## ht1
    cat(cr)
    cat("### Fitting HT1 - day+1|day+0 ###########",cr)

        cat("Model a0_b0", cr, "nllh", htJ$ht1.0$nllh, cr, "par", htJ$ht1.0$par, cr)
        mod1  <- htJ$ht1.info$ht1.1AB$mod
        nllh1 <- htJ$ht1.info$ht1.1AB$nllh
        cat(cr,'ht1 AB1: Order of importance for alpha and beta:',cr,mod1[sort(nllh1, index=T)$ix],cr)
        # cat(cr,'Order of importance for mu and sigma:',cr,s0MS.mod[sort(s0MS.nllh, index=T)$ix],cr,cr)

        if(!is.null(htJ$ht1.info$ht1.AB$mod)) {
            mod1  <- htJ$ht1.info$ht1.AB$mod
            nllh1 <- htJ$ht1.info$ht1.AB$nllh
            cat(cr,'ht1 AB: Order of importance for combined covaraites:',cr,mod1[sort(nllh1, index=T)$ix],cr,cr)
            cat("ht1.top", cr, "Model",mod1[1], cr, "nllh", htJ$ht1$nllh, cr, "par", htJ$ht1$par, cr, cr)
        } else {
            cat(cr,'AB: No significant covariates for ht1',cr,cr)
            cat("ht1.top", cr, "Model a0_b0", cr, "nllh", htJ$ht1$nllh, cr, "par", htJ$ht1$par, cr, cr)
        }
    ## htN
    cat("### Fitting HT - day+N|day+1 ###########",cr)

        cat("Model a0_b0", cr, "nllh", htJ$htN.0$nllh, cr, "par", htJ$htN.0$par, cr)
        mod1  <- htJ$htN.info$htN.1AB$mod
        nllh1 <- htJ$htN.info$htN.1AB$nllh
        cat(cr,'htN AB1: Order of importance for alpha and beta:',cr,mod1[sort(nllh1, index=T)$ix],cr)
        # cat(cr,'Order of importance for mu and sigma:',cr,s0MS.mod[sort(s0MS.nllh, index=T)$ix],cr,cr)

        if(!is.null(htJ$htN.info$htN.AB$mod)) {
            mod1  <- htJ$htN.info$htN.AB$mod # alredy sorted
            nllh1 <- htJ$htN.info$htN.AB$nllh
            cat(cr,'htN AB: Order of importance for combined covaraites:',cr,mod1,cr,cr)
            cat("htN.top", cr, "Model",mod1[1], cr, "nllh", htJ$htN$nllh, cr, "par", htJ$htN$par, cr, cr)
        } else {
            cat(cr,'AB: No significant covariates for htN',cr,cr)
            cat("htN.top", cr, "Model a0_b0", cr, "nllh", htJ$htN$nllh, cr, "par", htJ$htN$par, cr, cr)
        }

    ### write out to file
    ## ht1
    writeLines(capture.output(cat(cr)), con=conDiag)
    writeLines(capture.output(cat("### Fitting HT1 - day+1|day+0 ###########",cr)), con=conDiag)

        writeLines(capture.output(cat("Model a0_b0", cr, "nllh", htJ$ht1.0$nllh, cr, "par", htJ$ht1.0$par, cr)), con=conDiag)
        mod1  <- htJ$ht1.info$ht1.1AB$mod
        nllh1 <- htJ$ht1.info$ht1.1AB$nllh
        writeLines(capture.output(cat(cr,'ht1 AB1: Order of importance for alpha and beta:',cr,mod1[sort(nllh1, index=T)$ix],cr)), con=conDiag)
        # writeLines(capture.output(cat(cr,'Order of importance for mu and sigma:',cr,s0MS.mod[sort(s0MS.nllh, index=T)$ix],cr,cr)), con=conDiag)

        if(!is.null(htJ$ht1.info$ht1.AB$mod)) {
            mod1  <- htJ$ht1.info$ht1.AB$mod
            nllh1 <- htJ$ht1.info$ht1.AB$nllh
            writeLines(capture.output(cat(cr,'ht1 AB: Order of importance for combined covaraites:',cr,mod1[sort(nllh1, index=T)$ix],cr,cr)), con=conDiag)
            writeLines(capture.output(cat("ht1.top", cr, "Model",mod1[1], cr, "nllh", htJ$ht1$nllh, cr, "par", htJ$ht1$par, cr, cr)), con=conDiag)
        } else {
            writeLines(capture.output(cat(cr,'AB: No significant covariates for ht1',cr,cr)), con=conDiag)
            writeLines(capture.output(cat("ht1.top", cr, "Model a0_b0", cr, "nllh", htJ$ht1$nllh, cr, "par", htJ$ht1$par, cr, cr)), con=conDiag)
        }
    ## htN
    writeLines(capture.output(cat("### Fitting HT - day+N|day+1 ###########",cr)), con=conDiag)

        writeLines(capture.output(cat("Model a0_b0", cr, "nllh", htJ$htN.0$nllh, cr, "par", htJ$htN.0$par, cr)), con=conDiag)
        mod1  <- htJ$htN.info$htN.1AB$mod
        nllh1 <- htJ$htN.info$htN.1AB$nllh
        writeLines(capture.output(cat(cr,'htN AB1: Order of importance for alpha and beta:',cr,mod1[sort(nllh1, index=T)$ix],cr)), con=conDiag)
        # writeLines(capture.output(cat(cr,'Order of importance for mu and sigma:',cr,s0MS.mod[sort(s0MS.nllh, index=T)$ix],cr,cr)), con=conDiag)

        if(!is.null(htJ$htN.info$htN.AB$mod)) {
            mod1  <- htJ$htN.info$htN.AB$mod # alredy sorted
            nllh1 <- htJ$htN.info$htN.AB$nllh
            writeLines(capture.output(cat(cr,'htN AB: Order of importance for combined covaraites:',cr,mod1,cr,cr)), con=conDiag)
            writeLines(capture.output(cat("htN.top", cr, "Model",mod1[1], cr, "nllh", htJ$htN$nllh, cr, "par", htJ$htN$par, cr, cr)), con=conDiag)
        } else {
            writeLines(capture.output(cat(cr,'AB: No significant covariates for htN',cr,cr)), con=conDiag)
            writeLines(capture.output(cat("htN.top", cr, "Model a0_b0", cr, "nllh", htJ$htN$nllh, cr, "par", htJ$htN$par, cr, cr)), con=conDiag)
        }

    close(conDiag)

}

##############################################################################
fn_htJointPlot <- function(htJ, dodays=c(3,4,5,6), stplot=NULL){

    if(is.null(stplot)) x11() else pdf(file=stplot, 10, 10)

    # lag01 <- fn_lag_events01(events01, scale.doy)
    # list2env(lag01$om,env=parent.frame())


    ### ht1 ################################################################################
    # d0.l     <- om12[om.1E1,]
    d0.l     <- qlaplace(htJ$ht1$om)
    d5.l     <- rep(d0.l[,1],5)
    mod.th.l <- qlaplace(htJ$ht1$mod.lev.u)
    # cov.iobs <- om12.iobs[om.1E1,]
    # cov.cc   <-   om12.cc[om.1E1,]
    # cov.cosd  <-  om12.cosd[om.1E1,]
    # cov.age  <-  om12.age[om.1E1,]
    # cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.doy[,1], cov.age[,1])
    # st.c     <- c('I','C','D','A')
    # colnames(cov.ICDA) <- st.c

    if(length(htJ$ht1$par)==4) { # stationalry ht1 model

        sim0 <- sim.from.htHrc(htJ$ht1, d5.l, cov=htJ$ht1$model$cvte , nonoise=FALSE)

        up.2by2()
        igs  <- which(sim0 > ht.exclude.l)
        plot(d5.l[igs], sim0[igs],  pch=3,cex=.3,col=2, main='Ht1 Day1 vs Day2', xlab='Day1', ylab='Day2',ylim=c(-2,max(sim0)))
        points(  d0.l[,1],         d0.l[,2],    pch=20)

        ### QQ verification plots
        qqplot(d0.l[,2], sample(sim0[igs],length(d0.l[,2])),pch=3,cex=.3, main=paste('HT1 Day ',2,': obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
        if(is.null(stplot)) readline("Next?")

    } else { # non stationary ht1 model

        cvte   <- htJ$ht1$model$cvte
        n.cvte <- dimnames(cvte)
        covsim <- NULL
        for(i in 1:dim(cvte)[2]) {
            if(n.cvte[[2]][i]=='I') covsim <- cbind(covsim, c(0,0,1)) else covsim <- cbind(covsim, quantile(cvte[,i],c(0.05,.5,.95)))
        }
        dimnames(covsim) <- n.cvte

        sim0 <- sim.from.htHrc(htJ$ht1, d5.l, cov=rbind(cvte,cvte,cvte,cvte,cvte), nonoise=FALSE)
        simX <- NULL
        simN <- NULL
        for(i in 1:dim(cvte)[2]) {
            cov1      <- covsim[2,]
            cov1[i]   <- max(covsim[,i],na.rm=TRUE)
            simX[[i]] <- sim.from.htHrc(htJ$ht1, d5.l, cov=cov1 , nonoise=FALSE)
            cov1[i]   <- min(covsim[,i],na.rm=TRUE)
            simN[[i]] <- sim.from.htHrc(htJ$ht1, d5.l, cov=cov1 , nonoise=FALSE)
        }

        up.2by2()
        par(mfrow=c(2,1+ceiling(dim(cvte)[2]/2)))
        igs  <- which(sim0 > ht.exclude.l)
        plot(d5.l[igs],  sim0[igs],  pch=3,cex=.3,col=2, main='HT1 Day1 vs Day2', xlab='Day1', ylab='Day2',ylim=c(-2,max(sim0)))
        points(d0.l[,1],  d0.l[,2],    pch=20)

        ### QQ verification plots
        qqplot(d0.l[,2], sample(sim0[igs],length(d0.l[,2])),pch=3,cex=.3, main=paste('HT1 Day ',2,': obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
        if(is.null(stplot)) readline("Next?")

        ### covariate
        for(i in 1:dim(cvte)[2]) {
            if(grepl(n.cvte[[2]][i], htJ$ht1$nmod, fixed=TRUE)) {
                # qqplot(d0.l[,2], sample(sim0[igs],length(d0.l[,2])),pch=3,cex=.3, main=paste('HT1 Day ',2,': obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
                igs1   <- which(simN[[i]]         >mod.th.l)
                if (length(igs1)<200) igs1 <- seq_along(simN[[i]])
                igs2   <- which(simX[[i]]         >mod.th.l)
                if (length(igs2)<200) igs2 <- seq_along(simX[[i]])
                # r2 <- paste(round(covsim[c(1,3)],1),collapse=' ')
                qqplot(simN[[i]][igs1], simX[[i]][igs2], pch=3,cex=.3, main=paste('HT1 Effect of',n.cvte[[2]][i],'(max v min)'),xlab='Low',ylab='High')
            }
        }
        if(is.null(stplot)) readline("Next?")

    } # end non stationary ht1 model

    ### htN ################################################################################
    # d0.l     <- om12[om.1G1,]
    d0.l     <- qlaplace(htJ$htN$om)
    d5.l     <- rep(d0.l[,1],5)
    mod.th.l <- qlaplace(htJ$htN$mod.lev.u)
    # cov.iobs <- om12.iobs[om.1E1,]
    # cov.cc   <-   om12.cc[om.1E1,]
    # cov.cosd  <-  om12.cosd[om.1E1,]
    # cov.age  <-  om12.age[om.1E1,]
    # cov.ICDA <- cbind(cov.iobs[,1], cov.cc[,1], cov.doy[,1], cov.age[,1])
    # st.c     <- c('I','C','D','A')
    # colnames(cov.ICDA) <- st.c
    cvte   <- htJ$htN$model$cvte
    n.cvte <- dimnames(cvte)

    if(length(htJ$htN$par)==4) { # stationalry htN model

        sim0 <- sim.from.htHrc(htJ$htN, d5.l, cov=cvte , nonoise=FALSE)

        up.2by3()
        igs  <- which(sim0 > ht.exclude.l)
        plot(d5.l[igs], sim0[igs],  pch=3,cex=.3,col=2, main='HtN DayN vs DayN+1', xlab='DayN', ylab='DayN+1',ylim=c(-2,max(sim0)))
        points(  d0.l[,1],         d0.l[,2],    pch=20)

        ### QQ verification plots
        iA <- which(n.cvte[[2]]=='A')
        for(dd in dodays) {
            i0 <- which(cvte[,iA]==(dd-1)) # index for day before
            x0 <- cov1 <- NULL
            while(length(x0)<500) {
                x0   <- c(x0, d0.l[i0,1])
                cov1 <- rbind(cov1,cvte[i0,])
            }
            sim0 <- sim.from.htHrc(htJ$htN, x0, cov=cov1 , nonoise=FALSE)
            igs  <- which(sim0 > ht.exclude.l)
            qqplot(d0.l[i0,2], sample(sim0[igs],length(d0.l[i0,2])),pch=3,cex=.3, main=paste('HTN Day',dd,'obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
        }
        if(is.null(stplot)) readline("Next?")

    } else { # non stationary htN model

        up.2by3()

        sim0 <- sim.from.htHrc(htJ$htN, d5.l, cov=rbind(cvte,cvte,cvte,cvte,cvte), nonoise=FALSE)
        par(mfcol=c(2,1+ceiling(dim(cvte)[2]/2)))
        igs  <- which(sim0 > ht.exclude.l)
        plot(d5.l[igs],  sim0[igs],  pch=3,cex=.3,col=2, main='HtN DayN vs DayN+1', xlab='DayN', ylab='DayN+1',ylim=c(-2,max(sim0)))
        points(d0.l[,1],  d0.l[,2],    pch=20)

        # qqplot(d0.l[,2], sample(sim0[igs],length(d0.l[,2])),pch=3,cex=.3, main=paste('HTN DayN+1 : obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
        ### QQ verification plots
        iA <- which(n.cvte[[2]]=='A')
        for(dd in dodays) {
            i0 <- which(cvte[,iA]==(dd-1)) # index for day before
            x0 <- cov1 <- NULL
            while(length(x0)<500) {
                x0   <- c(x0, d0.l[i0,1])
                cov1 <- rbind(cov1,cvte[i0,])
            }
            sim0 <- sim.from.htHrc(htJ$htN, x0, cov=cov1 , nonoise=FALSE)
            igs  <- which(sim0 > ht.exclude.l)
            qqplot(d0.l[i0,2], sample(sim0[igs],length(d0.l[i0,2]),replace=TRUE),pch=3,cex=.3, main=paste('HTN Day',dd,'obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
        }
        if(is.null(stplot)) readline("Next?")

        ### covariate plots
        up.2by2()
        covsim <- NULL
        for(i in 1:dim(cvte)[2]) {
            if(n.cvte[[2]][i]=='I') covsim <- cbind(covsim, c(0,0,1)) else covsim <- cbind(covsim, quantile(cvte[,i],c(0.05,.5,.95)))
        }
        dimnames(covsim) <- n.cvte
        simX <- NULL
        simN <- NULL
        for(i in 1:dim(cvte)[2]) {
            cov1      <- covsim[2,]
            cov1[i]   <- max(covsim[,i],na.rm=TRUE)
            simX[[i]] <- sim.from.htHrc(htJ$htN, d5.l, cov=cov1 , nonoise=FALSE)
            cov1[i]   <- min(covsim[,i],na.rm=TRUE)
            simN[[i]] <- sim.from.htHrc(htJ$htN, d5.l, cov=cov1 , nonoise=FALSE)
        }

        for(i in 1:dim(cvte)[2]) {
            if(grepl(n.cvte[[2]][i], htJ$htN$nmod, fixed=TRUE)) {
                # qqplot(d0.l[,2], sample(sim0[igs],length(d0.l[,2])),pch=3,cex=.3, main=paste('HT1 Day ',2,': obs v sim'),xlab='Obs',ylab='Sim',ylim=c(2,max(sim0)),xlim=c(2,max(sim0)))
                igs1   <- which(simN[[i]]         >mod.th.l)
                if (length(igs1)<200) igs1 <- seq_along(simN[[i]])
                igs2   <- which(simX[[i]]         >mod.th.l)
                if (length(igs2)<200) igs2 <- seq_along(simX[[i]])
                # r2 <- paste(round(covsim[c(1,3)],1),collapse=' ')
                qqplot(simN[[i]][igs1], simX[[i]][igs2], pch=3,cex=.3, main=paste('HTN Effect of',n.cvte[[2]][i],'(max v min)'),xlab='Low',ylab='High')
            }
        }
        if(is.null(stplot)) readline("Next?")

    } # end non stationary htN model

    if(!is.null(stplot)) dev.off()

}

##############################################################################
fn_lag_events01 <- function(events01, scale.doy) {

    ### simpler lagging. unwrap the 2d chain array to 1d
    ## obs
    ndata    <- length(events01$obs$ch.st.l)
    o12      <- array(NA,dim=c(ndata,2))
    o12.age  <- array(NA,dim=c(ndata,2))
    o12.doy  <- array(NA,dim=c(ndata,2))
    o12.cc   <- array(NA,dim=c(ndata,2))
    # day t=t data
    o12[ ,1]     <- events01$obs$ch.st.l[1:ndata]  # unwrap the 2d array to 1d
    o12.age[ ,1] <-  events01$obs$ch.age[1:ndata]
    o12.doy[ ,1] <-  events01$obs$ch.doy[1:ndata]
    o12.cc[ ,1]  <- rep(events01$obs$ch.gmst, each=dim(events01$obs$ch.st.l)[1]) # need to do something different for d1.cc
    # day t=t+1 data
        o12[1:(ndata-1) ,2] <- events01$obs$ch.st.l[2:ndata]
    o12.age[1:(ndata-1) ,2] <-  events01$obs$ch.age[2:ndata]
    o12.doy[1:(ndata-1) ,2] <-  events01$obs$ch.doy[2:ndata]
     o12.cc[            ,2] <- o12.cc[ ,1]
    # filter removing cluster of points near zero, NAs and day before and after which have ages of 0
    # o.igood <- NULL
    o12.age[is.na(o12.age)] <- 0

    o.iday1E1  <- which(o12.age[,1]==1)
    o.iday1G1  <- which(o12.age[,1]>1)
        #   plot(o12[1:100,1],ty='l',ylim=c(0,12))
        # points(o12[1:100,1],pch=3,cex=.3)
        #  lines(o12[1:100,2],col=2)
        # points(o12[1:100,2],pch=3,cex=.3)
        # points(o.iday1G1,o12[o.iday1G1,1],col=4,pch=1,cex=1)
        # points(o.iday1E1,o12[o.iday1E1,1],col=3,pch=1,cex=1)

    ## model
    ndata    <- length(events01$mod$ch.st.l)
    m12      <- array(NA,dim=c(ndata,2))
    m12.age  <- array(NA,dim=c(ndata,2))
    m12.doy  <- array(NA,dim=c(ndata,2))
    m12.cc   <- array(NA,dim=c(ndata,2))
    # day t=t data
    m12[ ,1]     <- events01$mod$ch.st.l[1:ndata]  # unwrap the 2d array to 1d
    m12.age[ ,1] <-  events01$mod$ch.age[1:ndata]
    m12.doy[ ,1] <-  events01$mod$ch.doy[1:ndata]
    m12.cc[ ,1]  <- rep(events01$mod$ch.gmst, each=dim(events01$mod$ch.st.l)[1]) # need to do something different for d1.cc
    # day t=t+1 data
        m12[1:(ndata-1) ,2] <- events01$mod$ch.st.l[2:ndata]
    m12.age[1:(ndata-1) ,2] <-  events01$mod$ch.age[2:ndata]
    m12.doy[1:(ndata-1) ,2] <-  events01$mod$ch.doy[2:ndata]
        m12.cc[            ,2] <- m12.cc[ ,1]
    # filter removing cluster of points near zero, NAs and day before and after which have ages of 0
    m.igood <- NULL
    m12.age[is.na(m12.age)] <- 0

    m.iday1E1  <- which(m12.age[,1]==1)
    m.iday1G1  <- which(m12.age[,1]>1)
        #   plot(m12[1:100,1],ty='l',ylim=c(0,12))
        # points(m12[1:100,1],pch=3,cex=.3)
        #  lines(m12[1:100,2],col=2)
        # points(m12[1:100,2],pch=3,cex=.3)
        # points(m.iday1G1,m12[m.iday1G1,1],col=4,pch=1,cex=1)
        # points(m.iday1E1,m12[m.iday1E1,1],col=3,pch=1,cex=1)

    ## combine
    om12      <- rbind(o12,m12)
    om12.age  <- rbind(o12.age,m12.age)
    om12.sdoy <- rbind(o12.doy/scale.doy$obs, m12.doy/scale.doy$mod)
    ## scale doy via cos((doy - midsummer_doy)/yearlength)
    tmp1.o    <- (o12.doy-median(events01$obs$hotseas$hotSdoy))/scale.doy$obs
    tmp1.m    <- (m12.doy-median(events01$mod$hotseas$hotSdoy))/scale.doy$mod
    om12.cosd <- cos(2*pi*(rbind(tmp1.o,tmp1.m)))
    om12.cc   <- rbind(o12.cc ,m12.cc )
    om12.iobs <- rbind(array(1,dim=dim(o12.cc)), array(0,dim=dim(m12.cc)))
    om12.class<- factor(!om12.iobs[,1],labels=c('obs','mod')) # weired factor thing labels have to be the ranked order of the data hence !

    om.1E1 <- c(o.iday1E1, m.iday1E1 +length(o12[,1]))
    om.1G1 <- c(o.iday1G1, m.iday1G1 +length(o12[,1]))
        #     k1 <- 6500:6800
        #  plot(k1,om12[k1,1],ty='l',ylim=c(0,12))
        # lines(k1,om12[k1,2],col=2)
        # points(om.1E1,om12[om.1E1,1],col=3,pch=1,cex=1)
        # points(om.1G1,om12[om.1G1,1],col='grey70',pch=1,cex=1)
        # points(om.1E1,om12[om.1E1,2],col=4,pch=1,cex=1)

    lag01 <- list()
    lag01$obs <- list(o12 =o12,   o12.age=o12.age,      o12.doy=o12.doy,    o12.cc=o12.cc,                    o.iday1E1=o.iday1E1, o.iday1G1=o.iday1G1)
    lag01$mod <- list(m12 =m12,   m12.age=m12.age,      m12.doy=m12.doy,    m12.cc=m12.cc,                    m.iday1E1=m.iday1E1, m.iday1G1=m.iday1G1)
    lag01$om  <- list(om12=om12, om12.age=om12.age,   om12.sdoy=om12.sdoy, om12.cc=om12.cc, om12.iobs=om12.iobs,
                               om12.class=om12.class, om12.cosd=om12.cosd, om.1E1=om.1E1,       om.1G1=om.1G1)
    return(lag01)

}

##############################################################################
sim.from.ht <- function(ht, x, cov=NULL, nonoise=FALSE){
    # returns:  vector of X(t=1)
    # ht:       HT-model reurned by BveHTDep or BveHTDepCov
    # x:        vector if X(t=0) data
        origz.l <- qlaplace(ht$origz)
        z.samp  <- sample(x=origz.l, size=length(x), replace=TRUE)
        par     <- ht$par
        if(length(par)==4) {
            alpha      <-     par[1]
            beta       <-     par[2]
            if(nonoise) sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * mean(z.samp) else
                        sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * z.samp
        } else {
            alpha      <-     tanh(par[1]+(par[2]*cov))
            beta       <-     tanh(par[3]+(par[4]*cov))
            mu         <-          par[5] +par[6]*cov
            sqrt.sigma <- sqrt(exp(par[7] +par[8]*cov))
            if(nonoise) sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * (mu+sqrt.sigma*mean(z.samp)) else
                        sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * (mu+sqrt.sigma*z.samp)
        }
    sim0
}

get_abms_indexes <- function(ht){
    m       <- ht$model
    nparalp <- length(m$alpcovind)+1
    nparbet <- length(m$betcovind)+1
    nparmu  <- length(m$mucovind)+1
    nparsig <- length(m$sigcovind)+1
    ipa <-  1                        : nparalp      # which fitted parameters correspond to alpha beta mu & sigma
    ipb <- (nparalp+1)               :(nparalp+nparbet)
    ipm <- (nparalp+nparbet+1)       :(nparalp+nparbet+nparmu)
    ips <- (nparalp+nparbet+nparmu+1):(nparalp+nparbet+nparmu+nparsig)
    list(a=ipa,b=ipb,m=ipm,s=ips)
}

sim.from.htHrc <- function(ht, x, cov=NULL, nonoise=FALSE, RETURNPAR=FALSE, PARBOOT=TRUE, nboot){
    # returns:  vector of X(t=1)
    # ht:       HT-model reurned by BveHTDepCovHrc
    # x:        vector if X(t=0) data
    nx      <- length(x)
    origz.l <- qlaplace(ht$origz)
    z.samp  <- sample(x=origz.l, size=nx, replace=TRUE)
    pam     <- ht$par
    m       <- ht$model
    nparalp <- length(m$alpcovind)+1
    nparbet <- length(m$betcovind)+1
    nparmu  <- length(m$mucovind)+1
    nparsig <- length(m$sigcovind)+1
    totpar  <- sum(nparalp,nparbet,nparmu,nparsig)
    ipa <-  1                        : nparalp      # which fitted parameters correspond to alpha beta mu & sigma
    ipb <- (nparalp+1)               :(nparalp+nparbet)
    ipm <- (nparalp+nparbet+1)       :(nparalp+nparbet+nparmu)
    ips <- (nparalp+nparbet+nparmu+1):(nparalp+nparbet+nparmu+nparsig)
    ja <- m$alpcovind
    jb <- m$betcovind
    jm <- m$mucovind
    js <- m$sigcovind
    n.cov <- dim(m$cvte)[2]

    if(any(is.character(c(ja,jb,jm,js)))) {
        if(is.null(dim(cov))) cov.names <- names(cov) else cov.names <- dimnames(cov)[[2]]
        if(!is.null(ja)) ja <- apply(cbind(ja),1,FUN=function(x,covn) {which(covn==x)},covn=cov.names)
        if(!is.null(jb)) jb <- apply(cbind(jb),1,FUN=function(x,covn) {which(covn==x)},covn=cov.names)
        if(!is.null(jm)) jm <- apply(cbind(jm),1,FUN=function(x,covn) {which(covn==x)},covn=cov.names)
        if(!is.null(js)) js <- apply(cbind(js),1,FUN=function(x,covn) {which(covn==x)},covn=cov.names)
    }

    if(length(pam) == 4) { # stationary model
        alpha <- rep(m$alplink(sum(pam[ipa])), nx)
        beta  <- rep(m$betlink(sum(pam[ipb])), nx)
        mu    <- rep( m$mulink(sum(pam[ipm])), nx)
        sigma <- rep(m$siglink(sum(pam[ips])), nx)
    } else {    # covariates needed
       if(is.null(dim(cov))) {  # just one set of covaraites for all x
            if(length(cov)!=n.cov){
                    cat("sim.from.htHrc: wrong number of covariates",cr)
                    cat("Need",n.cov,cr)
                    return(NA)
            } else { # use this one set for every x
                cov1  <- c(1,cov)   # add 1s for the a0, b0, m0 & s0 terms
                ica   <- c(1,ja+1)  # alpha's index's for the correct covariate
                icb   <- c(1,jb+1)
                icm   <- c(1,jm+1)
                ics   <- c(1,js+1)
                alpha <- m$alplink(sum(pam[ipa]*cov1[ica]))
                beta  <- m$betlink(sum(pam[ipb]*cov1[icb]))
                mu    <-  m$mulink(sum(pam[ipm]*cov1[icm]))
                sigma <- m$siglink(sum(pam[ips]*cov1[ics]))
            }
       } else { # unique set of covariates for each x
            if(dim(cov)[2] != dim(m$cvte)[2]){
                cat("sim.from.htHrc: wrong number of covariates",cr)
                cat("Need",dim(m$cvte)[2],cr)
                return(NA)
            }
            if(dim(cov)[1]>1 & dim(cov)[1]!=nx){
                cat("sim.from.htHrc: length of x and cov different",cr)
                cat("Need",nx,cr)
                return(NA)
            }
            cov1  <- cbind(1,cov) # add column of 1s for the a0, b0, m0 & s0 terms
            ica   <- c(1,ja+1)  # alpha's index's for the correct covariate
            icb   <- c(1,jb+1)
            icm   <- c(1,jm+1)
            ics   <- c(1,js+1)
            alpha <- m$alplink(rowSums(t(pam[ipa]*t(cov1[,ica]))))
            beta  <- m$betlink(rowSums(t(pam[ipb]*t(cov1[,icb]))))
            mu    <-  m$mulink(rowSums(t(pam[ipm]*t(cov1[,icm]))))
            sigma <- m$siglink(rowSums(t(pam[ips]*t(cov1[,ics]))))

       }

    }

    sqrt.sigma <- sqrt(sigma)
    if(nonoise) sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * (mu+sqrt.sigma*mean(z.samp)) else
                sim0       <- alpha * pmax(x, 0) + (pmax(x,0)^beta) * (mu+sqrt.sigma*z.samp)

    if(RETURNPAR) {
        return(list(sim=sim0,par=list(alpha=alpha, beta=beta, mu=mu, sigma=sigma)))
    } else {
        return(sim0)
    }
}

boot.par.htHrc <- function(ht, nboot=100, ADDZNOISE=TRUE){
    mod.lev.l <- qlaplace(ht$mod.lev.u)
    origz.l   <- qlaplace(ht$origz)
    lam.pen   <- ht$lam.pen
    bl        <- c(-1, -10, -10^6, 1e-04)
    bu        <- c(1, 1, 10^6, 10^6)

    # boot fit for uncertainty
    nsimb     <- length(origz.l)
    b.par   <- array(0.0,dim=c(nboot,4+1))
    for(i in seq(nboot)){
        sdat         <- htsim(ht$par, origz.l, mod.lev.l, nsimb, zfac=0.025, ADDZNOISE=ADDZNOISE)
        fit1         <- optim(ht$par, cond.llik, method="L-BFGS-B", x=sdat$x,  y=sdat$y, lower=bl, upper=bu, hessian=FALSE)
        b.par[i,1:4] <- fit1$par
        b.par[i,5]   <- fit1$value
    }
    return(b.par)

}

ht.cond.llik <- function(p, x, y) {
    ll <- (length(x)/2) * (log(2 * pi)) + (length(x)/2) * log(p[4]) +
          p[2] * sum(log(x)) + (1/((2 * p[4]))) *
          sum(((y - p[1] * x - (p[3] * x^(p[2])))^2)/(x^(2 * p[2]))) +
          lam.pen*(p[3]^2)

    if(is.finite(ll)) return(ll) else return(10^6)
}

htsim <- function(par, origz.l, th.l, nsim0, zfac=0.025, ADDZNOISE=ADDZNOISE) {
    xn.l   <- th.l + rexp(nsim0)
    zn.l   <- sample(origz.l, nsim0, replace = T)
    if(ADDZNOISE) zn.l <- zn.l + rnorm(nsim0, mean=0,sd=zfac*sd(origz.l))
    yn.l   <- par[1] * xn.l + (xn.l^(par[2])) * zn.l
    cat("htsim is wrong?!!!",cr)
    return(NA)
    # return(list(x=xn.l, y=yn.l, z=zn.l))
}


### Termination models
### GLM
##############################################################################
fn_fitAutoMultiJointTermGLM <- function(events01, scale.doy, stdiag=NULL, stplot=NULL ){

    if(!is.null(stdiag)) conDiag    <- file(stdiag,'w')
    if(!is.null(stdiag)) writeLines(stdiag,con=conDiag)

    # if(!is.null(stdiag)) writeLines(capture.output(print()), con=conDiag)

    lag01 <- fn_lag_events01(events01, scale.doy)
    list2env(lag01$om,env=parent.frame())

    # remove missing values
    i0          <- which(om12.age[,1] > 0)
    iterm       <- which(om12.age[i0,2] == 0)
    fin         <- rep(0, length(i0))
    live        <- rep(1, length(i0))
     fin[iterm] <- 1
    live[iterm] <- 0
    st.c        <- c('I','C','D','A')

   ### day 1 term ################################################################
    i1        <- which(om12.age[i0,1]==1)
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], I=om12.iobs[i0[i1],1], C=om12.cc[i0[i1],1], D=om12.cosd[i0[i1],1], A=om12.age[i0[i1],1])

    # calculate all combinations of ICDA
    all.com <- NULL
    for(i in 1:3) all.com[[i]] <- combn(st.c[1:3],i)
    ncomb <- sum(sapply(all.com,dim)[2,])

    # termination model
    if(!is.null(stdiag)) writeLines("\n### Day 1 Termination model ###", con=conDiag)
    term.1 <- NULL
    fmla.1 <- NULL
    count0 <- 1
    term.1[[count0]]      <- glm(fin ~ tmp, data=term_vars,family="binomial")
    term.1[[count0]]$fmla <- 'fin ~ tmp'
    fmla.1[count0]        <- 'fin ~ tmp'
    count0 <- 2
    for(i in seq_along(all.com)) {
        for(j in seq_along(all.com[[i]][1,])) {
            fmla <- paste('fin ~ tmp',paste('+',all.com[[i]][,j],collapse=' ', sep=''))
            term.1[[count0]]      <- glm(fmla, data=term_vars, family="binomial")
            term.1[[count0]]$fmla <- fmla
            fmla.1[count0]        <- fmla
            cat(fmla, cr)
            # if(!is.null(stdiag)) writeLines(capture.output(print(fmla)), con=conDiag)
            if(!is.null(stdiag)) writeLines(fmla, con=conDiag)
            if(!is.null(stdiag)) writeLines(capture.output(print(round(summary(term.1[[count0]])$coeff,3))), con=conDiag)
            count0                <- count0 +1
        }
    }
    aic.1      <- sapply(term.1,"[[","aic")
    iterm.best <- which.min(aic.1)
    cat("fn_fitJointTerm: best term 1 model:", fmla.1[iterm.best],cr)
    if(!is.null(stdiag)) writeLines(capture.output(cat("fn_fitJointTerm: best term 1 model:", fmla.1[iterm.best],cr)), con=conDiag)

    # living model
    if(!is.null(stdiag)) writeLines("\n### Day 1 Living model ###", con=conDiag)
    live.1 <- NULL
    fmla.1 <- NULL
    count0 <- 1
    live.1[[count0]]      <- glm(live ~ tmp, data=term_vars,family="binomial")
    live.1[[count0]]$fmla <- 'live ~ tmp'
    fmla.1[count0]        <- 'live ~ tmp'
    count0 <- 2
    for(i in seq_along(all.com)) {
        for(j in seq_along(all.com[[i]][1,])) {
            fmla <- paste('live ~ tmp',paste('+',all.com[[i]][,j],collapse=' ', sep=''))
            live.1[[count0]]      <- glm(fmla, data=term_vars, family="binomial")
            live.1[[count0]]$fmla <- fmla
            fmla.1[count0]        <- fmla
            cat(fmla, cr)
            if(!is.null(stdiag)) writeLines(fmla, con=conDiag)
            if(!is.null(stdiag)) writeLines(capture.output(print(round(summary(live.1[[count0]])$coeff,3))), con=conDiag)
            count0                <- count0 +1
        }
    }
    aic.1      <- sapply(live.1,"[[","aic")
    ilive.best <- which.min(aic.1)
    cat("fn_fitJointTerm: best live 1 model:", fmla.1[ilive.best],cr)
    if(!is.null(stdiag)) writeLines(capture.output(cat("fn_fitJointTerm: best live 1 model:", fmla.1[ilive.best],cr)), con=conDiag)

    TermModel <- term.1[[iterm.best]]
    LiveModel <- live.1[[ilive.best]]
    day1 <- list(term=TermModel, live=LiveModel)

   # predict(day1$term, data.frame(tmp=3.2, D=0.5),type="response")

   ### day N term ################################################################
    i1        <- which(om12.age[i0,1]>1)
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], I=om12.iobs[i0[i1],1], C=om12.cc[i0[i1],1], D=om12.cosd[i0[i1],1], A=om12.age[i0[i1],1])

    # calculate all combinations of ICDA
    all.com <- NULL
    for(i in 1:4) all.com[[i]] <- combn(st.c,i)
    ncomb <- sum(sapply(all.com,dim)[2,])

    # termination model
    if(!is.null(stdiag)) writeLines("\n### Day N Termination model ###", con=conDiag)
    term.N <- NULL
    fmla.N <- NULL
    count0 <- 1
    term.N[[count0]]      <- glm(fin ~ tmp, data=term_vars,family="binomial")
    term.N[[count0]]$fmla <- 'fin ~ tmp'
    fmla.N[count0]        <- 'fin ~ tmp'
    count0 <- 2
    for(i in seq_along(all.com)) {
        for(j in seq_along(all.com[[i]][1,])) {
            fmla <- paste('fin ~ tmp',paste('+',all.com[[i]][,j],collapse=' ', sep=''))
            term.N[[count0]]      <- glm(fmla, data=term_vars, family="binomial")
            term.N[[count0]]$fmla <- fmla
            fmla.N[count0]        <- fmla
            if(!is.null(stdiag)) writeLines(fmla, con=conDiag)
            if(!is.null(stdiag)) writeLines(capture.output(print(round(summary(term.N[[count0]])$coeff,3))), con=conDiag)
            count0                <- count0 +1
            cat(fmla, cr)
        }
    }
    aic.N      <- sapply(term.N,"[[","aic")
    iterm.best <- which.min(aic.N)
    cat("fn_fitJointTerm: best term N model:", fmla.N[iterm.best],cr)
    if(!is.null(stdiag)) writeLines(capture.output(cat("fn_fitJointTerm: best term N model:", fmla.N[iterm.best],cr)), con=conDiag)
    # plot the effect of all covariates
        fn_termJointGLMPlot(term.N[[length(term.N)]])

    # living model
    if(!is.null(stdiag)) writeLines("\n### Day N Living model ###", con=conDiag)
    live.N <- NULL
    fmla.N <- NULL
    count0 <- 1
    live.N[[count0]]      <- glm(live ~ tmp, data=term_vars,family="binomial")
    live.N[[count0]]$fmla <- 'live ~ tmp'
    fmla.N[count0]        <- 'live ~ tmp'
    count0 <- 2
    for(i in seq_along(all.com)) {
        for(j in seq_along(all.com[[i]][1,])) {
            fmla <- paste('live ~ tmp',paste('+',all.com[[i]][,j],collapse=' ', sep=''))
            live.N[[count0]]      <- glm(fmla, data=term_vars, family="binomial")
            live.N[[count0]]$fmla <- fmla
            fmla.N[count0]        <- fmla
            if(!is.null(stdiag)) writeLines(fmla, con=conDiag)
            if(!is.null(stdiag)) writeLines(capture.output(print(round(summary(live.N[[count0]])$coeff,3))), con=conDiag)
            count0                <- count0 +1
            cat(fmla, cr)
        }
    }
    aic.N      <- sapply(live.N,"[[","aic")
    ilive.best <- which.min(aic.N)
    cat("fn_fitJointTerm: best live N model:", fmla.N[ilive.best],cr)
    if(!is.null(stdiag)) writeLines(capture.output(cat("fn_fitJointTerm: best live N model:", fmla.N[ilive.best],cr)), con=conDiag)

    TermModel <- term.N[[iterm.best]]
    LiveModel <- live.N[[ilive.best]]
    dayN <- list(term=TermModel, live=LiveModel)
   #
    close(conDiag)

    hw.Jterm <- list(day1=day1, dayN=dayN)

    return(hw.Jterm)

}

##############################################################################
fn_fitJointTermGLM <- function(events01, fmla1N, scale.doy, stdiag=NULL, stplot=NULL, DOPLOT=FALSE ){

    # fmla must be preselected to be whole year or for hotseas as required
    # fmla1N$day1 & fmla1N$dayN

    if(!is.null(stdiag)) conDiag    <- file(stdiag,'w')
    if(!is.null(stdiag)) writeLines(stdiag,con=conDiag)
    # if(!is.null(stdiag)) writeLines(capture.output(print()), con=conDiag)

    if(!is.null(stplot)) pdf(file=stplot, 10, 10)

    lag01 <- fn_lag_events01(events01, scale.doy)
    list2env(lag01$om,env=parent.frame())

    # remove non-relevant values
    i0          <- which(om12.age[,1] > 0)
    iterm       <- which(om12.age[i0,2] == 0)
    fin         <- rep(0, length(i0))
    live        <- rep(1, length(i0))
     fin[iterm] <- 1
    live[iterm] <- 0

    ### day 1 ################################################################
    if(!is.null(stdiag)) {
        writeLines("\n###################################",       con=conDiag)
        writeLines("### Day 1 Termination GLM model ###",       con=conDiag)
    }
    i1        <- which(om12.age[i0,1]==1)
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], class=om12.iobs[i0[i1],1], gmst=om12.cc[i0[i1],1], D=om12.sdoy[i0[i1],1], A=om12.age[i0[i1],1])

    # termination model
    fmla1  <- fmla1N$day1  ##  termGAM.fmla$hseas$day1 ##
    term.1 <- NULL
    if(class(fmla1N$day1)=='formula') {
        # term.1 <- gam(fmla1 , data=term_vars, family='gaussian')
        cat(paste(fmla1), cr)
        term.1      <- glm(fmla1, data=term_vars, family="binomial")
        term.1$fmla <- fmla1
        if(!is.null(stdiag)) {
            # writeLines("\n### Day 1 Termination GLM model ###",               con=conDiag)
            writeLines(capture.output(print(names(fmla1)) ),                     con=conDiag)
            writeLines(capture.output(print(round(summary(term.1)$coeff,3))), con=conDiag)
            writeLines("##########\n",       con=conDiag)
            # writeLines('### End Day 1 Termination GLM model ###\n\n',         con=conDiag)
        }
        if(DOPLOT) {
            up.3by2()
            par(oma=c(0,0,2,0))
            # plot(term.1)
            fn_termJointGLMPlot(term.1, main0="Term Day 1")
            mtext(paste("Day 1:",paste(capture.output(print(fmla1)),collapse='')),outer=TRUE,col=1,line=-0.2) # center justfied
            if(is.null(stplot)) readline("Continue to Day N?")
        }
    } else {
        for(i in seq_along(fmla1)) {
            term.1[[i]] <- glm(fmla1[[i]] , data=term_vars, family='binomial')

            # writeLines("\n### Day 1 Termination GLM model ###",                    con=conDiag)
            writeLines(capture.output(print(names(fmla1)[i]) ),                     con=conDiag)
            writeLines(capture.output(print(round(summary(term.1[[i]])$coeff,3))), con=conDiag)
            writeLines("##########\n",       con=conDiag)
            # writeLines('### End Day 1 Termination GLM model ###\n\n',              con=conDiag)

            if(DOPLOT) {
                up.3by2()
                par(oma=c(0,0,2,0))
                # plot(term.1[[i]])
                fn_termJointGLMPlot(term.1[[i]], main0="Term Day 1")
                mtext(paste("Day 1:",paste(capture.output(print(fmla1[[i]])),collapse='')),outer=TRUE,col=1,line=-0.2) # center justfied
                if(is.null(stplot)) readline("Continue to next/Day N?")
            }
        }
    }

    ### day N ################################################################
    if(!is.null(stdiag)) {
        writeLines("\n###################################",       con=conDiag)
        writeLines("### Day N Termination GLM model ###",       con=conDiag)
    }
    i1        <- which(om12.age[i0,1]>1)
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], class=om12.iobs[i0[i1],1], gmst=om12.cc[i0[i1],1], D=om12.cosd[i0[i1],1], A=om12.age[i0[i1],1])

    # termination model
    fmlaN  <- fmla1N$dayN  ##  termGAM.fmla$hseas$dayN  ##
    term.N <- NULL
    live.N <- NULL
    if(class(fmla1N$dayN)=='formula') {
        # term.N <- gam(fmlaN , data=term_vars, family='gaussian')
        cat(paste(fmlaN), cr)
        term.N      <- glm(fmlaN , data=term_vars, family='binomial')
        term.N$fmla <- fmlaN
        if(!is.null(stdiag)) {
            # writeLines("\n### Day N Termination GLM model ###",               con=conDiag)
            writeLines(capture.output(print(names(fmlaN)) ),                     con=conDiag)
            writeLines(capture.output(print(round(summary(term.N)$coeff,3))), con=conDiag)
            writeLines("##########\n",       con=conDiag)
            # writeLines('### End Day N Termination GLM model ###\n\n',         con=conDiag)
        }
        if(DOPLOT) {
            up.3by2()
            par(oma=c(0,0,2,0))
            # plot(term.N)
            fn_termJointGLMPlot(term.N, main0="Term Day N")
            mtext(paste("Day N:",paste(capture.output(print(fmlaN)),collapse='')),outer=TRUE,col=1,line=-0.2) # center justfied
        }
    } else {
        for(i in seq_along(fmlaN)) {
            term.N[[i]] <- glm(fmlaN[[i]] , data=term_vars, family='binomial')

            # writeLines("\n### Day N Termination GLM model ###",                    con=conDiag)
            writeLines(capture.output(print(names(fmlaN)[i]) ),                     con=conDiag)
            writeLines(capture.output(print(round(summary(term.N[[i]])$coeff,3))), con=conDiag)
            writeLines("##########\n",       con=conDiag)
            # writeLines('### End Day N Termination GLM model ###\n\n',              con=conDiag)

            if(DOPLOT) {
                up.3by2()
                par(oma=c(0,0,2,0))
                # plot(term.N[[i]])
                fn_termJointGLMPlot(term.N[[i]], main0="Term Day N")
                mtext(paste("Day N:",paste(capture.output(print(fmlaN[[i]])),collapse='')),outer=TRUE,col=1,line=-0.2) # center justfied
                if(is.null(stplot)) readline("Continue to next/end?")
            }
        }
    }
    writeLines("###################################",       con=conDiag)

    writeLines("\n\n### AIC/BIC Day 1 Termination GLM model ###",                     con=conDiag)
    for(i in seq_along(fmla1)) writeLines(capture.output(print( paste(names(fmla1[i]),'     ',myround(AIC(term.1[[i]])),myround(BIC(term.1[[i]])))) ),  con=conDiag)
    writeLines("\n\n### AIC/BIC Day N Termination GLM model ###",                     con=conDiag)
    for(i in seq_along(fmlaN)) writeLines(capture.output(print( paste(names(fmlaN[i]),'     ',myround(AIC(term.N[[i]])),myround(BIC(term.N[[i]])))) ),  con=conDiag)

    # live models are exact opposites of fin models    so give nothing new
     # for(i in seq_along(fmlaN)) {
     #     print(summary(live.N[[i]]))
     #     up.3by2()
     #     plot(live.N[[i]])
     #     cat(cr,'##########################',cr,cr)
     #     readline("next?")
     # }


    if(!is.null(stdiag)) close(conDiag)
    if(!is.null(stplot)) dev.off()

    names(term.1) <- names(fmla1N$day1)
    names(term.N) <- names(fmla1N$dayN)
    return(list(day1=term.1, dayN=term.N))
}

##############################################################################
fn_termJointGLMPlot <- function(hw.term, main0='Term Day=>?:'){

    #  DAY > 1  ##############################################################
    # hw.cov <- hw.term$data
    hw.cov <- hw.term$model
    cov.n  <- attr(hw.term$terms,"term.labels")
    cov.q  <- NULL
    for(i in seq_along(cov.n)) {
        if(cov.n[i]=='I') cov.q <- cbind(cov.q, c(0,1,1)) else
        if(cov.n[i]=='A') cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
        cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
    }
    cov.x  <- NULL
    for(i in seq_along(cov.n)) {
        if(cov.n[i]=='I') cov.x[[i]] <- c(0,1) else
        # if(cov.n[i]=='A') cov.x[[i]] <-  else
        cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
    }

    df1   <- NULL
    for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
    df1   <- data.frame(df1)
    df0   <- df1
    df0$I <- 0
    up.2by2()
    for(i in seq_along(cov.n)) {
        if(cov.n[i]=='I') {}  else {
        # if(cov.n[i]=='A')  else {
            df0b             <- df0
            df0b[[cov.n[i]]] <- cov.x[[i]]
            df1b             <- df1
            df1b[[cov.n[i]]] <- cov.x[[i]]
            newy0            <- predict(hw.term,df0b,type="response")
            newy1            <- predict(hw.term,df1b,type="response")
             plot(cov.x[[i]], newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main0,cov.n[i]), ylab='Prob terminate', xlab="Covaraite")
            lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
            grid()
        }
    }
    legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')

}

fn_termJointGLMPlot2 <- function(hw.term1N, st.pdf=NULL){

    if(!is.null(st.pdf)) pdf(st.pdf, width=7, height=7) else x11()

    main0 <- c('Term GLM Day=1:','Term GLM Day>1:')

    for(k in seq_along(hw.term1N)) {

        hw.term <- hw.term1N[[k]] # hw.term1N$day1
        main1   <- main0[k]

        if(class(hw.term)[1]=="glm") {
            main2 <- paste(hw.term$formula)

            hw.cov  <- hw.term$model
            cov.n   <- attr(hw.term$terms,"term.labels")
            cov.q   <- NULL
            for(i in seq_along(cov.n)) {
                if(cov.n[i]=='I') cov.q <- cbind(cov.q, c(0,1,1)) else
                if(cov.n[i]=='A') cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
                cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
            }
            cov.x  <- NULL
            for(i in seq_along(cov.n)) {
                if(cov.n[i]=='I') cov.x[[i]] <- c(0,1) else
                cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
            }

            df1   <- NULL
            for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
            df1   <- data.frame(df1)
            df0   <- df1
            df0$I <- 0
            up.2by2()
            par(mfrow=c(2,2))
            par(mar=c(5,5,2,1))
            for(i in seq_along(cov.n)) {
                if(cov.n[i]=='I') {}  else {
                    df0b             <- df0
                    df0b[[cov.n[i]]] <- cov.x[[i]]
                    df1b             <- df1
                    df1b[[cov.n[i]]] <- cov.x[[i]]
                    newy0            <- predict(hw.term,df0b,type="response")
                    newy1            <- predict(hw.term,df1b,type="response")
                    plot(cov.x[[i]], newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covaraite")
                    lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
                    grid()
                }
            }
            legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
            mtext(paste(main2[c(2,1,3)],collapse=' '), side=1,outer=T,adj=0.05,line=-2,cex=0.8)
        }else {  # list of GAMS

            for(j in seq_along(hw.term)) {

                hw.term2 <- hw.term[[j]]
                main2    <- paste(hw.term2$formula)
                hw.cov   <- hw.term2$model
                cov.n    <- attributes(hw.term2$terms)$term.labels
                cov.q    <- NULL
                for(i in seq_along(cov.n)) {
                    if(cov.n[i]=='I') cov.q <- cbind(cov.q, c(0,1,1)) else
                    if(cov.n[i]=='A') cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
                    cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
                }
                cov.x  <- NULL
                for(i in seq_along(cov.n)) {
                    if(cov.n[i]=='I') cov.x[[i]] <- c(0,1) else
                    cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
                }

                df1   <- NULL
                for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
                df1   <- data.frame(df1)
                df0   <- df1
                df0$I <- 0
                up.2by2()
                par(mfrow=c(2,2))
                par(mar=c(5,5,2,1))
                for(i in seq_along(cov.n)) {
                    if(cov.n[i]=='I') {}  else {
                        df0b             <- df0
                        df0b[[cov.n[i]]] <- cov.x[[i]]
                        df1b             <- df1
                        df1b[[cov.n[i]]] <- cov.x[[i]]
                        newy0            <- predict(hw.term2,df0b,type="response")
                        newy1            <- predict(hw.term2,df1b,type="response")
                        plot(cov.x[[i]], newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covaraite")
                        lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
                        grid()
                    }
                }
                legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
                mtext(paste(main2[c(2,1,3)],collapse=' '), side=1,outer=T,adj=0.05,line=-2,cex=0.8)
            }
        }
    }
    if(!is.null(st.pdf)) dev.off()

}


### GAM Termination models
##############################################################################
fn_fitJointTermGAM <- function(events01, fmla1N, scale.doy, stdiag=NULL, stplot=NULL, DOPLOT=FALSE  ){

    # fmla must be preselected to be whole year or for hotseas as required
    # fmla1N$day1 & fmla1N$dayN

    if(!is.null(stdiag)) conDiag    <- file(stdiag,'w')
    if(!is.null(stdiag)) writeLines(stdiag,con=conDiag)
    # if(!is.null(stdiag)) writeLines(capture.output(print()), con=conDiag)
    if(!is.null(stplot)) pdf(file=stplot, 10, 10)

    lag01 <- fn_lag_events01(events01, scale.doy)
    list2env(lag01$om,env=parent.frame())

    # remove non-relevant values
    i0          <- which(om12.age[,1] > 0)
    iterm       <- which(om12.age[i0,2] == 0)
    fin         <- rep(0, length(i0))
    live        <- rep(1, length(i0))
     fin[iterm] <- 1
    live[iterm] <- 0

    ### day 1 ################################################################
    if(!is.null(stdiag)) {
        writeLines("\n###################################",       con=conDiag)
        writeLines("### Day 1 Termination GAM model ###",       con=conDiag)
    }
    i1        <- which(om12.age[i0,1]==1)
    # D==DOY & A==Age
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], class=om12.class[i0[i1]], gmst=om12.cc[i0[i1],1], D=om12.sdoy[i0[i1],1], A=om12.age[i0[i1],1])

    # termination model
    fmla1  <- fmla1N$day1  ##  termGAM.fmla$hseas$day1 ##
    term.1 <- NULL
    if(class(fmla1N$day1)=='formula') {
        term.1        <- gam(fmla1 , data=term_vars, family='binomial')
        # names(term.1) <- names(fmla1)
        if(!is.null(stdiag)) {
            # writeLines("\n### Day 1 Termination GAM model ###",       con=conDiag)
            writeLines(capture.output(print(names(fmla1)) ),                  con=conDiag)
            writeLines(capture.output(print(summary(term.1))),        con=conDiag)
            writeLines("##########\n",       con=conDiag)
            # writeLines('### End Day 1 Termination GAM model ###\n\n', con=conDiag)
        }
        # if(DOPLOT) {
        #     up.3by2()
        #     par(oma=c(0,0,2,0))
        #     plot(term.1)
        #     mtext(paste("Day 1:",paste(capture.output(print(fmla1)),collapse='')),outer=TRUE,col=1,line=-0.2,cex=0.7) # center justfied
        #     if(is.null(stplot)) readline("Continue to Day N?")
        # }
    } else {
        for(i in seq_along(fmla1)) {
            term.1[[i]]        <- gam(fmla1[[i]] , data=term_vars, family='binomial')
            # names(term.1[[i]]) <- names(fmla1[[i]])

            if(!is.null(stdiag)) {
                # writeLines("\n### Day 1 Termination GAM model ###",                     con=conDiag)
                writeLines(capture.output(print(names(fmla1)[i]) ),                     con=conDiag)
                writeLines(capture.output(print(summary(term.1[[i]]))),                 con=conDiag)
                writeLines(capture.output(print('AIC    BIC') ),                        con=conDiag)
                writeLines(capture.output(print(c(AIC(term.1[[i]]),BIC(term.1[[i]]))) ),con=conDiag)
                writeLines("##########\n",       con=conDiag)
                # paste(myround(lapply(termGAM$day1, AIC)),myround(lapply(termGAM$day1, BIC)))
                # writeLines('### End Day 1 Termination GAM model ###\n\n',               con=conDiag)
            }
            if(DOPLOT) {
                up.3by2()
                par(oma=c(0,0,2,0))
                plot(term.1[[i]])
                mtext(paste("Day 1:",paste(capture.output(print(fmla1[[i]])),collapse='')),outer=TRUE,col=1,line=-0.2,cex=0.7) # center justfied
                if(is.null(stplot)) readline("Continue to next/Day N?")
            }
        }
    }

    ### day N ################################################################
    if(!is.null(stdiag)) {
        writeLines("\n###################################",       con=conDiag)
        writeLines("### Day N Termination GAM model ###",       con=conDiag)
    }
    i1        <- which(om12.age[i0,1]>1)
    # D==DOY & A==Age
    term_vars <- data.frame(fin=fin[i1], live=live[i1], tmp=om12[i0[i1],1], class=om12.class[i0[i1]], gmst=om12.cc[i0[i1],1], D=om12.sdoy[i0[i1],1], A=om12.age[i0[i1],1])

    # termination model
    fmlaN  <- fmla1N$dayN  ##  termGAM.fmla$hseas$dayN  ##
    term.N <- NULL
    live.N <- NULL
    if(class(fmla1N$dayN)=='formula') {
        term.N        <- gam(fmlaN , data=term_vars, family='binomial')
        # names(term.N) <- names(fmlaN)

        if(!is.null(stdiag)) {
            # writeLines("\n### Day N Termination GAM model ###",       con=conDiag)
            writeLines(capture.output(print(names(fmlaN)) ),          con=conDiag)
            writeLines(capture.output(print(summary(term.N))),        con=conDiag)
            writeLines("##########\n",       con=conDiag)
        # writeLines('### End Day N Termination GAM model ###\n\n', con=conDiag)
        }
        # if(DOPLOT) {
        #     up.3by2()
        #     par(oma=c(0,0,2,0))
        #     plot(term.N)
        #     mtext(paste("Day N:",paste(capture.output(print(fmlaN)),collapse='')),outer=TRUE,col=1,line=-0.2,cex=0.7) # center justfied
        # }
    } else {
        for(i in seq_along(fmlaN)) {
            term.N[[i]]        <- gam(fmlaN[[i]] , data=term_vars, family='binomial')
            # names(term.N[[i]]) <- names(fmlaN[[i]])

            if(!is.null(stdiag)) {
                # writeLines("\n### Day N Termination GAM model ###",       con=conDiag)
                writeLines(capture.output(print(names(fmlaN)[i]) ),                     con=conDiag)
                writeLines(capture.output(print(summary(term.N[[i]]))),                 con=conDiag)
                writeLines(capture.output(print('AIC    BIC') ),                        con=conDiag)
                writeLines(capture.output(print(c(AIC(term.N[[i]]),BIC(term.N[[i]]))) ),con=conDiag)
                writeLines("##########\n",       con=conDiag)
                # writeLines('### End Day N Termination GAM model ###\n\n', con=conDiag)
            }
            if(DOPLOT) {
                up.3by2()
                par(oma=c(0,0,2,0))
                plot(term.N[[i]])
                mtext(paste("Day N:",paste(capture.output(print(fmlaN[[i]])),collapse='')),outer=TRUE,col=1,line=-0.2,cex=0.7) # center justfied
                if(is.null(stplot)) readline("Continue to next/end?")
            }
        }
    }
    writeLines("###################################",       con=conDiag)

    writeLines("\n\n### AIC/BIC Day 1 Termination GAM model ###",                     con=conDiag)
    for(i in seq_along(fmla1)) writeLines(capture.output(print( paste(names(fmla1[i]),'     ',myround(AIC(term.1[[i]])),myround(BIC(term.1[[i]])))) ),  con=conDiag)
    writeLines("\n\n### AIC/BIC Day N Termination GAM model ###",                     con=conDiag)
    for(i in seq_along(fmlaN)) writeLines(capture.output(print( paste(names(fmlaN[i]),'     ',myround(AIC(term.N[[i]])),myround(BIC(term.N[[i]])))) ),  con=conDiag)
    # writeLines(capture.output(print() ),                        con=conDiag)
    # writeLines(capture.output(print() ),                        con=conDiag)

    # live models are exact opposites of fin models    so give nothing new
    # for(i in seq_along(fmlaN)) {
    #     print(summary(live.N[[i]]))
    #     up.3by2()
    #     plot(live.N[[i]])
    #     cat(cr,'##########################',cr,cr)
    #     readline("next?")
    # }

    # this all worksbut nor sure where it leaves us.  still rising term prob with high temperatures.
    # no other strong dependence

    if(!is.null(stdiag)) close(conDiag)
    if(!is.null(stplot)) dev.off()

    names(term.1) <- names(fmla1N$day1)
    names(term.N) <- names(fmla1N$dayN)
    return(list(day1=term.1, dayN=term.N))
}

##############################################################################
fn_termJointGAMPlot <- function(hw.term, main0='Term Day=>? ~', st.pdf=NULL){

    if(!is.null(st.pdf)) pdf(st.pdf, width=7, height=7) else x11()

    if(class(hw.term)[1]!="gam") {
        cat("fn_termJointGAMPlot: must be a GAM object",cr)
        return(NULL)
    }
    # hw.cov <- hw.term$data
    hw.cov <- hw.term$model
    cov.n  <- attributes(hw.term$terms)$term.labels
    if(any(grepl('class',cov.n))) {
        ISCLASS <- TRUE
        cov.n  <- cov.n[cov.n!='class']
    } else ISCLASS <- FALSE

    cov.q  <- NULL
    cov.x  <- NULL
    for(i in seq_along(cov.n)) {
        # if(cov.n[i]=='class') cov.q <- cbind(cov.q, c(0,1,1)) else
        if(cov.n[i]=='A')   cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
        cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
        # if(cov.n[i]=='class') cov.x[[i]] <- c(0,1) else
        # if(cov.n[i]=='A') cov.x[[i]] <-  else
        cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
    }

    df1   <- NULL
    for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
    df1   <- data.frame(df1)
    df1$A <- 5
    df0   <- df1
    up.2by3()
    if(ISCLASS) {
        df1$class <- 'obs'
        df0$class <- 'mod'
        for(i in seq_along(cov.n)) {
            # if(cov.n[i]=='class') {}  else {
            # if(cov.n[i]=='A')  else {
                df0b             <- df0
                df0b[[cov.n[i]]] <- cov.x[[i]]
                df1b             <- df1
                df1b[[cov.n[i]]] <- cov.x[[i]]
                newy0            <- predict(hw.term,df0b,type="response")
                newy1            <- predict(hw.term,df1b,type="response")
                plot(cov.x[[i]],  newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main0,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
                lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
                grid()
            # }
        }
        legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
    } else {
        df1$class <- 'obs'
        for(i in seq_along(cov.n)) {
            df1b             <- df1
            df1b[[cov.n[i]]] <- cov.x[[i]]
            newy1            <- predict(hw.term,df1b,type="response")
            plot(cov.x[[i]], newy1, ylim=c(0,1), ty='l', lwd=2, main=paste(main0,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
            grid()
        }
        legend('bottomleft',c('Obs'),col=1, pch=NA, lty=1, lwd=2, bty='n')
    }

    up.3by2()
    par(oma=c(0,0,2,0))
    plot(hw.term)
    # mtext(paste("Day 1:",paste(capture.output(print(fmla1)),collapse='')),outer=TRUE,col=1,line=-0.2,cex=0.7) # center justfied

    if(!is.null(st.pdf)) dev.off()

}

##############################################################################
fn_termJointGAMPlot2 <- function(hw.term1N, st.pdf=NULL){

    if(!is.null(st.pdf)) pdf(st.pdf, width=7, height=7) else x11()

    main0   <- c('Term GAM Day=1:','Term GAM Day>1:')

    for(k in seq_along(hw.term1N)) {

        hw.term <- hw.term1N[[k]]
        main1   <- main0[k]

        if(class(hw.term)[1]=="gam") {
            main2 <- paste(paste(paste(hw.term$formula)[1:2],collapse=''),paste(hw.term$formula)[3],sep='=')

            hw.cov <- hw.term$model
            cov.n  <- attributes(hw.term$terms)$term.labels
            if(any(grepl('class',cov.n))) {
                ISCLASS <- TRUE
                cov.n  <- cov.n[cov.n!='class']
            } else ISCLASS <- FALSE
            cov.q  <- NULL
            cov.x  <- NULL
            for(i in seq_along(cov.n)) {
                if(cov.n[i]=='A')   cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
                cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
                cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
            }
            df1   <- NULL
            for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
            df1   <- data.frame(df1)
            df1$A <- 5
            df0   <- df1
            up.2by3()
            par(mar=c(5,5,2,1))
            if(ISCLASS) {
                df1$class <- 'obs'
                df0$class <- 'mod'
                for(i in seq_along(cov.n)) {
                        df0b             <- df0
                        df0b[[cov.n[i]]] <- cov.x[[i]]
                        df1b             <- df1
                        df1b[[cov.n[i]]] <- cov.x[[i]]
                        newy0            <- predict(hw.term,df0b,type="response")
                        newy1            <- predict(hw.term,df1b,type="response")
                        plot(cov.x[[i]],  newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
                        lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
                        grid()
                }
                legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
            } else {
                df1$class <- 'obs'
                for(i in seq_along(cov.n)) {
                    df1b             <- df1
                    df1b[[cov.n[i]]] <- cov.x[[i]]
                    newy1            <- predict(hw.term,df1b,type="response")
                    plot(cov.x[[i]], newy1, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
                    grid()
                }
                legend('bottomleft',c('Obs'),col=1, pch=NA, lty=1, lwd=2, bty='n')
            }
            mtext(main2, side=1,outer=T,adj=0.2,line=-2,cex=0.8)
            up.3by2()
            par(mar=c(5,5,2,1))
            par(oma=c(0,0,2,0))
            plot(hw.term)
            mtext(main2, side=1,outer=T,adj=0.2,line=-2,cex=0.8)
        } else {  # list of GAMS
            for(j in seq_along(hw.term)) {

                hw.term2 <- hw.term[[j]]
                main2 <- paste(paste(paste(hw.term2$formula)[1:2],collapse=''),paste(hw.term2$formula)[3],sep='=')
                hw.cov <- hw.term2$model
                cov.n  <- attributes(hw.term2$terms)$term.labels
                if(any(grepl('class',cov.n))) {
                    ISCLASS <- TRUE
                    cov.n  <- cov.n[cov.n!='class']
                } else ISCLASS <- FALSE
                cov.q  <- NULL
                cov.x  <- NULL
                for(i in seq_along(cov.n)) {
                    if(cov.n[i]=='A')   cov.q <- cbind(cov.q, c(min(hw.cov[[cov.n[i]]],na.rm=T),trunc(max(hw.cov[[cov.n[i]]],na.rm=T)/2),max(hw.cov[[cov.n[i]]],na.rm=T))) else
                    cov.q <- cbind(cov.q, quantile( hw.cov[[cov.n[i]]] ,c(0.1,0.5, 0.9), na.rm=T))
                    cov.x[[i]] <- seq(min(hw.cov[[cov.n[i]]],na.rm=T), max(hw.cov[[cov.n[i]]],na.rm=T), length=100)
                }
                df1   <- NULL
                for(i in seq_along(cov.n)) df1[[cov.n[i]]] <- rep(cov.q[2,i],100)
                df1   <- data.frame(df1)
                df1$A <- 5
                df0   <- df1
                up.2by3()
                par(mar=c(5,5,2,1))
                if(ISCLASS) {
                    df1$class <- 'obs'
                    df0$class <- 'mod'
                    for(i in seq_along(cov.n)) {
                            df0b             <- df0
                            df0b[[cov.n[i]]] <- cov.x[[i]]
                            df1b             <- df1
                            df1b[[cov.n[i]]] <- cov.x[[i]]
                            newy0            <- predict(hw.term2,df0b,type="response")
                            newy1            <- predict(hw.term2,df1b,type="response")
                            plot(cov.x[[i]],  newy0, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
                            lines(cov.x[[i]], newy1, col=2, lwd=2, lty=3)
                            grid()
                    }
                    legend('bottomleft',c('Model','Obs'),col=1:2, pch=NA, lty=c(1,3), lwd=2, bty='n')
                } else {
                    df1$class <- 'obs'
                    for(i in seq_along(cov.n)) {
                        df1b             <- df1
                        df1b[[cov.n[i]]] <- cov.x[[i]]
                        newy1            <- predict(hw.term2,df1b,type="response")
                        plot(cov.x[[i]], newy1, ylim=c(0,1), ty='l', lwd=2, main=paste(main1,cov.n[i]), ylab='Prob terminate', xlab="Covariate")
                        grid()
                    }
                    legend('bottomleft',c('Obs'),col=1, pch=NA, lty=1, lwd=2, bty='n')
                }
                mtext(main2, side=1,outer=T,adj=0.2,line=-2,cex=0.8)
                up.3by2()
                par(mar=c(5,5,2,1))
                par(oma=c(0,0,2,0))
                plot(hw.term2)
                mtext(main2, side=1,outer=T,adj=0.2,line=-2,cex=0.8)

            }
        }

    }

    if(!is.null(st.pdf)) dev.off()

}

##############################################################################
fn_ <- function(){

}
