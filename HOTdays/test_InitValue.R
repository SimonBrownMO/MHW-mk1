# source("test_InitValue.R")

###  RESULTS
### 2025.03.21
### Both qqplots of all doy and quantile plots wrt doy both indicate that the IV model is working fine.

library(glue)
source("/home/h03/hadsx/extremes/R/Rutils/sjb_colours.R")
tblack <- transcol('black',percent=80)
tred   <- transcol('red',percent=50)

# MSref.file <-                  "/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/HadCM3/v2b/England/01/MakeStationary/MSref/ukgd_cpm85_5k_x2y18.MSref.2025-03-03.RData"
# Par.file <- gsub('YYY',simyear,"/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/HadCM3/v2b/England/01/SimHD/Parameters/ukgd_cpm85_5k_x2y18.YYY-YYY-112-295-0.94.RData")
st.pwd     <- system("pwd", intern=TRUE)


if(FALSE) {
    list2env(HDconfig$initV , envir = .GlobalEnv)

    load(st_events,   verb=TRUE)
    load(st_msdata01, verb=TRUE)

    gpd.th.u  <- initV$th.u

    ### ensure the doy basis fn is consistent with hotseason if present
    if(events01$info$USEHOTSEASON & events01$hotseas$IsHotSeason) {
        fmla.IVgpd <- initV$fmla.hseas.IVgpd
    } else {
        fmla.IVgpd <- initV$fmla.ally.IVgpd
    }

    allIVfits <- list()
    for(i in seq_along(fmla.IVgpd)) {
        allIVfits[[i]]      <- fn_fitInitVal(events01, gpd.th.u, fmla.IVgpd[[i]], IsJoint=TRUE)
        allIVfits[[i]]$name <- names(fmla.IVgpd)[i]
        cat("Fitted IV",i,names(fmla.IVgpd)[i],cr)
    }

    ### filter out bad fits
    isFitGood         <- NULL
    for(i in seq_along(fmla.IVgpd)) {
        x1 <- summary(allIVfits[[i]])
        pval <- c( x1[[1]]$logscale$'Pr(>|t|)', x1[[2]]$logscale$'Pr(>|t|)', x1[[2]]$shape$'Pr(>|t|)', x1[[1]]$shape$'Pr(>|t|)')
        if(all(pval < 1e-6)) {isFitGood[i] <- FALSE; cat(i,names(fmla.IVgpd)[i],pval,cr)} else isFitGood[i] <- TRUE
    }
    allIVfits  <-  allIVfits[isFitGood]
    fmla.IVgpd <- fmla.IVgpd[isFitGood]

    # ### save statistics of each fit
    # st.InitVal.stats <- paste(dirname(HDconfig$initV$st_InitValue),'diag',sub('.RData','.stats',basename(HDconfig$files$st_InitValue)),sep='/')
    # for(i in seq_along(allIVfits)) {
    #     if(i==1) fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=TRUE) else
    #                 fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=FALSE)
    # }

    ### choose best fit
    allIVfits.aic        <- lapply(allIVfits, AIC)
    names(allIVfits.aic) <- names(fmla.IVgpd)
    cat("AIC",cr)
    for(i in seq_along(fmla.IVgpd)) cat(names(fmla.IVgpd)[i],allIVfits.aic[[i]],cr)

    for(i in seq_along(fmla.IVgpd)) {
        cat(cr,'################',cr,names(fmla.IVgpd)[i],cr)
        print(summary(allIVfits[[i]]))
    }

    # I dont trust the AIC as it does not pick out the failed fits
    cat(cr,cr,"##############################",cr,"Best IV GPD fit",cr)
    ix1 <- which.min(allIVfits.aic)
    cat(names(fmla.IVgpd)[ix1],allIVfits.aic[[ix1]],cr)
    chosen.IVgpd.name <- names(fmla.IVgpd)[ix1]
    chosen.IVgpd      <- allIVfits[[ix1]]

    fn_plotInitVal(chosen.IVgpd, chosen.IVgpd$fit.data, xplot=1, retp=NULL, st.pdf=NULL)


    ##############################################################
    ### by hand
    ##############################################################

    d1.o <- events01$obs$ch.st.l[2,] # chains start the day before heatwave
    d1.m <- events01$mod$ch.st.l[2,] # chains start the day before heatwave

    sdoy.o <- events01$obs$ch.doy[2,]/scale.doy$obs
    sdoy.m <- events01$mod$ch.doy[2,]/scale.doy$mod

    gmst.o <- events01$obs$ch.gmst
    gmst.m <- events01$mod$ch.gmst

    class.o <- rep('obs',length(d1.o))
    class.m <- rep('mod',length(d1.m))

    # evdata01        <- data.frame(x=c(d1.o,d1.m), sdoy=c(sdoy.o,sdoy.m), gmst=c(gmst.o,gmst.m), class=c(class.o,class.m), iv.k=1, doy=1, gmst=1)
    evdata01        <- data.frame(x=c(d1.o,d1.m), sdoy=c(sdoy.o,sdoy.m), gmst=c(gmst.o,gmst.m), class=c(class.o,class.m), iv.k=1, gmst=1)
    evdata01$excess <- evdata01$x - gpd.th.l
    iexcess         <- which(evdata01$excess >0)

    evdata01.ex     <- subset(evdata01, excess > 0)

    plot(evdata01.ex$x, pch=3, cex=.3)
    plot(evdata01.ex$excess, pch=3, cex=.3)
    plot(evdata01.ex$sdoy, evdata01.ex$excess, pch=3, cex=.3)
    im <- which(evdata01.ex$class=='mod')
    plot(evdata01.ex$sdoy, evdata01.ex$excess, ty='n')
    points(evdata01.ex$sdoy[im],  evdata01.ex$excess[im], pch=3, cex=.3, col=2)
    points(evdata01.ex$sdoy[-im], evdata01.ex$excess[-im], pch=3, cex=.3)


    # fmla.IVgpd[[1]]
    # [[1]] excess ~ class + s(sdoy, bs = "tp", k = iv.k$sdoy, by = class)
    # [[2]] ~1

    doclass    <- factor(c(1,2),labels=c('obs','mod'))[2]
    data01.sim <- data.frame(sdoy=0.52, isobs=0, iv.k=1, class=doclass)
    IVgpd.par  <- predict(allIVfits[[1]], data01.sim, type="response")

    gpd.init.th.l  <- qlaplace(initV$th.u)
    hw.simthresh.l <- qlaplace(new.crit.lev.u)

    newscale <- IVgpd.par[1] + IVgpd.par[2]*(hw.simthresh.l-gpd.init.th.l) # compensate for a different critical threshold if present
    IVgpd3   <- unlist(c(newscale,  IVgpd.par[2], hw.simthresh.l))
    init.sim <- rgpd(1e4,loc=IVgpd3[3], scale=IVgpd3[1], shape=IVgpd3[2])

    qqplot(events01$mod$ch.st.l[2,], init.sim)

}

do.doy1        <-  0 # 200 # doy start : ==0 if using whole hot season
do.doy2        <-  0 # 200 # doy end   : ==0 if using whole hot season
nyears2sim     <-  1e5 # 100000  #000
simyear        <- 3950 # 5700  # 3950 7533
stversion      <- 'v2c'
new.crit.lev.u <- 0.94

MSref.file <- glue("/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/HadCM3/{stversion}/England/01/MakeStationary/MSref/ukgd_cpm85_5k_x2y18.MSref.2025-03-20.RData")
source(paste(st.pwd,"/../setup_all.R",sep=''))
source(paste(st.pwd,"/setup_HotDays.R",sep=''))

Par.file   <- glue("/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/HadCM3/{stversion}/England/01/SimHD/Parameters/ukgd_cpm85_5k_x2y18.{simyear}-{simyear}-112-295-0.94.RData")
load(Par.file, verb=TRUE)

load(st_events, verb=TRUE)

# if(do.doy1==0) do.doy <-sort(unique(events01$hotseas$hotSdoy$o)) else
# if(do.doy1==0) do.doy <-sort(unique(events01$mod$ch.doy[2,])) else
if(do.doy1==0) {
    do.doy.i <- sort(unique(events01$mod$ch.doy[2,]))
    i2       <- indexforcommon(do.doy.i, sim.par$do.doy )
    do.doy   <- do.doy.i[i2[,1]]
} else do.doy <- do.doy1:do.doy2

iy1 <- grep(simyear,names(sim.par))  # for par.sim

data01.sim <- sim.par[[iy1]]$data01.sim
gpdpar     <- sim.par[[iy1]]$gpdpar
th.gpd.sim <- sim.par[[iy1]]$th.gpd.sim
Ipr        <- unlist(sim.par[[iy1]]$Ipr)
IVgpd      <- array(unlist(sim.par[[iy1]]$IVgpd), dim=c(2,length(sim.par[[iy1]]$IVgpd)))

# nsimperdoy     <- nyears2sim * (1 - new.crit.lev.u)  ?wrong this is number of days that are above threshold not number of events STARTING
nsimperdoy     <- nyears2sim * Ipr
hw.simthresh.l <- qlaplace(new.crit.lev.u)
# gpd.init.th.l  <- qlaplace(chosen.IVgpd$gpd.th.u)
gpd.init.th.l  <- qlaplace(0.94) # sim.par[[paste('y',simyear[1],sep='')]]$gpd.init.th.l
delta.sdoy     <- diff(data01.sim$sdoy)[1]

init.sim  <- list()
init.doy  <- list()
for(ddoy in do.doy) {

    idoy       <- which(data01.sim$doy == ddoy)
    indat      <- data01.sim[idoy,]
    delta.sdoy <- diff(data01.sim$sdoy)[1]

    ###### number of events starting on this day #####################################
    # Ipr[idoy] <- predict(IE.model,data01.sim[idoy,],type="response") # prob of event starting this day
    nev       <- length(which(rbinom(nyears2sim,1,Ipr[idoy])==1))                  # total number of events starting this day over nyears2sim years
    if(nev>nyears2sim) stop("PR.INIT gone wrong, nev too large. Stopping!")

    ###### intialise each event with suitable temperature #####################################
    # IVgpd    <- predict(IVgpd.model, indat, type="response")
    newscale <- IVgpd[1,idoy] + IVgpd[2,idoy]*(hw.simthresh.l-gpd.init.th.l) # compensate for a different critical threshold if present
    IVgpd3   <- unlist(c(newscale,  IVgpd[2,idoy], hw.simthresh.l))
    init.sim[[idoy]] <- rgpd(nev,loc=IVgpd3[3], scale=IVgpd3[1], shape=IVgpd3[2])
    init.doy[[idoy]] <- rep(ddoy,nev)

} # end for each doy selected
cat("Completed Done-1 doy:",ddoy, idoy,cr)

init.sim.all <- unlist(init.sim)
init.doy.all <- unlist(init.doy)

## qqplot of all doy
i1 <- which(events01$mod$ch.st.l[2,]<=max(init.sim.all) & events01$mod$ch.doy[2,] %in% do.doy)
qqplot(events01$mod$ch.st.l[2,i1], sample(init.sim.all,1e4),xlab='model day 1 input',ylab='model day 1 sim', main='Laplace stationary margnins')

## plot wrt doy
    # add gam quantiles
    fmlai.Q1 <- list(x ~ s(doy, k=-1) , ~ s(doy))
    q0       <- c(0.2,0.8,0.95)
    inna <- which(!is.na(events01$mod$ch.st.l[2,i1]))
    if(length(inna)>1e5) inna <- sample(inna,1e5)
    in.q <- mqgam(fmlai.Q1, data.frame(x=events01$mod$ch.st.l[2,i1][inna], doy=events01$mod$ch.doy[2,i1][inna]), q0)
    in.q.l <- list()
    in.q.l <- qdo(in.q, q0, predict, newdata=data.frame(doy=do.doy))

    snna <- which(!is.na(init.sim.all))
    if(length(snna)>1e5) snna <- sample(snna,1e5)
    si.q <- mqgam(fmlai.Q1, data.frame(x=init.sim.all[snna], doy=init.doy.all[snna]), q0)
    si.q.l <- list()
    si.q.l <- qdo(si.q, q0, predict, newdata=data.frame(doy=do.doy))

plot(events01$mod$ch.doy[2,], events01$mod$ch.st.l[2,], pch=3, cex=.3, ty='n',xlim=c(90,300),xlab='DOY',ylab='day 1 ', main='Laplace stationary margins')
points(init.doy.all, init.sim.all, pch=2, cex=.3,col=2)
points(events01$mod$ch.doy[2,i1], events01$mod$ch.st.l[2,i1], pch=3, cex=.3)
    for(i in 1:3) lines(do.doy, in.q.l[[i]], col=5,lwd=2)
    for(i in 1:3) lines(do.doy, si.q.l[[i]], col='orange',lty=2,lwd=3)
legend('topright',c('In','In','sim','sim'),col=c(1,5,2,'orange'),lty=c(NA,1,NA,2),pch=c(3,NA,2,NA))









#
