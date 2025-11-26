# source("def_qgam_byyear.R")

### Conclusions
#
# 1) gamma smoothing is a bad idea as it affected the doy term too
# 2) for the comparision with Hobday use 
#           fmla.qgam  <- list(x ~ s(sdoy, bs="cc", k=12) +stime, ~ 1 ) 
#       ie linear term for time
# 
# 3) for interaction term tried
# f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, bs="cc", k=1, by=gmst),      ~ 1 )  # all about the same
# f.qgam.DLGi2 <- list(x ~ s(sdoy)                +s(sdoy, bs="cc", k=12,by=gmst),      ~ 1 )  # all about the same
# f.qgam.DLGi3 <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy,               by=gmst),      ~ 1 )  # all about the same
#   all give similar results, with k=12 the interation term is very wiggly, with k=1 it becomes a sinusoid
# not clear if this is real or just fitting to interannual noise.


# 4) 2025.11.17
#       - multi step approach seems to be the way to go
#           - fit and remove gmst effect 
#               f.qgam.DLG  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
#           - Use summer only
#           - STRAND 1 fit and remove year to year effect (using 90th quantile to capture summer variability)
#           - STRAND 2 fit and remove year to year effect (using 90th quantile to capture summer variability)

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
source("../../libs/fn_MakeStationary.R")
source("../../libs/fn_HotDays.R")



# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v3/MSref/ostia_cdr_nrt_regions.MSref.2025-11-04.RData"
st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-20.RData"
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"

load(st_msref, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)


#   argGam: A list of parameters to be passed to ‘mgcv::gam’. This list
#           can potentially include all the arguments listed in ‘?gam’,
#           with the exception of ‘formula’, ‘family’ and ‘data’.

data01$year <- trunc(data01$time)

### simple linear gmst effect as baseline

q3         <- c(0.1,0.5,0.9)
do.ptiles  <- q3
f.qgam.DLG  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
qgam.DLG    <- mqgam(f.qgam.DLG, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))
q.DLG       <- qdo(qgam.DLG, q3, predict, newdata=data01 )

# f.qgam.DLGi2  <- list(x ~ gmst +s(gmst, bs="cc", k=4,by=sdoy), ~ 1 )
# qgam.DLGi2    <- mqgam(f.qgam.DLGi2, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))

# data01$fdoy <- factor(data01$doy)
# f.qgam.fDLGi2  <- list(x ~ gmst +s(gmst, bs="tp", k=4,by=fdoy), ~ 1 )  ### fails

# data01$cdoy <- cos(2*pi*data01$sdoy)
# f.qgam.cDLGi2  <- list(x ~ gmst +s(gmst, bs="tp", k=4,by=cdoy), ~ 1 )  ### this is quite...BAD

plot(data01$time, data01$x, pch=20, cex=.3)
points(q.DLG[[3]] ~ data01$time, pch=20, cex=.3,col=2)
grid()
readline("Continue 1 ?")

# does very high k allow interannual variability to be captured?  Answer: NO
f.qgam.D64LG  <- list(x ~ s(sdoy, bs="cc", k=128) +gmst, ~ 1 )
qgam.D64LG    <- mqgam(f.qgam.D64LG, data=data01, qu=q3)
q.D64LG       <- qdo(qgam.D64LG, q3, predict, newdata=data01 )
plot(data01$time, data01$x, pch=20, cex=.3)
points(q.D64LG[[3]] ~ data01$time, pch=20, cex=.3,col=2)
# points(q.cDLGi2 ~ data01$time, pch=20, cex=.3,col=3)
grid()
readline("Continue 1b ?")

# ## trend with time
# nd <- data01
# nd$sdoy <- 0.6202186
# q.DLGi2       <- qdo(qgam.DLGi2, q3, predict, newdata=nd )
# plot(q.DLGi2[[2]] ~ data01$time, pch=20, cex=.3,ylim=c(8,14))
# nd$sdoy <- 0.1775956
# q.DLGi2       <- qdo(qgam.DLGi2, q3, predict, newdata=nd )
# points(q.DLGi2[[2]] ~ data01$time, pch=20, cex=.3,col=2)
# grid()
# readline("Continue 2 ?")

# # annual cycle
# i1y <- which(data01$time>=2022 & data01$time<2023)
# nd <- data01[i1y,]
# q.DLGi2       <- qdo(qgam.DLGi2, q3, predict, newdata=nd )
# plot(q.DLGi2[[2]] ~ data01$doy[i1y], pch=20, cex=.3,ylim=c(0,18))
# i1y <- which(data01$time>=1980 & data01$time<1981)
# nd <- data01[i1y,]
# q.DLGi2       <- qdo(qgam.DLGi2, q3, predict, newdata=nd )
# points(q.DLGi2[[2]] ~ data01$doy[i1y], pch=20, cex=.3,col=2)
# readline("Continue 3 ?")

# q.DLGi2       <- qdo(qgam.DLGi2, q3, predict, newdata=data01 )
# plot(q.DLGi2[[2]] ~ data01$time, pch=20, cex=.3,ylim=c(0,18))
# readline("Continue 4 ?")


### year to year variability model

    load(st_msdata01, verb=TRUE)

    ### need to add here byYear covariates
    y.all <- trunc(data01$time)
    y.u   <- unique(y.all)
    cov.y <- array(0, dim=c(length(data01$time), length(y.u)) ) 
    fac.y <- NULL
    for(i in seq_along(y.u)) {
        iy <- which(y.all==y.u[i])
        cov.y[iy,i] <- 1
        fac.y       <- c(fac.y, rep(i, length(iy)) )
    }

### define seasons from a year in the middle
    iy2  <- which(data01$time>=1999.9 & data01$time<2001.4)
    iy2a <- which(data01$time>=1999.9 & data01$time<2000.8)
    iy2b <- which(data01$time>=2000.4 & data01$time<2001.3)
    clim50  <- qdo(qgam.DLG, 0.5, predict, newdata=data01[iy2,] )
    clim50a <- qdo(qgam.DLG, 0.5, predict, newdata=data01[iy2a,] )
    clim50b <- qdo(qgam.DLG, 0.5, predict, newdata=data01[iy2b,] )
    plot(data01$time[iy2], data01$x[iy2], pch=20, cex=.3)
    lines(data01$time[iy2a], clim50a, col=2, lwd=2)
    lines(data01$time[iy2b], clim50b, col=4, lwd=2, lty=2)
    iclim50a.n <- which.min(clim50a)
    iclim50a.x <- which.max(clim50a)
    iclim50b.n <- which.min(clim50b)
    iclim50b.x <- which.max(clim50b)
    abline(v=data01$time[iy2a][iclim50a.n], col=2, lty=1)
    abline(v=data01$time[iy2a][iclim50a.x], col=2, lty=1)
    abline(v=data01$time[iy2b][iclim50b.n], col=4, lty=2)
    abline(v=data01$time[iy2b][iclim50b.x], col=4, lty=3)
    mdoy.winter <- data01$doy[ iy2a[iclim50a.n] ]                                                # 65 
    mdoy.summer <- data01$doy[ iy2a[iclim50a.x] ]                                                # 227
    mdoy.spring <- data01$doy[ iy2a[iclim50a.n] + round((iy2a[iclim50a.x]-iy2a[iclim50a.n])/2) ] # 146
    mdoy.autumn <- data01$doy[ iy2b[iclim50b.n] + round((iy2b[iclim50b.x]-iy2b[iclim50b.n])/2) ] # 329
    abline(v=data01$time[iy2a][which(data01$doy[iy2a]==mdoy.spring)], col=3, lty=1)
    abline(v=data01$time[iy2b][which(data01$doy[iy2b]==mdoy.autumn)], col=3, lty=1)
    abline(v=data01$time[iy2b][iclim50b.n], col=4, lty=2)

    diff.w2s <- (mdoy.summer       - mdoy.winter)
    diff.s2w <- (mdoy.winter + 365 - mdoy.summer)

    doy.spring <- (mdoy.spring -round(diff.w2s/4)) : ((mdoy.spring +round(diff.w2s/4)) )
    doy.summer <- (mdoy.summer -round(diff.w2s/4)) : ((mdoy.summer +round(diff.s2w/4)) )
    doy.autumn <- (mdoy.autumn -round(diff.s2w/4)) : ((mdoy.autumn +round(diff.s2w/4)) )
    ig366 <- which(doy.autumn>=366)
    if(length(ig366)>0) {
        doy.autumn[ig366] <- doy.autumn[ig366] - 366 +1
    }
    doy.winter <- (mdoy.winter -round(diff.s2w/4)) : ((mdoy.winter +round(diff.w2s/4)) )

    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter)] ], rep(12,length(doy.winter)), pch=20, col=3, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.summer)] ], rep(12,length(doy.summer)), pch=20, col=2, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.spring)] ], rep(12.1,length(doy.spring)), pch=20, col=4, cex=.3)
    points(data01$time[ iy2b[which(data01$doy[iy2b] %in% doy.autumn)] ], rep(12.1,length(doy.autumn)), pch=20, col=1, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter)]+365 ], rep(12,length(doy.winter)), pch=20, col=3, cex=.3)

    range(doy.spring)
    range(doy.summer)
    paste(doy.autumn[1], doy.autumn[length(doy.autumn)])
    range(doy.winter)

    ## trim so no overlap in seasons
    doy.autumn2 <- doy.autumn[ which( !(doy.autumn %in% doy.summer) & !(doy.autumn %in% doy.winter) ) ]          
    doy.summer2 <- doy.summer[ which( !(doy.summer %in% doy.spring) & !(doy.summer %in% doy.autumn2) ) ]
    doy.spring2 <- doy.spring[ which( !(doy.spring %in% doy.summer2) & !(doy.spring %in% doy.winter) ) ]
    doy.winter2 <- doy.winter[ which( !(doy.winter %in% doy.spring2) & !(doy.winter %in% doy.autumn2) ) ]

    range(doy.spring2)
    range(doy.summer2)
    paste(doy.autumn2[1], doy.autumn2[length(doy.autumn2)])
    range(doy.winter2)
### END define seasons from a year in the middle

# ### plan 1
      # fit just for gmst and remove
      # fit for year to get year to year variability and remove
      # fit for sdoy to get annual cycle +gmst interaction

      ### add years as factors to data01
      data01$year  <- factor(fac.y,labels=y.u) # weired factor thing labels have to be the ranked order of the data
      data01$ryear <- fac.y/max(fac.y)
      data01$x0    <- data01$x

      # gmst
      fmla.gmst       <- list(x0 ~ s(gmst, bs='tp', k=4), ~ 1 )
      fit.gmst        <- mqgam(fmla.gmst, data=data01, qu=c(0.5,0.9))
      fit.gmst$fmla   <- fmla.gmst
      fit.gmst$ptiles <- c(0.5,0.9)
      q50       <- qdo(fit.gmst, 0.5, predict)
      data01$x2 <- data01$x0 - q50
            plot(data01$time, data01$x0, pch=20, cex=.3)
          points(data01$time, data01$x2+mean(data01$x0), pch=20, cex=.3, col=4)
           lines(data01$time, q50, col=2)
readline("Continue 3a ?")

      # year as factor, using 90th quantile to represent summer annual variability
      fmla.year       <- list(x2 ~ year, ~ 1 ) 
      fit.year        <- mqgam(fmla.year, data=data01, qu=c(0.5,0.9))
      fit.year$fmla   <- fmla.year
      fit.year$ptiles <- c(0.5,0.9)
      q90y      <- qdo(fit.year, 0.9, predict)

      data01$x <- data01$x2 - q90y
            plot(data01$time, data01$x0, pch=20, cex=.3)
          points(data01$time, q90y+mean(data01$x0), pch=20, cex=.3, col=2)
readline("Continue 3b1 ?")
          plot(data01$time, data01$x,  pch=20, cex=.3, col=2)
          abline(h=0,col=4)
readline("Continue 3b2 ?")

      # annual cycle
      fmla.htdata <- list(x ~ s(sdoy, bs='cc',k=12) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy)) 
      fit.htdata  <- mqgam(fmla.htdata, data=data01, qu=do.ptiles)
      fit.htdata$fmla   <- fmla.htdata
      fit.htdata$ptiles <- do.ptiles
      q50ht             <- qdo(fit.htdata, 0.5, predict)

         plot(data01$time, data01$x, pch=20, cex=.3)
        lines(data01$time, q50ht, col=2)
        grid()

         plot(data01$time, data01$x-q50ht, pch=20, cex=.3)
        lines(data01$time, q50ht, col=2)
        grid()
still not that satisfactory


readline("stop Continue 3 ?")

##############################################################################################
### plan 2 all together ######################################################################
load(st_msdata01, verb=TRUE)

    ### add years as factors to data01
    data01$year  <- factor(fac.y,labels=y.u) # weired factor thing labels have to be the ranked order of the data
    data01$ryear <- fac.y/max(fac.y)
    data01$x0    <- data01$x

    do.ptiles <- c(0.5,0.9 )
    fmla.MSqgamY <- list(x ~ year +gmst +s(sdoy, bs='cc',k=12) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))  
        # ok but discontinuities with year factors and too wak year to year
    # fmla.MSqgamY <- list(x ~       gmst +s(sdoy, bs='cc',k=12, by=year) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))  
        # very very slow
        # and terrible fit.  Diagnostics suggest good fit but actually unusable
        # this could be due to to dominance of winter interannual variability
    # that 
        
#   DO WE NEED TO RESTRICT TO SUMMER AS STILL LOADS OF interannual variability IN WINTER ??

    fit.MSqgamY  <- mqgam(fmla.MSqgamY, data=data01, qu=do.ptiles, argGam=list(correlation=corAR1(form=~1|year))) # adding correlation makes no difference

    # fit.MSqgamY  <- mqgam(fmla.MSqgamY, data=data01, qu=do.ptiles)
    fit.MSqgamY$fmla   <- fmla.MSqgamY
    fit.MSqgamY$ptiles <- do.ptiles

    pr90y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.9, predict)
    pr50y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.5, predict)
    # pr10y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.1, predict)
        
plot(data01$time, data01$x, pch=20, cex=.3, ylim=c(5,19))
lines(data01$time, pr90y,  pch=20, cex=.3, col=2)
lines(data01$time, pr50y,  pch=20, cex=.3, col=3)
# lines(data01$time, pr10y,  pch=20, cex=.3, col=4)
readline("Continue 4 ?")

##############################################################################################
### plan3 have mean seasonas in the fit ######################################################

data01$seas <- rep(NA, length(data01$doy))
data01$seas[ which(data01$doy %in% doy.winter2) ] <- 1/4
data01$seas[ which(data01$doy %in% doy.spring2) ] <- 2/4
data01$seas[ which(data01$doy %in% doy.autumn2) ] <- 3/4
data01$seas[ which(data01$doy %in% doy.summer2) ] <- 4/4
data01$fseas <- factor(data01$seas, levels=c('winter','spring','autumn','summer') )

# q3          <- c(0.1,0.5,0.9)
# f.DsLG      <- list(x ~ s(sdoy,seas, bs=c("cc","tp"), k=c(12,12)) +gmst, ~ 1 )   # fail cc only 1d
# qgam.DsLG   <- mqgam(f.DsLG, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))
# q.DsLG      <- qdo(qgam.DsLG, q3, predict, newdata=data01 )
# not sure this will do anything useful - seasonal effect already in sdoy term


#