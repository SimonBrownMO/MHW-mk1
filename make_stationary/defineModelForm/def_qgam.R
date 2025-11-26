# source("def_qgam.R")

### Conclusions
#
# 1) gamma smoothing is a bad idea as it affected the doy term too
# 2) for the comparision with Hobday use 
#           fmla.qgam2  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )  
#       ie linear term for gmst which can be different for different quantiles
# 
# 3) for interaction term tried
f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, bs="cc", k=1,by=gmst),       ~ 1 )  # all about the same
f.qgam.DLGi2 <- list(x ~ s(sdoy) +s(sdoy, bs="cc", k=12,by=gmst),                      ~ 1 )  # all about the same
f.qgam.DLGi3 <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, by=gmst),                     ~ 1 )  # all about the same
#   all give similar results, with k=12 the interation term is very wiggly, with k=1 it becomes a sinusoid
# not clear if this is real or just fitting to interannual noise.

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
st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-10-29.RData"
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
f.qgam.DLG  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
qgam.DLG    <- mqgam(f.qgam.DLG, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))
q.DLG       <- qdo(qgam.DLG, q3, predict, newdata=data01 )

if(FALSE) {
  # FAIL f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst +ti(sdoy,gmst, bs=c('cc','tp'), k=c(12,4)), ~ 1 )
  # FAIL f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst +ti(sdoy,gmst, bs=c('cc','1'), k=c(12,4)), ~ 1 )
  # FAIL f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12,by=gmst) +gmst, ~ 1 )
  # FAIL f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, bs="cc", k=12,by=gmst) +gmst, ~ 1 )
  # FAIL f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12,by=gmst) ,       ~ 1 ) 
  # FAIL f.qgam.DLGi  <- list(x ~ gmst +s(sdoy, bs="cc", k=12,by=gmst) ,       ~ 1 ) 
  f.qgam.DLGi  <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, bs="cc", k=1,by=gmst),       ~ 1 )  # all about the same
  f.qgam.DLGi2 <- list(x ~ s(sdoy) +s(sdoy, bs="cc", k=12,by=gmst),                      ~ 1 )  # all about the same
  f.qgam.DLGi3 <- list(x ~ s(sdoy, bs="cc", k=12) +s(sdoy, by=gmst),                     ~ 1 )  # all about the same
  qgam.DLGi3    <- mqgam(f.qgam.DLGi3, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))
  q.DLGi3       <-     qdo(qgam.DLGi3, q3, predict, newdata=data01 )

  qgam.DLGi    <- mqgam(f.qgam.DLGi, data=data01, qu=q3, argGam=list(correlation=corAR1(form=~1|year)))
  q.DLGi       <-     qdo(qgam.DLGi, q3, predict, newdata=data01 )

  plot(  data01$time, data01$x,  pch=20, cex=.3)
  points(data01$time, q.DLG[[3]], pch=20, cex=.3, col=2)
  idoy <- which(data01$doy==220)    
  lines(range(data01$time),range(q.DLG[[3]][idoy]),lwd=1, lty=2)
  lines(range(data01$time),range(q.DLG[[2]][idoy]),lwd=1, lty=2)
  grid()
  # > diff(range(q.DLG[[1]][idoy])) [1] 1.344837
  # > diff(range(q.DLG[[2]][idoy])) [1] 1.222375
  # > diff(range(q.DLG[[3]][idoy])) [1] 1.089935
  lines(range(data01$time),range(q.DLGi[[3]][idoy]),lwd=1, lty=3)
  points(data01$time, q.DLGi[[3]], pch=20, cex=.3, col=3)

  readline("Continue -1 ?")

  i0      <- which(data01$time>=2000 & data01$time<2001)
  nd      <- data01[i0,]
  nd$gmst <- min(data01$gmst)
  p1      <- qdo(qgam.DLG,   q3[3], predict, newdata=nd )
  p1i     <- qdo(qgam.DLGi , q3[3], predict, newdata=nd )
  p1i2    <- qdo(qgam.DLGi2, q3[3], predict, newdata=nd )
  p1i3    <- qdo(qgam.DLGi3, q3[3], predict, newdata=nd )
  nd$gmst <- max(data01$gmst)
  p2      <- qdo(qgam.DLG,   q3[3], predict, newdata=nd )
  p2i     <- qdo(qgam.DLGi , q3[3], predict, newdata=nd )
  p2i2    <- qdo(qgam.DLGi2, q3[3], predict, newdata=nd )
  p2i3    <- qdo(qgam.DLGi3, q3[3], predict, newdata=nd )
  r1 <- range(c(p1,p2, nd$x))
  plot(  nd$time, nd$x,  pch=20, cex=.3, ylim=r1)
  points(nd$time, p1, pch=20, cex=.3, col=1)
  points(nd$time, p2, pch=20, cex=.3, col=1)
  points(nd$time, p1i,  pch=20, cex=.3, col=4)
  points(nd$time, p2i,  pch=20, cex=.3, col=4)
  points(nd$time, p1i2, pch=20, cex=.3, col=2)
  points(nd$time, p2i2, pch=20, cex=.3, col=2)
  points(nd$time, p1i3, pch=20, cex=.3, col=3)
  points(nd$time, p2i3, pch=20, cex=.3, col=3)

  r1 <- range(c(p2i3-p1i3,p2i2-p1i2, p2-p1))  # TBH this complex change in annual cycle is dubious
  plot(p2-p1, pch=20, cex=.3, ylim=r1)        # thus probably best to avoid interaction term
  points(p2i2-p1i2, pch=20, cex=.3, col=2)
  points(p2i3-p1i3, pch=20, cex=.3, col=3)

  r1 <- range(c(p2i3-p1i3,p2i2-p1i2, p2-p1))  # TBH this complex change in annual cycle is dubious
    plot(p2-p1, pch=20, cex=.3, ylim=r1)        # thus probably best to avoid interaction term
  points(p2i-p1i,   pch=20, cex=.3, col=4)
  points(p2i2-p1i2, pch=20, cex=.3, col=2)
  points(p2i3-p1i3, pch=20, cex=.3, col=3)

  readline("Continue -1b ?")

  i0      <- which(data01$doy==220)
  nd      <- data01[i0,]
  p1      <- qdo(qgam.DLG,   q3[3], predict, newdata=nd )
  p1i     <- qdo(qgam.DLGi , q3[3], predict, newdata=nd )
  p1i2    <- qdo(qgam.DLGi2, q3[3], predict, newdata=nd )
  p1i3    <- qdo(qgam.DLGi3, q3[3], predict, newdata=nd )

  r1 <- range(c(p1i3,p1i2,p1))
    plot(p1,     pch=20, cex=.3, ylim=r1, ty='l')
  lines(p1i2, pch=20, cex=.3, col=2)
  lines(p1i3, pch=20, cex=.3, col=3)

  plot(p1i3-p1, ty='l', ylim=c(-0.2,0.2)) 
  lines(p1i2-p1, col=2)                   
  lines(p1i2-p1i3, col=3)
  grid()

  ### look at anomalies
  p1      <- qdo(qgam.DLG,   q3[3], predict, newdata=data01 )
  p1i2    <- qdo(qgam.DLGi2, q3[3], predict, newdata=data01 )
  p1i3    <- qdo(qgam.DLGi3, q3[3], predict, newdata=data01 )
    plot(data01$time, data01$x - p1,   pch=20, cex=.3, ylim=c(0,2))
  points(data01$time, data01$x - p1i2, pch=20, cex=.3, col=2)
  points(data01$time, data01$x - p1i3, pch=20, cex=.3, col=3)
  grid()

    plot(data01$doy, data01$x - p1,   pch=20, cex=.3, ylim=c(0,2))
  points(data01$doy, data01$x - p1i2, pch=20, cex=.3, col=2)
  points(data01$doy, data01$x - p1i3, pch=20, cex=.3, col=3)
  grid()







  gamm1 <- mgcv::gamm(formula = x ~ s(sdoy, bs="cc", k=12) +s(gmst, bs="cr", k=5) +ti(sdoy,gmst), data=data01, correlation=corAR1(form=~1|year))                    
  
  f.qgam  <-               list(x ~ s(sdoy, bs="cc", k=12) +s(gmst, bs="tp", k=5), ~ 1 )
  f.qgam  <-               list(x ~ s(sdoy, bs="cc", k=12) +s(gmst, bs="cr", k=6), ~ 1 )
  qgam50.200 <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=200, correlation=corAR1(form=~1|year)))
  qgam50.50  <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=50,  correlation=corAR1(form=~1|year)))
  qgam50.80  <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=80,  correlation=corAR1(form=~1|year)))
  qgam50     <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(           correlation=corAR1(form=~1|year)))
  plot(qgam50)  
  readline("Continue 0 ?")
  
  nd <- data01
  nd$sdoy <- 0.6
  pr.50     <- predict(qgam50,    newdata=nd, type='response')
  pr.50.200 <- predict(qgam50.200,newdata=nd, type='response')
  pr.50.80  <- predict(qgam50.80, newdata=nd, type='response')
  pr.50.50  <- predict(qgam50.50, newdata=nd, type='response')
  plot(data01$time, data01$x, pch=20, cex=.3, ylim=c(14,20))
  points(data01$time, pr.50,     pch=20, cex=.3, col=2)
  points(data01$time, pr.50.200, pch=20, cex=.3, col=3)
  points(data01$time, pr.50.80,  pch=20, cex=.3, col=4)
  points(data01$time, pr.50.50,   pch=20, cex=.3, col=5)
  legend('topright', legend=c('data','qgam50','qgam50.200','qgam50.80','qgam50.50'), col=1:5, pch=20)      
  # lines(c(1980,2025),c(15.5,16.8),lwd=2)
  lines(range(data01$time),range(pr.50.200),lwd=1, lty=2)
  
  
  f.qgam2  <- list(x ~ s(sdoy, bs="cc", k=12) +s(gmst, bs="cr", k=4), ~ 1 )
  f.qgam2  <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
  qgam502  <- qgam(f.qgam2, data=data01, qu=0.5, argGam=list(correlation=corAR1(form=~1|year)))
  
  pr.502 <- predict(qgam502, newdata=nd, type='response')
  points(data01$time, pr.502,     pch=1, cex=.3, col='orange2')
  
  readline("Continue 1 ?")
  
  ### plot doy effect
  pr.50     <- predict(qgam50,    newdata=data01, type='response')
  pr.50.200 <- predict(qgam50.200,newdata=data01, type='response')
  pr.50.80  <- predict(qgam50.80, newdata=data01, type='response')
  pr.50.50  <- predict(qgam50.50, newdata=data01, type='response')
  plot(data01$time, data01$x, pch=20, cex=.3, ylim=c(14,20))
  points(data01$time, pr.50,     pch=20, cex=.3, col=2)
  points(data01$time, pr.50.200, pch=20, cex=.3, col=3)
  points(data01$time, pr.50.80,  pch=20, cex=.3, col=4)
  points(data01$time, pr.50.50,   pch=20, cex=.3, col=5)
  legend('topright', legend=c('data','qgam50','qgam50.200','qgam50.80','qgam50.50'), col=1:5, pch=20)      
  # lines(c(1980,2025),c(15.5,16.8),lwd=2)
  lines(range(data01$time),range(pr.50.200),lwd=1, lty=2)
  
  
  
  pr.502 <- predict(qgam502, newdata=nd, type='response')
  points(data01$time, pr.502,     pch=1, cex=.3, col='orange2')
  readline("Continue 1b ?")
  
  ### round 2 what gamma value to use to get very smooth gmst effect
  f.qgam  <-               list(x ~ s(sdoy, bs="cc", k=12) +s(gmst, bs="tp", k=20), ~ 1 )
  qgA.1600<- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=1600, correlation=corAR1(form=~1|year)))
  qgA.160 <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=160, correlation=corAR1(form=~1|year)))
  qgA.80  <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=80,  correlation=corAR1(form=~1|year)))
  qgA.40  <- qgam(f.qgam, data=data01, qu=0.5, argGam=list(gamma=40,  correlation=corAR1(form=~1|year)))
  
  nd       <- data01
  nd$sdoy  <- 0.6
  prA.1600 <- predict(qgA.1600,newdata=nd, type='response')
  prA.160  <- predict(qgA.160, newdata=nd, type='response')
  prA.80   <- predict(qgA.80,  newdata=nd, type='response')
  prA.40   <- predict(qgA.40,  newdata=nd, type='response')
  plot(data01$time, data01$x, pch=20, cex=.3, ylim=c(14,20))
  points(data01$time, prA.160, pch=20, cex=.3, col=2)
  points(data01$time, prA.80,  pch=20, cex=.3, col=3)
  points(data01$time, prA.40,  pch=20, cex=.3, col=4)
  # points(data01$time, prA.1600, pch=20, cex=.3, col=5)  # CAUTION introduces a bias - lower values
  legend('topright', legend=c('data','qgamA.160','qgamA.80','qgamA.40'), col=1:4, pch=20)      
  # lines(c(1980,2025),c(15.5,16.8),lwd=2)
  lines(range(data01$time),range(prA.160),lwd=2)      
  
  # whats happenign with doy
  i0       <- which(data01$time>1981 & data01$time<2020)
  nd       <- data01[i0,]
  prA2.1600 <- predict(qgA.1600,newdata=nd, type='response')
  prA2.40   <- predict(qgA.40,  newdata=nd, type='response')
  
  plot(  nd$time, nd$x, pch=20, cex=.3) #, ylim=c(14,20))
  points(nd$time, prA2.1600,   pch=20, cex=.3, col=2)
  points(nd$time, prA2.40,   pch=20, cex=.3, col=3)
  lines(range(nd$time),range(prA.160),lwd=2)      
  
  readline("Continue 2 ?")
} # skip

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

### plan 1
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
        # plot(data01$time, data01$x, pch=20, cex=.3)
        # points(data01$time, data01$x2+mean(data01$x), pch=20, cex=.3, col=4)

    # year as factor, using 90th quantile to represent summer annual variability
    fmla.year       <- list(x2 ~ year, ~ 1 ) 
    fit.year        <- mqgam(fmla.year, data=data01, qu=c(0.5,0.9))
    fit.year$fmla   <- fmla.year
    fit.year$ptiles <- c(0.5,0.9)
    q90y      <- qdo(fit.year, 0.9, predict)

    data01$x <- data01$x2 - q90y
        # plot(data01$time, data01$x, pch=20, cex=.3)
        # points(data01$time, q90y+mean(data01$x), pch=20, cex=.3, col=2)
        # plot(data01$time, data01$x3,  pch=20, cex=.3, col=2)
        # abline(h=0,col=4)

    # annual cycle
    fmla.htdata <- list(x ~ s(sdoy, bs='cc',k=12) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy)) 
    fit.htdata  <- mqgam(fmla.htdata, data=data01, qu=do.ptiles)
    fit.htdata$fmla   <- fmla.htdata
    fit.htdata$ptiles <- do.ptiles


readline("Continue 3 ?")

### plan 2 all together
load(st_msdata01, verb=TRUE)

    ### add years as factors to data01
    data01$year  <- factor(fac.y,labels=y.u) # weired factor thing labels have to be the ranked order of the data
    data01$ryear <- fac.y/max(fac.y)
    data01$x0    <- data01$x

    do.ptiles <- c(0.1,0.5,0.9 )
    fmla.MSqgamY <- list(x ~ year +gmst +s(sdoy, bs='cc',k=12) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))  
        # ok but discontinuities with year factors and too wak year to year
    fmla.MSqgamY <- list(x ~       gmst +s(sdoy, bs='cc',k=12, by=year) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy))  
        # very very slow
        
    # fit.MSqgamY  <- mqgam(fmla.MSqgamY, data=data01, qu=do.ptiles, argGam=list(correlation=corAR1(form=~1|year))) # adding correlation makes no difference
    fit.MSqgamY  <- mqgam(fmla.MSqgamY, data=data01, qu=do.ptiles)
    fit.MSqgamY$fmla   <- fmla.MSqgamY
    fit.MSqgamY$ptiles <- do.ptiles

    pr90y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.9, predict)
    pr50y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.5, predict)
    pr10y <- qdo(fit.MSqgamY, newdata=data01, type='response', qu=0.1, predict)
        
plot(data01$time, data01$x, pch=20, cex=.3, ylim=c(5,19))
lines(data01$time, pr90y,  pch=20, cex=.3, col=2)
lines(data01$time, pr50y,  pch=20, cex=.3, col=3)
lines(data01$time, pr10y,  pch=20, cex=.3, col=4)





#