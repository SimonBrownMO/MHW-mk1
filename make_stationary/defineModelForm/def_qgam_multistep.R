# source("def_qgam_multistep.R")

# parent:def_qgam_mm.R

### following Woods 371 more closely for mixed model quantile GAM fitting
#
## 1) fit gamm to get mean gmst term and baseline annual cycle
## 2) fit qgam to get quantiles of residuals from 1) as function of gmst and sdoy
## 3) combine to get full quantile estimates
# 4) adjust data to be stationary at desired quantile levels
# 5) fit final qgam & evgam to adjusted data to get stationary quantile estimates
# NB: this is similar to def_qgam_byyear but uses mixed model fitting for step 1

## also explore adding a seasonal mean model to the mixture to address interannual variability more completely
## hopefully keeping the cyclic seasonal component in the mixed model fitting contiuous over year-ends

source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
# library(evd)
# library(float)
source("../../libs/fn_MakeStationary.R")
source("../../libs/fn_HotDays.R")


# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v3/MSref/ostia_cdr_nrt_regions.MSref.2025-11-04.RData"
st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"

load(st_msref, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)
data01$year <- trunc(data01$time)
data01$x0   <- data01$x

### step 0: overit to capture interannual variability
fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=128), ~ 1 )
ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
 q.qgam0 <- predict(ft.qgam0, newdata=data01, type="response" )

# ft.qgam0_128 <- ft.qgam0
# ft.qgam0_256 <- ft.qgam0

#   plot(data01$time, data01$x0, pch=20, cex=.3, main=paste('Overfit',fm.qgam0))
#  lines(data01$time, q.qgam0, col=2, lwd=2)
#  grid()
### these plots indicate k==128 is insufficient
q.qgam0 <- predict(ft.qgam0_256, newdata=data01, type="response" )
plot(data01$time, data01$x0 - q.qgam0, pch=20, cex=.3, main=paste('Overfit',fm.qgam0))
plot(data01$doy,  data01$x0 - q.qgam0, pch=20, cex=.3, main=paste('Overfit',fm.qgam0))






### step 1: fit mixed model to get mean gmst term and baseline annual cycle
q3         <- c(0.1,0.5,0.9)
fm.gamm1    <- list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
ft.gamm1    <- gamm(   fm.gamm1[[1]],    data=data01 )
 q.gamm1    <- predict(ft.gamm1$gam,  newdata=data01 ) 
iy2000 <- which(data01$year==2000)
i0n    <- which.min(q.gamm1[ iy2000] )
i0x    <- which.max(q.gamm1[ iy2000] )
doyn   <- data01$doy[ iy2000][i0n]
doyx   <- data01$doy[ iy2000][i0x]
sdoyn  <- data01$sdoy[iy2000][i0n]
sdoyx  <- data01$sdoy[iy2000][i0x]
nd       <- data01
nd$sdoy  <- sdoyn
qn.gamm1 <- predict(ft.gamm1$gam,  newdata=nd )
nd$sdoy  <- sdoyx
qx.gamm1 <- predict(ft.gamm1$gam,  newdata=nd )

 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time, q.gamm1,  col=2)
lines(data01$time, qn.gamm1, col=4)
lines(data01$time, qx.gamm1, col=3)
grid()
readline("Stop2 1")

# remove fit1 from data
data01$q.gamm1 <- q.gamm1
data01$x1      <- data01$x0 - q.gamm1

### step 2: fit qgam to  residuals from step 1
fm.qgam2 <- list( x1 ~ s(stime, bs="tp", k=128), ~ 1 )
ft.qgam2 <- qgam( fm.qgam2[[1]], data=data01, qu=0.5 )
 q.qgam2 <- predict(ft.qgam2, newdata=data01, type="response" )

  plot(data01$time, data01$x1, pch=20, cex=.3, main=paste(fm.qgam2))
 lines(data01$time, q.qgam2, col=2)
 grid()

fm.qgam2b <- list( x1 ~ s(stime,doy, bs=c("tp","tp"), k=16), ~ 1 )
ft.qgam2b <- qgam(   fm.qgam2b[[1]], data=data01, qu=0.5 )
 q.qgam2b <- predict(ft.qgam2b,   newdata=data01, type="response" )

  plot(data01$time[1:3650], data01$x1[1:3650], pch=20, cex=.3, main=paste(fm.qgam2))
 lines(data01$time[1:3650],  q.qgam2b[1:3650], col=2)
 grid()

readline("Stop2 2")

### residuals between overfit and gamm1

plot(data01$time, q.qgam0-q.gamm1, pch=20, cex=.3, main='residuals between overfit and gamm1')


  plot(data01$time[3260:3310], q.qgam0[3260:3310], pch=20, cex=.3,,ylim=c(7,10))
points(data01$time[3260:3310], q.gamm1[3260:3310], pch=20, cex=.3,col=2)

## gamm 2 






#