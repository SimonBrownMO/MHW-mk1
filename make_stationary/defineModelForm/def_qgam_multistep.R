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

### Conclusions ###############################################################################
### Overfit model to capture interannual variability is best with k=128 tp smoother spline
    # although this leaves a bi-annual cycle in the residuals
    # k=256 is more wiggly but without bi-annual cycle in the residuals
    # HOWEVER simple smooth.spline with spar=0.01 is even better
      # > summary(data01$x0 - smsp.x0$y)
      #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
      # -1.389471 -0.233948 -0.005821  0.000000  0.214433  1.987574 

    # it is still unclear whether the simple spline fitted to all data is better than the qgam fitted to the first residuals

    # Not using the simple spline overfit
      # step1: list(x ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
      #   -> data01: x1, q.gamm1
      # step2: list( x1 ~ s(stime, bs="tp", k=128), ~ 1 )
      #   -> data01: x2, q.qgam2a

### END Conclusions ###########################################################################



source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
# library(evd)
# library(float)
source("../../libs/fn_MakeStationary.R")
source("../../libs/fn_HotDays.R")


# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v3/MSref/ostia_cdr_nrt_regions.MSref.2025-26-04.RData"
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"

st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-28.RData"

load(st_msref, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)
data01$year <- trunc(data01$time)
data01$x0   <- data01$x

###############################################################################################
### step 0: overfit to capture interannual variability
###############################################################################################

tit2 <- 'Overfit smooth.spline'
  # in fact simple smooth.spline the best
  smsp.x0 <- smooth.spline(data01$stime,data01$x0, spar=0.01)
    plot(data01$time, data01$x0, pch=20, cex=.3, main=tit2)
  lines(data01$time, smsp.x0$y, col=2, lwd=2)
  grid()
readline("Stop0 1")

  #
      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=512), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_512 <- ft.qgam0
      # q.qgam0.512 <- predict(ft.qgam0_512, newdata=data01, type="response" )
      # save(file="ft.qgam0_512.RData", ft.qgam0_512,  fm.qgam0)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=286), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_286 <- ft.qgam0
      # q.qgam0.286 <- predict(ft.qgam0_286, newdata=data01, type="response" )
      # save(file="ft.qgam0_286.RData", ft.qgam0_286,  fm.qgam0)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=200), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_200 <- ft.qgam0
      # q.qgam0.200 <- predict(ft.qgam0_200, newdata=data01, type="response" )
      # save(file="ft.qgam0_200.RData", ft.qgam0_200,  fm.qgam0)
      # lines(data01$time, q.qgam0.200-q.gamm1, col=3, lwd=2)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=128), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_128 <- ft.qgam0
      # q.qgam0.128 <- predict(ft.qgam0_128, newdata=data01, type="response" )
      # save(file="ft.qgam0_128.RData", ft.qgam0_128,  fm.qgam0)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="cs", k=128), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_128cs <- ft.qgam0
      # q.qgam0.128cs <- predict(ft.qgam0_128cs, newdata=data01, type="response" )
      # save(file="ft.qgam0_128cs.RData", ft.qgam0_128cs,  fm.qgam0)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=90), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_90 <- ft.qgam0
      # q.qgam0.90 <- predict(ft.qgam0_90, newdata=data01, type="response" )
      # save(file="ft.qgam0_90.RData", ft.qgam0_90,  fm.qgam0)

      # fm.qgam0 <- list( x0 ~ s(stime, bs="tp", k=110), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_110 <- ft.qgam0
      # q.qgam0.110 <- predict(ft.qgam0_110, newdata=data01, type="response" )
      # save(file="ft.qgam0_110.RData", ft.qgam0_110,  fm.qgam0)

      # # very smooth like GAMM gmst term
      # fm.qgam0 <- list( x0 ~ s(stime, bs="ds", k=128), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )
      # ft.qgam0_128ds <- ft.qgam0
      # q.qgam0.128ds <- predict(ft.qgam0_128ds, newdata=data01, type="response" )
      # save(file="ft.qgam0_128ds.RData", ft.qgam0_128ds,  fm.qgam0)

      # # not worth it
      # fm.qgam0 <- list( x0 ~ s(stime,doy, bs=c("ds","tp"), k=128), ~ 1 )
      # ft.qgam0 <- qgam(   fm.qgam0[[1]], data=data01, qu=0.5 )


      #   plot(data01$time, data01$x0, pch=20, cex=.3, main=paste('Overfit',fm.qgam0[[1]][3]))
      #  lines(data01$time, q.qgam0, col=2, lwd=2)
      #  grid()
  #

  plot(data01$time, data01$x0 - smsp.x0$y, pch=20, cex=.3, main=tit2)
  grid()
readline("Stop0 2")
  plot(data01$doy,  data01$x0 - smsp.x0$y, pch=20, cex=.3, main=tit2)
  grid()
readline("Stop0 3")
### end step 0 ###########################################################################

###############################################################################################
### step 1: fit mixed model to get mean gmst term and baseline annual cycle
###############################################################################################
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

plot(data01$time, data01$x0 - q.gamm1, pch=20, cex=.3, main=paste('MM :',fm.gamm1[[1]][3]))
grid()
readline("Stop2 1b")

plot(data01$doy,  data01$x0 - q.gamm1, pch=20, cex=.3, main=paste('MM: ',fm.gamm1[[1]][3]))
grid()
### -> good fit to seasonal cycle and gmst

# remove fit1 from data
data01$q.gamm1 <- q.gamm1
data01$x1      <- data01$x0 - q.gamm1
readline("Stop2 end")



###############################################################################################
### step 1b: determine interannual variabitily from over fit & mid season #####################
###############################################################################################
doy.seasons <- find_doy_for_seasons(data01, ft.gamm1)
readline("Stop2b end")
doy.midseas <- NULL
i2001       <- which(data01$year==2001)
name.seas   <- names(doy.seasons)

iseas                                   <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='winter')]])
doy.midseas[which(name.seas=='winter')] <- data01$doy[i2001][iseas][which.min(q.gamm1[i2001][iseas])]
iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='spring')]])
doy.midseas[which(name.seas=='spring')] <- round(median(data01$doy[i2001][iseas]))
iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='summer')]])
doy.midseas[which(name.seas=='summer')] <- data01$doy[i2001][iseas][which.max(q.gamm1[i2001][iseas])]
iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='autumn')]])
doy.midseas[which(name.seas=='autumn')] <- round(median(data01$doy[i2001][iseas]))

names(doy.midseas) <- names(doy.seasons)
print(doy.midseas)

# plot(data01$doy[i2001], q.gamm1[i2001])
# abline(v=doy.midseas)
intann.seas.anom <- list()
for(i in 1:4){
  i1                      <- which(data01$doy==doy.midseas[i])
  intann.seas.anom[[i]]   <- (smsp.x0$y-q.gamm1)[i1]
  intann.seas.anom[[i+4]] <- data01$time[i1]
}
names(intann.seas.anom) <- c(names(doy.seasons), paste0(names(doy.seasons),".years"))

plot(data01$time,      smsp.x0$y-q.gamm1, ty='n')
 for(i in c(2)){ # for(i in 1:4){
   lines(intann.seas.anom[[i+4]], intann.seas.anom[[i]], col=i+2, pch=20, cex=1)
  points(intann.seas.anom[[i+4]], intann.seas.anom[[i]], col=i+2, pch=20, cex=1)
}
readline("Stop1b 1")

### inter season dependence
r1 <- range(-1.5,1)
# names(doy.seasons) "spring" "summer" "autumn" "winter"
plot(intann.seas.anom[[4]],       intann.seas.anom[[1]], pch=20, cex=1, xlim=r1, ylim=r1, xlab='winter', ylab='spring') # correlation 0.5036951
readline("Stop1b 12")
plot(intann.seas.anom[[1]],       intann.seas.anom[[2]], pch=20, cex=1, xlim=r1, ylim=r1, xlab='spring', ylab='summer')# correlation 0.2241365
readline("Stop1b 13")
plot(intann.seas.anom[[2]][1:45], intann.seas.anom[[3]], pch=20, cex=1, xlim=r1, ylim=r1, xlab='summer', ylab='autumn') # correlation 0.2815403
readline("Stop1b 14")
plot(intann.seas.anom[[3]], intann.seas.anom[[4]][2:46], pch=20, cex=1, xlim=r1, ylim=r1, xlab='autumn', ylab='winter') # correlation 0.500645
readline("Stop1b 15")

plot(intann.seas.anom[[4]], intann.seas.anom[[2]], pch=20, cex=1, xlim=r1, ylim=r1, xlab='winter', ylab='summer') # correlation 0.3624361
readline("Stop1b 16")
# NONE are significant

### decide the next fit
 plot(data01$time, q.gamm1, ty='l')
lines(data01$time, smsp.x0$y,  col=2)
points(intann.seas.anom[[2+4]], intann.seas.anom[[2]]+12, col=3, pch=20, cex=1)
abline(h=12)
readline("Stop1b 17")

 plot(data01$time, data01$x1, pch=20, cex=.3, main='x1 residuals from step 1')
points(intann.seas.anom[[1+4]], intann.seas.anom[[1]], col=3, pch=20, cex=2)
points(intann.seas.anom[[2+4]], intann.seas.anom[[2]], col=2, pch=20, cex=2)
points(intann.seas.anom[[3+4]], intann.seas.anom[[3]], col=5, pch=20, cex=2)
points(intann.seas.anom[[4+4]], intann.seas.anom[[4]], col=4, pch=20, cex=2)

# difference too confusiog
 plot(data01$time, smsp.x0$y-q.gamm1, ty='l')

readline("Stop1b 18")
readline("WHAT TO DO?")
# WHAT TO DO?
















###############################################################################################
### step 2: fit qgam to  residuals from step 1 ################################################
###############################################################################################

fm.qgam2a <- list( x1 ~ s(stime, bs="tp", k=128), ~ 1 )
ft.qgam2a <- qgam(   fm.qgam2a[[1]], data=data01, qu=0.5 )
 q.qgam2a <- predict(ft.qgam2a,   newdata=data01, type="response" )

fm.qgam2b <- list( x1 ~ s(stime, bs="tp", k=96), ~ 1 )
ft.qgam2b <- qgam(   fm.qgam2b[[1]], data=data01, qu=0.5 )
 q.qgam2b <- predict(ft.qgam2b,   newdata=data01, type="response" )

fm.qgam2c <- list( x1 ~ s(stime, bs="tp", k=64), ~ 1 )
ft.qgam2c <- qgam(   fm.qgam2c[[1]], data=data01, qu=0.5 )
 q.qgam2c <- predict(ft.qgam2c,   newdata=data01, type="response" )

  smsp.x1 <- smooth.spline(data01$stime,data01$x1)
    plot(data01$time, data01$x1, pch=20, cex=.3, main='smooth.spline x1')
  lines(data01$time, smsp.x1$y, col=2, lwd=2)
  grid()
readline("Stop3 a")

  plot(data01$time, data01$x1, pch=20, cex=.3, main=paste(fm.qgam2))
#  lines(data01$time, q.qgam2a, col=2, lwd=2)
#  lines(data01$time, q.qgam2b, col=3, lwd=2)
 lines(data01$time, q.qgam2c, col=4, lwd=2)
 lines(data01$time, smsp.x1$y, col=5, lwd=4)
 grid()
readline("Stop3 b")

#   plot(data01$doy, data01$x1, pch=20, cex=.3, main=paste(fm.qgam2))
#  lines(data01$doy, q.qgam2a, col=2, lwd=2)
#  lines(data01$doy, q.qgam2b, col=3, lwd=2)
#  lines(data01$doy, q.qgam2c, col=4, lwd=2)
#  lines(data01$doy, smsp.x1$y, col=5, lwd=4)
#  grid()

## no good  fm.qgam2b <- list( x1 ~ s(stime,doy, bs=c("tp","tp"), k=16), ~ 1 )

### residuals between overfit and gamm1
# plot(data01$time, q.qgam0.128, ty='l', main='model comparison')
plot(data01$time, smsp.x0$y, ty='l', main='model comparison')
lines(data01$time, q.gamm1, col=2)
readline("Stop3 c")

 plot(data01$sdoy, smsp.x0$y, ty='l', main='model comparison')
lines(data01$sdoy, q.gamm1, col=2)
readline("Stop3 d")

# plot(data01$time, q.qgam0.128-q.gamm1, pch=20, cex=.3, main='residuals between overfit and gamm1')
plot(data01$time, smsp.x0$y-q.gamm1, pch=20, cex=.3, main='residuals between overfit and gamm1')
lines(data01$time, q.qgam2c, col=2)
## good agreement here means we dont need  ft.qgam0
grid()
readline("Stop3 d")

plot(data01$time, smsp.x0$y-q.gamm1, pch=20, cex=.3, main='residuals between overfit and gamm1')
lines(data01$time, 0.2*(smsp.x0$y-mean(smsp.x0$y)), col=2, lwd=1)
grid()
readline("Stop3 e")

rm.x0 <- runmed(data01$x0, k=91)
rm.x1 <- runmed(data01$x1, k=91)
plot(data01$time, rm.x0, pch=20, cex=.3, main='running median')


## see second degree residuals
data01$q.qgam2a <-  q.qgam2a
data01$x2 <-  data01$x1 - q.qgam2a

plot(data01$time, data01$x2, pch=20, cex=.3, main='second degree residuals')
grid()  
readline("Stop3 e")
plot(data01$doy, data01$x2, pch=20, cex=.3, main='second degree residuals')
grid()  
readline("Stop3 f")

### compare with overfit residuals
plot(data01$x2, data01$x0 - smsp.x0$y, pch=20, cex=.3, main='overfit residuals')
abline(0,1)
grid()
readline("Stop3 g")

plot(data01$time, data01$x2 - (data01$x0 - smsp.x0$y), pch=20, cex=.3, main='overfit residuals')
abline(0,1)
grid()
readline("Stop3 h")
plot(data01$doy, data01$x2 - (data01$x0 - smsp.x0$y), pch=20, cex=.3, main='overfit residuals')
abline(0,1)
grid()
readline("Stop3 i")


### review
### check all
plot(data01$time, data01$x0, pch=20, cex=.3)
plot(data01$time, data01$x1, pch=20, cex=.3)
plot(data01$time, data01$x2, pch=20, cex=.3)
plot(data01$doy, data01$x0, pch=20, cex=.3)
plot(data01$doy, data01$x1, pch=20, cex=.3)
plot(data01$doy, data01$x2, pch=20, cex=.3)

### check seasons
doy.seasons <- find_doy_for_seasons()
for(s1 in 1:4) {
  i1 <- which(data01$doy %in% doy.seasons[[s1]])
  d2 <- data01[i1,]
  plot(d2$time, d2$x0, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
  plot(d2$time, d2$x1, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
  plot(d2$time, d2$x2, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
  plot(d2$doy,  d2$x0, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
  plot(d2$doy,  d2$x1, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
  plot(d2$doy,  d2$x2, pch=20, cex=.3, main=names(doy.seasons)[s1])
  readline("continue?")
}

### now need to see what structure remains in x1 & x2 for all quantiles
# do.ptiles <- seq(from=1, to=99, by=1)/100.0     
do.ptiles <- seq(from=10, to=90, by=10)/100.0     
ms.k$doy  <- 12
ms.k$gmst <- 4

### x1 first
fm3.qgam.x1     <- list(x1 ~ s(sdoy, bs="cc", k=ms.k$doy) + s(gmst, bs='tp', k=ms.k$gmst) +    ti(sdoy, gmst, bs=c("cc", "tp"))  , ~ s(sdoy))
# fm3.qgam.x1   <- list(x1 ~ gmst  +s(sdoy, bs='cc',k=ms.k$doy)                                                                  , ~ s(sdoy))

ft3.qgam.x1        <- mqgam(fm3.qgam.x1, data=data01, qu=do.ptiles)
ft3.qgam.x1$ptiles <- do.ptiles
ft3.qgam.x1$fmla   <- fm3.qgam.x1
cat('Chosen MS QGAM:',cr)
print(qdo(ft3.qgam.x1, 0.9, summary))

q3 <- list()
for(p1 in 1:length(do.ptiles)) {
  q3[[p1]] <- qdo(ft3.qgam.x1, do.ptiles[p1], predict )
}
names(q3) <- paste0('q',do.ptiles*100)  

plot(data01$time, data01$x1, pch=20, cex=.3, main='QGAM fits to x1')
for(p1 in 1:length(do.ptiles)) {
  lines(data01$time, q3[[p1]], col=p1+1, lwd=2)
}
grid()
readline("Stop4 a")
### this looks nuts

### x2 second
fm3.qgam.x2     <- list(x2 ~ s(sdoy, bs="cc", k=ms.k$doy) + s(gmst, bs='tp', k=ms.k$gmst) +    ti(sdoy, gmst, bs=c("cc", "tp"))  , ~ s(sdoy))
# fm3.qgam.x2   <- list(x2 ~ gmst  +s(sdoy, bs='cc',k=ms.k$doy)                                                                  , ~ s(sdoy))

ft3.qgam.x2        <- mqgam(fm3.qgam.x2, data=data01, qu=do.ptiles)
ft3.qgam.x2$ptiles <- do.ptiles
ft3.qgam.x2$fmla   <- fm3.qgam.x2
cat('Chosen MS QGAM:',cr)
print(qdo(ft3.qgam.x2, 0.9, summary))

q3 <- list()
for(p1 in 1:length(do.ptiles)) {
  q3[[p1]] <- qdo(ft3.qgam.x2, do.ptiles[p1], predict )
}
names(q3) <- paste0('q',do.ptiles*100)  

plot(data01$time, data01$x2, pch=20, cex=.3, main='QGAM fits to x2')
for(p1 in 1:length(do.ptiles)) {
  lines(data01$time, q3[[p1]], col=p1+1, lwd=2)
}
grid()
readline("Stop4 a")



### define seasons from a leap year in the middle
find_doy_for_seasons <- function(DOPLOT=TRUE) {
    iy2  <- which(data01$time>=1999.9 & data01$time<2001.4)
    iy2a <- which(data01$time>=1999.9 & data01$time<2000.8)
    iy2b <- which(data01$time>=2000.4 & data01$time<2001.3)
    clim50  <- predict(ft.gamm1$gam, newdata=data01[iy2,] )
    clim50a <- predict(ft.gamm1$gam, newdata=data01[iy2a,] )
    clim50b <- predict(ft.gamm1$gam, newdata=data01[iy2b,] )

    iclim50a.n <- which.min(clim50a)
    iclim50a.x <- which.max(clim50a)
    iclim50b.n <- which.min(clim50b)
    iclim50b.x <- which.max(clim50b)

    mdoy.winter <- data01$doy[ iy2a[iclim50a.n] ]                                                # 65 
    mdoy.summer <- data01$doy[ iy2a[iclim50a.x] ]                                                # 227
    mdoy.spring <- data01$doy[ iy2a[iclim50a.n] + round((iy2a[iclim50a.x]-iy2a[iclim50a.n])/2) ] # 146
    mdoy.autumn <- data01$doy[ iy2b[iclim50b.n] + round((iy2b[iclim50b.x]-iy2b[iclim50b.n])/2) ] # 329


    diff.w2s <- (mdoy.summer       - mdoy.winter)
    diff.s2w <- (mdoy.winter + 365 - mdoy.summer)

    doy.spring <- (mdoy.spring -round(diff.w2s/4)) : ((mdoy.spring +round(diff.w2s/4)) )
    doy.summer <- (mdoy.summer -round(diff.w2s/4)) : ((mdoy.summer +round(diff.s2w/4)) )
    doy.autumn <- (mdoy.autumn -round(diff.s2w/4)) : ((mdoy.autumn +round(diff.s2w/4)) )
    ig366 <- which(doy.autumn>366) # need to keep leap years
    if(length(ig366)>0) {
        doy.autumn[ig366] <- doy.autumn[ig366] - 366
    }
    doy.winter <- (mdoy.winter -round(diff.s2w/4)) : ((mdoy.winter +round(diff.w2s/4)) )

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

  if(DOPLOT) {
    plot(data01$time[iy2], data01$x[iy2], pch=20, cex=.3, main="Defining seasons from mean annual cycle")
    lines(data01$time[iy2a], clim50a, col=2, lwd=4)
    lines(data01$time[iy2b], clim50b, col=4, lwd=8, lty=2)

    abline(v=data01$time[iy2a][iclim50a.n], col=2, lty=1, lwd=4)
    abline(v=data01$time[iy2a][iclim50a.x], col=2, lty=1, lwd=4)
    abline(v=data01$time[iy2b][iclim50b.n], col=4, lty=2, lwd=8)
    abline(v=data01$time[iy2b][iclim50b.x], col=4, lty=3, lwd=8)

    abline(v=data01$time[iy2a][which(data01$doy[iy2a]==mdoy.spring)], col=3, lty=1)
    abline(v=data01$time[iy2b][which(data01$doy[iy2b]==mdoy.autumn)], col=3, lty=1)
    abline(v=data01$time[iy2b][iclim50b.n], col=1, lty=2, lwd=2)

    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter2)] ],     rep(13,  length(doy.winter2)), pch=20, col=3, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.summer2)] ],     rep(13,  length(doy.summer2)), pch=20, col=2, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.spring2)] ],     rep(13.1,length(doy.spring2)), pch=20, col=4, cex=.3)
    points(data01$time[ iy2b[which(data01$doy[iy2b] %in% doy.autumn2)] ],     rep(13.1,length(doy.autumn2)), pch=20, col=1, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter2)]+366 ], rep(13,  length(doy.winter2)), pch=20, col=3, cex=.3) # 2000 is a leap year
  }

  doy.seasons <- list( spring=doy.spring2, summer=doy.summer2, autumn=doy.autumn2, winter=doy.winter2 )
  return(doy.seasons)
}

#