 # source("def_qgam_mm.R")

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
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"
# st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"

load(st_msref, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)
data01$year <- trunc(data01$time)

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

## gamm 2 where season mean is introduced as random effect

### define seasons from a leap year in the middle
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

    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter2)] ],     rep(13,  length(doy.winter2)), pch=20, col=3, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.summer2)] ],     rep(13,  length(doy.summer2)), pch=20, col=2, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.spring2)] ],     rep(13.1,length(doy.spring2)), pch=20, col=4, cex=.3)
    points(data01$time[ iy2b[which(data01$doy[iy2b] %in% doy.autumn2)] ],     rep(13.1,length(doy.autumn2)), pch=20, col=1, cex=.3)
    points(data01$time[ iy2a[which(data01$doy[iy2a] %in% doy.winter2)]+366 ], rep(13,  length(doy.winter2)), pch=20, col=3, cex=.3) # 2000 is a leap year


### END define seasons from a year in the middle
data01$seasmean <- rep(NA, nrow(data01))
i0 <- which(data01$doy %in% doy.winter2)
data01$seasmean[i0] <- NA
i1 <- which(data01$doy %in% doy.spring2)
data01$seasmean[i1] <- NA
i2 <- which(data01$doy %in% doy.summer2)
data01$seasmean[i2] <- NA
i3 <- which(data01$doy %in% doy.autumn2)
data01$seasmean[i3] <- NA

for (y in unique(data01$year)) {
  iy <- which(data01$year==y)
  data01$seasmean[iy] <- mean(data01$x[iy])
}

## creat a factor for each season in each year
data01$season <- rep(NA, nrow(data01))
for (y in sort(unique(data01$year))) {
  iy   <- which(data01$year==y)
#   iym1 <- which(data01$year==(y-1)) # need for autumn spanning end of year
  i1   <- which(data01$doy[iy] %in% doy.winter2)
  i2   <- which(data01$doy[iy] %in% doy.spring2)
  i3   <- which(data01$doy[iy] %in% doy.summer2)
  i4a  <- sort(which(data01$doy[iy] %in% doy.autumn2[doy.autumn2>260]))
  i4b  <- sort(which(data01$doy[iy] %in% doy.autumn2[doy.autumn2<90]))
  data01$season[iy[i1]]  <- y+0     # paste0("Wi",y)
  data01$season[iy[i2]]  <- y+.25   # paste0("Sp",y)
  data01$season[iy[i3]]  <- y+.5    # paste0("Su",y)
  data01$season[iy[i4a]] <- y+.75   # paste0("Au",y)
  data01$season[iy[i4b]] <- y-1+.75 # paste0("Au",y-1)
}

data01$fseason <- as.factor(data01$season)

# as random effect
ft.gamm2    <- gamm(   fm.gamm1[[1]],    data=data01, random=list(fseason=~1) )
 q.gamm2    <- predict(ft.gamm2$gam,  newdata=data01 )  
plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM fseason random: ',fm.gamm1[[1]][3]))
lines(data01$time, q.gamm2,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()
# seemingly identical to ft.gamm1
readline("Stop2 2")

# as fixed effect fseason
fm.gamm3    <- list(x ~ s(sdoy, bs="cc", k=12) +gmst +fseason, ~ 1 )
ft.gamm3    <- gamm(   fm.gamm3[[1]],    data=data01 )
 q.gamm3    <- predict(ft.gamm3$gam,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM fseason fixed: ',fm.gamm3[[1]][3]))
lines(data01$time, q.gamm3,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  
readline("Stop2 3")

# as fixed effect season BAD no real impact
fm.gamm4    <- list(x ~ s(sdoy, bs="cc", k=12) +gmst +season, ~ 1 )
ft.gamm4    <- gamm(   fm.gamm4[[1]],    data=data01 )
 q.gamm4    <- predict(ft.gamm4$gam,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM season fixed: ',fm.gamm4[[1]][3]))
lines(data01$time, q.gamm4,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  
# seemingly identical to ft.gamm1
readline("Stop2 4")

# as fixed effect seasmean Good but season discontinuities
fm.gamm5    <- list(x ~ s(sdoy, bs="cc", k=12) +gmst +seasmean, ~ 1 )
ft.gamm5    <- gamm(   fm.gamm5[[1]],    data=data01 )
 q.gamm5    <- predict(ft.gamm5$gam,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM seasmean fixed: ',fm.gamm5[[1]][3]))
lines(data01$time, q.gamm5,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  
readline("Stop2 5")

# as theo effect seasmean BAD
data01$aseasmean <- data01$seasmean - mean(data01$seasmean)
fm.gamm6    <- list(x ~ s(sdoy, bs="cc", k=12, by=aseasmean) +gmst, ~ 1 )
ft.gamm6    <- gamm(   fm.gamm6[[1]],    data=data01 )
 q.gamm6    <- predict(ft.gamm6$gam,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM by=aseasmean: ',fm.gamm6[[1]][3]))
lines(data01$time, q.gamm6,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  # TOTAL FAIL
readline("Stop2 6")

# doy-seasmean interaction, Prob best of this set but still season discontinuities
fm.gamm7    <- list(x ~ s(sdoy, bs="cc", k=12) +gmst + te(doy,aseasmean,bs=c('cc','cc')), ~ 1 )
ft.gamm7    <- gamm(   fm.gamm7[[1]],    data=data01 )
 q.gamm7    <- predict(ft.gamm7$gam,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM te(doy,aseasmean): ',fm.gamm7[[1]][3]))
lines(data01$time, q.gamm7,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  
readline("Stop2 3d")

# # theo factored seasonal mean
# # VERY SLOW
# fm.gamm4    <- list(x ~ s(sdoy, bs="cc", k=12, by=fseason) +gmst, ~ 1 )
# ft.gamm4    <- gamm(   fm.gamm4[[1]],    data=data01 )
#  q.gamm4    <- predict(ft.gamm4$gam,  newdata=data01 )  
#  plot(data01$time, data01$x, pch=20, cex=.3)
# lines(data01$time, q.gamm4,  col=2)
# grid()  
# readline("Stop2 4")

iy <- which(data01$year>=2006 & data01$year<=2008)
plot(data01$time[iy], data01$x[iy], pch=20, cex=.3)
lines(data01$time[iy], q.gamm1[iy], col=2)
lines(data01$time[iy], q.gamm2[iy], col=3)
lines(data01$time[iy], q.gamm3[iy], col=4)
# lines(data01$time[iy], q.gamm4[iy], col=4)
grid()
readline("Stop2 5")

iy <- which(data01$year>=2006 & data01$year<=2008)



data01$seasmean <- rep(NA, nrow(data01))
i0 <- which(data01$doy %in% doy.winter2)
data01$seasmean[i0] <- NA
i1 <- which(data01$doy %in% doy.spring2)
data01$seasmean[i1] <- NA
i2 <- which(data01$doy %in% doy.summer2)
data01$seasmean[i2] <- NA
i3 <- which(data01$doy %in% doy.autumn2)
data01$seasmean[i3] <- NA

for (y in unique(data01$year)) {
  iy <- which(data01$year==y)
  data01$seasmean[iy] <- mean(data01$x[iy])
}


### 2026.02.04 Theo's new random effect
library(data.table)
data01 <- data.table(data01)
# introduce years
data01[,year := trunc(time,0)]
# factor year
data01[,fYear := factor(year)]

ft.theo1 <- gam(x ~ s(sdoy,bs="cc") # seasonal cycle
                  + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")),
                    knots=list(sdoy=c(0,1)),data=data01
)
q.theo1    <- predict(ft.theo1,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM te(doy,aseasmean): ',ft.theo1$formula[3]))
lines(data01$time, q.theo1,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  
readline("Stop2 2")

ft.theo2 <- gam(x ~ s(sdoy,bs="cc") # seasonal cycle
                  + s(fYear,bs="sz",xt=list(bs="cc")),
                    knots=list(sdoy=c(0,1)),data=data01
)
ft.theo2b <- gam(x ~ s(sdoy,bs="cc") # seasonal cycle
                  + s(fYear,bs="sz"),data=data01
)
q.theo2    <- predict(ft.theo2,  newdata=data01 )  
q.theo2b   <- predict(ft.theo2b, newdata=data01 )  
q.theo2c   <- predict(ft.theo2c, newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM te(doy,aseasmean): ',ft.theo2$formula[3]))
lines(data01$time, q.theo2,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()  

 plot(data01$doy, data01$x, pch=20, cex=.3, main=paste('MM te(doy,aseasmean): ',ft.theo2$formula[3]))
lines(data01$doy, q.theo2,  col=2)

ft.regam1   <- gam(x ~ s(sdoy, bs="cc", k=12) + gmst + s(fYear, bs="re"), data=data01 )
 q.regam1    <- predict(ft.regam1,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM fseason random: ',fm.gamm1[[1]][3]))
lines(data01$time, q.regam1,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()

ft.regam2   <- gam(x ~ s(sdoy, bs="cc", k=12) + gmst + s(fseason, bs="re"), data=data01 )
 q.regam2    <- predict(ft.regam2,  newdata=data01 )  
 plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM fseason random: ',fm.gamm1[[1]][3]))
lines(data01$time, q.regam1,  col=4)
lines(data01$time, q.regam2,  col=2)
lines(data01$time, qx.gamm1, col=3)
grid()

readline("Stop2 2")
 plot(data01$doy, data01$x, pch=20, cex=.3, main=paste('MM fseason random: ',fm.gamm1[[1]][3]))
lines(data01$doy, q.regam1,  col=4)
# lines(data01$doy, q.regam2,  col=2)
grid()


# 
plot(data01$doy, data01$x, pch=20, cex=.3, main=paste('MM te(doy,aseasmean): ',ft.theo2$formula[3]), ty='n')
lines(data01$doy, q.theo2,  col=2)

iy2 <- which(data01$year %in% 2000:2001)
y   <- seq_along(iy2)
plot(y, data01[iy2,x],pch=46, xlim=c(0,720), ylim=range(data01$x))
for(y1 in 1980:2000){
    iy2 <- which(data01$year %in% y1:(y1+1))
    y   <- seq_along(iy2)
    points(y, data01[iy2,x],pch=46)
    lines(y,q.theo2[iy2],col="red")
}
for(y1 in 1980:2000){
    iy2 <- which(data01$year %in% y1:(y1+1))
    y   <- seq_along(iy2)
    lines(y,q.theo2c[iy2],col=3)
}

### notes
# q.theo2 & q.theo2b very similar to q.regam1
# both show discontinuities at year ends









ft.theo1  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo2  <- gam(x ~ s(sdoy,bs="cc")              + s(fYear,     bs="sz",     xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo2b <- gam(x ~ s(sdoy,bs="cc")              + s(fYear,     bs="sz"),                                               data=data01)
ft.theo2c <- gam(x ~ s(sdoy,bs="cc")              + s(fYear,     bs="re"),                                               data=data01)
ft.regam1 <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(fYear,     bs="re"),                                               data=data01)

