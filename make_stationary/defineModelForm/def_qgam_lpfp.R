# source("def_qgam_lpfp.R")

# 2025-12-19
# try using the previous n days as covariate los pass filter predictor "lpfp"
# need to know the temperature anomaly with respect to the climatological previous n days
# eg
# pre-sst-cov = filter(actual(previous n days)) - filter(climatological(previous n days))

kern          <- rep(1/kern.len, kern.len)



source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
source("../../libs/fn_MakeStationary.R")
source("../../libs/fn_HotDays.R")

st_msref  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-28.RData"


load(st_msref, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)
data01$year <- trunc(data01$time)
data01$x0   <- data01$x

do.season <- "summer"  # "spring" "summer" "autumn" "winter"
do.ptiles         <- seq(from=1, to=99, by=1)/100.0     
do.ptiles         <- seq(from=10, to=90, by=10)/100.0     

###############################################################################################
### step 1: fit mixed model to get mean gmst term and baseline annual cycle
###############################################################################################
# this does not work ft.gmst    <- gamm(   x ~ gmst,    data=data01 , random=list(stime=~1))
ft.gmst0 <- lm(x0 ~ gmst,        data=data01)
 q.gmst0 <- predict(ft.gmst0, newdata=data01) 

fm.gamm1    <- list(x0 ~ s(sdoy, bs="cc", k=12) +gmst, ~ 1 )
ft.gamm1    <- gamm(   fm.gamm1[[1]],    data=data01 )
 q.gamm1    <- predict(ft.gamm1$gam,  newdata=data01 ) 
     plot(data01$time, data01$x0, pch=20, cex=.3, main='x0 s(doy,gmst)')
    lines(data01$time, q.gamm1, col=2)
    readline("Stop 0 ")
data01$x1   <- data01$x0 - q.gamm1


### calculate low pass filter of sst
i0 <- 1:(365*7)
# kern          <- rep(1/kern.len, kern.len)
kern.len <- seq(15,180, by=15)
lp.q1 <- list()
lp.x0 <- list()
lp.dq1 <- list()
plot(data01$time[i0], data01$x0[i0], pch=20, cex=.3, main='sst & lpfp')
for(i in seq_along(kern.len)){
  kern        <- rep(1/kern.len[i], kern.len[i])
  lp.q1[[i]]  <- stats::filter(q.gamm1,   filter=kern, sides=1)
  lp.x0[[i]]  <- stats::filter(data01$x0, filter=kern, sides=1)
  lp.dq1[[i]] <- lp.x0[[i]] - lp.q1[[i]]
  points(data01$time[i0], lp.dq1[[i]][i0]+12,col=i, pch=3, cex=.4)
} 

# i2 <- i0 + floor(kern.len/2)
# i1 <- i2     + floor(kern.len/2) +1

HERE 23.12.2025

data01$lpfp <- lp.dq1[[which(kern.len==30)]]

fm.gamm2    <- list(x0 ~ s(sdoy, bs="cc", k=12) +gmst +s(lpfp), ~ 1 )
ft.gamm2    <- gamm(   fm.gamm2[[1]],    data=data01 )
 q.gamm2    <- predict(ft.gamm2$gam,  newdata=data01 ) 
 plot(data01$time, data01$x0, pch=20, cex=.3, main='sst & lpfp')
lines(data01$time,   q.gamm2, col=2)

q.gamm2  <- list()
fm.gamm2 <- list(x0 ~ s(sdoy, bs="cc", k=12) +gmst +s(lpfp), ~ 1 )
 plot(data01$time, data01$x0, pch=20, cex=.3, main='sst & lpfp')
for(i in seq_along(kern.len)){
  data01$lpfp <- lp.dq1[[i]]
  ft.gamm2    <- gamm(   fm.gamm2[[1]],    data=data01 )
  q.gamm2[[i]]<- predict(ft.gamm2$gam,  newdata=data01 ) 
  lines(data01$time,   q.gamm2[[i]], col=i+1)
}

 plot(data01$time, data01$x0, pch=20, cex=.3, main='sst & lpfp')
lines(data01$time,   q.gamm1, col=2)
lines(data01$time,   q.gamm2[[1]], col=3)
lines(data01$time,   q.gamm2[[4]], col=4)
lines(data01$time,   q.gamm2[[8]], col='orange3')


###############################################################################################
### define seasons
    doy.seasons  <- find_doy_for_seasons(data01, ft.gamm1)
    doy.midseas  <- NULL
    sdoy.midseas <- NULL
    i2001        <- which(data01$year==2001)
    name.seas    <- names(doy.seasons)
    iwhich.seas  <- which(name.seas==do.season)

    iseas                                    <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='winter')]])
    doy.midseas[which(name.seas=='winter')] <-  data01$doy[i2001][iseas][which.min(q.gamm1[i2001][iseas])]
    sdoy.midseas[which(name.seas=='winter')] <- data01$sdoy[i2001][iseas][which.min(q.gamm1[i2001][iseas])]
    iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='spring')]])
    doy.midseas[which(name.seas=='spring')] <- round(median(data01$doy[i2001][iseas]))
    sdoy.midseas[which(name.seas=='spring')] <- median(data01$sdoy[i2001][iseas])
    iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='summer')]])
    doy.midseas[which(name.seas=='summer')] <-  data01$doy[i2001][iseas][which.max(q.gamm1[i2001][iseas])]
    sdoy.midseas[which(name.seas=='summer')] <- data01$sdoy[i2001][iseas][which.max(q.gamm1[i2001][iseas])]
    iseas       <- which(data01$doy[i2001] %in% doy.seasons[[which(name.seas=='autumn')]])
    doy.midseas[which(name.seas=='autumn')] <- round(median(data01$doy[i2001][iseas]))
    sdoy.midseas[which(name.seas=='autumn')] <- median(data01$sdoy[i2001][iseas])

    names(doy.midseas)  <- names(doy.seasons)
    names(sdoy.midseas) <- names(doy.seasons)
    print("mid-seasons DOY:")
    print(doy.midseas)
    print(sdoy.midseas)

    # iseas       <- which(data01$doy %in% doy.seasons[[do.season]])
    #  plot(data01$time[iseas],data01$x0[iseas], pch=20, cex=.3, main=paste('seasonal fit:',do.season))
    # lines(data01$time[iseas],q.gamm1[iseas], col=2)

###############################################################################################

### step 1b: get seasonal inter-annual variability term
# xsdoy         <- sdoy.midseas[iwhich.seas]
# q.gamm1.xsdoy <- predict(ft.gamm1$gam,  newdata=data.frame(gmst=data01$gmst, sdoy=rep(xsdoy,length(data01$sdoy)) ))

    #  plot(data01$time, data01$x0, pch=20, cex=.3, main='gmst')
    # lines(data01$time, q.gmst0, col=2)
    # lines(data01$time, q.gamm1.xsdoy, col=3)
    # readline("Stop1")

    #  plot(q.gmst0 -q.gamm1.xsdoy, pch=20, cex=.3, main='gmst')
    # readline("Stop1")

# data01$x1   <- data01$x0 - q.gamm1.xsdoy
# OR


###############################################################################################
### step 2: fit to seasonal anomalies all do.ptiles
###############################################################################################
iseas   <- which(data01$doy %in% doy.seasons[[do.season]])
seasd01 <- data01[iseas,]
plot(seasd01$time, seasd01$x1, pch=20, cex=.3, main='x1: x0 - q.gamm1')
readline("Stop 1 1")
plot(seasd01$doy,  seasd01$x1, pch=20, cex=.3, main='x1: x0 - q.gamm1')
readline("Stop 1 2")

## 2a fit stright to seasonal anomalies from q.gamm1
fm.x1 <- list(x1 ~ s(sdoy, bs="tp", k=12) , ~ 1 )
ft.x1 <- mqgam(   fm.x1[[1]],    data=seasd01, qu=do.ptiles)
 q.x1 <- qdo(ft.x1, qu=do.ptiles,  predict ) 

     plot(seasd01$time,  seasd01$x1, pch=20, cex=.3, main='x1')
    lines(seasd01$time, q.x1[[1]], col=2)
    lines(seasd01$time, q.x1[[2]], col=2)
readline("Stop 1 3")

     plot(seasd01$doy,  seasd01$x1, pch=20, cex=.3, main='x1')
    lines(seasd01$doy, q.x1[[1]], col=2)
    lines(seasd01$doy, q.x1[[2]], col=2)
readline("Stop 1 4")


## 2b fit to x1 relative to seas mean 
## subtract seas mean for each year -> x1b
x1b    <- seasd01$x1
seas.y <- trunc(intann.seas.anom[[iwhich.seas+4]])
for(y in seas.y){
    iy      <- which(seasd01$year==y)
    x1b[iy] <- seasd01$x1[iy] - mean(seasd01$x1[iy])
}
    plot(seasd01$time, x1b, pch=20, cex=.3, main='x1b')
readline("Stop 2 1")
    plot(seasd01$doy,  x1b, pch=20, cex=.3, main='x1b')
readline("Stop 2 2")
seasd01$x1b <- x1b
seasd01$x2 <- x1b

fm.x1b <- list(x1b ~ s(sdoy, bs="tp", k=12) , ~ 1 )
ft.x1b <- mqgam(   fm.x1b[[1]],    data=seasd01, qu=do.ptiles)
 q.x1b <- qdo(ft.x1b, qu=do.ptiles,  predict ) 

     plot(seasd01$time,  seasd01$x1b, pch=20, cex=.3, main='x1b: x1 - seas_mean')
    # for(i in seq(10,90,by=10)) lines(seasd01$time, q.x1b[[i]], col=2)
    for(i in seq(1,9,by=1)) lines(seasd01$time, q.x1b[[i]], col=2)
readline("Stop 2 3")

     plot(seasd01$doy,  seasd01$x1b, pch=20, cex=.3, main='x1')
    pdoy <- sort(unique(seasd01$doy))
    # for(i in seq(10,90,by=10)) lines(pdoy, q.x1b[[i]][seq_along(pdoy)], col=2)
    for(i in seq(1,9,by=1)) lines(pdoy, q.x1b[[i]][seq_along(pdoy)], col=2)
readline("Stop 2 4")


###############################################################################################
### step EVgam: fit to x1 & x1b
###############################################################################################
p.rl50 <- 1/nrow(seasd01) # 1/(50*length(doy.seasons[[do.season]]))
# calc MSgpd
fm.gpdS0G0 <- list(excess ~ 1 ,                     ~ 1 )
fm.gpdSDG0 <- list(excess ~ s(sdoy, bs='tp',k=12) , ~ 1 )

th.x1.gpd      <- q.x1[[9]] # q.x1[[90]]
seasd01$excess <- seasd01$x1 - th.x1.gpd
seasd01.gpd    <- subset(seasd01, excess > 0)
      plot(seasd01$time,     seasd01$x1, pch=20, cex=.3, main='x1 gpd excess')
    points(seasd01$time,     th.x1.gpd, col=3, pch=46)
    points(seasd01.gpd$time, seasd01.gpd$x1, col=2, cex=.3)
ft.x1gpdS0G0     <- evgam(fm.gpdS0G0, seasd01.gpd, family="gpd", trace=2)
ft.x1gpdSDG0     <- evgam(fm.gpdSDG0, seasd01.gpd, family="gpd", trace=2)
      plot(seasd01$time,     seasd01$x1, pch=20, cex=.3, main='x1 gpd excess')
    points(seasd01$time,     th.x1.gpd, col=3, pch=46)
    gpdpar   <- predict(ft.x1gpdSDG0, newdata=seasd01, type="response")
    qgpdqs01 <- th.x1.gpd + gpdpar$scale * (p.rl50^(-gpdpar$shape) -1)/gpdpar$shape
    points(seasd01$time,     qgpdqs01, col=4, pch=20, cex=.3)

      plot(seasd01$doy,     seasd01$x1, pch=20, cex=.3, main='x1 gpd excess')
    points(seasd01$doy,     qgpdqs01, col=4, pch=20, cex=.3)

      plot(seasd01$time,     seasd01$x1-qgpdqs01, pch=20, cex=.3, main='x1 gpd excess')
      abline(h=0)

th.x1b.gpd     <- q.x1b[[9]] # q.x1b[[90]]
seasd01$excess <- seasd01$x1b - th.x1b.gpd
seasd01.gpd    <- subset(seasd01, excess > 0)
      plot(seasd01$time,     seasd01$x1b, pch=20, cex=.3, main='x1 gpd excess')
    points(seasd01$time,     th.x1b.gpd, col=3, pch=46)
    points(seasd01.gpd$time, seasd01.gpd$x1b, col=2, cex=.3)
ft.x1bgpdS0G0  <- evgam(fm.gpdS0G0, seasd01.gpd, family="gpd", trace=2)
ft.x1bgpdSDG0  <- evgam(fm.gpdSDG0, seasd01.gpd, family="gpd", trace=2)
      plot(seasd01$time,     seasd01$x1b, pch=20, cex=.3, main='x1b gpd excess')
    points(seasd01$time,     th.x1b.gpd, col=3, pch=46)
    gpdpar   <- predict(ft.x1bgpdSDG0, newdata=seasd01, type="response")
    qgpdqs01 <- th.x1b.gpd + gpdpar$scale * (p.rl50^(-gpdpar$shape) -1)/gpdpar$shape
    points(seasd01$time,     qgpdqs01, col=4, pch=20, cex=.3)

      plot(seasd01$doy,     seasd01$x1b, pch=20, cex=.3, main='x1b gpd excess', ylim=range(c(seasd01$x1b,qgpdqs01)))
    points(seasd01$doy,     qgpdqs01, col=4, pch=20, cex=.3)

      plot(seasd01$time,     seasd01$x1b-qgpdqs01, pch=20, cex=.3, main='x1b gpd excess')
      abline(h=0, col=2)


### so I think we have a reasonable approach now to capture seasonal inter-annual variability


### plot roundup
st_pdf  <- paste0("def_qgam_season_01.pdf")
if(SAVEPLOTS) pdf(st_pdf, width=10, height=7)

 plot(data01$time, data01$x0, pch=20, cex=.3, main='NWS x0 gamm(doy,gmst)')
lines(data01$time, q.gamm1, col=2)
grid()

 plot(data01$time[iseas], data01$x0[iseas], pch=20, cex=.3, main='NWS x0 gamm(doy,gmst) summer')
lines(data01$time[iseas],   q.gamm1[iseas], col=2)
grid()

plot(seasd01$time, seasd01$x1, pch=20, cex=.3, main='NWS x1: x0 - gamm(doy,gmst) summer')
abline(h=0)
grid()

plot(seasd01$time,  seasd01$x2, pch=20, cex=.3, main='NWS x2: x1 - seas_mean(y)')
abline(h=0)
grid()

plot(seasd01$doy,  seasd01$x2, pch=20, cex=.3, main='NWS x2: x1 - seas_mean(y)')
abline(h=0)
grid()

dev.off()


#