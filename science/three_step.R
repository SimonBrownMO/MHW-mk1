# source("three_step.R")

library(mgcv)
library(data.table)

# load the data
load("data01_NWS.RData")
data01 <- data.table(data01)

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

    ft1.gam <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst, data=data01)
    q.gam1  <- predict(ft1.gam) 
    nd$sdoy <- sdoyx
    q.gam1x <- predict(ft1.gam,  newdata=nd) 

    # these are identical
###

### step 2: fit gam based on time to residuals to isolate greater than seasonal variability
data01$x1 <- data01$x - q.gam1
# ft1.gamb <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(stime,     bs="re"), data=data01)
# q.gam1b  <- predict(ft1.gamb) 

# ft2.gam <- gam(x1 ~ s(sdoy,bs="cc", k=12) +s(time, k=128) , data=data01)
ft2.gam <- gam(x1 ~ s(stime, k=128) , data=data01)
sp.gam2 <- smooth.spline(data01$stime, data01$x1, spar=0.3)
q.gam2  <- predict(ft2.gam) 
  plot(data01$time, data01$x1, pch=20, cex=.3, main="x1")
points(data01$time, data01$sdoy-0.5,   col=5, pch=20, cex=.3)
 lines(data01$time, sp.gam2$y, col=3, lwd=4)
 lines(data01$time, q.gam2,    col=2, lwd=4)
abline(h=0)

    ### try fixing knots
    i1 <- which(data01$doy==doyn | data01$doy==doyx)
    knots1 =list(time=data01$time[i1])
    ft2b.gam <- gam(x1 ~ s(time, k=128), knots=knots1 , data=data01)

### step 3: qgam to residuals to isolate low frequency interannual variability
data01$x2 <- data01$x1 - q.gam2
ft3.qgam   <- mqgam(list(x2 ~ s(sdoy, bs='cc', k=12)  +gmst +ti(sdoy,gmst, bs=c('cc','tp')), ~ s(sdoy)), data=data01, qu=c(0.1,0.5,0.9))
q.qgam3    <- qdo(ft3.qgam, 0.5, predict)
q10.qgam3  <- qdo(ft3.qgam, 0.1, predict)
q90.qgam3  <- qdo(ft3.qgam, 0.9, predict)
 plot(data01$time, data01$x2, pch=20, cex=.3, main="x2")
lines(data01$time, q.qgam3,     col=4, lwd=1)
lines(data01$time, q10.qgam3,   col=2, lwd=1)
lines(data01$time, q90.qgam3,   col=2, lwd=1)
abline(h=0) 



### hardcopy
pdf("three_step.pdf", width=14, height=10)

# fit1
 plot(data01$time, data01$x, pch=20, cex=.3, main='x: gam(x ~ s(sdoy,bs=cc, k=12) + gmst')
lines(data01$time, q.gam1,    col=2, lwd=2)
lines(data01$time, q.gam1x,   col=3, lwd=2)
abline(h=mean(data01$x), col=1, lwd=2)
grid()

# fit2
 plot(data01$time, data01$x1, pch=20, cex=.3, main="x1: x1 ~ s(time, k=128)")
lines(data01$time, q.gam2,    col=2, lwd=2)
points(data01$time, data01$sdoy-2,   col=5, pch=20, cex=.1)
abline(h=0)
grid()  

# fit3
 plot(data01$time, data01$x2, pch=20, cex=.3, main="x2")
lines(data01$time, q.qgam3,     col=4, lwd=1)
lines(data01$time, q10.qgam3,   col=2, lwd=1)
lines(data01$time, q90.qgam3,   col=2, lwd=1)
abline(h=0) 
grid()
legend('topleft', legend=c("90%ile", "50%ile", "10%ile"), col=c(2,4,2), lwd=2)

dev.off()


# 1 2 and 3 in one go

ft1.gam  <- gam(       x  ~ s(sdoy,bs="cc", k=12) + gmst, data=data01)
ft2.gam  <- gam(       x1 ~ s(time, k=128) , data=data01)
ft3.qgam <- mqgam(list(x2 ~ s(sdoy, bs='cc', k=12)  +gmst +ti(sdoy,gmst, bs=c('cc','tp')), ~ s(sdoy)), data=data01, qu=c(0.1,0.5,0.9))




ft123.qgam <- mqgam(list(x ~ s(sdoy, bs='cc', k=12) +s(gmst,k=4) +ti(sdoy,gmst, bs=c('cc','tp'), k=c(12,4))+s(stime, k=128), ~ s(sdoy)), data=data01, qu=c(0.1,0.5,0.9))

qdo(ft123.qgam, 0.5, plot, page=1)
  q.123.qgam  <- qdo(ft123.qgam, 0.5, predict)
q10.123.qgam  <- qdo(ft123.qgam, 0.1, predict)
q90.123.qgam  <- qdo(ft123.qgam, 0.9, predict)

    nd       <- data01
    nd$doy   <- doyx
    nd$sdoy  <- sdoyx
    nd$time  <- data01$time[which(data01$doy==doyx)[1]]

  qx.123.qgam  <- qdo(ft123.qgam, 0.5, predict, newdata=nd)
q10x.123.qgam  <- qdo(ft123.qgam, 0.1, predict, newdata=nd)
q90x.123.qgam  <- qdo(ft123.qgam, 0.9, predict, newdata=nd)

 plot(data01$time, data01$x, pch=20, cex=.3, main="x")
lines(data01$time,   q.123.qgam,   col=4, lwd=1)
# lines(data01$time, q10.123.qgam,   col=2, lwd=1)
# lines(data01$time, q90.123.qgam,   col=2, lwd=1)
lines(data01$time,    qx.123.qgam, col=3, lwd=2)
lines(data01$time,  q90x.123.qgam, col=2, lwd=2)
lines(data01$time,  q10x.123.qgam, col=2, lwd=2)
abline(h=0) 
grid()
legend('topleft', legend=c("90%ile", "50%ile", "10%ile"), col=c(2,4,2), lwd=2)



    nd       <- data01
    nd$doy   <- doyx
    nd$sdoy  <- sdoyx
    # nd$time  <- data01$time[which(data01$doy==doyx)[1]]
    nd$gmst  <- tail(data01$gmst[which(data01$doy==doyx)],1)

  qx.123.qgam  <- qdo(ft123.qgam, 0.5, predict, newdata=nd)
plot(data01$time,    qx.123.qgam, col=3, lwd=2)
crazy




### theo's approach
ft.regam2 <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
