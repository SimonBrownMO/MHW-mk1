# source("theo_annual_random_effect_qgam.R")

# trying to bring random effects to qgam


source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
c21 <- distinct21colours() 
source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")
library(evgam,lib="/home/users/simon.brown/code/R/libs/R-4.4.1-2024_12_04/lib/R/library")
library(mgcv)
library(qgam)
library(data.table)
source("../libs/fn_MakeStationary.R")
source("../libs/fn_HotDays.R")

st_msref <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"

load(st_msref, verb=TRUE)
list2env(MSconfig, envir = .GlobalEnv)
list2env(MSconfig$files, envir = .GlobalEnv)

load(st_preproc, verb=TRUE)
load(st_msdata01, verb=TRUE)
data01 <- data.table(data01)
data01[,year := trunc(time,0)]
# factor year
data01[,fYear := factor(year)]

### step 1: fit mixed model to get mean gmst term and baseline annual cycle
    q3         <- c(0.1,0.5,0.9)
    fm.gamm1    <- list(x ~ s(sdoy, bs="cc", k=12) + s(gmst,k=4), ~ 1 )
    ft.gamm1    <- gamm(   fm.gamm1[[1]],    data=data01 )
        # Approximate significance of smooth terms:
        #           edf Ref.df     F p-value    
        # s(sdoy) 9.965 10.000 49123  <2e-16 ***
        # s(gmst) 2.993  2.993  2288  <2e-16 ***
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
readline("Stop 1")
###

ft.theo1  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo2  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",     xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo3  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1                 ), knots=list(sdoy=c(0,1)),data=data01)
ft.theo4  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")),                         data=data01)
ft.theo5  <- gam(x ~ s(sdoy)                      + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
q.theo1  <- predict(ft.theo1 )
q.theo2  <- predict(ft.theo2 )
q.theo3  <- predict(ft.theo3 )
q.theo4  <- predict(ft.theo4 )
q.theo5  <- predict(ft.theo5 )
summary(diff(q.theo1))

iy <- which(data01$year %in% 1992:1999)
c21 <- distinct21colours() 

 plot(data01$time[iy], data01$x[iy], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time[iy],  q.theo1[iy], col=c21[1],lwd=2, lty=1)
lines(data01$time[iy],  q.theo2[iy], col=c21[2],lwd=4, lty=3)
lines(data01$time[iy],  q.theo3[iy], col=c21[3],lwd=2, lty=1)
lines(data01$time[iy],  q.theo4[iy], col=c21[4],lwd=4, lty=3)
lines(data01$time[iy],  q.theo5[iy], col=c21[5],lwd=4, lty=1)
legend("topleft", legend=c("theo1","theo2","theo3","theo4","theo5"), col=c21[1:5], lwd=2)
grid()  

# NB 2024.doy 123 is missing
  plot(diff(q.theo1), pch=1, cex=.7, ylim=c(-1,1), main='diff(q.theoN)')
points(diff(q.theo2), pch=2, cex=.7, col=c21[2])
points(diff(q.theo3), pch=3, cex=.7, col=c21[3])
points(diff(q.theo4), pch=4, cex=.7, col=c21[4])
points(diff(q.theo5), pch=5, cex=.7, col=c21[5])
legend("topleft", legend=c("theo1","theo2","theo3","theo4","theo5"), col=c21[1:5], lwd=2)
grid()

# find minimum step jump at doy=[365,366]
i2 <- which(abs(diff(data01$doy)) > 2)
# points(i2,diff(q.theo1)[i2],col=2)
                                 #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
summary(abs(diff(q.theo1))[i2])  # 0.00869  0.14672  0.23499  0.28272  0.41018  0.95348
summary(abs(diff(q.theo2))[i2])  # 0.001102 0.140721 0.230678 0.279734 0.414310 0.946467
summary(abs(diff(q.theo3))[i2])  # 0.007724 0.161290 0.231071 0.288168 0.386632 1.086348
summary(abs(diff(q.theo4))[i2])  # 0.007344 0.117452 0.221716 0.282699 0.439797 0.853300
summary(abs(diff(q.theo5))[i2])  # 0.0203   0.1100   0.2391   0.2841   0.4201   0.8271
# not much in it but theo3 probably the best

### refine theo3 -> cant
ft.theo3   <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1), knots=list(sdoy=c(0,1)),data=data01)
ft.theo3b  <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz")     , knots=list(sdoy=c(0,1)),data=data01)
ft.theo3c  <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1),                         data=data01)
## all differences tiny. 3b slow
ft.theo3d  <- gam(x ~ s(sdoy,bs="cc",k=6) + s(sdoy,fYear,bs="sz",id=1),       knots=list(sdoy=c(0,1)),data=data01)
ft.theo3d  <- gam(x ~ s(sdoy,bs="cc",k=24) + s(sdoy,fYear,bs="sz",id=1, k=6), knots=list(sdoy=c(0,1)),data=data01)
## no amount of fiddling with k values seems to improve anything
q.theo3  <- predict(ft.theo3 )
q.theo3b <- predict(ft.theo3b)
q.theo3c <- predict(ft.theo3c)
q.theo3d <- predict(ft.theo3d)
plot(q.theo3 - q.theo3d, pch=20, cex=.3, main="theo3 - theo3d")
 plot(data01$time[iy],  data01$x[iy], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time[iy],   q.theo3[iy], col=c21[3],lwd=2, lty=1)
# lines(data01$time[iy],  q.theo3b[iy], col=c21[4],lwd=2, lty=1)
# lines(data01$time[iy],  q.theo3c[iy], col=c21[5],lwd=2, lty=1)
lines(data01$time[iy],  q.theo3d[iy], col=c21[5],lwd=2, lty=1)
legend("topleft", legend=c("theo3","theo3d"), col=c(c21[3],c21[5]), lwd=2)

### add gmst to theo3
k.gmst <- sort(c(range(data01$gmst),quantile(data01$gmst, c(.33,.66))))
ft.theo4a <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1) + gmst       , knots=list(sdoy=c(0,1),gmst=k.gmst),data=data01)
ft.theo4b <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1) + s(gmst,k=4), knots=list(sdoy=c(0,1),gmst=k.gmst),data=data01)
# ft.theo4c - not really a viable solution although seemingly has solved the negative gmst issue
ft.theo4c <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1) + ti(sdoy,gmst,bs=c('cc','tp')), knots=list(sdoy=c(0,1),gmst=k.gmst),data=data01)
# ft.theo4d - just screwy
ft.theo4d <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1) + te(sdoy,gmst,bs=c('cc','tp')), knots=list(sdoy=c(0,1),gmst=k.gmst),data=data01)
q.theo4a  <- predict(ft.theo4a )
q.theo4b  <- predict(ft.theo4b )
q.theo4c  <- predict(ft.theo4c )
q.theo4d  <- predict(ft.theo4d )

nd        <- data.table(sdoy=data01$sdoy, gmst=data01$gmst, fYear="2000")
q.theo4a.fixy <- predict(ft.theo4a,   newdata=nd)
q.theo4b.fixy <- predict(ft.theo4b,   newdata=nd)
q.theo4c.fixy <- predict(ft.theo4c,   newdata=nd)
q.theo4d.fixy <- predict(ft.theo4d,   newdata=nd)
ix        <- which(data01$doy==doyx)
qx.theo4a.fixy <- predict(ft.theo4a,  newdata=nd[ix,])
qx.theo4b.fixy <- predict(ft.theo4b,  newdata=nd[ix,])
qx.theo4c.fixy <- predict(ft.theo4c,  newdata=nd[ix,])
qx.theo4d.fixy <- predict(ft.theo4d,  newdata=nd[ix,])

nd             <- data.table(sdoy=data01$sdoy, gmst=k.gmst[1], fYear=data01$fYear)
q.theo4a.fixg0 <- predict(ft.theo4a,   newdata=nd)
q.theo4b.fixg0 <- predict(ft.theo4b,   newdata=nd)
q.theo4c.fixg0 <- predict(ft.theo4c,   newdata=nd)
q.theo4d.fixg0 <- predict(ft.theo4d,   newdata=nd)
qx.theo4a.fixg0 <- predict(ft.theo4a,   newdata=nd[ix,])
qx.theo4b.fixg0 <- predict(ft.theo4b,   newdata=nd[ix,])
qx.theo4c.fixg0 <- predict(ft.theo4c,   newdata=nd[ix,])
qx.theo4d.fixg0 <- predict(ft.theo4d,   newdata=nd[ix,])
nd             <- data.table(sdoy=data01$sdoy, gmst=k.gmst[4], fYear=data01$fYear)
q.theo4a.fixg1 <- predict(ft.theo4a,   newdata=nd)
q.theo4b.fixg1 <- predict(ft.theo4b,   newdata=nd)
q.theo4c.fixg1 <- predict(ft.theo4c,   newdata=nd)
q.theo4d.fixg1 <- predict(ft.theo4d,   newdata=nd)

# to all data01
 plot(data01$time,       data01$x, pch=20, cex=.3, main=paste(ft.theo3$formula))
lines(data01$time,        q.theo3, col=c21[3],lwd=2, lty=1)
lines(data01$time,       q.theo4a, col=c21[5],lwd=2, lty=1)
lines(data01$time,       q.theo4c, col=c21[6],lwd=2, lty=1)
lines(data01$time,       q.theo4d, col=c21[7],lwd=2, lty=1)


 plot(data01$time,       data01$x,      pch=20, cex=.3, main=paste(ft.theo4b$formula[3]), ylim=range(c(q.theo4b.fixg0, q.theo4b.fixy)), ty='n')
lines(data01$time,      q.theo4b.fixy,  col=c21[6],lwd=2, lty=1)
lines(data01$time,      q.theo4b.fixg0, col=c21[8],lwd=2, lty=1)
grid()
legend("topleft", legend=c("q.theo4b.fixy","q.theo4b.fixg0"), col=c21[c(6,8)], lwd=2)

plot(q.theo3 - q.theo4a, pch=20, cex=.3, main="theo3 - theo4a")
plot(q.theo4a.g1 - q.theo4a.g0, pch=20, cex=.3, main="theo3 - theo4a")

 plot(data01$time,       data01$x,      pch=20, cex=.3, main=paste(ft.theo4c$formula[3]))
lines(data01$time,      q.theo4c.fixy,  col=c21[6],lwd=2, lty=1)
lines(data01$time,      q.theo4c.fixg0, col=c21[8],lwd=2, lty=1)
grid()
legend("topleft", legend=c("q.theo4c.fixy","q.theo4c.fixg0"), col=c21[c(6,8)], lwd=2)

plot(data01$time,       data01$x,      pch=20, cex=.3, main=paste(ft.theo4c$formula[3]), ty='n')
plot(data01$time,      q.theo4c-q.theo4c.fixy,  col=c21[6],lwd=2, lty=1, ty='l')

# ADDING THE HIDDEN PROCESS CAUSES GMST TERM TO CHANGE SIGN

pdf(file='gmst_models.pdf', width=12,height=9)
plot(ft.gamm1$gam, page=1, main=paste(fm.gamm1[[1]][3]))
plot(ft.theo3,     page=1, main=paste(ft.theo3$formula[3]))
plot(ft.theo4bx,   page=1, main=paste(ft.theo4bx$formula[3]))

     plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
    lines(data01$time, q.gamm1,  col=2)
    lines(data01$time, qn.gamm1, col=4)
    lines(data01$time, qx.gamm1, col=3)
    grid()

 plot(data01$time,       data01$x, pch=20, cex=.3, main=paste(ft.theo4b$formula[3]))
lines(data01$time,        q.theo3, col=c21[3],lwd=2, lty=1)
lines(data01$time,       q.theo4a, col=c21[5],lwd=2, lty=1)


 plot(data01$time,       data01$x,      pch=20, cex=.3, main=paste(ft.theo4b$formula[3]), ylim=range(c(q.theo4b.fixg0, q.theo4b.fixy)), ty='n')
lines(data01$time,      q.theo4b.fixy,  col=c21[6],lwd=2, lty=1)
lines(data01$time,      q.theo4a.fixg0, col=c21[8],lwd=2, lty=1)
grid()
legend("topleft", legend=c("fix fYear to 2000","fix gmst to first value"), col=c21[c(6,8)], lwd=2)

 plot(data01$time,       data01$x,      pch=20, cex=.3, main=paste('DOY fixed:',ft.theo4b$formula[3]), ylim=range(c(qx.theo4b.fixg0, qx.theo4b.fixy)), ty='n')
lines(data01$time[ix],      qx.theo4b.fixy,  col=c21[6],lwd=2, lty=1)
lines(data01$time[ix],      qx.theo4a.fixg0, col=c21[8],lwd=2, lty=1)
grid()
legend("topleft", legend=c("fix fYear to 2000","fix gmst to first value"), col=c21[c(6,8)], lwd=2)

plot( (qx.theo4b.fixg0 + qx.theo4b.fixy)/2,  col=c21[6],lwd=2, lty=1,ty='l', main="mean of fixg0 and fixy", ylim=c(17, 20))
grid() # the warming in this approach is about half of that with the GAMM
dev.off()

CURRENT STATUS: 2026.02.10 - cant currently get sensible model combining random effects and gmst. 

########################################################################################
OLD BELOW HERE
########################################################################################


iy <- which(data01$year == 2025)
 plot(data01$time[iy], data01$x[iy], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
  plot(data01$time[iy],  q.theo4a.g0[iy], col=c21[3],lwd=2, lty=1, ylim=range(q.theo4a.g0[iy], q.theo4a.g1[iy]))
points(data01$time[iy],  q.theo4a.g1[iy], col=c21[3],lwd=2, lty=1)
  plot(q.theo4a.g1[iy]- q.theo4a.g0[iy], col=c21[3])


iy <- which(data01$year !=2025)
ft.theo3x  <- gam(x ~ s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1), knots=list(sdoy=c(0,1)),data=data01[iy,])
ft.theo4bx <- gam(x ~ s(gmst,k=4) + s(sdoy,bs="cc") + s(sdoy,fYear,bs="sz",id=1), knots=list(sdoy=c(0,1),gmst=k.gmst),data=data01[iy,])











ft.theo1  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo1  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)
ft.theo1  <- gam(x ~ s(sdoy,bs="cc")              + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)




ft.theo2b <- gam(x ~ s(sdoy,bs="cc")              + s(fYear,     bs="sz"),                                               data=data01)
ft.theo2c <- gam(x ~ s(sdoy,bs="cc")              + s(fYear,     bs="re"),                                               data=data01)
ft.regam1 <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(fYear,     bs="re"),                                               data=data01)
ft.regam2 <- gam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")), knots=list(sdoy=c(0,1)),data=data01)


q.theo2b <- predict(ft.theo2b)
q.theo2c <- predict(ft.theo2c)
q.regam1 <- predict(ft.regam1)
q.regam2 <- predict(ft.regam2)

iy <- which(data01$year %in% 1997:1999)
 plot(data01$time[iy], data01$x[iy], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time[iy],  q.theo1[iy], col=2)
lines(data01$time[iy],  q.theo2[iy], col=3)
lines(data01$time[iy], q.theo2b[iy], col=4)
lines(data01$time[iy], q.theo2c[iy], col=5)
lines(data01$time[iy], q.regam1[iy], col=6)
lines(data01$time[iy], q.regam2[iy], col=7)
grid()


y   <- seq_along(iy2)
plot(y, data01[iy2,x],pch=20, cex=.5)
lines(y,SC[iy2],col="red")
grid()















    nd       <- data01
    nd$sdoy   <- sdoyn
    qn.regam2 <- predict(ft.regam2,  newdata=nd )
    nd$sdoy   <- sdoyx
    qx.regam2 <- predict(ft.regam2,  newdata=nd )
    q.regam2  <- predict(ft.regam2 )
     plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
    lines(data01$time, q.regam2,  col=2)
    lines(data01$time, qn.regam2, col=4)
    lines(data01$time, qx.regam2, col=3)
    grid()
readline("Stop 2")

     plot(data01$time, data01$x-q.regam2, pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))



readline("Stop 3")
# fmla.qgam  <- list(x ~ class +s(sdoy, bs='cc', k=qgam.k$doy, by=class) +s(gmst, bs=bs_cc, k=qgam.k$gmst) +ti(sdoy, gmst,  bs=c('cc','tp')), ~ s(sdoy) +class)

# ft.regam1 <- qgam(x ~ s(sdoy,bs="cc", k=12) + gmst + s(fYear,     bs="re"),                                               data=data01)






### straight Theo 2c
fm.1       <- list(x ~ s(sdoy, bs='cc', k=12) +gmst +s(fYear,     bs="re"), ~ s(sdoy))
ft.reqgam1 <- mqgam(fm.1, data=data01, qu=c(0.5,0.9))
q90.1      <- qdo(ft.reqgam1, 0.9, predict)
lines(data01$time, q90.1,  col=c21[6], lwd=1)

### running smooth of temperature as covariate
kwin <- 91
# rm_x <- runmed(data01$x, kwin) # does not work at season turnovers - goes flat 
k1 <- rep(1,kwin)/kwin
# k1   <- kernel("daniell", kwin)

rm_x <- zapsmall(convolve(data01$x, rev(k1), type="o")) # padds begining and end. need to remove
rm_x <- rm_x[ (kwin-1)/2 + 1 : (nrow(data01)) ]
# and need to fill the start and end with something sensible
for(i in 1:((kwin-1)/2)){
    rm_x[i] <- mean(data01$x[1: (i+(i-1)) ])
    # cat(1, i, (i+(i-1)), cr)
}   
for(i in nrow(data01):(nrow(data01)-((kwin-1)/2))){
    rm_x[i] <- mean(data01$x[(nrow(data01)-(nrow(data01)-i)*2): nrow(data01) ])
    # cat((nrow(data01)-(nrow(data01)-i)*2), i, nrow(data01),  cr)
}   
p1   <- 1:nrow(data01) # 1:(365*5) #1:200 # (nrow(data01)-190):(nrow(data01))
 plot(data01$time[p1], data01$x[p1], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time[p1],  q.gamm1[p1], col=2)
lines(data01$time[p1],     rm_x[p1], col=7, lwd=2)
# points(data01$time[p1],     rm_x[p1], pch=3, cex=1,col=6, lwd=1)
grid()

data01[,rm_x := rm_x]
fm.2       <- list(x ~ s(sdoy, bs='cc', k=12) +gmst +s(rm_x, bs="re"), ~ s(sdoy))
fm.2       <- list(x ~ s(sdoy, bs='cc', k=12) +gmst +s(rm_x), ~ s(sdoy))
ft.reqgam2 <- mqgam(fm.2, data=data01, qu=c(0.5,0.9))
q50.2      <- qdo(ft.reqgam2, 0.5, predict)
 plot(data01$time[p1], data01$x[p1], pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
lines(data01$time[p1],  q.gamm1[p1], col=2)
lines(data01$time,            q50.2, col=c21[12], lwd=1)

################################################################
### creat ditribution of rm_x ~ sdoy from which to simulate from
s.sdoy <- seq(0,1,length=365)
plot(data01$sdoy, data01$rm_x, pch=20, cex=.3)
ft.rm_x <- gam(rm_x ~ s(sdoy, bs='cc', k=12), data=data01)
nd.rm_x <- data.table(sdoy = s.sdoy)
nd.rm_x[,rm_x := predict(ft.rm_x, newdata=nd.rm_x)]
lines(nd.rm_x$sdoy, nd.rm_x$rm_x, col=2,       lwd=6)

## Now see about integrating out the year-specific seasonal cycle
## by simulating from the posterior of the coefficients.
# simulate from posterior
n.sims   <- 200
b_sims   <- rmvn(n.sims, coef(ft.rm_x), ft.rm_x$Vp )
# coefficients index relating to just the doy-year interaction
DoYindex <- grep("sdoy",names(coef(ft.rm_x)))
# grid of years-doy 
newd        <- data.table( expand.grid(seq(0,1,length=365), 1)) 
names(newd) <- c("sdoy","B")
# newd        <- data.table( sdoy=seq(0,1,length=365))
# names(newd) <- c("sdoy")
X           <- predict(ft.rm_x, newd, type="lpmatrix")
SC_sims     <- tcrossprod( X[,DoYindex],  b_sims[,DoYindex] )
# now put these simulations into a table with DoY samples
# irrespective of year
SC_sims <- data.table(SC_sims)
# include the "observed" DoYs from the data
SC_sims[,sdoy := newd$sdoy]
# list of 366 matrices
SC_list <- split(SC_sims[, -"sdoy"], SC_sims[["sdoy"]], drop = TRUE)
# list of 366 vactors
SC_list <- lapply(SC_list, function(x){unlist(x)})
# put in a table
SC_table <- do.call(rbind,SC_list)

pred0 <- tcrossprod( X[,-DoYindex] ,b_sims[,-DoYindex])
pred1 <- SC_table[,sample(1:ncol(SC_table),n.sims)]
# Now add n.sims of the sdoy simulations
pred <- pred0 + pred1

plot(s.sdoy,apply(pred,1,mean),type="n",ylim=range(pred),xlab="Day of Year",ylab="Seasonal Cycle")
for(i in seq_along(pred[1,])) lines(s.sdoy, pred[,i],col=rgb(0,0,1,0.1))
lines(s.sdoy,apply(pred,1,mean),lwd=2,col="red4")

plot(data01$sdoy, data01$rm_x, pch=20, cex=.3)
lines(nd.rm_x$sdoy, nd.rm_x$rm_x, col=2,       lwd=6)




p343.Vp <- X %*% t(rmvn(n.sims, coef(ft.rm_x), ft.rm_x$Vp ))
p343.Ve <- X %*% t(rmvn(n.sims, coef(ft.rm_x), ft.rm_x$Ve ))



plot(1:366,apply(pred,1,mean),type="n",ylim=range(pred),xlab="Day of Year",ylab="Seasonal Cycle")
for(i in seq_along(pred[1,])) lines(1:366, pred[,i],col=rgb(0,0,1,0.1))









readline("Stop 1")
