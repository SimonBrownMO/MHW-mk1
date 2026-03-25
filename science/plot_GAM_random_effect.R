# source("plot_GAM_random_effect.R")

SAVEPLOT <- TRUE
load(file="./GAM_random_effect.RData", verbose=TRUE)
# ft1.regamH1d <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype), 
#                         data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor )
# min_sp_floor[1:10] 0 0 0.5 0.5 0 0 1 1 1 1
# s(sdoy):ftypeobs        14.286     22 11.733  <2e-16 ***
# s(sdoy):ftypemod        21.064     22 17.508  <2e-16 ***
# s(gmst)                  2.503      9  0.554  <2e-16 ***
# ti(sdoy,gmst)           12.000     12  0.062  <2e-16 ***
# s(sdoy,fYear):ftypeobs 185.320    459 19.157  <2e-16 ***
# s(sdoy,fYear):ftypemod 401.094    995 26.798  <2e-16 ***

library(gratia)
pdf(file="./GAM_random_effect.pdf", width=12, height=9)

# model fit
plot(ft1.regamH1d, page=1)

appraise(ft1.regamH1d)

sm.names <-smooths(ft1.regamH1d)
up.2by2()
for (smi in sm.names[1:4]) {
  cat("Plotting",smi,"\n")
  p <- draw(ft1.regamH1d, select=smi, main=smi)
  print(p)
  if(!SAVEPLOT) readline("continue?")
}


basis.H1d <- lapply(sm.names, function(smi){
  basis(ft1.regamH1d, sm=smi)
})
names(basis.H1d) <- sm.names

basis.H1d <- list()
for (smi in sm.names) {
    basis.H1d[[smi]] <- basis(ft1.regamH1d, sm=smi)
}


basis.H1d <- basis(ft1.regamH1d)
draw(basis.H1d)

# prediction
    q1     <- predict(ft1.regamH1d) 
    iy2001 <- which(data01$year==2001 & data01$isobs==1)
    i0x    <- which.max(q1[ iy2001] )
    i0n    <- which.min(q1[ iy2001] )
    doyn   <- data01$doy[ iy2001][i0n]
    sdoyn  <- data01$sdoy[iy2001][i0n]
    doyx   <- data01$doy[ iy2001][i0x]
    sdoyx  <- data01$sdoy[iy2001][i0x]
    nd       <- data01
    nd$sdoy  <- sdoyn
    q1.n     <- predict(ft1.regamH1d,  newdata=nd )
    nd$sdoy  <- sdoyx
    q1.x     <- predict(ft1.regamH1d,  newdata=nd )

    # remove interannual variability
    nd       <- data01
    nd$fYear <- "2003"
    q1.2003  <- predict(ft1.regamH1d, newdata=nd) 

# all data + fit  
pdf(file="./GAM_random_effect_regamH1d.pdf", width=12, height=9)
up.1()
plot(data01$x, pch=20, cex=.3, main=paste("NWS ~ ",ft1.regamH1d$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
lines(q1,   col=2, lwd=2)
lines(q1.n, col=4, lwd=2)
lines(q1.x, col=3, lwd=2)
ix <- which(data01$doy==doyx)
# lines(ix, q1.x[ix], col=6, lwd=2)
lines(ix, q1.2003[ix], col=6, lwd=2)
legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Summer max", "Summer max @2003 IAV"), col=c(1,2,4,3,6), lwd=c(NA,2,2,2,2), pch=c(20,NA,NA,NA,20), bty="n", cex=1.2)
grid()
# dev.off()
#   if(!SAVEPLOT) readline("continue?")

up.1()
ix <- 16000:18000-185
plot(data01$x[ix], pch=20, cex=.3, main=paste("NWS ~ ",ft1.regamH1d$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
lines(q1[ix],   col=2, lwd=2)
lines(q1.n[ix], col=4, lwd=2)
lines(q1.x[ix], col=3, lwd=2)
iy <- which(data01$doy[ix]==doyx)
lines(iy, q1.2003[ix][iy], col=6, lwd=2); points(iy, q1.2003[ix][iy], col=6, pch=20, cex=2.0)
legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Summer max", "Summer max @2003 IAV"), col=c(1,2,4,3,6), lwd=c(NA,2,2,2,2), pch=c(20,NA,NA,NA,20), bty="n", cex=1.2)
grid()

dev.off()
  if(!SAVEPLOT) readline("continue?")

### focus on 1997-1998
# for obs
iy2 <- which(data01$year %in% 1987:1999 & data01$isobs==1)
y   <- seq_along(iy2)
plot(y, data01[iy2,'x'],pch=20, cex=.5, main="Obs 1987:1999")
lines(y,q1[iy2],col="red")
grid()
  if(!SAVEPLOT) readline("continue?")
# for model
iy2 <- which(data01$year %in% 1987:1999 & data01$isobs==0)
y   <- seq_along(iy2)
plot(y, data01[iy2,'x'],pch=20, cex=.5, main="Model 1987:1999")
lines(y,q1[iy2],col="red")
grid()
  if(!SAVEPLOT) readline("continue?")

# summer and winter limits
i1 <- which(data01$doy==doyn)
i2 <- which(data01$doy==doyx)
 plot(data01$time[c(i1,i2)], data01$x[c(i1,i2)], pch=20, cex=1.3, main=ft1.regamH1d$formula, xlab="Time", ylab="SST")
lines(data01$time[i1],       q1.n[i1], lwd=2, col=4)
lines(data01$time[i2],       q1.x[i2], lwd=2, col=3)  
legend("topleft", legend=c("Winter min", "Summer max"), col=c(4,3), lwd=c(2,2), pch=NA, bty="n")          
  if(!SAVEPLOT) readline("continue?")

# residuals
plot(data01$x-q1, pch=20, cex=.3, main="Residuals", xlab="index", ylab="Residual")
  if(!SAVEPLOT) readline("continue?")




  if(SAVEPLOT) dev.off()



#