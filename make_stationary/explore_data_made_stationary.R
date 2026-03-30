# source("")

source("../setup_all.R")

load("/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/MSref/ostia_cdr_nrt_regions.MSref.2026-03-26.RData", verb=T)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_msdata01,verb=T)
d01 <- data01

plot(qlaplace(d01$uqgam),qlaplace(d01$u), xlim=c(1.5,5.0), ylim=c(1.5,5.0), pch=3, cex=.3)
grid()
abline(0,1)
max(qlaplace(d01$u))

plot(qlaplace(d01$uqgam),qlaplace(d01$u), xlim=c(1.5,5.0), ylim=c(1.5,5.0), pch=3, cex=.3)

ina1 <- which(is.na(d01$uqgam))
ina2 <- which(is.na(d01$u))
plot(d01$residStep1[ina1], pch=3, cex=.3) # qgam fials for high and low residuals
plot(d01$residStep1[ina2], pch=3, cex=.3) # evgam fails for low residuals (as no lower tial model)

plot(d01$residStep1[-ina2],qlaplace(d01$u)[-ina2], xlim=c(0,ceiling(max(d01$residStep1[-ina2]))), ylim=c(0,ceiling(max(qlaplace(d01$u[-ina2])))), pch=3, cex=.3)
grid()

plot(d01$time, d01$x, ylim=c(0,12), pch=20, cex=.3)
plot(d01$x, qlaplace(d01$u), ylim=c(0,12), pch=3, cex=.3)

x.l <- qlaplace(d01$u)
i1  <- which(x.l>qlaplace(1-(10/52693))) # ~ 10 points above this threshold

plot(d01$x, pch=20, cex=.3)
points(seq_along(x.l)[i1], d01$x[i1], pch=20, cex=2, col=2)

plot(d01$residStep1, pch=20, cex=.3)
points(seq_along(x.l)[i1], rnorm(length(i1))/10+d01$residStep1[i1], pch=20, cex=2, col=2)

### 2026.03.30 - getting the right number of exceedances, interesting that the top x.l have quite a range of residStep1 values.
### so everything at this point is looking OK










#