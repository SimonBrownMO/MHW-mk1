# source("plot_Hobday.R")   

st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))
# source("../libs/fn_MakeStationary.R")
# source("../libs/fn_HotDays.R")

# Hobday
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
tit1       <- "Hobday"

# Hobday with time
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ time"

# Hobday with GMST
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# tit1       <- "Hobday ~ GMST"

SAVEPLIOTS <- TRUE

# > str1(events01$obs)
# List of 14
#  $ ch.pk       :List of 6
#  $ ch.st.o     : num [1:186, 1:120] -0.0824 0.0058 0.0718 0.3195 0.7945 ...
#  $ ch.age      : num [1:186, 1:120] 0 1 2 3 4 5 6 7 8 9 ...
#  $ ch.i        : int [1:186, 1:120] 3772 3773 3774 3775 3776 3777 3778 3779 3780 3781 ...
#  $ ch.time     :List of 120
#  $ ch.sev      : num [1:120] 19.7 11.6 15 20.5 19 ...
#  $ ch.mean.sev : num [1:120] 0.705 0.895 0.683 0.66 0.949 ...
#  $ ch.pkval    : num [1:120] 1.61 1.49 1.44 1.35 1.33 ...
#  $ ch.pkval.day: num [1:120] 8 7 9 13 12 14 11 29 6 92 ...
#  $ ch.duration : num [1:120] 28 13 22 31 20 19 23 45 9 116 ...
#  $ th.x1       : num 0
#  $ hotseas     :List of 3
#  $ time        : num [1:16723] 1980 1980 1980 1980 1980 ...
#  $ hob_cat     : num [1:186, 1:120] 0 1 1 1 1 2 2 2 2 2 ...





st_pdf     <- sub('.RData','.pdf', st_events)

load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_msdata01, verb=TRUE)
load(st_events,   verb=TRUE)

## select the all year hobday events
events01 <- events01$events01_ally

# load(st_qgam, verb=TRUE)
# th90 <- qdo(fit.qgam, newdata=data01, type='response', qu=0.9, predict)

###  plots  #####################################################
###  wrt time 
c21      <- distinct21colours()

if(SAVEPLIOTS) pdf(st_pdf, width=11, height=8.5)

# all years observed values
isum <- seq_along(data01$doy)
# r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
r1   <- range(data01$x, na.rm=TRUE )   # range(data01$x[isum] )
 plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Temperature (°C)')   
lines(data01$time[isum], data01$x[isum], col=1, lwd=1)   
grid()
abline(h=17.4, lty=2)

# with fitted th90
plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Temperature (°C)')   
lines(data01$time[isum], th5090$th90[isum], col='red', lwd=1)
grid()
abline(h=17.4, lty=2)

# with fitted th90 and events
plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Temperature (°C)')   
lines(data01$time[isum], th5090$th90[isum], col='red', lwd=1)
for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
abline(h=17.4, lty=2)

if(!SAVEPLIOTS) readline("Continue?")


# 1994 2003 observed values
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2026)
r1   <- range(data01$x[isum] ) # range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # 
plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Temperature (°C)' )    

lines(data01$time[isum], th5090$th90[isum], col='red', lwd=1)

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
if(!SAVEPLIOTS) readline("Continue?")


# all year anomalies
isum <- which(data01$doy >0 )
r1   <- range(data01$anomaly[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01$anomaly[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Anomalies (°C)' )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
if(!SAVEPLIOTS) readline("Continue?")

# 1994 2003 anomalies
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2026)
r1   <- range(data01$anomaly[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01$anomaly[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,do.region), xlab='Year', ylab='Anomalies (°C)' )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
if(!SAVEPLIOTS) readline("Continue?")

# NEED to
#     - PLOT classes
#     - compare MS probabiliites with Hobday categories
#     - show the difference a simple trend climC term makes to event extraction
#         - and severity classification

# colour by class
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2024)
r1   <- range(data01$anomaly[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01$anomaly[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,'Class',do.region), xlab='Year', ylab='Anomalies (°C)' )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    c.hd <- finite(events01$obs$hob_cat[,k])
    points(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()
if(!SAVEPLIOTS) readline("Continue?")

# colour by class  all
isum <- which(data01$doy >0 & data01$time>=1980 & data01$time<2026)
r1   <- range(data01$anomaly[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01$anomaly[isum], pch=20, cex=.3, ylim=r1, main=paste(tit1,'Class',do.region), xlab='Year', ylab='Anomalies (°C)'  )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    c.hd <- finite(events01$obs$hob_cat[,k])
    points(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01$anomaly[finite(events01$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()
if(!SAVEPLIOTS) readline("Continue?")

# probabilities on laplace, colour by class  all
isum <- which(data01$doy >0 & data01$time>=1980 & data01$time<2026)
r1   <- range(qlaplace(data01$u[finite(events01$obs$ch.i)]), na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], qlaplace(data01$u[isum]), pch=20, cex=.3, ylim=r1, main=paste(tit1,', Hobday classes:',do.region), xlab='Year', ylab='Laplace'  )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    c.hd <- finite(events01$obs$hob_cat[,k])
    # if(c.hd[1]==1) cex1 <- 1.0 else if(c.hd[1]==2) cex1 <- 2.2 else if(c.hd[1]==2) cex1 <- 2.8 else cex1 <- 0.3
    cex1 <- c.hd * 1.2
    points(ch.t, qlaplace(data01$u[finite(events01$obs$ch.i[,k])]), col=c.hd, pch=20, cex=cex1)
    #  lines(ch.t, qlaplace(data01$u[finite(events01$obs$ch.i[,k])]), col='grey60',lwd=2)
}
grid()

if(!SAVEPLIOTS) readline("Continue?")

if(SAVEPLIOTS) dev.off()





#