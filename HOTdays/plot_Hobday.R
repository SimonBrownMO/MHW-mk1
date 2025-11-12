# source("plot_Hobday.R")   

st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))

MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_events  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9_Hobday.RData"


load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_msdata01, verb=TRUE)
load(st_events,   verb=TRUE)

## select the all year hobday events
events01 <- events01.hob$events01_ally

# load(st_qgam, verb=TRUE)
# th90 <- qdo(fit.qgam, newdata=data01, type='response', qu=0.9, predict)

###  plots  #####################################################
###  wrt time 

c21      <- distinct21colours()

# all years observed values
isum <- seq_along(data01$doy)
r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    

lines(data01$time[isum], th.hob$th90[isum], col='red', lwd=1)

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
abline(h=17.4, lty=2)

readline("Continue?")


# 1994 2003 observed values
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2026)
r1   <- range(data01$x[isum] ) # range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # 
plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    

lines(data01$time[isum], th.hob$th90[isum], col='red', lwd=1)

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()
readline("Continue?")


# all year anomalies
isum <- which(data01$doy >0 )
r1   <- range(data01_an_hob[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()

readline("Continue?")


# 1994 2003 anomalies
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2026)
r1   <- range(data01_an_hob[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    points(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()

NEED to
    - PLOT classes
    - compare MS probabiliites with Hobday categories
    - show the difference a simple trend climC term makes to event extraction
        - and severity classification


# colour by class
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2024)
r1   <- range(data01_an_hob[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    c.hd <- finite(events01$obs$hob_cat[,k])
    points(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()

# colour by class  all
isum <- which(data01$doy >0 & data01$time>=1980 & data01$time<2026)
r1   <- range(data01_an_hob[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(events01$obs$ch.sev)) {
    ch.t <- events01$obs$ch.time[[k]]
    c.hd <- finite(events01$obs$hob_cat[,k])
    points(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01_an_hob[finite(events01$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()

readline("Continue?")






























# ###  wrt time  #####################################################
# # summer half year:  day 150 to day 315
# isum <- which(data01$doy >= 150 & data01$doy <= 315 & data01$time >= 1994 & data01$time <= 2004)
# r1 <- c(13,19)
# plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    

# lines(data01$time[isum], th90[isum], col='red', lwd=1)


# ch.lev.l <- qlaplace(events01$obs$ch.pk$thresh.u)
# c21      <- distinct21colours()
# for(k in 1:dim(events01$obs$ch.doy)[2]) {
#     ch.t <- events01$obs$ch.time[[k]]
#     points(ch.t, events01$obs$ch.st.l[,k][seq_along(ch.t)],col=c21[k %% 21 +1])
#      lines(ch.t, events01$obs$ch.st.l[,k][seq_along(ch.t)],col=c21[k %% 21 +1])
# }
# abline(h=ch.lev.l, lty=2, col='grey70')     
# grid()
# readline("Continue?")








#