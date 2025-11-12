# source("plot_compair_Hobday_qgam.R")   

st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))

MSref.file    <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_events_hd  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9_Hobday.RData"
st_events     <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"


load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_msdata01, verb=TRUE)

## select the all year 
load(st_events_hd,   verb=TRUE)
ev.hd <- events01.hob$events01_ally
load(st_events,   verb=TRUE)
ev.qg <- events01$events01_ally

# load(st_qgam, verb=TRUE)
# th90 <- qdo(fit.qgam, newdata=data01, type='response', qu=0.9, predict)

###  plots  #####################################################
c21      <- distinct21colours()

###  wrt time 
# all years observed values

# all year Hob anomalies
isum <- which(data01$doy >0 )
r1   <- range(data01_an_hob[finite(ev.hd$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    
for(k in seq_along(ev.hd$obs$ch.sev)) {
    ch.t <- ev.hd$obs$ch.time[[k]]
    points(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()

readline("Continue?")

# all year laplace temperatures
isum <- which(data01$doy >0 )
r1   <- c(qlaplace(0.817),max(qlaplace(data01$u)[isum], na.rm=TRUE ) )  # range(data01$x[isum] )
plot(data01$time[isum], qlaplace(data01$u)[isum], pch=20, cex=.5, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    
for(k in seq_along(ev.qg$obs$ch.sev)) {
    ch.t <- ev.qg$obs$ch.time[[k]]
    points(ch.t, qlaplace(data01$u)[finite(ev.qg$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, qlaplace(data01$u)[finite(ev.qg$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()

readline("Continue?")












# 1994 2003 anomalies
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2026)
r1   <- range(data01_an_hob[finite(ev.hd$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(ev.hd$obs$ch.sev)) {
    ch.t <- ev.hd$obs$ch.time[[k]]
    points(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c21[k %% 21 +1])
     lines(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
}
grid()

# NEED to
#     - PLOT classes
#     - compare MS probabiliites with Hobday categories
#     - show the difference a simple trend climC term makes to event extraction
#         - and severity classification


# colour by class
# isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
isum <- which(data01$doy >0 & data01$time>=2015 & data01$time<2024)
r1   <- range(data01_an_hob[finite(ev.hd$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(ev.hd$obs$ch.sev)) {
    ch.t <- ev.hd$obs$ch.time[[k]]
    c.hd <- finite(ev.hd$obs$hob_cat[,k])
    points(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()

# colour by class  all
isum <- which(data01$doy >0 & data01$time>=1980 & data01$time<2026)
r1   <- range(data01_an_hob[finite(ev.hd$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
plot(data01$time[isum], data01_an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Hobday anomalies & events RRR') )    

for(k in seq_along(ev.hd$obs$ch.sev)) {
    ch.t <- ev.hd$obs$ch.time[[k]]
    c.hd <- finite(ev.hd$obs$hob_cat[,k])
    points(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col=c.hd, pch=20, cex=2.6)
     lines(ch.t, data01_an_hob[finite(ev.hd$obs$ch.i[,k])], col='grey60',lwd=2)
}
grid()

readline("Continue?")








#