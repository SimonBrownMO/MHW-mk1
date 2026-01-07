# source("plot_data_explore.R")   

### make exploratory plots
### derived from plot_KatieHodge.R

# convert -trim -density 150x150 -resize 800x600 -quality 95 KatieHodge.pdf png/KatieHodge.png


st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))
# source("../libs/fn_MakeStationary.R")
# source("../libs/fn_HotDays.R")

SAVEPLOTS <- TRUE
st_pdf    <- 'data_explore_01.pdf'
c21       <- distinct21colours()
c2 <- c21


# # Hobday
# MSref.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
# st_ev.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# load(MSref.vh, verb=TRUE)
# list2env(MSconfig , envir = .GlobalEnv)
# list2env(MSconfig$files , envir = .GlobalEnv)
# load(st_ev.vh,   verb=TRUE,  temp_env <- new.env())
# l.vh <- as.list(temp_env)
# vh.d01 <- l.vh$data01
# vh.e  <- l.vh$events01$events01_ally$obs

# Hobday with time
MSref.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-26.RData"
load(MSref.v1, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)   
st_ev.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
load(st_ev.v1,   verb=TRUE,  temp_env <- new.env())
l.v1 <- as.list(temp_env)        
# v1.d01 <- l.v1$data01
v1.e   <- l.v1$events01$events01_ally$obs
load(st_q2p,   verb=TRUE,  temp_env <- new.env())
q2p.v1 <- as.list(temp_env)        
load(st_p2q,   verb=TRUE,  temp_env <- new.env())
p2q.v1 <- as.list(temp_env)        
load(st_msdata01,   verb=TRUE,  temp_env <- new.env())
data01.v1 <- as.list(temp_env)        
v1.d01 <- data01.v1$data01

# # Hobday with GMST
# MSref.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
# st_ev.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"
# load(MSref.v2, verb=TRUE)
# list2env(MSconfig , envir = .GlobalEnv)
# list2env(MSconfig$files , envir = .GlobalEnv)   
# load(st_ev.v2,   verb=TRUE,  temp_env <- new.env())
# l.v2   <- as.list(temp_env)
# v2.d01 <- l.v2$data01
# v2.e   <- l.v2$events01$events01_ally$obs


if(SAVEPLOTS) pdf(st_pdf, width=10, height=7)

iplt <- which(v1.d01$doy >0 ) #& v1.d01$time>=2022.4 & v1.d01$time<2024.4)

### Plot 1: all data with summer/winter peak trend line
r1   <- range( c(v1.d01$x[iplt], l.v1$th5090$th90[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.4, ylim=r1, main='North West Shelf', xlab='Year', ylab='Temperature (°C)' )   
maxdoy <- v1.d01$doy[iplt][which.max(l.v1$th5090$th50[iplt])]
ixdoy  <- which(v1.d01$doy==maxdoy)
lines(v1.d01$time[ixdoy], l.v1$th5090$th50[ixdoy], col='red3', lwd=2, lty=2)
mindoy <- v1.d01$doy[iplt][which.min(l.v1$th5090$th50[iplt])]
indoy  <- which(v1.d01$doy==mindoy)
lines(v1.d01$time[indoy], l.v1$th5090$th50[indoy], col='blue3', lwd=2, lty=2)
grid()
legend('bottomright', legend=c('Observation','Peak summer trend','Peak winter trend'), col=c('black','red3','blue3'), 
                lwd=2, lty=c(NA,1,1), pch=c(20,NA,NA), bty='n', cex=1.0)

### Plot 2: all data with NS 50%ile
r1   <- range( c(v1.d01$x[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.3, ylim=r1, main='North West Shelf', xlab='Year', ylab='Temperature (°C)' )   
lines(v1.d01$time[iplt], l.v1$th5090$th50[iplt], col='red3', lwd=1, lty=1)
lines(v1.d01$time[ixdoy], l.v1$th5090$th50[ixdoy], col='red3', lwd=2, lty=2)

grid()
legend('bottomright', legend=c('Observation','median'), col=c('black','red3'), 
                lwd=2, lty=c(NA,1), pch=c(20,NA), bty='n', cex=1.5)


### Plot 3: all data minus NS 50%ile
x1 <- v1.d01$x[iplt]-l.v1$th5090$th50[iplt]
r1 <- range( x1 ) 
plot(v1.d01$time[iplt], x1, pch=20, cex=.3, ylim=r1, main='North West Shelf', xlab='Year', ylab='Anomaly (°C)' )   
# lines(v1.d01$time[iplt], , col='red3', lwd=1, lty=1)
lines(smooth.spline(x1 ~ v1.d01$time[iplt], df=20), col='blue3', lwd=2)
lines(smooth.spline(x1 ~ v1.d01$time[iplt], df=70), col='red2', lwd=2)
grid()
abline(h=0)
legend('bottomright', legend=c('Observation','splines'), col=1, lwd=2, lty=c(NA,1), pch=c(20,NA), bty='n', cex=1.5)


### Plot 4: all data actuals and anomalys
x1 <- v1.d01$x[iplt]-l.v1$th5090$th50[iplt]
x2 <- 2.0*(l.v1$th5090$th50[iplt]-mean(l.v1$th5090$th50[iplt],na.rm=TRUE))/diff(range(l.v1$th5090$th50[iplt],na.rm=TRUE))
r1 <- range( x1 ) 
plot(v1.d01$time[iplt], x1, pch=20, cex=.3, ylim=r1, main='North West Shelf', xlab='Year', ylab='Anomaly (°C)' )   
lines(v1.d01$time[iplt],x2, col='red3', lwd=1, lty=1)
grid()
abline(h=0)

legend('bottomright', legend=c('Observation','median'), col=c('black','red3'), 
                lwd=2, lty=c(NA,1), pch=c(20,NA), bty='n', cex=1.5)


plot(l.v1$th5090$th50[iplt], x1, pch=20, cex=.3, ylim=r1, main='North West Shelf', xlab='NS annual cycle', ylab='Anomaly (°C)' )   
grid()
abline(h=0)

if(SAVEPLOTS) dev.off()





readline("Stop now 2")
######################################################################
######################################################################
## Hobday+time

### Hobdy+time
r1   <- range( c(v1.d01$x[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.4, ylim=r1, main='Trend North West Shelf', xlab='Year', ylab='Temperature (°C)' )   
lines(v1.d01$time[iplt], l.v1$th5090$th90[iplt], col='red3', lwd=1, lty=1)

# lines(v1.d01$time[iplt], l.vh$th5090$th90[iplt], col='red', lwd=1)
# lines(v1.d01$time[iplt], l.vh$th5090$th50[iplt], col=1, lwd=1)
# points(v1.d01$time[iplt], l.v1$th5090$th50[iplt], col=2, pch=20, cex=.3)
for(y in unique(floor(v1.d01$time[iplt]))) {
    iy <- which(floor(v1.d01$time[iplt])==y)
    lines(v1.d01$time[iplt][iy], l.v1$th5090$th50[iplt][iy], col=2, lwd=2)
}   
grid()
legend('topright', legend=c('Observation','daily 90th percentile'), col=c('black','red'), lty=c(NA,1), pch=c(20,NA), bty='n')

readline("Stop now")

### Hobdy+time
# pdf(width=12, height=8)
r1   <- range( c(v1.d01$x[iplt], l.v1$th5090$th90[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.3, ylim=r1, main=do.region, xlab='Year', ylab='Temperature (°C)' )    
lines(v1.d01$time[iplt], l.v1$th5090$th90[iplt], col='red3', lwd=2, lty=2)

grid()
legend('topright', legend=c('Obs','Hobday trend 90%ile','MHW fixed','MHW trend'), 
                col=c(1,'red','red3',1,1), lty=c(NA,1,2,1,2), lwd=c(1,1,2,2,3), pch=c(20,NA,NA,NA,NA),bty='n')


### method explainers
### all data +50%ile
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
iplt <- which(v1.d01$doy >0 )
r1   <- range( c(v1.d01$x[iplt], l.v1$th5090$th90[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.5, ylim=r1, main=paste(do.region,'SST'), xlab='Year', ylab='Temperature (°C)' )    
lines(v1.d01$time[iplt], l.v1$th5090$th50[iplt], col='red', lwd=2, lty=1)
grid()
legend('topleft', legend=c('Obs','Hobday trend 50%ile'), col=c(1,'red'), lty=c(NA,1), lwd=c(1,2), pch=c(20,NA),bty='n')

### Alll %iles for years 1980 & 2024

iplt1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)
iplt2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)

qs      <- c(1:9*10)/100 # c(1,1:9*10,99)/100
qs.1980 <- qdo(qgam.v1$fit.qgam, newdata=v1.d01[iplt1,], type='response', qu=qs, predict)
qs.2024 <- qdo(qgam.v1$fit.qgam, newdata=v1.d01[iplt2,], type='response', qu=qs, predict)
r1      <- range( c(unlist(qs.1980), unlist(qs.2024)) )

                         plot(v1.d01$time[iplt1], v1.d01$x[iplt1], ty='n', ylim=r1, main=paste('Annual cycle of percentiles',do.region), xlab='Year', ylab='Temperature (°C)' )    
for(i in seq_along(qs)) lines(v1.d01$time[iplt1], qs.1980[[i]], col=1, lwd=2, lty=1)
tim2 <- v1.d01$time[iplt2] - (2024 - 1980)
for(i in seq_along(qs)) lines(tim2, qs.2024[[i]], col='red', lwd=2, lty=2)
grid()
legend('topleft', legend=c('1980 Percentiles 10%, 20%...90%','2024 Percentiles 10%, 20%...90%'), col=c(1,'red'), lty=c(1,2), lwd=c(1,2), pch=NA,bty='n')


### cdf & icdf fns
iplt1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)
iplt2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)
doy.n <-  v1.d01$doy[iplt1][which.min(l.v1$th5090$th90[iplt1])]
doy.x <-  v1.d01$doy[iplt1][which.max(l.v1$th5090$th90[iplt1])]

iplt1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)[seq(from=doy.n, length=13, by=14)]
iplt2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)[seq(from=doy.n, length=13, by=14)]

qs.1980 <- qs.2024 <- NULL
ps.0 <- 5:95/100
for(i in seq_along(iplt1)) qs.1980[[i]] <- p2q.v1$qgam.p2q.fn[[iplt1[i]]](ps.0)
for(i in seq_along(iplt2)) qs.2024[[i]] <- p2q.v1$qgam.p2q.fn[[iplt2[i]]](ps.0)
r1      <- range( c(unlist(qs.1980), unlist(qs.2024)) )

plot(qs.1980[[1]], ps.0, ty='n', xlim=r1, ylim=c(0,1), ylab='Probability', xlab='Temperature (°C)', main='Cumulative distribution fns, fortnightly, winter to summer' )
# for(i in seq_along(iplt1)) qs.1980[[i]] <- p2q.v1$qgam.p2q.fn[[i]](ps.0)
for(i in seq_along(iplt1)) lines(qs.1980[[i]], ps.0, col=1, lwd=2, lty=1)
for(i in seq_along(iplt2)) lines(qs.2024[[i]], ps.0, col='red', lwd=2, lty=2)
grid()
legend('topleft', legend=c('1980 CDFs','2024 CDFs'), col=c(1,'red'), lty=c(1,2), lwd=c(1,2), pch=NA,bty='n')    




### compare actuals with laplace events

iplt <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.8)
r1   <- range( c(v1.d01$x[iplt], l.v1$th5090$th90[iplt]) ) 
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.8, ylim=r1, main=paste('Obs vs 90%ile',do.region), xlab='Year', ylab='Temperature (°C)' )    
# lines(v1.d01$time[iplt], l.vh$th5090$th90[iplt], col='red', lwd=1)
# lines(v1.d01$time[iplt], l.vh$th5090$th50[iplt], col=1, lwd=1)
lines(v1.d01$time[iplt], l.v1$th5090$th90[iplt], col='red3', lwd=2, lty=2)
grid()

par(mfcol=c(1,2))
par(mar=c(13.1, 4.1, 4.1, 2.1))

# anomaly plot
v1.anom <- v1.d01$x - l.v1$th5090$th90
r1   <- c(-1.0, max(v1.anom[iplt])) # range( v1.anom[iplt] ) 
plot(v1.d01$time[iplt], v1.anom[iplt], pch=20, cex=.8, ylim=r1, main=paste('Anomaly event definition',do.region), xlab='Year', ylab='Temperature (°C)' )    
grid()
abline(h=0, col=1, lwd=2, lty=2)
for(k in seq_along(v1.e$ch.sev)) {
    if(v1.e$ch.time[[k]][1] >= 2022.4 & v1.e$ch.time[[k]][1] < 2024.4) { 
        ch.t <- v1.e$ch.time[[k]]
         lines(ch.t, v1.anom[finite(v1.e$ch.i[,k])], col=c2[k],lwd=4,lty=1, lend=2)
    }
}

# laplace plot
v1.l <- qlaplace(v1.d01$u)
r1   <- c(-2.0, 7) # max(v1.l[iplt])) # range( v1.anom[iplt] ) 
plot(v1.d01$time[iplt], v1.l[iplt], pch=20, cex=.8, ylim=r1, main=paste('Probability event definition',do.region), xlab='Year', ylab='Probability (Laplace scale)' )    
grid()
abline(h=qlaplace(0.9), col=1, lwd=2, lty=2)
for(k in seq_along(v1.e$ch.sev)) {
    if(v1.e$ch.time[[k]][1] >= 2022.4 & v1.e$ch.time[[k]][1] < 2024.4) { 
        ch.t <- v1.e$ch.time[[k]]
         lines(ch.t, v1.l[finite(v1.e$ch.i[,k])], col=c2[k],lwd=4,lty=1, lend=2)
    }
}

abline(h=qlaplace(0.9), col=1, lty=2, lwd=2)
abline(h=qlaplace(1 - (1/(1*365))),  col='red3',  lty=2, lwd=2)
abline(h=qlaplace(1 - (1/(5*365))), col='green3', lty=2, lwd=2)
legend('bottomright', legend=c('90%ile','1 year','5 year'), col=c(1,'red3','green3'), lty=2, bty='n', lwd=2, cex=1.2)



par(mfcol=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))





if(SAVEPLOTS) dev.off()



readline("Stop now")

#############################################################################################
#############################################################################################
#############################################################################################






### compare with Laplace transform
set.seed(1)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
iplt <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
i1 <- which(v1.d01$anomaly>0)
r1 <- qlaplace(c(0.9, 1 - (1/(50*365)) ))
plot(v1.d01$anomaly[i1], qlaplace(v1.d01$u[i1]), pch=20, cex=.3, 
        xlab='Hobday Anomaly (C)', ylab='Hobday Probability (Laplace Transformed)',main=paste('Hobday+time',do.region) ,ylim=r1)
for(k in seq_along(v1.e$ch.sev)) {
    c.hd <- finite(v1.e$hob_cat[,k])
    points(v1.d01$anomaly[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=20, cex=0.6)
    if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
        c.hd <- finite(v1.e$hob_cat[,k])
        points(v1.d01$anomaly[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=20, cex=2.6)
    }
}
grid()
legend('bottomright', legend=c('90%ile','1 year','10 year','50 year'), col=c(1,'red3','green3','blue'), lty=2, bty='n', lwd=2, cex=1.2)


## classify based on Laplace data
v1.d01$lap_cat                             <- rep(0, length(v1.d01$anomaly))
v1.d01$lap_cat[ which( v1.d01$u > 0.90 )             ] <- 1
v1.d01$lap_cat[ which( v1.d01$u > 1 - (1/(1*365)) )  ] <- 2
v1.d01$lap_cat[ which( v1.d01$u > 1 - (1/(10*365)) ) ] <- 3
v1.d01$lap_cat[ which( v1.d01$u > 1 - (1/(50*365)) ) ] <- 4
abline(h=qlaplace(0.9), col=1, lty=2)
abline(h=qlaplace(1 - (1/(1*365))),  col='red3',   lty=2)
abline(h=qlaplace(1 - (1/(10*365))), col='green3', lty=2)
abline(h=qlaplace(1 - (1/(50*365))), col='blue',   lty=2)
# for(k in seq_along(v1.e$ch.sev)) {
#     # c.hd <- finite(v1.e$lap_cat[,k])
#     c.hd <- finite(v1.d01$lap_cat[finite(v1.e$ch.i[,k])])
#     points(v1.d01$anomaly[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=1, cex=1.2, lwd=2)
#     if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
#         c.hd <- finite(v1.e$lap_cat[,k])
#         # points(v1.d01$anomaly[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=1, cex=2.8, lwd=2)
#     }
# }


### Hobdy+time
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
iplt <- which(v1.d01$doy >0 & v1.d01$time>=2023.45 & v1.d01$time<2023.73)
# r1   <- range( c(v1.d01$x[iplt], l.v1$th5090$th90[iplt]) ) 
r1   <- c(14,18.1)
plot(v1.d01$time[iplt], v1.d01$x[iplt], pch=20, cex=.3, ylim=r1, main=paste('Hobday+time',do.region), xlab='Year', ylab='Temperature (°C)' )    
# lines(v1.d01$time[iplt], l.vh$th5090$th90[iplt], col='red', lwd=1)
# lines(v1.d01$time[iplt], l.vh$th5090$th50[iplt], col=1, lwd=1)
lines(v1.d01$time[iplt], l.v1$th5090$th90[iplt], col='red3', lwd=2, lty=2)

for(k in seq_along(v1.e$ch.sev)) {
    if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
        ch.t <- v1.e$ch.time[[k]]
        c.hd <- finite(v1.e$hob_cat[,k])
        lines(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c2[k],lwd=4)
        points(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c.hd, pch=3, cex=1,lwd=2)
    }
}
grid()
legend('topright', legend=c('Hobday 90%ile','Hobday 50%ile','Hobday+Time 90%ile'), col=c('red','black','red3'), lty=c(1,1,2), lwd=c(1,1,2), bty='n')
for(k in seq_along(v1.e$ch.sev)) {
    if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
        ch.t <- v1.e$ch.time[[k]]
        c.hd <- finite(v1.d01$lap_cat[finite(v1.e$ch.i[,k])]) # finite(v1.e$hob_cat[,k])
        # lines(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c2[k],lwd=4)
        points(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c.hd, pch=1, cex=2, lwd=2.5)
    }
}
legend('topleft',  legend=c('Hobday Cat','Probability Cat'), col=1, pch=c(3,1), lwd=2, lty=NA, cex=1.6, bty='n')

dev.off()

