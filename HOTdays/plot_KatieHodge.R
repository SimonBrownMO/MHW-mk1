# source("plot_KatieHodge.R")   

### compare different methods at the event level. Specifically 2023

# convert -trim -density 150x150 -resize 800x600 -quality 95 KatieHodge.pdf png/KatieHodge.png


st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))
# source("../libs/fn_MakeStationary.R")
# source("../libs/fn_HotDays.R")



# Hobday
MSref.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
st_ev.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

# Hobday with time
MSref.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
st_ev.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

# Hobday with GMST
MSref.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-21.RData"
st_ev.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

SAVEPLOTS<- TRUE

st_pdf   <- 'KatieHodge.pdf'
c21      <- distinct21colours()


load(MSref.vh, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

load(st_ev.vh,   verb=TRUE,  temp_env <- new.env())
l.vh <- as.list(temp_env)

load(MSref.v1, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)   
load(st_ev.v1,   verb=TRUE,  temp_env <- new.env())
l.v1 <- as.list(temp_env)        
load(st_q2p,   verb=TRUE,  temp_env <- new.env())
q2p.v1 <- as.list(temp_env)        
load(st_p2q,   verb=TRUE,  temp_env <- new.env())
p2q.v1 <- as.list(temp_env)        

load(MSref.v2, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)   
load(st_ev.v2,   verb=TRUE,  temp_env <- new.env())
l.v2 <- as.list(temp_env)

vh.d01 <- l.vh$data01
v1.d01 <- l.v1$data01
v2.d01 <- l.v2$data01

vh.e  <- l.vh$events01$events01_ally$obs
v1.e  <- l.v1$events01$events01_ally$obs
v2.e  <- l.v2$events01$events01_ally$obs


if(SAVEPLOTS) pdf(st_pdf, width=12, height=8)

######################################################################
######################################################################
## Hobday+time

### Hobdy+time
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
isum <- which(v1.d01$doy >=180 & v1.d01$doy <=280)
r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.4, ylim=r1, main='Trend North West Shelf', xlab='Year', ylab='Temperature (°C)' )   

# lines(v1.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
# lines(v1.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
# lines(v1.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)
# points(v1.d01$time[isum], l.v1$th5090$th50[isum], col=2, pch=20, cex=.3)
for(y in unique(floor(v1.d01$time[isum]))) {
    iy <- which(floor(v1.d01$time[isum])==y)
    lines(v1.d01$time[isum][iy], l.v1$th5090$th50[isum][iy], col=2, lwd=2)
}   
grid()
legend('topright', legend=c('Observation','daily 90th percentile'), col=c('black','red'), lty=c(NA,1), pch=c(20,NA), bty='n')

## mk2
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
# isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
isum   <- which(v1.d01$doy >=180 & v1.d01$doy <=280)
maxdoy <- v1.d01$doy[isum][which.max(l.v1$th5090$th50[isum])]
ixdoy  <- which(v1.d01$doy==maxdoy)

r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.4, ylim=r1, main='Trend North West Shelf', xlab='Year', ylab='Temperature (°C)' )   
lines(v1.d01$time[ixdoy], l.v1$th5090$th50[ixdoy], col='red3', lwd=2, lty=2)
grid()
legend('topright', legend=c('Observation','Summer trend'), col=c('black','red3'), lty=c(NA,1), pch=c(20,NA), bty='n')


### Hobdy+time
# pdf(width=12, height=8)
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
isum <- which(v1.d01$doy >0 & v1.d01$time>=2022.4 & v1.d01$time<2024.4)
r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.3, ylim=r1, main=paste('Hobday with trend',do.region), xlab='Year', ylab='Temperature (°C)' )    
lines(v1.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
# lines(v1.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
lines(v1.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)

for(k in seq_along(vh.e$ch.sev)) {
    if(vh.e$ch.time[[k]][1] >= 2022.4 & vh.e$ch.time[[k]][1] < 2024.4) { 
        ch.t <- vh.e$ch.time[[k]]
        c.hd <- finite(vh.e$hob_cat[,k])
        lines(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c2[k],lwd=2,lty=1)
        # points(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c.hd, pch=3, cex=1,lwd=2)
    }
}
for(k in seq_along(v1.e$ch.sev)) {
    if(v1.e$ch.time[[k]][1] >= 2022.4 & v1.e$ch.time[[k]][1] < 2024.4) { 
        ch.t <- v1.e$ch.time[[k]]
        c.hd <- finite(v1.e$hob_cat[,k])
         lines(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c2[k],lwd=4,lty=2, lend=2)
        # points(ch.t, v1.d01$x[finite(v1.e$ch.i[,k])], col=c2[k], pch=3, cex=.5,lwd=1)
    }
}
grid()
legend('topright', legend=c('Obs','Hobday fixed 90%ile','Hobday trend 90%ile','MHW fixed','MHW trend'), 
                col=c(1,'red','red3',1,1), lty=c(NA,1,2,1,2), lwd=c(1,1,2,2,3), pch=c(20,NA,NA,NA,NA),bty='n')


### method explainers
### all data +50%ile
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
isum <- which(v1.d01$doy >0 )
r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.5, ylim=r1, main=paste(do.region,'SST'), xlab='Year', ylab='Temperature (°C)' )    
lines(v1.d01$time[isum], l.v1$th5090$th50[isum], col='red', lwd=2, lty=1)
grid()
legend('topleft', legend=c('Obs','Hobday trend 50%ile'), col=c(1,'red'), lty=c(NA,1), lwd=c(1,2), pch=c(20,NA),bty='n')

### Alll %iles for years 1980 & 2024

isum1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)
isum2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)

qs      <- c(1:9*10)/100 # c(1,1:9*10,99)/100
qs.1980 <- qdo(qgam.v1$fit.qgam, newdata=v1.d01[isum1,], type='response', qu=qs, predict)
qs.2024 <- qdo(qgam.v1$fit.qgam, newdata=v1.d01[isum2,], type='response', qu=qs, predict)
r1      <- range( c(unlist(qs.1980), unlist(qs.2024)) )

                         plot(v1.d01$time[isum1], v1.d01$x[isum1], ty='n', ylim=r1, main=paste('Annual cycle of percentiles',do.region), xlab='Year', ylab='Temperature (°C)' )    
for(i in seq_along(qs)) lines(v1.d01$time[isum1], qs.1980[[i]], col=1, lwd=2, lty=1)
tim2 <- v1.d01$time[isum2] - (2024 - 1980)
for(i in seq_along(qs)) lines(tim2, qs.2024[[i]], col='red', lwd=2, lty=2)
grid()
legend('topleft', legend=c('1980 Percentiles 10%, 20%...90%','2024 Percentiles 10%, 20%...90%'), col=c(1,'red'), lty=c(1,2), lwd=c(1,2), pch=NA,bty='n')


### cdf & icdf fns
isum1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)
isum2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)
doy.n <-  v1.d01$doy[isum1][which.min(l.v1$th5090$th90[isum1])]
doy.x <-  v1.d01$doy[isum1][which.max(l.v1$th5090$th90[isum1])]

isum1 <- which(v1.d01$doy >0  & v1.d01$time>=1980.0 & v1.d01$time<1981.0)[seq(from=doy.n, length=13, by=14)]
isum2 <- which(v1.d01$doy >0  & v1.d01$time>=2024.0 & v1.d01$time<2025.0)[seq(from=doy.n, length=13, by=14)]

qs.1980 <- qs.2024 <- NULL
ps.0 <- 5:95/100
for(i in seq_along(isum1)) qs.1980[[i]] <- p2q.v1$qgam.p2q.fn[[isum1[i]]](ps.0)
for(i in seq_along(isum2)) qs.2024[[i]] <- p2q.v1$qgam.p2q.fn[[isum2[i]]](ps.0)
r1      <- range( c(unlist(qs.1980), unlist(qs.2024)) )

plot(qs.1980[[1]], ps.0, ty='n', xlim=r1, ylim=c(0,1), ylab='Probability', xlab='Temperature (°C)', main='Cumulative distribution fns, fortnightly, winter to summer' )
# for(i in seq_along(isum1)) qs.1980[[i]] <- p2q.v1$qgam.p2q.fn[[i]](ps.0)
for(i in seq_along(isum1)) lines(qs.1980[[i]], ps.0, col=1, lwd=2, lty=1)
for(i in seq_along(isum2)) lines(qs.2024[[i]], ps.0, col='red', lwd=2, lty=2)
grid()
legend('topleft', legend=c('1980 CDFs','2024 CDFs'), col=c(1,'red'), lty=c(1,2), lwd=c(1,2), pch=NA,bty='n')    




### compare actuals with laplace events

isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.8)
r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.8, ylim=r1, main=paste('Obs vs 90%ile',do.region), xlab='Year', ylab='Temperature (°C)' )    
# lines(v1.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
# lines(v1.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
lines(v1.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)
grid()

par(mfcol=c(1,2))
par(mar=c(13.1, 4.1, 4.1, 2.1))

# anomaly plot
v1.anom <- v1.d01$x - l.v1$th5090$th90
r1   <- c(-1.0, max(v1.anom[isum])) # range( v1.anom[isum] ) 
plot(v1.d01$time[isum], v1.anom[isum], pch=20, cex=.8, ylim=r1, main=paste('Anomaly event definition',do.region), xlab='Year', ylab='Temperature (°C)' )    
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
r1   <- c(-2.0, 7) # max(v1.l[isum])) # range( v1.anom[isum] ) 
plot(v1.d01$time[isum], v1.l[isum], pch=20, cex=.8, ylim=r1, main=paste('Probability event definition',do.region), xlab='Year', ylab='Probability (Laplace scale)' )    
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
isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
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
isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.45 & v1.d01$time<2023.73)
# r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
r1   <- c(14,18.1)
plot(v1.d01$time[isum], v1.d01$x[isum], pch=20, cex=.3, ylim=r1, main=paste('Hobday+time',do.region), xlab='Year', ylab='Temperature (°C)' )    
lines(v1.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
lines(v1.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
lines(v1.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)

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

