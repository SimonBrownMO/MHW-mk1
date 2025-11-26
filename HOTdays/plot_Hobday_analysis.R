# source("plot_Hobday_analysis.R")   

### compare different methods at the event level. Specifically 2023


st.pwd    <- system("pwd", intern=TRUE)
source(paste(st.pwd,"/../setup_all.R",sep=''))
# source("../libs/fn_MakeStationary.R")
# source("../libs/fn_HotDays.R")



# Hobday
MSref.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_ev.vh  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

# Hobday with time
MSref.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_ev.v1  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

# Hobday with GMST
MSref.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_ev.v2  <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/Events/ostia_cdr_nrt_regions_EventsTh0.9.RData"

SAVEPLIOTS <- TRUE

st_pdf   <- 'hobday_event_analysis_1.pdf'
c21      <- distinct21colours()

if(SAVEPLIOTS) pdf(st_pdf, width=11, height=8.5)

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

### Hobdy 1
set.seed(1)
c2 <- sample(c21, size=length(vh.e$ch.sev), replace=TRUE)
isum <- which(vh.d01$doy >0 & vh.d01$time>=2023.4 & vh.d01$time<2023.85)
r1   <- range( c(vh.d01$x[isum], l.v2$th5090$th90[isum]) ) # range(vh.d01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # 
plot(vh.d01$time[isum], vh.d01$x[isum], pch=20, cex=.3, ylim=r1, main=paste('Hobday',do.region), xlab='Year', ylab='Temperature (°C)' )    
lines(vh.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
lines(vh.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
# lines(vh.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)
# lines(vh.d01$time[isum], l.v2$th5090$th90[isum], col='red4', lwd=2, lty=2)
for(k in seq_along(vh.e$ch.sev)) {
    if(vh.e$ch.time[[k]][1] >= 2023.3 & vh.e$ch.time[[k]][1] < 2024) { 
        ch.t <- vh.e$ch.time[[k]]
        c.hd <- finite(vh.e$hob_cat[,k])
        lines(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c2[k],lwd=4)
        points(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c.hd, pch=3, cex=1,lwd=2)
    }
}
grid()
legend('topright', legend=c('Hobday 90%ile','Hobday 50%ile'), col=c('red','black'), lty=c(1,1), lwd=c(1,1), bty='n')

### Hobdy 1 + Laplace classification
set.seed(1)
c2 <- sample(c21, size=length(vh.e$ch.sev), replace=TRUE)
isum <- which(vh.d01$doy >0 & vh.d01$time>=2023.45 & vh.d01$time<2023.73)
r1   <- range( c(vh.d01$x[isum], l.v2$th5090$th90[isum]) ) # range(vh.d01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # 
r1   <- c(14,18.1)
plot(vh.d01$time[isum], vh.d01$x[isum], pch=20, cex=.3, ylim=r1, main=paste('Hobday',do.region), xlab='Year', ylab='Temperature (°C)' )    
lines(vh.d01$time[isum], l.vh$th5090$th90[isum], col='red', lwd=1)
lines(vh.d01$time[isum], l.vh$th5090$th50[isum], col=1, lwd=1)
# lines(vh.d01$time[isum], l.v1$th5090$th90[isum], col='red3', lwd=2, lty=2)
# lines(vh.d01$time[isum], l.v2$th5090$th90[isum], col='red4', lwd=2, lty=2)
for(k in seq_along(vh.e$ch.sev)) {
    if(vh.e$ch.time[[k]][1] >= 2023.3 & vh.e$ch.time[[k]][1] < 2024) { 
        ch.t <- vh.e$ch.time[[k]]
        c.hd <- finite(vh.e$hob_cat[,k])
        lines(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c2[k],lwd=4)
        points(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c.hd, pch=3, cex=1,lwd=2)
    }
}
grid()
legend('topright', legend=c('Hobday 90%ile','Hobday 50%ile'), col=c('red','black'), lty=c(1,1), lwd=c(1,1), bty='n')
legend('topleft',  legend=c('Hobday Cat','Probability Cat'), col=1, pch=c(3,1), lwd=2, lty=NA, cex=1.6, bty='n')

## classify based on Laplace data
vh.d01$lap_cat                             <- rep(0, length(vh.d01$an_hob))
vh.d01$lap_cat[ which( vh.d01$u > 0.90 )             ] <- 1
vh.d01$lap_cat[ which( vh.d01$u > 1 - (1/(1*365)) )  ] <- 2
vh.d01$lap_cat[ which( vh.d01$u > 1 - (1/(10*365)) ) ] <- 3
vh.d01$lap_cat[ which( vh.d01$u > 1 - (1/(50*365)) ) ] <- 4
# abline(h=qlaplace(0.9), col=1, lty=2)
# abline(h=qlaplace(1 - (1/(1*365))),  col='red3',   lty=2)
# abline(h=qlaplace(1 - (1/(10*365))), col='green3', lty=2)
# abline(h=qlaplace(1 - (1/(50*365))), col='blue',   lty=2)
# for(k in seq_along(vh.e$ch.sev)) {
#     c.hd <- finite(vh.d01$lap_cat[finite(vh.e$ch.i[,k])])
#     points(vh.d01$an_hobwt[finite(vh.e$ch.i[,k])], qlaplace(vh.d01$u[finite(vh.e$ch.i[,k])]), col=c.hd, pch=1, cex=1.2, lwd=2)
#     if(vh.e$ch.time[[k]][1] >= 2023.3 & vh.e$ch.time[[k]][1] < 2024) { 
#         c.hd <- finite(vh.e$lap_cat[,k])
#         # points(vh.d01$an_hobwt[finite(vh.e$ch.i[,k])], qlaplace(vh.d01$u[finite(vh.e$ch.i[,k])]), col=c.hd, pch=1, cex=2.8, lwd=2)
#     }
# }
for(k in seq_along(vh.e$ch.sev)) {
    if(vh.e$ch.time[[k]][1] >= 2023.3 & vh.e$ch.time[[k]][1] < 2024) { 
        ch.t <- vh.e$ch.time[[k]]
        c.hd <- finite(vh.d01$lap_cat[finite(vh.e$ch.i[,k])]) # finite(vh.e$hob_cat[,k])
        # lines(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c2[k],lwd=4)
        points(ch.t, vh.d01$x[finite(vh.e$ch.i[,k])], col=c.hd, pch=1, cex=2, lwd=2.5)
    }
}

### compare with Laplace transform
## Hobday
set.seed(1)
c2 <- sample(c21, size=length(vh.e$ch.sev), replace=TRUE)
isum <- which(vh.d01$doy >0 & vh.d01$time>=2023.4 & vh.d01$time<2023.85)
i1 <- which(vh.d01$an_hob>0)
plot(vh.d01$an_hob[i1], qlaplace(vh.d01$u[i1]), pch=20, cex=.3, xlab='Hobday Anomaly (C)', ylab='Probability (Laplace Transformed)',main=paste('Hobday',do.region) )
for(k in seq_along(vh.e$ch.sev)) {
    c.hd <- finite(vh.e$hob_cat[,k])
    points(vh.d01$an_hob[finite(vh.e$ch.i[,k])], qlaplace(vh.d01$u[finite(vh.e$ch.i[,k])]), col=c.hd, pch=20, cex=0.6)
    if(vh.e$ch.time[[k]][1] >= 2023.3 & vh.e$ch.time[[k]][1] < 2024) { 
        c.hd <- finite(vh.e$hob_cat[,k])
        points(vh.d01$an_hob[finite(vh.e$ch.i[,k])], qlaplace(vh.d01$u[finite(vh.e$ch.i[,k])]), col=c.hd, pch=20, cex=2.6)
    }
}
grid()
abline(h=qlaplace(0.9), col=1, lty=2)
abline(h=qlaplace(1 - (1/(1*365))),  col='red3',   lty=2)
abline(h=qlaplace(1 - (1/(10*365))), col='green3', lty=2)
abline(h=qlaplace(1 - (1/(50*365))), col='blue',   lty=2)
legend('bottomright', legend=c('90%ile','1 year','10 year','50 year'), col=c(1,'red3','green3','blue'), lty=2, bty='n', lwd=2, cex=1.2)


######################################################################
######################################################################
## Hobday+time

### Hobdy+time
set.seed(3)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
r1   <- range( c(v1.d01$x[isum], l.v2$th5090$th90[isum]) ) 
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


### compare with Laplace transform
set.seed(1)
c2 <- sample(c21, size=length(v1.e$ch.sev), replace=TRUE)
isum <- which(v1.d01$doy >0 & v1.d01$time>=2023.4 & v1.d01$time<2023.85)
i1 <- which(v1.d01$an_hobwt>0)
r1 <- qlaplace(c(0.9, 1 - (1/(50*365)) ))
plot(v1.d01$an_hobwt[i1], qlaplace(v1.d01$u[i1]), pch=20, cex=.3, 
        xlab='Hobday Anomaly (C)', ylab='Hobday Probability (Laplace Transformed)',main=paste('Hobday+time',do.region) ,ylim=r1)
for(k in seq_along(v1.e$ch.sev)) {
    c.hd <- finite(v1.e$hob_cat[,k])
    points(v1.d01$an_hobwt[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=20, cex=0.6)
    if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
        c.hd <- finite(v1.e$hob_cat[,k])
        points(v1.d01$an_hobwt[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=20, cex=2.6)
    }
}
grid()
legend('bottomright', legend=c('90%ile','1 year','10 year','50 year'), col=c(1,'red3','green3','blue'), lty=2, bty='n', lwd=2, cex=1.2)


## classify based on Laplace data
v1.d01$lap_cat                             <- rep(0, length(v1.d01$an_hobwt))
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
#     points(v1.d01$an_hobwt[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=1, cex=1.2, lwd=2)
#     if(v1.e$ch.time[[k]][1] >= 2023.3 & v1.e$ch.time[[k]][1] < 2024) { 
#         c.hd <- finite(v1.e$lap_cat[,k])
#         # points(v1.d01$an_hobwt[finite(v1.e$ch.i[,k])], qlaplace(v1.d01$u[finite(v1.e$ch.i[,k])]), col=c.hd, pch=1, cex=2.8, lwd=2)
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

