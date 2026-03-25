# source("GAM_random_effect.R")

library(mgcv)
library(data.table)
library(PCICt)

# load the data
load("data01_NWS.RData", verb=TRUE)
data01 <- data01[ ,c("time", "x", "gmst", "sdoy", "doy") ]
data01 <- data.table(data01)
data01[,isobs := 1]

### add in future modelled data
iregion <- 25  # which region and mask from digest_mass_files.R
load("../DATA/cpm_gmst.RData", verb=TRUE)  # cpm_gmst$m001$gmst
st.model <- "../DATA/mass_dump/../1980-2100-Region_mean_timeseries_Daily_r001i1p00000.RData"
load(st.model, verb=TRUE)
m.sdoy <- m.sst$doy / 360 # can do this as model has 360 day calendar
# CPM data only goes to 2080-11-30 12:00:00 so need to crop m.sst to match
i1       <- which(m.sst$date         %in% cpm_gmst$m001$time) 
i2       <- which(cpm_gmst$m001$time %in% m.sst$date[i1])
data01.m <- data.frame(time=m.sst$time[i1], x=m.sst$sst[iregion,i1], gmst=cpm_gmst$m001$gmst[i2], sdoy=m.sdoy[i1], doy=m.sst$doy[i1] )
data01.m <- data.table(data01.m)
data01.m[,isobs := 0]

data01 <- rbind(data01, data01.m)

data01[,year := trunc(time,0)]
data01[,fYear := factor(year)]  
data01[,ftype := factor(isobs, levels=c(1,0), labels=c("obs","mod"))]

### global temperatures
gmst_ref_period <- 1981:2000  # anomalies normalised to this as per UKCP18
iob  <- which(data01$isobs==1)
imo  <- which(data01$isobs==0)
iob2 <- which(data01$isobs==1 & data01$year %in% gmst_ref_period)
imo2 <- which(data01$isobs==0 & data01$year %in% gmst_ref_period)
data01$gmst[iob] <- data01$gmst[iob] - mean(data01$gmst[iob2])
data01$gmst[imo] <- data01$gmst[imo] - mean(data01$gmst[imo2])

plot(data01$time, data01$gmst, ty='n', main="gmst" )
points(data01$time[iob], data01$gmst[iob], pch=20, cex=.3, col=1)
points(data01$time[imo], data01$gmst[imo], pch=20, cex=.3, col=2)

if(FALSE) {
    ### step 1: fit mixed model to get mean gmst term and baseline annual cycle
        fm.gamm1    <- list(x ~ s(sdoy, bs="cc", k=12) + s(gmst,k=4), ~ 1 )
        fm.gamm1    <- list(x ~ s(sdoy, bs="cc", k=12, by=ftype) + s(gmst,k=4), ~ 1 )
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

        plot(data01$time, data01$x, pch=20, cex=.3, main=paste('MM : ',fm.gamm1[[1]][3]))
        lines(data01$time, q.gamm1,  col=2)
        lines(data01$time, qn.gamm1, col=4)
        lines(data01$time, qx.gamm1, col=3)
        grid()
        readline("continue? 1")


        ft1.gam <- gam(x ~ ftype +s(sdoy,bs="cc", k=12, by=ftype) + s(gmst,k=4) + ti(sdoy,gmst,bs=c("cc","tp")), data=data01)
        q.gam1  <- predict(ft1.gam) 
          plot(data01$x, pch=20, cex=.3)
        points(q.gam1,   pch=20, cex=.3,  col=2)
        grid()
        readline("continue? 2")

        sm1 <- smoothCon(s(sdoy, bs = 'cc', k = 12), data = data01, knots = NULL)

        ### fix knots. cc needs 12 but 0 and 1 are the same point so need to remove one of them
        nk <- 13
        do.knots <- list(sdoy=seq(0,1,length=nk)[-nk], gmst=seq(min(data01$gmst), max(data01$gmst), length=4))
        ft1.regam <- gam(x ~ ftype +s(sdoy,bs="cc", k=12, by=ftype) +s(gmst,k=4) +ti(sdoy,gmst,bs=c("cc","tp"),k=c(12,4)) +s(sdoy,fYear,bs="sz",id=1), knots=do.knots,data=data01)
        # poor fit

    # ft1.regamA <- gam(x ~ ftype +s(sdoy,bs="cc",by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), method='REML', select=TRUE, data=data01)
    # # vv slow

    # try gamma to penalise more and get smoother fit
    ft1.regamB <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst,k=4) +ti(sdoy,gmst,bs=c("cc","tp")) 
                                +s(sdoy,fYear,bs="sz",id=1), method="GCV.Cp", gamma=1.4, data=data01)

    # try sp in random effect smooth
    ft1.regamC <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst,k=4) +ti(sdoy,gmst,bs=c("cc","tp")) 
                                +s(sdoy,fYear,bs="sz",id=1, sp=3), method="GCV.Cp", data=data01)

    # try select=TRUE
    ft1.regamD <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst,k=4) +ti(sdoy,gmst,bs=c("cc","tp")) 
                                +s(sdoy,fYear,bs="sz",id=1), method="GCV.Cp", select=TRUE, data=data01)

    # try min.sp to force more smoothing in the random effect
    ft1.regamE <- gam(x ~ ftype +s(sdoy,bs="cc",by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) 
                                +s(sdoy,fYear,bs="sz",id=1), method="GCV.Cp", gamma=1.5, min.sp =c(NA,NA,NA,NA,2), select=TRUE, data=data01)


    # discretization only available with fREML
    # min.sp not supported with fast REML computation
    ft1.regamF <- bam(x ~ ftype +s(sdoy,bs="cc",by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,k=10), 
                          data=data01, method="fREML", gamma=1.3, select=TRUE, discrete=TRUE)
    # better

    # Identify which smoothing parameters belong to the sz term
    ft0 <- gam(x ~ ftype +s(sdoy,bs="cc",by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", gamma=1.3, select=TRUE)
    sp_names            <- names(ft0$sp)
    is_sz               <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
    min_sp_floor        <- rep(0, length(ft0$sp))     # initialize min.sp vector with zeros (no penalty)
    min_sp_floor[is_sz] <- 2                          # choose your floor (e.g., 2)
    ft1.regamG <- gam(x ~ ftype +s(sdoy,bs="cc",by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", gamma=1.3, select=TRUE, min.sp=min_sp_floor)
    # this looks like it has potential!!!!
    # s(sdoy):ftypeobs   7.841      8  545.000  <2e-16 ***
    # s(sdoy):ftypemod   7.954      8 2003.147  <2e-16 ***
    # s(gmst)            8.640      9   61.468  <2e-16 ***
    # ti(sdoy,gmst)     11.773     12    1.243  <2e-16 ***
    # s(sdoy,fYear)    373.726    999   18.285  <2e-16 ***
    # gmst too wiggly

    # add k to sdoy_1, constrain wigglyness of gmst
    is_gmst               <- grepl("\\(gmst\\)", sp_names)  
    is_sz                 <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
    min_sp_floor          <- rep(0, length(ft0$sp))     # initialize min.sp vector with zeros (no penalty)
    min_sp_floor[is_sz]   <- 2                          # choose your floor (e.g., 2)
    min_sp_floor[is_gmst] <- 2
    ft1.regamH <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", gamma=1.3, select=TRUE, min.sp=min_sp_floor)
    # s(sdoy):ftypeobs   9.546     10  388.007  <2e-16 ***
    # s(sdoy):ftypemod   9.845     10 1523.869  <2e-16 ***
    # s(gmst)            2.689      9   25.228  <2e-16 ***
    # ti(sdoy,gmst)     11.849     12    1.059  <2e-16 ***
    # s(sdoy,fYear)    377.451   1000   18.783  <2e-16 ***

    ft1.regamH0 <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", min.sp=min_sp_floor)
    # s(sdoy):ftypeobs   9.489  10.000  760.920  < 2e-16 ***
    # s(sdoy):ftypemod   9.597  10.000 2450.307  < 2e-16 ***
    # s(gmst)            1.695   2.124   11.021 0.000309 ***
    # ti(sdoy,gmst)      4.937  12.000    9.138  < 2e-16 ***
    # s(sdoy,fYear)    988.301 999.242   49.068  < 2e-16 ***
    # RE smooth not too large but too wiggly

    ft1.regamH1 <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
    # s(sdoy):ftypeobs   9.646     10  305.639  <2e-16 ***
    # s(sdoy):ftypemod   9.878     10 1389.559  <2e-16 ***
    # s(gmst)            2.689      9   25.218  <2e-16 ***
    # ti(sdoy,gmst)     11.946     12    0.962  <2e-16 ***
    # s(sdoy,fYear)    377.397   1000   18.783  <2e-16 ***

    ft1.regamH2 <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", gamma=1.3, min.sp=min_sp_floor)
    # s(sdoy):ftypeobs   9.419  10.000  899.32  < 2e-16 ***
    # s(sdoy):ftypemod   9.575  10.000 2655.41  < 2e-16 ***
    # s(gmst)            1.696   2.125   10.98 0.000326 ***
    # ti(sdoy,gmst)      5.046  12.000    9.14  < 2e-16 ***
    # s(sdoy,fYear)    984.061 998.763   48.97  < 2e-16 ***
    # RE smooth not too large but too wiggly

    ft1.regamH3 <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", select=TRUE)
    # s(sdoy):ftypeobs  10.000     10  320.77  <2e-16 ***
    # s(sdoy):ftypemod  10.000     10 1565.43  <2e-16 ***
    # s(gmst)            5.801      8  117.72  <2e-16 ***
    # ti(sdoy,gmst)     12.000     12   24.83  <2e-16 ***
    # s(sdoy,fYear)    998.986   1000   51.37  <2e-16 ***
    # poor. RE smooth is large and opposite to gmst smooth
    # THUS have to have min.sp constraint on RE smooth

    is_sz                 <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
    min_sp_floor          <- rep(0, length(ft0$sp))     # initialize min.sp vector with zeros (no penalty)
    min_sp_floor[is_sz]   <- 2                          # choose your floor (e.g., 2)
    ft1.regamH4 <- gam(x ~ ftype +s(sdoy,bs="cc",k=12,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                          data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
    # s(sdoy):ftypeobs   9.674     10  413.507  <2e-16 ***
    # s(sdoy):ftypemod   9.887     10 1571.615  <2e-16 ***
    # s(gmst)            8.789      9   61.950  <2e-16 ***
    # ti(sdoy,gmst)     11.806     12    1.213  <2e-16 ***
    # s(sdoy,fYear)    373.580    999   18.354  <2e-16 ***
    # gmst too wiggly, RE poss too smooth?


    # H0 and H2 have very busy RE smooths compared to H1 suggesting select=TRUE is the driver for calming down the RE smooth.
    # H0 and H2 have more uncertain evolution of gmst smooth compared to H1
    # H4 gsmt smooth too wiggly THUS need select=TRUE and min.sp for both gmst and RE
} # end if(FALSE)

ft0 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                        data=data01, method="GCV.Cp", select=TRUE, sp=3)
sp_names              <- names(ft0$sp)
is_gmst               <- grepl("\\(gmst\\)", sp_names)  
is_sz                 <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
min_sp_floor          <- rep(0, length(ft0$sp))     # initialize min.sp vector with zeros (no penalty)
min_sp_floor[is_sz]   <- 1                          # choose your floor (e.g., 2)
min_sp_floor[is_gmst] <- 0.5  # 0.2 leaves some inter-decadal variability in the gmst smooth
ft1.regamH1b <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                        data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
# s(sdoy):ftypeobs  14.273     22 111.184  <2e-16 ***
# s(sdoy):ftypemod  20.796     22 520.592  <2e-16 ***
# s(gmst)            3.555      9  35.628  <2e-16 ***
# ti(sdoy,gmst)     12.000     12   1.116  <2e-16 ***
# s(sdoy,fYear)    430.019   1000  21.405  <2e-16 ***
q1.regamH1b <- predict(ft1.regamH1b)
plot(data01$x, pch=20, cex=.3)
lines(q1.regamH1b ,   col=2, lwd=2)
grid()
plot(data01$x-q1.regamH1b, pch=20, cex=.3)
grid()

# Theos suggestions - quite slow

# no intervention
fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",by=ftype)
ft1.regamH1b0 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype)
ft1.regamH1b1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1)
ft1.regamH1b2 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)




fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype)
ft1.regamH1b0 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1)
ft0      <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
sp_names                                      <- names(ft0$sp)                   
min_sp_floor                                  <- rep(0,      length(sp_names))      # initialize min.sp vector with zeros (no penalty)
min_sp_floor[grepl("fYear",      sp_names)]   <- 1    # choose your floor (e.g., 2)
min_sp_floor[grepl("\\(gmst\\)", sp_names)]   <- 0.5  # 0.2 leaves some inter-decadal variability in the gmst smooth
ft1.regamH1b1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)


ft1.regamH1b2 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",by=ftype), 
                        data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)




# ft0 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz"), data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
sp_names              <- c() # names(ft0$sp)                   
is_gmst               <- grepl("\\(gmst\\)", sp_names) 
is_sz                 <- grepl("fYear",      sp_names)      
min_sp_floor          <- rep(0,      length(sp_names))      # initialize min.sp vector with zeros (no penalty)
min_sp_floor[is_sz]   <- 1    # choose your floor (e.g., 2)
min_sp_floor[is_gmst] <- 0.5  # 0.2 leaves some inter-decadal variability in the gmst smooth
ft1.regamH1b2 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz"), 
                        data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)


# vary min.sp for gmst and RE smooths
is_gmst               <- grepl("\\(gmst\\)", sp_names)  
is_sz                 <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
min_sp_floor          <- rep(0, length(ft0$sp))     # initialize min.sp vector with zeros (no penalty)
min_sp_floor[is_sz]   <- 0.2                        # choose your floor (e.g., 2)
min_sp_floor[is_gmst] <- 0.5  # 0.2 leaves some inter-decadal variability in the gmst smooth
ft1.regamH1c <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1), 
                        data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
# reducing RE below 1 starts to alter the gmst shape to be more complicated and more uncertain.  However with RE min.sp=1 there is 
# some structure in the residuals for obs only that looks like internal variability.  climate data seems structure free
q1.regamH1c <- predict(ft1.regamH1c)
plot(data01$x-q1.regamH1c, pch=20, cex=.3)
grid()

plot(ft1.regamH1c, page=1, main="H1c")


ft10 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype), 
                        data=data01, method="GCV.Cp", select=TRUE)
sp_names              <- names(ft10$sp)
is_gmst               <- grepl("\\(gmst\\)", sp_names)  
is_sz                 <- grepl("fYear", sp_names)   # pattern matches s(sdoy, fYear, bs="sz", ...)
min_sp_floor          <- rep(0, length(ft10$sp))     # initialize min.sp vector with zeros (no penalty)
min_sp_floor[is_sz]   <- 1.0                        # choose your floor (e.g., 2)
min_sp_floor[is_gmst] <- 0.5  # 0.2 leaves some inter-decadal variability in the gmst smooth
min_sp_floor          <- c(min_sp_floor,1.0,1.0) # fix bug in copilots code
# min_sp_floor[1:10] 0 0 0.5 0.5 0 0 1 1 1 1
# min_sp_floor          <- c(0,0, 0.5,0.5, 0,0, 1,1,1,1)


cat("Fitting ft1.regamH1d",cr)
min_sp_floor <- c(0,0, 0.5,0.5, 0,0, 1,1,1,1)
ft1.regamH1d <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype), 
                        data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor )
# s(sdoy):ftypeobs        14.286     22 11.733  <2e-16 ***
# s(sdoy):ftypemod        21.064     22 17.508  <2e-16 ***
# s(gmst)                  2.503      9  0.554  <2e-16 ***
# ti(sdoy,gmst)           12.000     12  0.062  <2e-16 ***
# s(sdoy,fYear):ftypeobs 185.320    459 19.157  <2e-16 ***
# s(sdoy,fYear):ftypemod 401.094    995 26.798  <2e-16 ***
q1.regamH1d <- predict(ft1.regamH1d)
plot(data01$x-q1.regamH1d, pch=20, cex=.3)
grid()

plot(ft1.regamH1d, page=1, main="H1d")

### odd thing happening with 2025 obs.  try excluding
i1 <- which(data01$year==2025 & data01$isobs==1)
ft1.regamH1d2 <- gam(x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype), 
                        data=data01[-i1,], method="GCV.Cp", select=TRUE, min.sp=min_sp_floor )

library(gratia)
draw(ft1.regamH1d, residuals=TRUE, main="H1d")


smooths(ft1.regamH1d)
[1] "s(sdoy):ftypeobs"       "s(sdoy):ftypemod"       "s(gmst)"               
[4] "ti(sdoy,gmst)"          "s(sdoy,fYear):ftypeobs" "s(sdoy,fYear):ftypemod"

sm1 <- smooth_estimates(ft1.regamH1d)
sm2 <- smooth_estimates(ft1.regamH1d, select="s(sdoy,fYear):ftypeobs")
iy2 <- which(data01$year %in% 1999 & data01$isobs==1)
sm2b <- smooth_estimates(ft1.regamH1d, select="s(sdoy,fYear):ftypeobs", data=data01[iy2,])

iobs <- which(data01$isobs==1)
range(data01$year[iobs])
pdf(file="./plot_RE_byYear_obs.pdf", width=12, height=9)
p <- draw(sm2)
print(p)
for( i in 1980:2025 ) {
  # iy2  <- which(data01$year %in% i & data01$isobs==1)
  iy2  <- which(data01$year == i & data01$isobs==1)
  cat(data01$year[iy2[1]], "\n")
  sm2b <- smooth_estimates(ft1.regamH1d, select="s(sdoy,fYear):ftypeobs", data=data01[iy2,])
  p <- draw(sm2b)
  print(p)
}
dev.off()

### sdoy, 1980+2024, obs and model
iy2  <- which(data01$year %in% c(1981,2024) ) #  & data01$isobs==1)
plot(data01$sdoy[iy2], data01$x[iy2], pch=20, cex=.3)
sm2b <- smooth_estimates(ft1.regamH1d, select="s(sdoy):ftypeobs", data=data01[iy2,])
draw(sm2b)
sm2b <- smooth_estimates(ft1.regamH1d, select="s(sdoy):ftypemod", data=data01[iy2,])
draw(sm2b)

pr.0 <- predict(ft1.regamH1d)
plot(data01$x, pch=20, cex=.3)
lines(pr.0, col=2, lwd=2)

iy2  <- which(data01$year %in% c(1981,2024) ) #  & data01$isobs==1)
tr.1 <- predict(ft1.regamH1d, newdata=data01[iy2,], type="terms")
pr.1 <- predict(ft1.regamH1d, newdata=data01[iy2,])



library("ggplot2")
library("dplyr")
sm2 |>
  add_confint() |>
  ggplot(aes(y = .estimate, x = sdoy)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.2, fill = "forestgreen"
  ) +
  geom_line(colour = "forestgreen", linewidth = 1.5) +
  labs(
    y = "Partial effect",
    # title = expression("Partial effect of" ~ f(x[2])),
    # x = expression(x[2])
    title = expression("Partial effect of" ~ s(sdoy,fYear):ftypeobs),
    x = expression(~ s(sdoy,fYear):ftypeobs)
  )









save(file="./GAM_random_effect.RData", ft1.regamH1d, data01)
