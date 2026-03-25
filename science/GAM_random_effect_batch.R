# source("GAM_random_effect_batch.R")

library(mgcv)
library(data.table)
library(PCICt)

plot_regam <- function(ft1,stpdf=NULL) {

    q1       <- predict(ft1) 
    iy2001   <- which(data01$year==2001 & data01$isobs==1)
    i0x      <- which.max(q1[ iy2001] )
    i0n      <- which.min(q1[ iy2001] )
    doyn     <- data01$doy[ iy2001][i0n]
    sdoyn    <- data01$sdoy[iy2001][i0n]
    doyx     <- data01$doy[ iy2001][i0x]
    sdoyx    <- data01$sdoy[iy2001][i0x]
    nd       <- data01
    nd$sdoy  <- sdoyn
    q1.n     <- predict(ft1,  newdata=nd )
    nd$sdoy  <- sdoyx
    q1.x     <- predict(ft1,  newdata=nd )
    # remove interannual variability
    nd       <- data01
    nd$fYear <- "2003"
    q1.2003  <- predict(ft1, newdata=nd) 
    if(!is.null(stpdf)) pdf(file=stpdf, width=12, height=9)

    up.1()
    plot(data01$x, pch=20, cex=.3, main=paste("NWS ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
    lines(q1,   col=2, lwd=2)
    lines(q1.n, col=4, lwd=2)
    lines(q1.x, col=3, lwd=2)
    ix <- which(data01$doy==doyx)
    lines(ix, q1.2003[ix], col=6, lwd=2)
    im <- which(data01$doy==doyn)
    lines(im, q1.2003[im], col=6, lwd=2)
    legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Winter min @2003 IAV", "Summer max", "Summer max @2003 IAV"), 
                         col=c(1,2,4,6,3,6), lwd=c(NA,2,2,2,2,2), pch=c(20,NA,NA,NA,NA,NA), bty="n", cex=1.2)
    grid()
    if(is.null(stpdf)) readline("continue?")

    up.1()
    ix <- 16000:18000-185
    plot(data01$x[ix], pch=20, cex=.3, main=paste("NWS ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
    lines(q1[ix],   col=2, lwd=2)
    lines(q1.n[ix], col=4, lwd=2)
    lines(q1.x[ix], col=3, lwd=2)
    iy <- which(data01$doy[ix]==doyx)
    lines(iy, q1.2003[ix][iy], col=6, lwd=2); points(iy, q1.2003[ix][iy], col=6, pch=20, cex=2.0)
    iz <- which(data01$doy[ix]==doyn)
    lines(iz, q1.2003[ix][iz], col=6, lwd=2); points(iz, q1.2003[ix][iz], col=6, pch=20, cex=2.0)
    legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Winter min @2003 IAV", "Summer max", "Summer max @2003 IAV"), 
                         col=c(1,2,4,6,3,6), lwd=c(NA,2,2,2,2,2), pch=c(20,NA,NA,20,NA,20), bty="n", cex=1.2)
    grid()
    if(!is.null(stpdf)) dev.off()
} ######################################################################################


plot_regam_RE <- function(ft1,stpdf=NULL, idx=1:4, do.years=1980:2025) {
    
    library(gratia)
    library("patchwork")

    if(!is.null(stpdf)) pdf(file=stpdf, width=12, height=9)

    sm1       <- smooth_estimates(ft1)
    names_ft1 <- unique(sm1[[1]])
    p1 <- draw(smooth_estimates(ft1, select=names_ft1[idx]))
    print(p1)
    if(is.null(stpdf)) readline("continue?")

    # sm.re.o <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypeobs", data=data01[which(data01$isobs==1),])
    # draw(sm.re.o)
    # sm.re.m <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypemod", data=data01[which(data01$isobs==0),])
    # draw(sm.re.m)

    ncol    <- 6
    nrow    <- 4
    newpage <- TRUE
    for( i in do.years ) {
      iy2  <- which(data01$year == i & data01$isobs==1)
      smo <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypeobs", data=data01[iy2,])
      if(newpage) {
        pp <- draw(smo) 
        newpage <- FALSE
      } else {
        pp <- pp + draw(smo)
      }
      iy2  <- which(data01$year == i & data01$isobs==0)
      smm <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypemod", data=data01[iy2,])
      pp <- pp + draw(smm)
      # print(p)
      # p1 + p2 +plot_layout(ncol = 2)
      # cat(data01$year[iy2[1]],length(pp), newpage, "\n")
      if(length(pp)==(ncol*nrow)) {
        print(pp + plot_layout(ncol=ncol, nrow=nrow))
        newpage <- TRUE
        readline("Continue?")
      }
    }
    print(pp + plot_layout(ncol=ncol, nrow=nrow))

    if(is.null(stpdf)) readline("continue?")
    if(!is.null(stpdf)) dev.off()
} ######################################################################################


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


# mk1 - very slow and large memory requirements
fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",by=ftype)
ft1.regamH1b1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
plot_regam(ft1.regamH1b1, stpdf="./GAM_random_effect_regamH1b1.pdf")
save(ft1.regamH1b1, data01, file="ft1_regamH1b_testing1.RData")




  # # plot residuals
  #   up.1()
  #   plot(ft1$re, pch=20, cex=.3, main=paste("RE: ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="RE")
  #   grid()
