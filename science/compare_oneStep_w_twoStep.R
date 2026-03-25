# source("compare_oneStep_w_twoStep.R")

library(mgcv)
library(data.table)
library(PCICt)
library(gratia)
source("copilot_fn.R")

DOONESTEP <- FALSE
DOTWOSTEP <- FALSE

######################################################################################
### onestep follows J4o
######################################################################################
if(DOONESTEP) {
    # load("RData/ft0-prelim-J4s.RData",verb=TRUE)
    # is_sz        <- which(grepl("sdoy,fYear,ftype", names(ft0$sp))) 
    # rm(ft0)
    # sp.sz        <- rep(20, length(is_sz))                  
    iob        <- which(data01$ftype=="obs")
    fYear2       <- factor(c(paste("o",data01$year[iob],sep=""), paste("m",data01$year[-iob],sep="")))
    data01[,fYear2 := fYear2]
    # fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
    # fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, ftype, fYear, bs="sz", sp=sp.sz)
    fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear2, bs="sz")
    ft0          <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
    sp.sz        <- rep(1, length(which(grepl("sdoy,fYear2", names(ft0$sp)))))                  
    rm(ft0)
    fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear2, bs="sz", sp=sp.sz)
    ft1.regamJ4o3 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 

    summary(ft1.regamJ4o3)
        # (Intercept)  12.1223     0.1014  119.57   <2e-16 ***
        # ftypemod      3.1917     0.1479   21.59   <2e-16 ***
        # s(sdoy):ftypeobs  14.84     26  2870.11  <2e-16 ***
        # s(sdoy):ftypemod  24.29     26 13918.92  <2e-16 ***
        # s(gmst)            4.00      4    76.24  <2e-16 ***
        # ti(sdoy,gmst)     12.00     12    58.42  <2e-16 ***
        # s(sdoy,fYear2)   584.95   1459    29.25  <2e-16 ***
        # R-sq.(adj) =  0.983   Deviance explained = 98.3%
        # GCV = 0.32873  Scale est. = 0.32473   n = 52693

    plot_regam(ft1.regamJ4o3, stpdf="plots/GAM_random_effect_regamJ4o3.pdf")
    save(ft1.regamJ4o3, data01,              file="RData/ft1_regamJ4o3.RData")
} else {
    load("RData/ft1_regamJ4o3.RData",verb=TRUE)
}
    #           plot(ft1.regamJ4o3,page=1)
    # plot_gam_terms(ft1.regamJ4o3, terms="s(sdoy,fYear2)", data=data01[-iob,], sort_by="time", col="red",  LINE=FALSE, SE=FALSE, main="s(sdoy,fYear2)")
    # plot_gam_terms(ft1.regamJ4o3, terms="s(sdoy,fYear2)", data=data01[ iob,], sort_by="time", col="blue", LINE=FALSE, SE=FALSE, add=TRUE)
    # grid()


######################################################################################
### two step flollows J4o but omits the RE effect by fYear2.  Insetead fits a second gam to the residuals of 
### the first step to capture multi-year interannual variability.  
######################################################################################
if(DOTWOSTEP) {
    iobs <- which(data01$ftype=="obs")
    ny.om <- length(unique(data01$fYear[iobs])) + length(unique(data01$fYear[-iobs]))
    ## fit 1
    fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=7, m=2) +ti(sdoy,gmst,bs=c("cc","ts")) 
    ft1.gam      <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 

    ## fit 2 A
    # can do multiple fits with different k to sepparate out different timescales
    # however, the final residuals are very similar, so unsure if multiple steps gain anything
        #     data01$resid1 <- resid(ft1.gam)
        #     fmla          <- resid1 ~ s(stime, bs='ts', k=floor(ny.om/4), by=ftype)
        #     ft2.gam       <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)  

        #     data01$resid2 <- resid(ft2.gam)
        #     fmla          <- resid2 ~ s(stime, bs='ts', k=ny.om, by=ftype)
        #     ft3.gam       <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)    
            
        #     data01$resid3 <- resid(ft3.gam)
        #     fmla          <- resid1 ~ s(stime, bs='ts', k=2*ny.om,by=ftype)
        #     ft4.gam       <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)    
            
        #     pr.twoA       <- predict(ft1.gam, se.fit=FALSE)
        #     pr.twoB       <- predict(ft2.gam, se.fit=FALSE)
        #     pr.ft3       <- predict(ft3.gam, se.fit=FALSE)
        #     pr.ft4       <- predict(ft4.gam, se.fit=FALSE)
        
        #       up.3()
        #     plot(ft1.gam$resid, pch=20, cex=0.3)
        #     lines(pr.twoB, col="red", lwd=2)

        #     plot(ft2.gam$resid, pch=20, cex=0.3)
        #     lines(pr.ft3, col="red", lwd=2)

        #     plot(ft3.gam$resid, pch=20, cex=0.3)
        #     lines(pr.ft4, col="red", lwd=2)

        #     plot(ft4.gam$resid, pch=20, cex=0.3)
    ## fit 2 B
    data01$resid1 <- resid(ft1.gam)
    fmla          <- resid1 ~ s(stime, bs='ts', k=2*ny.om, by=ftype)
    ft2.gam1      <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
    # # ts is identical to sz
    # fmla          <- resid1 ~ s(stime, bs='sz', k=2*ny.om, by=ftype)
    # ft2.gam2      <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
    # # sz is identical to ts
    # fmla          <- resid1 ~ s(stime, bs='gp', k=4*ny.om, by=ftype)
    # ft2.gam3      <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
        # pr.twoB1 <- predict(ft2.gam1, se.fit=FALSE)
        # pr.twoB2 <- predict(ft2.gam2, se.fit=FALSE)
        # pr.twoB3 <- predict(ft2.gam3, se.fit=FALSE)
        # plot(ft1.gam$resid, pch=20, cex=0.3)
        # lines(pr.twoB1, col="red",   lwd=2)
        # lines(pr.twoB2, col="green", lwd=2)
        # lines(pr.twoB3, col="blue",  lwd=2)
    ### cannot capture the shorter period variability with gp
    ft2.gam <- ft2.gam1

    summary(ft1.gam)
    summary(ft2.gam)

    plot_regam(ft1.gam, stpdf="plots/GAM_twostep_gam1.pdf")
    plot_regam(ft2.gam, stpdf="plots/GAM_twostep_gam2.pdf")
    save(ft1.gam, ft2.gam, data01, file="RData/ft1_twostep.RData")
} else {
    load("RData/ft1_twostep.RData",verb=TRUE)
}

# readline("Stop")
######################################################################################
### plotting
######################################################################################

pdf(file="plots/compare_oneStep_w_twoStep.pdf", width=9, height=6)

### one step #################################################
#
    #           plot(ft1.regamJ4o3,page=1)
    # plot_gam_terms(ft1.regamJ4o3, terms="s(sdoy,fYear2)", data=data01[-iob,], sort_by="time", col="red",  LINE=FALSE, SE=FALSE, main="s(sdoy,fYear2)")
    # plot_gam_terms(ft1.regamJ4o3, terms="s(sdoy,fYear2)", data=data01[ iob,], sort_by="time", col="blue", LINE=FALSE, SE=FALSE, add=TRUE)
    # grid()
    #
    # plot(ft2.gam$resid, pch=20, cex=0.3)
    # plot(data01$sdoy[ iobs], ft2.gam$resid[ iobs], pch=20, cex=0.3)
    # plot(data01$sdoy[-iobs], ft2.gam$resid[-iobs], pch=20, cex=0.3)
#
    pr.one  <- predict(ft1.regamJ4o3, se.fit=FALSE)
    
    plot(data01$x, pch=20, cex=0.3, main="One step")
    lines(pr.one,  pch=20, cex=0.3, col="red")

    term.one  <- predict(ft1.regamJ4o3, type = "terms", se.fit=FALSE)
    up.3()
    plot(term.one[,1], pch=20, cex=0.3, main="One step")
    plot(term.one[,2], pch=20, cex=0.3)
    plot(term.one[,3], pch=20, cex=0.3)
    up.3()
    plot(term.one[,4], pch=20, cex=0.3, main="One step")
    plot(term.one[,5], pch=20, cex=0.3)
    plot(term.one[,6], pch=20, cex=0.3)


### two step #################################################

    pr.twoA <- predict(ft1.gam, se.fit=FALSE)
    pr.twoB <- predict(ft2.gam, se.fit=FALSE)

    up.1()
    plot(data01$x, pch=20, cex=0.3, main="Two step A")
    lines(pr.twoA, pch=20, cex=0.3, col="red")

    plot(data01$resid1, pch=20, cex=0.3, main="Two step B")
    lines(pr.twoB, pch=20, cex=0.3, col="red",lwd=2)

    term.twoA  <- predict(ft1.gam, type = "terms", se.fit=FALSE)
    up.3()
    plot(term.twoA[,1], pch=20, cex=0.3, main="Two step A")
    plot(term.twoA[,2], pch=20, cex=0.3)
    plot(term.twoA[,3], pch=20, cex=0.3)
    up.3()
    plot(term.twoA[,4], pch=20, cex=0.3, main="Two step A")
    plot(term.twoA[,5], pch=20, cex=0.3)
    
    term.twoB  <- predict(ft2.gam, type = "terms", se.fit=FALSE)
    up.2()
    plot(term.twoB[,1], pch=20, cex=0.3, main="Two step B")
    plot(term.twoB[,2], pch=20, cex=0.3)

### compare one and two step #################################
    # pdf(file="plots/compare_10years.pdf", width=9, height=6)
    up.1()
    for(j in 1:5) {
        yrs <- 1980:1989 + (j-1)*10
        iy <- which(data01$fYear %in% yrs)
        plot(data01$x[iy], pch=20, cex=0.3, main=paste0("OBS/Model Years: ", yrs[1], "-", yrs[length(yrs)]))
        # lines(pr.twoA[iy], col="red", lwd=1)
        lines(pr.one[iy], col="red", lwd=1)
        lines(pr.twoA[iy]+pr.twoB[iy], col=rgb(0.15,0.7,0.8,1), lwd=1)
        grid()
        legend("bottomright", legend=c("One step","Two step"), col=c("red", rgb(0.15,0.7,0.8,1)), lwd=1, bty="n")
        if(STEP) readline("Continue?")
    }
    # dev.off()

    # pdf(file="plots/compare_10years_residuals.pdf", width=9, height=6)
    for(j in 1:5) {
        yrs <- 1980:1989 + (j-1)*10
        iy <- which(data01$fYear %in% yrs)
        # plot(data01$x[iy]-pr.twoA[iy], pch=20, cex=0.3, col="red", main=paste0("OBS/Model residuals Years: ", yrs[1], "-", yrs[length(yrs)]))
          plot(data01$x[iy]-pr.one[iy], pch=20, cex=0.3, col="red", main=paste0("OBS/Model residuals Years: ", yrs[1], "-", yrs[length(yrs)]))
        points(data01$x[iy]-(pr.twoA[iy]+pr.twoB[iy]), pch=20, cex=0.3, col=rgb(0.15,0.7,0.8,1))
        lines(1+0.1*pr.twoA[iy], col=1, lwd=1)
        grid()
        legend("bottomright", legend=c("One step","Two step"), col=c("red", rgb(0.15,0.7,0.8,1)), lwd=1, bty="n")
        if(STEP) readline("Continue?")
    }
    # dev.off()

    ### compare gmst term
      plot(term.one[ ,4], pch=20, cex=0.3, col='red', main="GMST term")
    points(term.twoA[,4], pch=20, cex=0.3, col=rgb(0.15,0.7,0.8,1))
    legend("bottomright", legend=c("One step","Two step"), col=c("red", rgb(0.15,0.7,0.8,1)), lwd=2, bty="n")

    ### compare doy-gmst interaction term
     plot(term.one[ ,5], ty='l', lwd=2, col='red', main="DOY-GMST interaction term")
    lines(term.twoA[,5],         lwd=1, col=rgb(0.15,0.7,0.8,1))
    legend("topleft", legend=c("One step","Two step"), col=c("red", rgb(0.15,0.7,0.8,1)), lwd=2, bty="n")


    ### compare doy + interaction term
     plot(rowSums(term.one[ ,c(2,3,5)]), ty='l', lwd=2, col='red', main="DOY-GMST interaction term")
    lines(rowSums(term.twoA[,c(2,3,5)]),         lwd=1, col=rgb(0.15,0.7,0.8,1))
    legend("topleft", legend=c("One step","Two step"), col=c("red", rgb(0.15,0.7,0.8,1)), lwd=2, bty="n")

    ### compare diff doy + interaction term
     plot(rowSums(term.one[ ,c(2,3,5)])-rowSums(term.twoA[,c(2,3,5)]), ty='l', lwd=2, col=1, main="One-Two: DOY-GMST interaction term")
    grid()

    ### compare diff doy + interaction term
     plot(rowSums(term.one[ ,2:5])-rowSums(term.twoA[,2:5]), ty='l', lwd=2, col=1, main="One-Two: DOY-GMST interaction term")
    grid()

dev.off()













#