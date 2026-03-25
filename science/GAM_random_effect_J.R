# source("GAM_random_effect_J.R")

library(mgcv)
library(data.table)
library(PCICt)
library(gratia)
source("copilot_fn.R")

# load the data
load("RData/data01_NWS.RData", verb=TRUE)
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


### continuous random effects
data01$stime <- (data01$time - 2000)/100
plot(data01$time, data01$stime, ty='l', main="stime" )

readline("Stop1?")

### No RE baseline model
fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) 
ft1.regamJ0 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
# plot(ft1.regamJ0, pages=1, shade=TRUE)
plot_regam(ft1.regamJ0, stpdf="./GAM_random_effect_regamJ0.pdf")
save(ft1.regamJ0, data01, file="ft1_regamJ0.RData")


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp')       +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="ts",by=ftype, k=140*8, sp=c(40,10))
ft1.regamJ1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ1)
# plot(ft1.regamJ1, pages=1)
plot_regam(ft1.regamJ1, stpdf="./GAM_random_effect_regamJ1.pdf")
save(ft1.regamJ1, data01, file="ft1_regamJ1.RData")
## seems to mess with the gmst term; trend in the stime term which we want to be stationary


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, k=140*8, sp=c(40,10))
ft1.regamJ2 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2)
# plot(ft1.regamJ2, pages=1)
plot_regam(ft1.regamJ2, stpdf="./GAM_random_effect_regamJ2.pdf")
save(ft1.regamJ2, data01, file="ft1_regamJ2.RData")

### m=1
fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(40,10))
ft1.regamJ2m1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m1)
# plot(ft1.regamJ2m1, pages=1)
plot_regam(ft1.regamJ2m1, stpdf="./GAM_random_effect_regamJ2m1.pdf")
save(ft1.regamJ2m1, data01, file="ft1_regamJ2m1.RData")

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(20,5))
ft1.regamJ2m1b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m1b)
# plot(ft1.regamJ2m1b, pages=1)
plot_regam(ft1.regamJ2m1b, stpdf="./GAM_random_effect_regamJ2m1b.pdf")
save(ft1.regamJ2m1b, data01, file="ft1_regamJ2m1b.RData")

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(10,2))
ft1.regamJ2m1c <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m1c)
# plot(ft1.regamJ2m1c, pages=1)
plot_regam(ft1.regamJ2m1c, stpdf="./GAM_random_effect_regamJ2m1c.pdf")
save(ft1.regamJ2m1c, data01, file="ft1_regamJ2m1c.RData")

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(1,0.02))
ft1.regamJ2m1d <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m1d)
# plot(ft1.regamJ2m1d, pages=1)
plot_regam(ft1.regamJ2m1d, stpdf="./GAM_random_effect_regamJ2m1d.pdf")
save(ft1.regamJ2m1d, data01, file="ft1_regamJ2m1d.RData")


### m=2
fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=2, k=140*8, sp=c(40,10))
ft1.regamJ2m2 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m2)
# plot(ft1.regamJ2m2, pages=1)
plot_regam(ft1.regamJ2m2, stpdf="./GAM_random_effect_regamJ2m2.pdf")
save(ft1.regamJ2m2, data01, file="ft1_regamJ2m2.RData")

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=2, k=140*8, sp=c(20,5))
ft1.regamJ2m2b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m2b)
plot_regam(ft1.regamJ2m2b, stpdf="./GAM_random_effect_regamJ2m2b.pdf")
save(ft1.regamJ2m2b, data01, file="ft1_regamJ2m2b.RData")


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=2, k=140*8, sp=c(10,2))
ft1.regamJ2m2c <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m2c)
plot_regam(ft1.regamJ2m2c, stpdf="./GAM_random_effect_regamJ2m2c.pdf")
save(ft1.regamJ2m2c, data01, file="ft1_regamJ2m2c.RData")






# NOT SURE THE PLOTTING OF GP random effects IS WORKING
# ft1 <- ft1.regamJ2m1
# iobs <- which(data01$isobs==1)
# smo <- smooth_estimates(ft1, select="s(stime):ftypeobs", data=data01[iobs,])
# pp <- draw(smo)
# print(pp)

# print(draw(smooth_estimates(ft1.regamJ2m1, select="s(stime):ftypeobs", data=data01[iobs,])))

# plot(ft1.regamJ2m1$residuals,  pch=20,cex=.3, main="ft1.regamJ2m1 residuals", xlab="index", ylab="residuals")
# plot(ft1.regamJ2m1b$residuals, pch=20,cex=.3, main="ft1.regamJ2m1 residuals", xlab="index", ylab="residuals")



# J3 aka bs='sz'
fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="sz",by=ftype, k=140*8, sp=c(40,10))
ft1.regamJ3 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ3)
# plot(ft1.regamJ3, pages=1)
plot_regam(ft1.regamJ3, stpdf="./GAM_random_effect_regamJ3.pdf")
save(ft1.regamJ3, data01, file="ft1_regamJ3.RData")
## identical to ft1.regamJ1

fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="sz",by=ftype, k=140*8, sp=c(20,5))
ft1.regamJ3b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ3b)
plot_regam(ft1.regamJ3b, stpdf="./GAM_random_effect_regamJ3b.pdf")
save(ft1.regamJ3b, data01, file="ft1_regamJ3b.RData")

fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="sz",by=ftype, k=140*8, sp=c(10,2))
ft1.regamJ3c <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ3c)
plot_regam(ft1.regamJ3c, stpdf="./GAM_random_effect_regamJ3c.pdf")
save(ft1.regamJ3c, data01, file="ft1_regamJ3c.RData")


# p <- draw(smooth_estimates(ft1.regamJ2m1))
# print(p)




# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy, fYear, bs="fs") +s(fYear, bs="re") 
fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
ft1.regamJ4a <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ4a)
plot_regam(ft1.regamJ4a, stpdf="./GAM_random_effect_regamJ4a.pdf")
save(ft1.regamJ4a, data01, file="ft1_regamJ4a.RData")

fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
ft1.regamJ4b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ4b)
plot_regam(ft1.regamJ4b, stpdf="./GAM_random_effect_regamJ4b.pdf")
save(ft1.regamJ4b, data01, file="ft1_regamJ4b.RData")

source("copilot_fn.R")
# ftype
# s(sdoy):ftypeobs
# s(sdoy):ftypemod
# s(gmst)
# ti(sdoy,gmst)
# s(sdoy,fYear)
st.terms <- c("ftype",smooths(ft1.regamJ4a))
iob      <- which(data01$isobs==1)

plot_gam_terms(ft1.regamJ4a, terms=st.terms[2:5], main='All terms',data=data01[-iob,], sort_by="time", col="red",  SE=FALSE, xlim=range(data01$time))
plot_gam_terms(ft1.regamJ4a, terms=st.terms[2:5], main='All terms',data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()

plot_gam_terms(ft1.regamJ4a, terms="s(gmst)", data=data01[iob,],  sort_by="time", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$time) )
plot_gam_terms(ft1.regamJ4a, terms="s(gmst)", data=data01[-iob,], sort_by="time", add=TRUE,  col="red")
grid()

plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)",    data=data01[iob,],  sort_by="time", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$time) )
plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)",    data=data01[-iob,], sort_by="time", add=TRUE,  col="red")
plot_gam_terms(ft1.regamJ4a, terms="s(sdoy):ftypemod", data=data01[-iob,], sort_by="time", add=TRUE, LINES=FALSE, SE=FALSE,  col=rgb(0.7,0,0,0.2))
grid()

plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)", data=data01[iob,],  sort_by="sdoy", add=FALSE, LINES=FALSE, SE=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$sdoy) )
plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)", data=data01[-iob,], sort_by="sdoy", add=TRUE, LINES=FALSE,  SE=FALSE, col="red")
grid()

plot_gam_terms(ft1.regamJ4a, terms="s(sdoy):ftypeobs", data=data01[iob,],  sort_by="sdoy", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$sdoy) )
plot_gam_terms(ft1.regamJ4a, terms="s(sdoy):ftypemod", data=data01[-iob,], sort_by="sdoy", add=TRUE,  col="red")
i1 <- which(data01$isobs==1 & data01$year %in% c(1981,2024,2079))
i2 <- which(data01$isobs==0 & data01$year %in% c(1981,2024,2079))
plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)", data=data01[i1,], sort_by="sdoy", LINES=FALSE, SE=FALSE, add=TRUE, col=rgb(0,0.6,0.6,0.2))
plot_gam_terms(ft1.regamJ4a, terms="ti(sdoy,gmst)", data=data01[i2,], sort_by="sdoy", LINES=FALSE, SE=FALSE, add=TRUE, col=rgb(0.7,0,0,0.2) )
grid()


plot_gam_terms(ft1.regamJ4a, terms=c("s(gmst)","ti(sdoy,gmst)"), data=data01[-iob,], sort_by="time", add=FALSE, col="red", xlim=range(data01$time) )
plot_gam_terms(ft1.regamJ4a, terms=c("s(gmst)","ti(sdoy,gmst)"), data=data01[iob,],  sort_by="time", add=TRUE,  col="blue")
grid()


plot_gam_terms(ft1.regamJ4a, terms=st.terms[-1], data=data01[-iob,], sort_by="time", col="red", xlim=range(data01$time), SE=FALSE)
plot_gam_terms(ft1.regamJ4a, terms=st.terms[-1], data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()

plot_gam_terms(ft1.regamJ4a, terms=st.terms[6], data=data01[-iob,], sort_by="time", col="red", xlim=range(data01$time), SE=FALSE)
plot_gam_terms(ft1.regamJ4a, terms=st.terms[6], data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()



























# ### plo diagnostics
# q1.regamJ2m1 <- predict(ft1.regamJ2m1)
# plot(data01$x, pch=20, cex=.3, main="ft1.regamJ2m1", xlab="index", ylab="Temperature")
# points(q1.regamJ2m1, col=2, pch=20, cex=.3)
# grid()
# plot(data01$x-q1.regamJ2m1, pch=20, cex=.3, main="ft1.regamJ2m1", xlab="index", ylab="Residuals")

# ### try Theo's sim
# model <- ft1.regamJ2m1
# DoYindex   <- grep("sdoy", names(coef(model)))
# Stimeindex <- grep("stime",names(coef(model)))
# # model matrix
# X <- predict(model,type="lpmatrix")  # [1:366, 1:414] [knots, coefs]
# # GAM coefficients
# b <- coef(model)                     # [414]    
# # estimated seasonal cycles
# SC <- X[,DoYindex] %*%  b[DoYindex]  # [1:16723, 1]
# # add the intercept to put it on the scale of the data
# SC <- SC + b[1] + b[2]
# # Check the seasonal cycle fit
# plot(data01$x,pch=20, cex=.3, main="s(sdoy,bs=cc)+s(sdoy,fYear,bs=sz)")
# lines(SC,col="red",lwd=2)
# readline("Continue?2")
# plot(data01$time, data01$x-SC,pch=20, cex=.3, main="x - s(sdoy,bs=cc)+s(sdoy,fYear,bs=sz)")
# grid()
# readline("Continue?2b")




  # plot residuals
    # up.1()
    # plot(ft1$re, pch=20, cex=.3, main=paste("RE: ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="RE")
    # grid()

