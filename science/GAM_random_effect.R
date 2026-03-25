# source("GAM_random_effect.R")

library(mgcv)
library(data.table)
library(PCICt)
library(gratia)

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

# plot(data01$time, data01$gmst, ty='n', main="gmst" )
# points(data01$time[iob], data01$gmst[iob], pch=20, cex=.3, col=1)
# points(data01$time[imo], data01$gmst[imo], pch=20, cex=.3, col=2)

### original
# fmla <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1)
# ft0  <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
# sp_names                                      <- names(ft0$sp)                   
# min_sp_floor                                  <- rep(0,      length(sp_names)) # initialize min.sp vector with zeros (no penalty)
# min_sp_floor[grepl("fYear",      sp_names)]   <- 1                             # choose your floor (e.g., 2)
# min_sp_floor[grepl("\\(gmst\\)", sp_names)]   <- 0.5                           # 0.2 leaves some inter-decadal variability in the gmst smooth
# ft1.regamH1b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
# q1.regamH1b  <- predict(ft1.regamH1b)
# plot(data01$x, pch=20, cex=.3)
# lines(q1.regamH1b ,   col=2, lwd=2)
# grid()
# plot(data01$x-q1.regamH1b, pch=20, cex=.3)
# grid()


# continuous random effects
data01$stime <- (data01$time - 2000)/100
plot(data01$time, data01$stime, ty='l', main="stime" )

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) 
ft1.regamJ0 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
plot(ft1.regamJ0, pages=1, shade=TRUE)

fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp')       +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="ts",by=ftype, k=140*8, sp=c(40,10))
ft1.regamJ1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ1)
plot(ft1.regamJ1, pages=1)
## seems to mess with the gmst term; trend in the stime term which we want to be stationary


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(40,10))
ft1.regamJ2m1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ2m1)
plot(ft1.regamJ2m1, pages=1)


fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(stime,bs="sz",by=ftype, k=140*8, sp=c(40,10))
ft1.regamJ3 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ3)
plot(ft1.regamJ3, pages=1)
## very similar to ft1.regamJ1


p <- draw(smooth_estimates(ft1.regamJ2m1))
print(p)

### add some of the missing bits to theos model
  fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc"),by=ftype)
  # constrained
  sp_names <- c( "s(sdoy):ftypeobs", "s(sdoy):ftypemod", "s(gmst)1"        
    , "s(gmst)2"        , "ti(sdoy,gmst)1"  , "ti(sdoy,gmst)2"  
    , "s(sdoy,fYear)1"  , "s(sdoy,fYear)2"  , "s(sdoy,fYear)3"  
    , "s(sdoy,fYear)4"  , "s(sdoy,fYear)5"  , "s(sdoy,fYear)6"  
    , "s(sdoy,fYear)7"  , "s(sdoy,fYear)8"  , "s(sdoy,fYear)9"  
    , "s(sdoy,fYear)10" , "s(sdoy,fYear)11" , "s(sdoy,fYear)12" 
    , "s(sdoy,fYear)13" , "s(sdoy,fYear)14" , "s(sdoy,fYear)15" 
    , "s(sdoy,fYear)16" , "s(sdoy,fYear)17" , "s(sdoy,fYear)18" 
    , "s(sdoy,fYear)19" , "s(sdoy,fYear)20" , "s(sdoy,fYear)21" 
    , "s(sdoy,fYear)22" , "s(sdoy,fYear)23" , "s(sdoy,fYear)24" 
    , "s(sdoy,fYear)25" , "s(sdoy,fYear)26" , "s(sdoy,fYear)27" 
    , "s(sdoy,fYear)28" , "s(sdoy,fYear)29" , "s(sdoy,fYear)30" 
    , "s(sdoy,fYear)31" , "s(sdoy,fYear)32" , "s(sdoy,fYear)33" 
    , "s(sdoy,fYear)34" , "s(sdoy,fYear)35" , "s(sdoy,fYear)36" 
    , "s(sdoy,fYear)37" , "s(sdoy,fYear)38" , "s(sdoy,fYear)39" 
    , "s(sdoy,fYear)40" , "s(sdoy,fYear)41" , "s(sdoy,fYear)42" 
    , "s(sdoy,fYear)43" , "s(sdoy,fYear)44" , "s(sdoy,fYear)45" 
    , "s(sdoy,fYear)46" , "s(sdoy,fYear)47" , "s(sdoy,fYear)48" 
    , "s(sdoy,fYear)49" , "s(sdoy,fYear)50" , "s(sdoy,fYear)51" 
    , "s(sdoy,fYear)52" , "s(sdoy,fYear)53" , "s(sdoy,fYear)54" 
    , "s(sdoy,fYear)55" , "s(sdoy,fYear)56" , "s(sdoy,fYear)57" 
    , "s(sdoy,fYear)58" , "s(sdoy,fYear)59" , "s(sdoy,fYear)60" 
    , "s(sdoy,fYear)61" , "s(sdoy,fYear)62" , "s(sdoy,fYear)63" 
    , "s(sdoy,fYear)64" , "s(sdoy,fYear)65" , "s(sdoy,fYear)66" 
    , "s(sdoy,fYear)67" , "s(sdoy,fYear)68" , "s(sdoy,fYear)69" 
    , "s(sdoy,fYear)70" , "s(sdoy,fYear)71" , "s(sdoy,fYear)72" 
    , "s(sdoy,fYear)73" , "s(sdoy,fYear)74" , "s(sdoy,fYear)75" 
    , "s(sdoy,fYear)76" , "s(sdoy,fYear)77" , "s(sdoy,fYear)78" 
    , "s(sdoy,fYear)79" , "s(sdoy,fYear)80" , "s(sdoy,fYear)81" 
    , "s(sdoy,fYear)82" , "s(sdoy,fYear)83" , "s(sdoy,fYear)84" 
    , "s(sdoy,fYear)85" , "s(sdoy,fYear)86" , "s(sdoy,fYear)87" 
    , "s(sdoy,fYear)88" , "s(sdoy,fYear)89" , "s(sdoy,fYear)90" 
    , "s(sdoy,fYear)91" , "s(sdoy,fYear)92" , "s(sdoy,fYear)93" 
    , "s(sdoy,fYear)94" , "s(sdoy,fYear)95" , "s(sdoy,fYear)96" 
    , "s(sdoy,fYear)97" , "s(sdoy,fYear)98" , "s(sdoy,fYear)99" 
  , "s(sdoy,fYear)100", "s(sdoy,fYear)101", "s(sdoy,fYear)102")
  min_sp_floor                                  <- rep(0,      length(sp_names)) # initialize min.sp vector with zeros (no penalty)
  min_sp_floor[grepl("fYear",      sp_names)]   <- 1                             # choose your floor (e.g., 2)
  min_sp_floor[grepl("\\(gmst\\)", sp_names)]   <- 0.5                           # 0.2 leaves some inter-decadal variability in the gmst smooth
  ft1.regamH1b3c2 <- gam(fmla, data=data01, method="GCV.Cp", knots=list(sdoy=c(0,1)),select=TRUE, min.sp=min_sp_floor)
  plot_regam(ft1.regamH1b3c2, stpdf="./GAM_random_effect_regamH1b3c2.pdf")
### this is quite a bit worse than ft1.regamH1b3c and does not solve the year to year discontinuity issue



plot_regam(ft1.regamJ0, stpdf="./GAM_random_effect_regamJ0.pdf")
save(ft1.regamJ0, data01, file="ft1_regamJ0.RData")

# Theos suggestions - quite slow
if(FALSE) {
  # mk0 
  # Dont think we ever want this as for the overlap there will be two independent instances of the RE for each year.
  # fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz")
  # ft1.regamH1b0 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
  # plot_regam(ft1.regamH1b0, stpdf="./GAM_random_effect_regamH1b0.pdf")
  # save(ft1.regamH1b0, data01, file="ft1_regamH1b_testing0.RData")

  # mk1 - very slow and large memory requirements
  fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",by=ftype)
  ft1.regamH1b1 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
  plot_regam(ft1.regamH1b1, stpdf="./GAM_random_effect_regamH1b1.pdf")
  save(ft1.regamH1b1, data01, file="ft1_regamH1b_testing1.RData")

  # mk2 
  # Dont think we ever want this as for the overlap there will be two independent instances of the RE for each year.
  # fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1)
  # ft1.regamH1b2 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
  # plot_regam(ft1.regamH1b2, stpdf="./GAM_random_effect_regamH1b2.pdf")
  # save(ft1.regamH1b2, data01, file="ft1_regamH1b_testing2.RData")

  # mk3
  # ideal if we can get away with it as is faster than mk1 and lower memory requirement
  fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype)
  ft1.regamH1b3 <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE)
  plot_regam(ft1.regamH1b3, stpdf="./GAM_random_effect_regamH1b3.pdf")
  # constrained
  sp_names <- c( "s(sdoy):ftypeobs", "s(sdoy):ftypemod", "s(gmst)1"        
    , "s(gmst)2"        , "ti(sdoy,gmst)1"  , "ti(sdoy,gmst)2"  
    , "s(sdoy,fYear)1"  , "s(sdoy,fYear)2"  , "s(sdoy,fYear)3"  
    , "s(sdoy,fYear)4"  , "s(sdoy,fYear)5"  , "s(sdoy,fYear)6"  
    , "s(sdoy,fYear)7"  , "s(sdoy,fYear)8"  , "s(sdoy,fYear)9"  
    , "s(sdoy,fYear)10" , "s(sdoy,fYear)11" , "s(sdoy,fYear)12" 
    , "s(sdoy,fYear)13" , "s(sdoy,fYear)14" , "s(sdoy,fYear)15" 
    , "s(sdoy,fYear)16" , "s(sdoy,fYear)17" , "s(sdoy,fYear)18" 
    , "s(sdoy,fYear)19" , "s(sdoy,fYear)20" , "s(sdoy,fYear)21" 
    , "s(sdoy,fYear)22" , "s(sdoy,fYear)23" , "s(sdoy,fYear)24" 
    , "s(sdoy,fYear)25" , "s(sdoy,fYear)26" , "s(sdoy,fYear)27" 
    , "s(sdoy,fYear)28" , "s(sdoy,fYear)29" , "s(sdoy,fYear)30" 
    , "s(sdoy,fYear)31" , "s(sdoy,fYear)32" , "s(sdoy,fYear)33" 
    , "s(sdoy,fYear)34" , "s(sdoy,fYear)35" , "s(sdoy,fYear)36" 
    , "s(sdoy,fYear)37" , "s(sdoy,fYear)38" , "s(sdoy,fYear)39" 
    , "s(sdoy,fYear)40" , "s(sdoy,fYear)41" , "s(sdoy,fYear)42" 
    , "s(sdoy,fYear)43" , "s(sdoy,fYear)44" , "s(sdoy,fYear)45" 
    , "s(sdoy,fYear)46" , "s(sdoy,fYear)47" , "s(sdoy,fYear)48" 
    , "s(sdoy,fYear)49" , "s(sdoy,fYear)50" , "s(sdoy,fYear)51" 
    , "s(sdoy,fYear)52" , "s(sdoy,fYear)53" , "s(sdoy,fYear)54" 
    , "s(sdoy,fYear)55" , "s(sdoy,fYear)56" , "s(sdoy,fYear)57" 
    , "s(sdoy,fYear)58" , "s(sdoy,fYear)59" , "s(sdoy,fYear)60" 
    , "s(sdoy,fYear)61" , "s(sdoy,fYear)62" , "s(sdoy,fYear)63" 
    , "s(sdoy,fYear)64" , "s(sdoy,fYear)65" , "s(sdoy,fYear)66" 
    , "s(sdoy,fYear)67" , "s(sdoy,fYear)68" , "s(sdoy,fYear)69" 
    , "s(sdoy,fYear)70" , "s(sdoy,fYear)71" , "s(sdoy,fYear)72" 
    , "s(sdoy,fYear)73" , "s(sdoy,fYear)74" , "s(sdoy,fYear)75" 
    , "s(sdoy,fYear)76" , "s(sdoy,fYear)77" , "s(sdoy,fYear)78" 
    , "s(sdoy,fYear)79" , "s(sdoy,fYear)80" , "s(sdoy,fYear)81" 
    , "s(sdoy,fYear)82" , "s(sdoy,fYear)83" , "s(sdoy,fYear)84" 
    , "s(sdoy,fYear)85" , "s(sdoy,fYear)86" , "s(sdoy,fYear)87" 
    , "s(sdoy,fYear)88" , "s(sdoy,fYear)89" , "s(sdoy,fYear)90" 
    , "s(sdoy,fYear)91" , "s(sdoy,fYear)92" , "s(sdoy,fYear)93" 
    , "s(sdoy,fYear)94" , "s(sdoy,fYear)95" , "s(sdoy,fYear)96" 
    , "s(sdoy,fYear)97" , "s(sdoy,fYear)98" , "s(sdoy,fYear)99" 
  , "s(sdoy,fYear)100", "s(sdoy,fYear)101", "s(sdoy,fYear)102")
  min_sp_floor                                  <- rep(0,      length(sp_names)) # initialize min.sp vector with zeros (no penalty)
  min_sp_floor[grepl("fYear",      sp_names)]   <- 1                             # choose your floor (e.g., 2)
  min_sp_floor[grepl("\\(gmst\\)", sp_names)]   <- 0.5                           # 0.2 leaves some inter-decadal variability in the gmst smooth
  ft1.regamH1b3c <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE, min.sp=min_sp_floor)
  plot_regam(ft1.regamH1b3c, stpdf="./GAM_random_effect_regamH1b3c.pdf")

  # fails
  # fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs=c("cc","sc"),id=1,by=ftype)
  # fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sc",id=1,by=ftype)

  save(ft1.regamH1b3, ft1.regamH1b3c, data01, file="ft1_regamH1b_testing3.RData")


  # ft1.regamH1b0
  # x ~ ftype +s(sdoy, bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz")
  #                      edf Ref.df       F p-value    
  # s(sdoy):ftypeobs  13.405     22 114.875  <2e-16 ***
  # s(sdoy):ftypemod  20.802     22 530.477  <2e-16 ***
  # s(gmst)            3.495      9  35.483  <2e-16 ***
  # ti(sdoy,gmst)     12.000     12   0.946  <2e-16 ***
  # s(sdoy,fYear)    447.713   1000  21.373  <2e-16 ***
  # AIC 102031.9
  # BIC 106471.7
  # names(ft1.regamH1b2$sp)
  #   [1] "s(sdoy):ftypeobs" "s(sdoy):ftypemod" "s(gmst)1"        
  #   [4] "s(gmst)2"         "ti(sdoy,gmst)1"   "ti(sdoy,gmst)2"  
  #   [7] "s(sdoy,fYear)1"   "s(sdoy,fYear)2"   "s(sdoy,fYear)3"  
    #  [10] "s(sdoy,fYear)4"   "s(sdoy,fYear)5"   "s(sdoy,fYear)6"  
    #  [13] "s(sdoy,fYear)7"   "s(sdoy,fYear)8"   "s(sdoy,fYear)9"  
    #  [16] "s(sdoy,fYear)10"  "s(sdoy,fYear)11"  "s(sdoy,fYear)12" 
    #  [19] "s(sdoy,fYear)13"  "s(sdoy,fYear)14"  "s(sdoy,fYear)15" 
    #  [22] "s(sdoy,fYear)16"  "s(sdoy,fYear)17"  "s(sdoy,fYear)18" 
    #  [25] "s(sdoy,fYear)19"  "s(sdoy,fYear)20"  "s(sdoy,fYear)21" 
    #  [28] "s(sdoy,fYear)22"  "s(sdoy,fYear)23"  "s(sdoy,fYear)24" 
    #  [31] "s(sdoy,fYear)25"  "s(sdoy,fYear)26"  "s(sdoy,fYear)27" 
    #  [34] "s(sdoy,fYear)28"  "s(sdoy,fYear)29"  "s(sdoy,fYear)30" 
    #  [37] "s(sdoy,fYear)31"  "s(sdoy,fYear)32"  "s(sdoy,fYear)33" 
    #  [40] "s(sdoy,fYear)34"  "s(sdoy,fYear)35"  "s(sdoy,fYear)36" 
    #  [43] "s(sdoy,fYear)37"  "s(sdoy,fYear)38"  "s(sdoy,fYear)39" 
    #  [46] "s(sdoy,fYear)40"  "s(sdoy,fYear)41"  "s(sdoy,fYear)42" 
    #  [49] "s(sdoy,fYear)43"  "s(sdoy,fYear)44"  "s(sdoy,fYear)45" 
    #  [52] "s(sdoy,fYear)46"  "s(sdoy,fYear)47"  "s(sdoy,fYear)48" 
    #  [55] "s(sdoy,fYear)49"  "s(sdoy,fYear)50"  "s(sdoy,fYear)51" 
    #  [58] "s(sdoy,fYear)52"  "s(sdoy,fYear)53"  "s(sdoy,fYear)54" 
    #  [61] "s(sdoy,fYear)55"  "s(sdoy,fYear)56"  "s(sdoy,fYear)57" 
    #  [64] "s(sdoy,fYear)58"  "s(sdoy,fYear)59"  "s(sdoy,fYear)60" 
    #  [67] "s(sdoy,fYear)61"  "s(sdoy,fYear)62"  "s(sdoy,fYear)63" 
    #  [70] "s(sdoy,fYear)64"  "s(sdoy,fYear)65"  "s(sdoy,fYear)66" 
    #  [73] "s(sdoy,fYear)67"  "s(sdoy,fYear)68"  "s(sdoy,fYear)69" 
    #  [76] "s(sdoy,fYear)70"  "s(sdoy,fYear)71"  "s(sdoy,fYear)72" 
    #  [79] "s(sdoy,fYear)73"  "s(sdoy,fYear)74"  "s(sdoy,fYear)75" 
    #  [82] "s(sdoy,fYear)76"  "s(sdoy,fYear)77"  "s(sdoy,fYear)78" 
    #  [85] "s(sdoy,fYear)79"  "s(sdoy,fYear)80"  "s(sdoy,fYear)81" 
    #  [88] "s(sdoy,fYear)82"  "s(sdoy,fYear)83"  "s(sdoy,fYear)84" 
    #  [91] "s(sdoy,fYear)85"  "s(sdoy,fYear)86"  "s(sdoy,fYear)87" 
    #  [94] "s(sdoy,fYear)88"  "s(sdoy,fYear)89"  "s(sdoy,fYear)90" 
    #  [97] "s(sdoy,fYear)91"  "s(sdoy,fYear)92"  "s(sdoy,fYear)93" 
    # [100] "s(sdoy,fYear)94"  "s(sdoy,fYear)95"  "s(sdoy,fYear)96" 
    # [103] "s(sdoy,fYear)97"  "s(sdoy,fYear)98"  "s(sdoy,fYear)99" 
    # [106] "s(sdoy,fYear)100" "s(sdoy,fYear)101" "s(sdoy,fYear)102"
  #


  # ft1.regamH1b2
  # x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1)
  #                      edf Ref.df      F p-value    
  # s(sdoy):ftypeobs  19.546     22 142.68  <2e-16 ***
  # s(sdoy):ftypemod  21.916     22 696.62  <2e-16 ***
  # s(gmst)            6.023      9 104.58  <2e-16 ***
  # ti(sdoy,gmst)     12.000     12  26.07  <2e-16 ***
  # s(sdoy,fYear)    998.897   1000  51.66  <2e-16 ***
  # R-sq.(adj) =  0.983   Deviance explained = 98.4%
  # GCV = 0.31962  Scale est. = 0.31319   n = 52693
  # AIC 89413.89
  # BIC 98830.72
  # names(ft1.regamH1b2$sp)
  # [1] "s(sdoy):ftypeobs" "s(sdoy):ftypemod" 
  # [3] "s(gmst)1"         "s(gmst)2"        
  # [5] "ti(sdoy,gmst)1"   "ti(sdoy,gmst)2"   
  # [7]"s(sdoy,fYear)1"   "s(sdoy,fYear)2"  



  # ft1.regamH1b3
  # x ~ ftype + s(sdoy,bs="cc",k=24,by=ftype) + s(gmst) + ti(sdoy,gmst,bs=c("cc","tp")) + s(sdoy,fYear, bs="sz",id=1,by=ftype)
  #                              edf Ref.df      F  p-value    
  # s(sdoy):ftypeobs         9.49528     22  9.522  < 2e-16 ***
  # s(sdoy):ftypemod        16.28942     22 36.053  < 2e-16 ***
  # s(gmst)                  0.05457      9  0.000 5.49e-05 ***
  # ti(sdoy,gmst)            5.13316     12 19.119  < 2e-16 ***
  # s(sdoy,fYear):ftypeobs 458.83669    459 49.661  < 2e-16 ***
  # s(sdoy,fYear):ftypemod 999.27507   1000 88.404  < 2e-16 ***
  # R-sq.(adj) =   0.99   Deviance explained =   99%
  # GCV = 0.20079  Scale est. = 0.19511   n = 52693
  # AIC 64896.15
  # BIC 78134.28
  # names(ft1.regamH1b3$sp)
  # [1] "s(sdoy):ftypeobs"        "s(sdoy):ftypemod"       
  # [3] "s(gmst)1"                "s(gmst)2"               
  # [5] "ti(sdoy,gmst)1"          "ti(sdoy,gmst)2"         
  # [7] "s(sdoy,fYear):ftypeobs1" "s(sdoy,fYear):ftypeobs2"


  # ft1.regamH1b3c - constrained version of ft1.regamH1b3 min.sp=min_sp_floor
  # fmla     <- x ~ ftype +s(sdoy,bs="cc",k=24,by=ftype) +s(gmst) +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy,fYear,bs="sz",id=1,by=ftype)
  #                            edf Ref.df      F p-value    
  # s(sdoy):ftypeobs        15.699     22  8.642  <2e-16 ***
  # s(sdoy):ftypemod        21.063     22 17.757  <2e-16 ***
  # s(gmst)                  2.503      9  0.554  <2e-16 ***
  # ti(sdoy,gmst)           12.000     12  0.062  <2e-16 ***
  # s(sdoy,fYear):ftypeobs 185.313    459 19.175  <2e-16 ***
  # s(sdoy,fYear):ftypemod 401.094    995 26.798  <2e-16 ***
  # AIC 90924.56
  # BIC 96608.74
  # > names(ft1.regamH1b3c$sp)
  # [1] "s(sdoy):ftypeobs"        "s(sdoy):ftypemod"       
  # [3] "s(gmst)1"                "s(gmst)2"               
  # [5] "ti(sdoy,gmst)1"          "ti(sdoy,gmst)2"         
  # [7] "s(sdoy,fYear):ftypeobs1" "s(sdoy,fYear):ftypeobs2"
}

### plo diagnostics
q1.regamJ2m1 <- predict(ft1.regamJ2m1)
plot(data01$x, pch=20, cex=.3, main="ft1.regamJ2m1", xlab="index", ylab="Temperature")
points(q1.regamJ2m1, col=2, pch=20, cex=.3)
grid()
plot(data01$x-q1.regamJ2m1, pch=20, cex=.3, main="ft1.regamJ2m1", xlab="index", ylab="Residuals")

### try Theo's sim
model <- ft1.regamJ2m1
DoYindex   <- grep("sdoy", names(coef(model)))
Stimeindex <- grep("stime",names(coef(model)))
# model matrix
X <- predict(model,type="lpmatrix")  # [1:366, 1:414] [knots, coefs]
# GAM coefficients
b <- coef(model)                     # [414]    
# estimated seasonal cycles
SC <- X[,DoYindex] %*%  b[DoYindex]  # [1:16723, 1]
# add the intercept to put it on the scale of the data
SC <- SC + b[1] + b[2]
# Check the seasonal cycle fit
plot(data01$x,pch=20, cex=.3, main="s(sdoy,bs=cc)+s(sdoy,fYear,bs=sz)")
lines(SC,col="red",lwd=2)
readline("Continue?2")
plot(data01$time, data01$x-SC,pch=20, cex=.3, main="x - s(sdoy,bs=cc)+s(sdoy,fYear,bs=sz)")
grid()
readline("Continue?2b")


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




  # plot residuals
    # up.1()
    # plot(ft1$re, pch=20, cex=.3, main=paste("RE: ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="RE")
    # grid()


######################################################################################
simulate_gam <- function(model, newdata, nsim = 1000,
                         component = c("full", "smooth", "s(time)"),
                         response = TRUE) {
  
  component <- match.arg(component)
  
  # 1. LPMATRIX FOR NEWDATA
  X <- predict(model, newdata, type = "lpmatrix")
  
  # 2. POSTERIOR COVARIANCE OF COEFFICIENTS
  V <- vcov(model, unconditional = TRUE)
  beta_hat <- coef(model)
  
  # 3. EXTRACT COMPONENT IF REQUESTED
  if (component != "full") {
    
    # select columns belonging to the target smooth
    target_pattern <- if (component == "smooth") "s(" else component
    idx <- grepl(target_pattern, colnames(X), fixed = FALSE)
    
    if (!any(idx)) 
      stop("No columns found for component: ", component)
    
    # zero out other components
    X[, !idx] <- 0
  }
  
  # 4. SIMULATE COEFFICIENTS FROM MULTIVARIATE NORMAL
  library(MASS)
  beta_sim <- mvrnorm(nsim, mu = beta_hat, Sigma = V)
  
  # 5. SIMULATE PREDICTIONS ON LINK SCALE
  eta_sim <- beta_sim %*% t(X)   # matrix of nsim × nrow(newdata)
  
  # 6. OPTIONAL: TRANSFORM TO RESPONSE SCALE
  if (response) {
    eta_sim <- model$family$linkinv(eta_sim)
  }
  
  return(eta_sim)
} ######################################################################################