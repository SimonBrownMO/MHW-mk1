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



# # fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='tp') +ti(sdoy,gmst,bs=c("cc","tp")) +s(sdoy, fYear, bs="fs") +s(fYear, bs="re") 
# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
# ft1.regamJ4a <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4a)
# plot_regam(ft1.regamJ4a, stpdf="./GAM_random_effect_regamJ4a.pdf")
# save(ft1.regamJ4a, data01, file="ft1_regamJ4a.RData")
# # bad but getting there

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
# ft1.regamJ4b <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4b)
# plot_regam(ft1.regamJ4b, stpdf="./GAM_random_effect_regamJ4b.pdf")
# save(ft1.regamJ4b, data01, file="ft1_regamJ4b.RData")
# # so, so bad

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", m=1)
# ft1.regamJ4c <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4c)
# plot_regam(ft1.regamJ4c, stpdf="./GAM_random_effect_regamJ4c.pdf")
# save(ft1.regamJ4c, data01, file="ft1_regamJ4c.RData")
# # so, so bad

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", m=2)
# ft1.regamJ4d <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4d)
# plot_regam(ft1.regamJ4d, stpdf="./GAM_random_effect_regamJ4d.pdf")
# save(ft1.regamJ4d, data01, file="ft1_regamJ4d.RData")
# # so, so bad

# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=c(10,2))
# ft1.regamJ4e <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4e)
# plot_regam(ft1.regamJ4e, stpdf="./GAM_random_effect_regamJ4e.pdf")
# save(ft1.regamJ4e, data01, file="ft1_regamJ4e.RData")
# # so, so bad

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", m=1)
# ft1.regamJ4f <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4f)
# plot_regam(ft1.regamJ4f, stpdf="./GAM_random_effect_regamJ4f.pdf")
# save(ft1.regamJ4f, data01, file="ft1_regamJ4f.RData")
# worse than ft1.regamJ4a

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", m=2)
# ft1.regamJ4g <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4g)
# plot_regam(ft1.regamJ4g, stpdf="./GAM_random_effect_regamJ4g.pdf")
# save(ft1.regamJ4g, data01, file="ft1_regamJ4g.RData")
# worse than ft1.regamJ4a


# is_sz        <- which(grepl("fYear", names(ft1.regamJ4a$sp)))
# sp.sz        <- rep(10, length(is_sz))
# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
# ft1.regamJ4h <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4h)
# plot_regam(ft1.regamJ4h, stpdf="./GAM_random_effect_regamJ4h.pdf")
# save(ft1.regamJ4h, data01, file="ft1_regamJ4h.RData")
## looking really good.  Too much smoothing in gmst and RE terms

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
# ft1.regamJ4i <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4i)
# plot_regam(ft1.regamJ4i, stpdf="./GAM_random_effect_regamJ4i.pdf")
# save(ft1.regamJ4i, data01, file="ft1_regamJ4i.RData")
# # alot worse than ft1.regamJ4a

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
# ft1.regamJ4j <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4j)
# plot_regam(ft1.regamJ4j, stpdf="./GAM_random_effect_regamJ4j.pdf")
# save(ft1.regamJ4j, data01, file="ft1_regamJ4j.RData")
# so, so bad

# fmla        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', cp=10) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz")
# ft1.regamJ4k <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4k)
# plot_regam(ft1.regamJ4k, stpdf="./GAM_random_effect_regamJ4k.pdf")
# save(ft1.regamJ4k, data01, file="ft1_regamJ4k.RData")



### penailse both the gmst and the RE terms
# load("ft1_regamJ4a.RData",verb=T)
# is_sz        <- which(grepl("fYear", names(ft1.regamJ4a$sp)))
# sp.sz        <- rep(20, length(is_sz))
# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4, m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
# ft1.regamJ4l <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4l)
# plot_regam(ft1.regamJ4l, stpdf="./GAM_random_effect_regamJ4l.pdf")
# save(ft1.regamJ4l, data01, file="ft1_regamJ4l.RData")
# # same as ft1.regamJ4h too smooth both gmst & RE terms.  bit more wiggle early on in gmst

# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
# ft1.regamJ4m <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4m)
# plot_regam(ft1.regamJ4m, stpdf="./GAM_random_effect_regamJ4m.pdf")
# save(ft1.regamJ4m, data01, file="ft1_regamJ4m.RData")
# # same as ft1.regamJ4h

# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=4, sp=10) +ti(sdoy,gmst,bs=c("cc","ts")) +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
# ft1.regamJ4n <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4n)
# plot_regam(ft1.regamJ4n, stpdf="./GAM_random_effect_regamJ4n.pdf")
# save(ft1.regamJ4n, data01, file="ft1_regamJ4n.RData")
# # same as ft1.regamJ4h


# sp.sz        <- rep(20, length(is_sz))
# fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
# ft1.regamJ4o <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
# summary(ft1.regamJ4o)
# plot_regam(ft1.regamJ4o, stpdf="./GAM_random_effect_regamJ4o.pdf")
# save(ft1.regamJ4o, data01, file="ft1_regamJ4o.RData")
    # (Intercept) 12.18473    0.04109  296.51   <2e-16 ***
    # ftypemod     3.12388    0.04888   63.91   <2e-16 ***
    # s(sdoy):ftypeobs     14.67     26  5893.95  <2e-16 ***
    # s(sdoy):ftypemod     23.55     26 33708.72  <2e-16 ***
    # s(gmst)               7.00      7  1924.37  <2e-16 ***
    # ti(sdoy,gmst)        11.90     12   118.31  <2e-16 ***
    # s(sdoy,fYear,ftype) 253.37   1000    15.87  <2e-16 ***
    # R-sq.(adj) =  0.976   Deviance explained = 97.6%
    # GCV = 0.4562  Scale est. = 0.4535    n = 52693
# # very good.  GMST slight wiggles RE more variability needed

sp.sz        <- rep(5, length(is_sz))
fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=6, m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
ft1.regamJ4p <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ4p)
plot_regam(ft1.regamJ4p, stpdf="./GAM_random_effect_regamJ4p.pdf")
save(ft1.regamJ4p, data01, file="ft1_regamJ4p.RData")
# bigger low freq wiggles in gmst compared to J4o, just as much structure in the residuals - so worse than J4o

sp.sz        <- rep(5, length(is_sz))
fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=6, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
ft1.regamJ4q <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ4q)
plot_regam(ft1.regamJ4q, stpdf="./GAM_random_effect_regamJ4q.pdf")
save(ft1.regamJ4q, data01, file="ft1_regamJ4q.RData")
# like J4p but lower freq. - so worse than J4o

sp.sz        <- rep(1, length(is_sz))
fmla         <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYear, ftype, bs="sz", sp=sp.sz)
ft1.regamJ4r <- gam(fmla, data=data01, method="GCV.Cp", select=TRUE) 
summary(ft1.regamJ4r)
plot_regam(ft1.regamJ4r, stpdf="./GAM_random_effect_regamJ4r.pdf")
save(ft1.regamJ4r, data01, file="ft1_regamJ4r.RData")
# more RE variability than J4o as expected but also more wiggles in gmst, residuals look identical - so worse than J4o

### add stime term tp J4o


# fmla0        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts")) +s(sdoy, fYear, ftype, bs="sz") +s(stime,bs="gp",by=ftype)
# ft0          <- gam(fmla0, data=data01, method="GCV.Cp", select=TRUE) 
# # smooths(ft0)
# # [1] "s(sdoy):ftypeobs"    "s(sdoy):ftypemod"    "s(gmst)"            
# # [4] "ti(sdoy,gmst)"       "s(sdoy,fYear,ftype)" "s(stime):ftypeobs"  
# # [7] "s(stime):ftypemod" 
# save(ft0, data01, file="ft0-prelim-J4s.RData")
library(mgcv)
library(data.table)
library(PCICt)
library(gratia)
source("copilot_fn.R")
load("ft0-prelim-J4s.RData",verb=T)

is_sz        <- which(grepl("sdoy,fYear,ftype", names(ft0$sp))) 
is_gp        <- which(grepl("stime",            names(ft0$sp)))
sp.sz        <- rep(20, length(is_sz))                  # 203
sp.gp        <- rep(20, length(is_gp))                  # 4
# ft0b          <- gam(fmla1, data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
# length_min_sp <- sum(sapply(ft0b$smooth, function(x) length(x$S)))
# min.sp        <- rep(0,length_min_sp)    # c(rep(NA,5), rep(20,length(is_sz)), rep(20,length(is_gp)))

min.sp <- c(rep(0,5), rep(20,length(is_sz)), rep(1,length(is_gp)))  # 20,20  20,1

ny <- length(c(unique(data01$year[which(data01$isobs==1)]),unique(data01$year[which(data01$isobs==0)])))

fmla1        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts")) +s(sdoy, fYear, ftype, bs="sz") +s(stime,bs="gp",by=ftype, k=ny)
ft1.regamJ4s <- gam(fmla1, data=data01, method="GCV.Cp", select=TRUE, min.sp=min.sp) 
summary(ft1.regamJ4s)
plot_regam(ft1.regamJ4s, stpdf="./GAM_random_effect_regamJ4s.pdf")
save(ft1.regamJ4s, data01, min.sp, file="ft1_regamJ4s.RData")

fmla1        <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=8, m=1) +ti(sdoy,gmst,bs=c("cc","ts")) +s(sdoy, fYear, ftype, bs="sz")           +s(stime,bs="sz",by=ftype, k=ny)
ft1.regamJ4t <- gam(fmla1, data=data01, method="GCV.Cp", select=TRUE, min.sp=min.sp) 
summary(ft1.regamJ4t)
plot_regam(ft1.regamJ4t, stpdf="./GAM_random_effect_regamJ4t.pdf")
save(ft1.regamJ4t, data01, min.sp, file="ft1_regamJ4t.RData")

### 2026.03.10  - no evidence that this approach, using a second RE component to capture 1-5 year variability works
###             - seems like ft1.regamJ4o is the way to go for now