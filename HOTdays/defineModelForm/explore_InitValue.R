# source("explore_InitValue.R")

### Use this file to determine the form of Initialising an event model
st.pwd <- system("pwd", intern=TRUE)
# source(paste(st.pwd,"/../../setup_all.R",sep=''))
source("setup_dev.R")


    MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/MSref/ostia_cdr_nrt_regions.MSref.2026-03-26.RData"
    load(MSref.file, verb=TRUE)

    list2env(MSconfig ,       envir = .GlobalEnv)
    list2env(MSconfig$files , envir = .GlobalEnv)

    load(st_msdata01, verb=TRUE)
    i0 <- which(data01$isobs==1)

    source(paste(st.pwd,"/../../make_stationary/setup_MakeStationary.R",sep=''))

    st_events  <- list.files(glue(dirname(MSref.file),"/../Events"), pattern='_EventsTh', full.names=TRUE)
    st_events  <- st_events[grep('.RData', st_events)]
    load(st_events, verb=TRUE)
    list2env(events01$info , envir = .GlobalEnv)

    source(paste(st.pwd,"/../setup_HotDays.R",sep=''))

readline("Stop1")

############################################################################################################
### doFitInitValue.R #######################################################################################

list2env(HDconfig$initV , envir = .GlobalEnv)

gpd.th.u  <- initV$th.u
scale.doy <- list(obs=(366+1), mod=(360+1)) # bit of a bodge as assumes all obs are from leap years

####### evgam fomula to fit
#### ally
fmla.IVgpd             <- list()
fmla.IVgpd$S0.G0       <- list(excess ~ 1                                                                                                            , ~ 1 )
fmla.IVgpd$SB.G0       <- list(excess ~ class                                                                                                        , ~ 1 )
fmla.IVgpd$SBD.G0      <- list(excess ~ class +s(sdoy, bs="cc",k=iv.k$sdoy,by=class)                                                                 , ~ 1 )
fmla.IVgpd$SBDT.G0     <- list(excess ~ class +s(sdoy, bs="cc",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ 1 )
fmla.IVgpd$SBDTi.G0    <- list(excess ~ class +s(sdoy, bs="cc",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("cc","ts")) , ~ 1 )
fmla.IVgpd$Si.G0       <- list(excess ~                                                                              ti(sdoy, gmst, bs=c("cc","ts")) , ~ 1 )
fmla.IVgpd$SBi.G0      <- list(excess ~ class                                                                       +ti(sdoy, gmst, bs=c("cc","ts")) , ~ 1 )
fmla.IVgpd$SDTi.G0     <- list(excess ~        s(sdoy, bs="cc",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("cc","ts")) , ~ 1 )
fmla.IVgpd$ST.G0       <- list(excess ~                                               s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ 1 )
fmla.IVgpd$SD.G0       <- list(excess ~        s(sdoy, bs="cc",k=iv.k$sdoy,by=class)                                                                 , ~ 1 )

fmla.IVgpd$S0.GB       <- list(excess ~ 1                                                                                                            , ~ class )
fmla.IVgpd$SB.GB       <- list(excess ~ class                                                                                                        , ~ class )
fmla.IVgpd$SBD.GB      <- list(excess ~ class +s(sdoy, bs="cc",k=iv.k$sdoy,by=class)                                                                 , ~ class )
fmla.IVgpd$S0D.GB      <- list(excess ~        s(sdoy, bs="cc",k=iv.k$sdoy,by=class)                                                                 , ~ class )
fmla.IVgpd$S0T.GB      <- list(excess ~                                               s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ class )
fmla.IVgpd$S0DT.GB     <- list(excess ~        s(sdoy, bs="cc",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ class )
fmla.IVgpd$SBDTi.GB    <- list(excess ~ class +s(sdoy, bs="cc",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("cc","ts")) , ~ class )
initV$fmla.ally.IVgpd <- fmla.IVgpd

#### hotseas
fmla.IVgpd             <- list()
fmla.IVgpd$S0.G0       <- list(excess ~ 1                                                                                                            , ~ 1 )
fmla.IVgpd$SB.G0       <- list(excess ~ class                                                                                                        , ~ 1 )
fmla.IVgpd$SBD.G0      <- list(excess ~ class +s(sdoy, bs="ts",k=iv.k$sdoy,by=class)                                                                 , ~ 1 )
fmla.IVgpd$SBDT.G0     <- list(excess ~ class +s(sdoy, bs="ts",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ 1 )
fmla.IVgpd$SBDTi.G0    <- list(excess ~ class +s(sdoy, bs="ts",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("ts","ts")) , ~ 1 )
fmla.IVgpd$Si.G0       <- list(excess ~                                                                              ti(sdoy, gmst, bs=c("ts","ts")) , ~ 1 )
fmla.IVgpd$SBi.G0      <- list(excess ~ class                                                                       +ti(sdoy, gmst, bs=c("ts","ts")) , ~ 1 )
fmla.IVgpd$SDTi.G0     <- list(excess ~        s(sdoy, bs="ts",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("ts","ts")) , ~ 1 )
fmla.IVgpd$ST.G0       <- list(excess ~                                               s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ 1 )
fmla.IVgpd$SD.G0       <- list(excess ~        s(sdoy, bs="ts",k=iv.k$sdoy,by=class)                                                                 , ~ 1 )

fmla.IVgpd$S0.GB       <- list(excess ~ 1                                                                                                            , ~ class )
fmla.IVgpd$SB.GB       <- list(excess ~ class                                                                                                        , ~ class )
fmla.IVgpd$SBD.GB      <- list(excess ~ class +s(sdoy, bs="ts",k=iv.k$sdoy,by=class)                                                                 , ~ class )
fmla.IVgpd$S0D.GB      <- list(excess ~        s(sdoy, bs="ts",k=iv.k$sdoy,by=class)                                                                 , ~ class )
fmla.IVgpd$S0T.GB      <- list(excess ~                                               s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ class )
fmla.IVgpd$S0DT.GB     <- list(excess ~        s(sdoy, bs="ts",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst)                                  , ~ class )
fmla.IVgpd$SBDTi.GB    <- list(excess ~ class +s(sdoy, bs="ts",k=iv.k$sdoy,by=class) +s(gmst, bs="ts", k=iv.k$gmst) +ti(sdoy, gmst, bs=c("ts","ts")) , ~ class )
initV$fmla.hseas.IVgpd <- fmla.IVgpd






### choose use hotseason or not
## DO HOTSEASON
# do.events        <- events01
# fmla.IVgpd       <- initV$fmla.hseas.IVgpd
# st.InitVal.stats <- paste(dirname(HDconfig$initV$st_InitValue),'diag',sub('.RData','.hotseas.stats',basename(HDconfig$files$st_InitValue)),sep='/')
# st.InitVal.stats <- basename(st.InitVal.stats)
## DO ALL YEAR
do.events        <- events01$events01_ally
fmla.IVgpd       <- initV$fmla.ally.IVgpd
st.InitVal.stats <- paste(dirname(HDconfig$initV$st_InitValue),'diag',sub('.RData','.ally.stats',basename(HDconfig$files$st_InitValue)),sep='/')
st.InitVal.stats <- basename(st.InitVal.stats)

allIVfits <- list()
for(i in seq_along(fmla.IVgpd)) {
    allIVfits[[i]]      <- fn_jointFitInitVal(do.events, gpd.th.u, scale.doy, fmla.IVgpd[[i]], IsJoint=TRUE)
    allIVfits[[i]]$name <- names(fmla.IVgpd)[i]
    cat("Fitted IV",i,names(fmla.IVgpd)[i],cr)
}

### filter out bad fits
isFitGood         <- NULL
for(i in seq_along(fmla.IVgpd)) {
    x1 <- summary(allIVfits[[i]])
    pval <- c( x1[[1]]$logscale$'Pr(>|t|)', x1[[2]]$logscale$'Pr(>|t|)', x1[[2]]$shape$'Pr(>|t|)', x1[[1]]$shape$'Pr(>|t|)')
    if(all(pval < 1e-6)) {isFitGood[i] <- FALSE; cat(i,names(fmla.IVgpd)[i],pval,cr)} else isFitGood[i] <- TRUE
}
allIVfits  <-  allIVfits[isFitGood]
fmla.IVgpd <- fmla.IVgpd[isFitGood]

### save statistics of each fit

for(i in seq_along(allIVfits)) {
    if(i==1) fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=TRUE) else
             fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=FALSE)
}

for(i in seq_along(fmla.IVgpd)) {
    cat(cr,'################',cr,names(fmla.IVgpd)[i],cr)
    print(summary(allIVfits[[i]]))
}

### choose best fit
allIVfits.aic        <- lapply(allIVfits, AIC)
names(allIVfits.aic) <- names(fmla.IVgpd)
allIVfits.bic        <- lapply(allIVfits, BIC)
names(allIVfits.bic) <- names(fmla.IVgpd)
cat("BIC",cr)
isbic <- sort(unlist(allIVfits.bic), index.return=TRUE)$ix
for(i in isbic) cat(names(fmla.IVgpd)[i],allIVfits.bic[[i]],cr)

readline("Stop2")

cat(cr,cr,"##############################",cr,"IV GPD fit with lowest BIC",cr)
cat(names(fmla.IVgpd)[isbic[1]],allIVfits.bic[[isbic[1]]],cr)
summary(allIVfits[[isbic[1]]])

cat(cr,cr,"##############################",cr,"next",cr)
cat(names(fmla.IVgpd)[isbic[2]],allIVfits.bic[[isbic[2]]],cr)
summary(allIVfits[[isbic[2]]])

cat(cr,cr,"##############################",cr,"next",cr)
cat(names(fmla.IVgpd)[isbic[3]],allIVfits.bic[[isbic[3]]],cr)
summary(allIVfits[[isbic[3]]])

cat(cr,cr,"##############################",cr,"SD.G0",cr)
i1 <- which(names(fmla.IVgpd)[isbic]=='SD.G0')
cat(names(fmla.IVgpd)[isbic[i1]],allIVfits.bic[[isbic[i1]]],cr)
summary(allIVfits[[isbic[i1]]])

cat(cr,cr,"##############################",cr,"SBD.GB",cr)
i1 <- which(names(fmla.IVgpd)[isbic]=='SBD.GB')
cat(names(fmla.IVgpd)[isbic[i1]],allIVfits.bic[[isbic[i1]]],cr)
summary(allIVfits[[isbic[i1]]])

# ### 2026.04.01 - for whole year SBD.GB seems the way to go
# allyear (hotseason similar)
# summary(allIVfits[[13]]) # SBD.GB
# logscale
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -1.53       0.09  -16.46   <2e-16
# classmod        0.21       0.11    1.93   0.0269

# shape
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)     0.04       0.07    0.53    0.299
# classmod       -0.22       0.08   -2.74  0.00304

# ** Smooth terms **

# logscale
#                   edf max.df Chi.sq Pr(>|t|)
# s(sdoy):classobs 2.32      8  10.66   0.0088
# s(sdoy):classmod 6.33      8  40.27 4.88e-07
