# source("setup_HotDays.R")

### load MakeStationary reference
reload_MS_ref(MSref.file)
HDconfig <- list()

### create required names and directories (used later in this file so needs to be here)
hd_version     <- 'vHD2d'
HDSAVEDIR      <- paste(MSSAVEDIR,'HotDays/',hd_version,sep='')
# now set in setup_all.R # datestamp      <- paste(format(Sys.time(), "%Y-%m-%d-%H%M%S"),sep='_')
subdir         <- c('HDmeta','HDref','Events/plots','Hotseason/plots','InitEvent/plots','InitEvent/diag','InitValue/plots','InitValue/diag','HT/plots', 'Term/plots', 'Term/diag')
for(d1 in subdir) system(paste("mkdir -p",paste(HDSAVEDIR,d1,sep='/') ) )
st_metadata    <- sub('DDD',datestamp, sub('.RData','.HDmeta.DDD.RData',paste(HDSAVEDIR,subdir[1],basename(st_base),sep='/')))
st_HDconfig    <- sub('DDD',datestamp, sub('.RData','.HDref.DDD.RData', paste(HDSAVEDIR,subdir[2],basename(st_base),sep='/')))
# deg0C          <- 273.15
HDconfig$general             <- list()
HDconfig$general$cr          <- cr
HDconfig$general$deg0C       <- deg0C
HDconfig$general$hd_version  <- hd_version
HDconfig$general$HDSAVEDIR   <- HDSAVEDIR
HDconfig$general$datestamp   <- datestamp
HDconfig$general$st_metadata <- st_metadata

### DO DIAGNOSTIC PLOTS OR NOT
DOINITEVENTPLOT <- FALSE 
DOINITVALUEPLOT <- FALSE
DOHTPLOT        <- FALSE 
PARAMETERPLOT   <- FALSE 
PARPLTNEW       <- FALSE 
DOTERMPLOT      <- FALSE 

############################################################################################################
### fitted models to use ###################################################################################
# NB NB NB This bit is very important as it will decide what responses to clim change are allowed


# HDconfig$events       <- events

############################################################################################################
### doFitInitEvent.R #########################################################################################
SAVEINITEVENTDIAG    <- TRUE
ievent.th.u          <- events01$info$event.th.u # 0.90 # the GAM.fit can fail if this is too low with errors like: Error in eigen(hess1, symmetric = TRUE) : 0 x 0 matrix
ievent.k             <- list(doy=-1, cc=-1)  # doy=6 or 12, cc=5 or 3
ievent.covaraite     <- 'gmst'
# if(FUTURE_PROJECTION) ievent.k$cc <- -1 else ievent.k$cc <- -1     # GAM smoother (k) for climate change effect on heatwave initiation probability
st_InitEvent <- sub('.RData','_IEventThXXX.RData',basename(st_base))
st_InitEvent <- paste(HDSAVEDIR, 'InitEvent', sub('XXX',ievent.th.u,st_InitEvent),sep='/')

initI                   <- list()
initI$st_InitEvent      <- st_InitEvent
initI$SAVEINITEVENTDIAG <- SAVEINITEVENTDIAG
initI$DOINITEVENTPLOT   <- DOINITEVENTPLOT
initI$ievent.th.u       <- ievent.th.u
initI$ievent.k          <- ievent.k
initI$ievent.covaraite  <- ievent.covaraite
HDconfig$initI          <- initI


############################################################################################################
### doFitInitValue.R #########################################################################################
initV                   <- list()
initV$SAVEINITVALUEDIAG <- TRUE
initV$DOINITVALUEPLOT   <- DOINITVALUEPLOT
initV$th.u              <- events01$info$event.th.u # 0.90   ### NEED TO tune for best fit to data.  Has to be below the max of do.ptiles
initV$iv.k              <- list(sdoy=-1, gmst=-1)  # doy=6 or 12, cc=5 or 3
st_InitValue            <- sub('.RData','_IValueThXXX.RData',basename(st_base))
st_InitValue            <- paste(HDSAVEDIR, 'InitValue', sub('XXX', initV$th.u ,st_InitValue),sep='/')
initV$st_InitValue      <- st_InitValue

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
initV$fmla.ally.IVgpd$chosen <- NULL # "SBDTi.G0"  # OR NULL if want to dynamically select best fit.

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

### Either SELECT YOUR CHOSEN MODEL OR CHOOSE DYNAMIC selection ON AIC
### recommend using defineModelForm/explore_InitValue.R to explore the model fits and select the best one before setting this.
### If previously chosen
fmla.ally.IVgpd.chosen   <- "SD.G0"
fmla.hseas.IVgpd.chosen  <- "SD.G0"
### IF DYNAMIC
    # DO NOTHING
### ELSE
    # fmla.ally.IVgpd.chosen  <- "SBDTi.G0" # "SBD.G0"
    # fmla.hseas.IVgpd.chosen <- "SBDTi.G0" # "SBD.G0"
    # initV$fmla.ally.IVgpd   <- initV$fmla.ally.IVgpd[fmla.ally.IVgpd.chosen]
    # initV$fmla.hseas.IVgpd  <- initV$fmla.hseas.IVgpd[fmla.hseas.IVgpd.chosen]

HDconfig$initV         <- initV


############################################################################################################
### doFitHt.R #####################################################################################
SAVEHTDIAG    <- TRUE
ht.th.u       <- events01$info$event.th.u # 0.90
cr.th.u       <- events01$info$event.th.u # 0.90
ht.cov.sig.th <- 0.9    # likelihood ratio test needs to be above this for covariate to be used
lam.pen0      <- 0.0  # Needs to be zero unless you are ABSOLUTELY sure. penalty to reduce covaiate impact on mu & sigma

FORCEALPHASTRUCTURE <- FALSE # rather than let the testing decide
if(FORCEALPHASTRUCTURE) a.covs.force <- c(1,2,3)

ALLOW_MU_SIGMA_DEPENDENCE <- FALSE  # are Mu and Sigma allowed to depend on covariates?  Generally a bad idea.

### NB 1 == t+1|t+0; 2 == t+2|(t+1,t+0)
HTmaxorder   <- 1      ### 2019.10.04 Only use k=1 HT model and k=2 takes v long time
if(HTmaxorder>1) stop("HT models of order>1 are not proven to work")
lagnumber    <- HTmaxorder # HTmaxorder +1
nsim1        <- 2e4  # does not need to be large as most of the simulating is done later
ht.exclude.l <- -2 # 0.7  ### exclude temps below this as they do not represent synoptic conditinos associated with heatwave generating processes (based on tempterature)
                     ### used to be 0.001 to remove the clustering around zero
st_Ht <- sub('.RData','_HtThXXX.RData',basename(st_base))
st_Ht <- paste(HDSAVEDIR, 'HT', sub('XXX', ht.th.u ,st_Ht),sep='/')

fitht                           <- list()
fitht$st_Ht                     <- st_Ht
fitht$SAVEHTDIAG                <- SAVEHTDIAG
fitht$DOHTPLOT                  <- DOHTPLOT
# fitht$PLTFIT                    <- PLTFIT
fitht$PARAMETERPLOT             <- PARAMETERPLOT
fitht$PARPLTNEW                 <- PARPLTNEW
fitht$ht.th.u                   <- ht.th.u
fitht$cr.th.u                   <- cr.th.u
fitht$ht.cov.sig.th             <- ht.cov.sig.th
fitht$lam.pen0                  <- lam.pen0
fitht$FORCEALPHASTRUCTURE       <- FORCEALPHASTRUCTURE
# fitht$a.covs.force              <- a.covs.force
fitht$ALLOW_MU_SIGMA_DEPENDENCE <- ALLOW_MU_SIGMA_DEPENDENCE
fitht$HTmaxorder                <- HTmaxorder
fitht$lagnumber                 <- lagnumber
fitht$nsim1                     <- nsim1
fitht$ht.exclude.l              <- ht.exclude.l
# fitht$st_Ht                  <- st_Ht
HDconfig$fitht                  <- fitht

############################################################################################################
### doFitTerm.R #########################################################################################
SAVETERMDIAG    <- TRUE
term.th.u       <- events01$info$event.th.u # 0.90
durationmax     <- events01$info$event.length # 40 # HW terminated after this number of days.  This is not the value by which the covariate has been scaled, which is daymax
scaleOldTerm    <- 0.8  # long duration HW are pooled to calc a fixed term prob.
                           # In reality Pterm falls with age. This factor reduces the pooled probability to account for this
# term.tdac  <- gam(fin ~ s(tmp,k=3) +s(doy,bs=form.doy,k=4) +s(age,k=3) +s(cc,k=3),  data=termination_vars[ix1,],family="binomial")
tr.k            <- list(tmp=-1, doy=4, D=4, gmst=6, age=6, cc=6)
st_Term <- sub('.RData','_Term_ThXXX.RData',basename(st_base))
st_Term <- paste(HDSAVEDIR, 'Term', sub('XXX', term.th.u ,st_Term),sep='/')

### GAM
termGAM.fmla <- list()
# hotseas
 #day1
 day1          <- list()
 day1$T        <- fin ~  s(tmp, k=4)
 day1$TI       <- fin ~  s(tmp, k=4) +class
 day1$TIC      <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst)
 day1$TID      <- fin ~  s(tmp, k=4) +class                      +s(D,k=tr.k$D,by=class)
 day1$TICD     <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class)
 day1$TICDiCD  <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class) +ti(gmst,D)
 day1$TCDiCD   <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)          +ti(gmst,D)
 day1$TCD      <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)
 day1$TC       <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst)
 day1$TD       <- fin ~  s(tmp, k=4)                             +s(D,k=tr.k$D)

 #dayN
 dayN          <- list()
 dayN$T        <- fin ~  s(tmp)
 dayN$TI       <- fin ~  s(tmp) +class
 dayN$TIC      <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst)
 dayN$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class)
 dayN$TICDA    <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class) +s(A, k=5)
 dayN$TICDAiCD <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class) +s(A, k=5) +ti(gmst,D)
 dayN$TCDAiCD  <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)          +s(A, k=5) +ti(gmst,D)
 dayN$TCDA     <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)          +s(A, k=5)
 dayN$TCA      <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst)                         +s(A, k=5)
 dayN$TDA      <- fin ~  s(tmp)                             +s(D,k=tr.k$D)          +s(A, k=5)
 dayN$TIDA     <- fin ~  s(tmp) +class                      +s(D,k=tr.k$D)          +s(A, k=5)
 dayN$TA       <- fin ~  s(tmp)                                                     +s(A, k=5)
termGAM.fmla$hseas <- list(day1=day1, dayN=dayN)

# all year
 #day1
 day1          <- list()
 day1$T        <- fin ~  s(tmp, k=4)
 day1$TI       <- fin ~  s(tmp, k=4) +class
 day1$TIC      <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst)
 day1$TID      <- fin ~  s(tmp, k=4) +class                      +s(D,k=tr.k$D,bs='cc',by=class)
 day1$TICD     <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class)
 day1$TICDiCD  <- fin ~  s(tmp, k=4) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +ti(gmst,D,bs=c("ts",'cc'))
 day1$TCDiCD   <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +ti(gmst,D)
 day1$TCD      <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')
 day1$TC       <- fin ~  s(tmp, k=4)        +s(gmst,k=tr.k$gmst)
 day1$TD       <- fin ~  s(tmp, k=4)                             +s(D,k=tr.k$D,bs='cc')
 #dayN
 dayN          <- list()
 dayN$T        <- fin ~  s(tmp)
 dayN$TI       <- fin ~  s(tmp) +class
 dayN$TIC      <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst)
 dayN$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class)
 dayN$TICDA    <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +s(A, k=5)
 dayN$TA       <- fin ~  s(tmp)                                                             +s(A, k=5)
 dayN$TICDAiCD <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +s(A, k=5) +ti(gmst,D,bs=c("ts",'cc'))
 dayN$TCDAiCD  <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +s(A, k=5) +ti(gmst,D)
 dayN$TCDA     <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +s(A, k=5)
 dayN$TCA      <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst)                                 +s(A, k=5)
 dayN$TDA      <- fin ~  s(tmp)                             +s(D,k=tr.k$D,bs='cc')          +s(A, k=5)
 dayN$TIDA     <- fin ~  s(tmp) +class                      +s(D,k=tr.k$D,bs='cc')          +s(A, k=5)
 dayN$TA       <- fin ~  s(tmp)                                                             +s(A, k=5)

termGAM.fmla$ally   <- list(day1=day1, dayN=dayN)

### GLM
### this permutates by hand - more flexible
   # calculate all combinations of ICDA
    st.n      <- c('class','C','D','A')
    st.c      <- c('class','gmst','D','A') # c('class','gmst','DOY','Age')
    # day 1
    alld1.com <- NULL
    for(i in 1:3) alld1.com[[i]] <- combn(st.c[1:3],i)
    alld1n.com <- NULL
    for(i in 1:3) alld1n.com[[i]] <- combn(st.n[1:3],i)
    fm.glm.d1      <- list()
    fm.glm.d1n     <- NULL
    fm.glm.d1[[1]] <- 'fin ~ tmp'
    fm.glm.d1n[1]  <- 'T'
    count0         <- 2
    for(i in seq_along(alld1.com)) {
        for(j in seq_along(alld1.com[[i]][1,])) {
            fm.glm.d1[[count0]]        <- as.formula(paste('fin ~ tmp',paste('+',alld1.com[[i]][,j],collapse=' ', sep='')))
            # class(fm.glm.d1[[count0]]) <- "formula"
            fm.glm.d1n[count0]         <- paste('T',paste(alld1n.com[[i]][,j],collapse='', sep=''), sep='')
            count0                     <- count0 +1
        }
    }
    names(fm.glm.d1) <- c('T','TI','TC','TD','TIC','TID','TCD','TICD')
    # day N
    alldN.com <- NULL
    for(i in 1:4) alldN.com[[i]] <- combn(st.c,i)
    alldNn.com <- NULL
    for(i in 1:4) alldNn.com[[i]] <- combn(st.n,i)
    fm.glm.dN      <- list()
    fm.glm.dNn     <- NULL
    fm.glm.dN[[1]] <- 'fin ~ tmp'
    fm.glm.dNn[1]  <- 'T'
    count0         <- 2
    for(i in seq_along(alldN.com)) {
        for(j in seq_along(alldN.com[[i]][1,])) {
            fm.glm.dN[[count0]]       <- as.formula(paste('fin ~ tmp',paste('+',alldN.com[[i]][,j],collapse=' ', sep='')))
            # class(fm.glm.dN[[count0]]) <- "formula"
            fm.glm.dNn[count0]         <- paste('T',paste(alldNn.com[[i]][,j],collapse='', sep=''), sep='')
            count0                     <- count0 +1
        }
    }
    names(fm.glm.dN) <- c('T','TI','TC','TD','TA','TIC','TID','TIA','TCD','TCA','TDA','TICD','TICA','TIDA','TCDA','TICDA')

termGLM.fmla        <- list()
termGLM.fmla$hseas  <- list(day1=fm.glm.d1, dayN=fm.glm.dN)
termGLM.fmla$ally   <- list(day1=fm.glm.d1, dayN=fm.glm.dN) # both the same but to be consistent with GAM structure

term                     <- list()
term$st_Term             <- st_Term
term$SAVETERMDIAG        <- SAVETERMDIAG
term$DOTERMPLOT          <- DOTERMPLOT
# term$FORCETERMPROBFALL   <- FORCETERMPROBFALL
term$term.th.u           <- term.th.u
# term$term.prob.fall.tol  <- term.prob.fall.tol
# term$minprob.day1        <- minprob.day1
# term$minprob.dayN        <- minprob.dayN
# term$topout.tmp.l        <- topout.tmp.l
term$tr.k                <- tr.k
# term$mintermprob         <- mintermprob
term$fmGAM               <- termGAM.fmla
term$fmGLM               <- termGLM.fmla
### This next line is very important.  The use of GAMs for the termination probability is erratic.
### The GAM can be put off what looks physicaly sensible where there is little data, such as long duration HW
### or very high temperatures.  Unless you specificaly need/want GAMs for termination AND have proven the
### GAM is performing as expected I would STRONGLY advise the use of GLMs
### GMST for the UK does not seem to play an important part in the termination probability
term$chosen              <- list(day1=c('GLM','TD'), dayN=c('GLM','TDA'))
HDconfig$term            <- term

############################################################################################################
### No edits needed below here
############################################################################################################

# all setup variables
files <- list(
    # st_infile=paste(MSSAVEDIR,st_infile,sep='/'),
    msfiles=files,
    datestamp=datestamp,
    st_metadata=st_metadata,
    st_HDconfig=st_HDconfig,
    st_base=st_base,
    st_events=st_events,

    st_InitEvent=st_InitEvent,
    st_InitValue=initV$st_InitValue,
    st_Ht=st_Ht,
    st_Term=st_Term
    )
HDconfig$files  <- files
cat("File definitions are:",cr)
print(str(files))

save(file=files$st_HDconfig, HDconfig)

cat(cr)

cat("###############################################################################",cr)
sysdate   <- Sys.time()
stsysdate <- print(sysdate)
cat("All Done setup_HotDays.R",cr)
cat("###############################################################################",cr,cr)
