# source("explore_TermFit.R")

### Use this file to determine the form of Termination model
### and the covariates required
### Choices are:
        # GAM or GLM.  Use GLM if GAM produces rising probabilities of termination with temperatures
        # Covariates:
            # bias between obs & model (I or class)
            # climate change (CC)
            # day of year (D)
            # A (A)
        # interaction terms between:
            # CC & D

### This programme will try all permutations of covariates for both GLM & GAM
### produce plots and diagnostics

### Experience suggests that it is difficutl to get the GAM to work for terminating heatwaves
### as it always wants to go to 1.0 for the highest temperatures.
### Re D tge GAM struggles to find anything convincing which is a relief as D as a covaraite in GLM is tricky
### Re A GAM has to be calmed with low K otherwise goes mad.

st.pwd <- system("pwd", intern=TRUE)
# source(paste(st.pwd,"/../../setup_all.R",sep=''))
source("setup_dev.R")

DOHOTSEAS <- FALSE

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

############################################################################################################
############################################################################################################


# ### libraries
# library(Rcpp, lib="/home/h03/hadsx/extremes/R/packages")
# library(evgam, lib="/home/h03/hadsx/extremes/R/packages")
# library(mgcv)
# library(qgam)
# library(ncdf4)
# library(PCICt)
# source("../libs/fn_JointMakeStationary.R")
# source("../libs/fn_JointHotDays.R")
# source("../libs/lib_HotDay.R")
# source("/home/h03/hadsx/extremes/R/Rutils/sjb_colours.R")
# source("/home/h03/hadsx/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")

# # DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/"
# DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/CPM5km/v_1/England/01/"
# reload_MS_ref(paste0(DATA,"MakeStationary/MSref/ukgd_cpm85_5k_x146y42.MSref.2024-09-26-115727.RData"))
# reload_HD_ref(paste0(DATA,"HotDays/vHD2c/HDref/ukgd_cpm85_5k_x146y42.HDref.2024-10-07-151753.RData"))

# load(st_msdata01, verb=TRUE)
# load(st_events,verb=TRUE)





scale.doy <- list(obs=(366+1), mod=(360+1)) # bit of a bodge as assumes all obs are from leap years
stcovs    <- c('Index', 'CC', 'D', 'A')

###################################################################################
### GAM
###################################################################################
tr.k         <- list(D=4, gmst=6, A=6)
termGAM.fmla <- list()
# hotseas
 #day1
 day1          <- list()
 day1$T        <- fin ~  s(tmp)
 day1$TI       <- fin ~  s(tmp) +class
 day1$TIC      <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst)
 day1$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class)
 day1$TICDiCD  <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,by=class) +ti(gmst,D)
 day1$TCDiCD   <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)          +ti(gmst,D)
 day1$TCD      <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D)
 day1$TC       <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst)
 day1$TD       <- fin ~  s(tmp)                             +s(D,k=tr.k$D)

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
 dayN$TA       <- fin ~  s(tmp)                                                     +s(A, k=5)
termGAM.fmla$hseas <- list(day1=day1, dayN=dayN)

# all year
 #day1
 day1          <- list()
 day1$T        <- fin ~  s(tmp)
 day1$TI       <- fin ~  s(tmp) +class
 day1$TIC      <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst)
 day1$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class)
 day1$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class)
 day1$TICDiCD  <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +ti(gmst,D,bs=c('tp','cc'))
 day1$TCDiCD   <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +ti(gmst,D)
 day1$TCD      <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')
 day1$TC       <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst)
 day1$TD       <- fin ~  s(tmp)                             +s(D,k=tr.k$D,bs='cc')
 #dayN
 dayN          <- list()
 dayN$T        <- fin ~  s(tmp)
 dayN$TI       <- fin ~  s(tmp) +class
 dayN$TIC      <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst)
 dayN$TICD     <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class)
 dayN$TICDA    <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +s(A, k=5)
 dayN$TA       <- fin ~  s(tmp)                                                             +s(A, k=5)
 dayN$TICDAiCD <- fin ~  s(tmp) +class +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc',by=class) +s(A, k=5) +ti(gmst,D,bs=c('tp','cc'))
 dayN$TCDAiCD  <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +s(A, k=5) +ti(gmst,D)
 dayN$TCDA     <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst) +s(D,k=tr.k$D,bs='cc')          +s(A, k=5)
 dayN$TCA      <- fin ~  s(tmp)        +s(gmst,k=tr.k$gmst)                                 +s(A, k=5)
 dayN$TDA      <- fin ~  s(tmp)                             +s(D,k=tr.k$D,bs='cc')          +s(A, k=5)
 dayN$TA       <- fin ~  s(tmp)                                                             +s(A, k=5)
termGAM.fmla$ally <- list(day1=day1, dayN=dayN)

##############################################################################################
##############################################################################################
    # GAMfmla <- list(day1=termGAM.fmla$hseas$day1$TICDiCD, dayN=termGAM.fmla$hseas$dayN$TICDAiCD)
    GAMfmla    <- termGAM.fmla$hseas
    termJgam.s <- fn_fitJointTermGAM(events01, GAMfmla, scale.doy, stdiag='explore-term-diag.hotseas.gam.txt', stplot='explore-term-diag.hotseas.gam.pdf', DOPLOT=TRUE)
    names(termJgam.s$day1) <- names(GAMfmla$day1)
    names(termJgam.s$dayN) <- names(GAMfmla$dayN)

    GAMfmla    <- termGAM.fmla$ally
    termJgam.a <- fn_fitJointTermGAM(events01$events01_ally, GAMfmla, scale.doy, stdiag='explore-term-diag.ally.gam.txt', stplot='explore-term-diag.ally.gam.pdf', DOPLOT=TRUE)
    names(termJgam.a$day1) <- names(GAMfmla$day1)
    names(termJgam.a$dayN) <- names(GAMfmla$dayN)

    ### test single fit
    GAMfmla     <- termGAM.fmla$hseas
    i1          <- "TICD"   # grep("TICD",  names(GAMfmla$day1))
    i2          <- "TICDA"  # grep("TICDA", names(GAMfmla$dayN))
    GAMfmla1    <- list(day1=GAMfmla$day1[i1], dayN=GAMfmla$dayN[i2])
    termJgam.s1 <- fn_fitJointTermGAM(events01, GAMfmla1, scale.doy, stdiag='explore-term-diag1.hotseas.gam.txt', stplot='explore-term-diag1.hotseas.gam.pdf', DOPLOT=TRUE)

    GAMfmla     <- termGAM.fmla$ally
    i1          <- "TICD"   # grep("TICD",  names(GAMfmla$day1))
    i2          <- "TICDA"  # grep("TICDA", names(GAMfmla$dayN))
    GAMfmla1    <- list(day1=GAMfmla$day1[i1], dayN=GAMfmla$dayN[i2])
    termJgam.a1 <- fn_fitJointTermGAM(events01, GAMfmla1, scale.doy, stdiag='explore-term-diag1.ally.gam.txt', stplot='explore-term-diag1.ally.gam.pdf', DOPLOT=TRUE)


### plot the effect of covariates on probability of a HW terminating
# fn_termJointGAMPlot(termJgam$dayN$TICDA,main='Day N ~')

# 2026.04.02
    # problems with rising probs with high temp seem to have gone away
    # ally
    # - dat-1: T, C, important covariates but D is v weak, use of I questionable.
    # - day-N: T, C, A are important covariates but D is weaker and only seen for Obs, making the use of I questionable.
    # - however, for the weak D effect the GAM is producing very small parameters so wont have much impact on the probabilities
    # - suggest day-1:TCD and day-N:TCDA as the best GAM models.
    # hotseas TCD TCDA
    # pragmatic general approach would be TICD for day-1 and TICDA for day-N, but check smooths are not screwy
readline("Continue to GLM?")


# not sure we need GLM as tmp response seem s ol with GAM
# ###################################################################################
# ### GLM
# ###################################################################################
# GLMfmla    <- list(day1=c(termGLM.fmla$hseas$day1$TCD,termGLM.fmla$hseas$day1$TICD), dayN=c(termGLM.fmla$hseas$dayN$TCDA,termGLM.fmla$hseas$dayN$TICDA))
# termJglm.s <- fn_fitJointTermGLM(events01, GLMfmla, scale.doy, stdiag='explore-term-diag.hotseas.glm.txt', stplot='explore-term-diag.hotseas.glm.pdf', DOPLOT=TRUE)

# GLMfmla    <- list(day1=c(termGLM.fmla$ally$day1$TCD,termGLM.fmla$ally$day1$TICD), dayN=c(termGLM.fmla$ally$dayN$TCDA,termGLM.fmla$ally$dayN$TICDA))
# termJglm.a <- fn_fitJointTermGLM(events01$events01_ally, GLMfmla, scale.doy, stdiag='explore-term-diag.ally.glm.txt', stplot='explore-term-diag.ally.glm.pdf', DOPLOT=TRUE)

# ### this automaticaly permutates the covariates and returns the best model.
# termJmulti  <- fn_fitAutoMultiJointTermGLM(events01, scale.doy, stdiag='explore-term-diag.multiglm.txt')
# termJmulti2 <- fn_fitJointTermGLM(events01, termGLM.fmla$hseas, scale.doy, stdiag='explore-term-diag.multiglm.txt', stplot='explore-term-diag.multiglm.pdf', DOPLOT=TRUE)
# ### plot, for best model the effect of covariates on probability of a HW terminating
# fn_termJointGLMPlot(termJmulti$dayN$term, main0='Day N:')

# ### this permutates by hand - more flexible
#    # calculate all combinations of ICDA
#     st.n      <- c('I','C','D','A')
#     st.c      <- c('I','gmst','D','A')
#     # day 1
#     alld1.com <- NULL
#     for(i in 1:3) alld1.com[[i]] <- combn(st.c[1:3],i)
#     alld1n.com <- NULL
#     for(i in 1:3) alld1n.com[[i]] <- combn(st.n[1:3],i)
#     fm.glm.d1     <- list()
#     fm.glm.d1n    <- NULL
#     fm.glm.d1[[1]]  <- 'fin ~ tmp'
#     fm.glm.d1n[1] <- 'T'
#     count0        <- 2
#     for(i in seq_along(alld1.com)) {
#         for(j in seq_along(alld1.com[[i]][1,])) {
#             fm.glm.d1[[count0]] <- as.formula(paste('fin ~ tmp',paste('+',alld1.com[[i]][,j],collapse=' ', sep='')))

#             fm.glm.d1n[count0]  <- paste('T',paste(alld1n.com[[i]][,j],collapse='', sep=''), sep='')
#             count0              <- count0 +1
#         }
#     }
#     names(fm.glm.d1) <- c('T','TI','TC','TD','TIC','TID','TCD','TICD')
#     # day N
#     alldN.com <- NULL
#     for(i in 1:4) alldN.com[[i]] <- combn(st.c,i)
#     alldNn.com <- NULL
#     for(i in 1:4) alldNn.com[[i]] <- combn(st.n,i)
#     fm.glm.dN    <- list()
#     fm.glm.dNn   <- NULL
#     fm.glm.dN[[1]]  <- 'fin ~ tmp'
#     fm.glm.dNn[1] <- 'T'
#     count0       <- 2
#     for(i in seq_along(alldN.com)) {
#         for(j in seq_along(alldN.com[[i]][1,])) {
#             fm.glm.dN[[count0]] <- as.formula(paste('fin ~ tmp',paste('+',alldN.com[[i]][,j],collapse=' ', sep='')))
#             fm.glm.dNn[count0]  <- paste('T',paste(alldNn.com[[i]][,j],collapse='', sep=''), sep='')
#             count0              <- count0 +1
#         }
#     }
#     names(fm.glm.dN) <- c('T','TI','TC','TD','TA','TIC','TID','TIA','TCD','TCA','TDA','TICD','TICA','TIDA','TCDA','TICDA')
# fm.glm <- list(day1=fm.glm.d1, dayN=fm.glm.dN)
# # fm.glm <- list(day1='fin ~ tmp', dayN='fin ~ tmp')

# termJglm <- fn_fitJointTermGLM(events01, fm.glm, scale.doy, stdiag='explore-term-diag.glm.txt', stplot='explore-term-diag.glm.pdf',DOPLOT=TRUE)
# names(termJglm$day1) <- fm.glm.d1n
# names(termJglm$dayN) <- fm.glm.dNn
# fn_termJointGLMPlot(termJglm$dayN$TICDA, main0='Day N:')
