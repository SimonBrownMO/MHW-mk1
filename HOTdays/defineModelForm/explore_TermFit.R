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

### libraries
library(Rcpp, lib="/home/h03/hadsx/extremes/R/packages")
library(evgam, lib="/home/h03/hadsx/extremes/R/packages")
library(mgcv)
library(qgam)
library(ncdf4)
library(PCICt)
source("../libs/fn_JointMakeStationary.R")
source("../libs/fn_JointHotDays.R")
source("../libs/lib_HotDay.R")
source("/home/h03/hadsx/extremes/R/Rutils/sjb_colours.R")
source("/home/h03/hadsx/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")

# DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2/CPM/"
DATA <- "/data/users/hadsx/extremes/heatwaves/HadUKGrid/dur-clim/CPM5km/v_1/England/01/"
reload_MS_ref(paste0(DATA,"MakeStationary/MSref/ukgd_cpm85_5k_x146y42.MSref.2024-09-26-115727.RData"))
reload_HD_ref(paste0(DATA,"HotDays/vHD2c/HDref/ukgd_cpm85_5k_x146y42.HDref.2024-10-07-151753.RData"))

load(st_msdata01, verb=TRUE)
load(st_events,verb=TRUE)

scale.doy <- list(obs=data01.std.param$doy$o_max, mod=data01.std.param$doy$m_max)
stcovs    <- c('Index', 'CC', 'D', 'A')



###################################################################################
### GAM
###################################################################################
tr.k <- list(D=4, gmst=6, A=6)
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

# fmla1N <- list(day1=termGAM.fmla$hseas$day1$TICDiCD, dayN=termGAM.fmla$hseas$dayN$TICDAiCD)
fmla1N <- termGAM.fmla$hseas
termJgam <- fn_fitJointTermGAM(events01, fmla1N, scale.doy, stdiag='explore-term-diag.gam.txt', stplot='explore-term-diag.gam.pdf', DOPLOT=TRUE)
names(termJgam$day1) <- names(fmla1N$day1)
names(termJgam$dayN) <- names(fmla1N$dayN)
   # > str1(termJgam)
    # List of 2
    #  $ day1:List of 6
    #  $ dayN:List of 8
    #     > str1(termJgam$dayN)
    #     List of 8
    #     $ T       :List of 56
    #     $ TI      :List of 58
    #     $ TIC     :List of 58
    #     $ TID     :List of 58
    #     $ TIA     :List of 58
    #     $ TICD    :List of 58
    #     $ TICDA   :List of 58
    #     $ TICDAiCD:List of 58
### Now go through the stdiag file and see which covariates and model is most appropriate

### plot the effect of covariates on probability of a HW terminating
# fn_termJointGAMPlot(termJgam$dayN$TICDAiCD,main='Day N ~')




###################################################################################
### GLM
###################################################################################
fmla1N       <- list(day1=termGLM.fmla$hseas$day1$TICD, dayN=termGLM.fmla$hseas$dayN$TICDA)

### this automaticaly permutates the covariates and returns the best model.
termJmulti  <- fn_fitAutoMultiJointTermGLM(events01, scale.doy, stdiag='explore-term-diag.multiglm.txt')
termJmulti2 <- fn_fitJointTermGLM(events01, termGLM.fmla$hseas, scale.doy, stdiag='explore-term-diag.multiglm.txt', stplot='explore-term-diag.multiglm.pdf', DOPLOT=TRUE)
### plot, for best model the effect of covariates on probability of a HW terminating
fn_termJointGLMPlot(termJmulti$dayN$term, main0='Day N:')

### this permutates by hand - more flexible
   # calculate all combinations of ICDA
    st.n      <- c('I','C','D','A')
    st.c      <- c('I','gmst','D','A')
    # day 1
    alld1.com <- NULL
    for(i in 1:3) alld1.com[[i]] <- combn(st.c[1:3],i)
    alld1n.com <- NULL
    for(i in 1:3) alld1n.com[[i]] <- combn(st.n[1:3],i)
    fm.glm.d1     <- list()
    fm.glm.d1n    <- NULL
    fm.glm.d1[[1]]  <- 'fin ~ tmp'
    fm.glm.d1n[1] <- 'T'
    count0        <- 2
    for(i in seq_along(alld1.com)) {
        for(j in seq_along(alld1.com[[i]][1,])) {
            fm.glm.d1[[count0]] <- as.formula(paste('fin ~ tmp',paste('+',alld1.com[[i]][,j],collapse=' ', sep='')))

            fm.glm.d1n[count0]  <- paste('T',paste(alld1n.com[[i]][,j],collapse='', sep=''), sep='')
            count0              <- count0 +1
        }
    }
    names(fm.glm.d1) <- c('T','TI','TC','TD','TIC','TID','TCD','TICD')
    # day N
    alldN.com <- NULL
    for(i in 1:4) alldN.com[[i]] <- combn(st.c,i)
    alldNn.com <- NULL
    for(i in 1:4) alldNn.com[[i]] <- combn(st.n,i)
    fm.glm.dN    <- list()
    fm.glm.dNn   <- NULL
    fm.glm.dN[[1]]  <- 'fin ~ tmp'
    fm.glm.dNn[1] <- 'T'
    count0       <- 2
    for(i in seq_along(alldN.com)) {
        for(j in seq_along(alldN.com[[i]][1,])) {
            fm.glm.dN[[count0]] <- as.formula(paste('fin ~ tmp',paste('+',alldN.com[[i]][,j],collapse=' ', sep='')))
            fm.glm.dNn[count0]  <- paste('T',paste(alldNn.com[[i]][,j],collapse='', sep=''), sep='')
            count0              <- count0 +1
        }
    }
    names(fm.glm.dN) <- c('T','TI','TC','TD','TA','TIC','TID','TIA','TCD','TCA','TDA','TICD','TICA','TIDA','TCDA','TICDA')
fm.glm <- list(day1=fm.glm.d1, dayN=fm.glm.dN)
# fm.glm <- list(day1='fin ~ tmp', dayN='fin ~ tmp')

termJglm <- fn_fitJointTermGLM(events01, fm.glm, scale.doy, stdiag='explore-term-diag.glm.txt', stplot='explore-term-diag.glm.pdf',DOPLOT=TRUE)
names(termJglm$day1) <- fm.glm.d1n
names(termJglm$dayN) <- fm.glm.dNn
fn_termJointGLMPlot(termJglm$dayN$TICDA, main0='Day N:')
